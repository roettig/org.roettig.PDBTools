package org.roettig.PDBTools;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import org.biojava.bio.seq.Sequence;
import java.util.*;
import java.util.logging.Logger;
import org.biojava.bio.*;
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.StructureIOFile;
import org.roettig.SequenceTools.PairwiseAlignment;
import org.roettig.SequenceTools.SeqTools;
import org.biojava.bio.symbol.SimpleSymbolList;

public class PDBTools
{
    
    public static Logger logger = Logger.getLogger("org.roettig.PDBTools");
    
    public static Vector<Sequence> getSequences()
    {
        Vector<Sequence> ret = new Vector<Sequence>();
        return ret;
    }
    
    public static boolean hasGroup(Structure struc, ResidueLocator rl)
    {
	List<Chain> chains    = struc.getChains(0);
        boolean found = false;
        
        for(Chain ch: chains)
        {
          if(!ch.getName().equals(rl.chain))
            continue;
          List<Group> agr = ch.getAtomGroups();
          for(Group g: agr)
          {
             if( g.getPDBName().equals(rl.name) && Integer.parseInt(g.getPDBCode())==rl.idx)
             {
                found = true;
                logger.info("found requested group "+rl);
             }
          }
        }
        if(!found)
        {
            logger.info("did not find requested group "+rl);
        }
      return found;
    }
    
    private static String getElement(Atom at)
    {
	if(at.getName().contains("C"))
	    return "C";
	if(at.getName().contains("N"))
	    return "N";
	if(at.getName().contains("O"))
	    return "O";
	if(at.getName().contains("H"))
	    return "H";
	if(at.getName().contains("P"))
	    return "P";
	if(at.getName().contains("S"))
	    return "S";
	return " ";
    }
    
    public static void writeAtomData(PrintWriter out, List<Group> groups)
    {
	
	long idx = 1;
	
	for(Group g: groups)
	{
	    String Record_ID = null;
	    if(g instanceof AminoAcid)
		Record_ID = "ATOM";
	    else
		Record_ID = "HETATM";
	    for(Atom at: g.getAtoms())
	    {
		//System.out.println(String.format("%5d",idx));
		//                                         "%5d %-4.4s%c%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4.4s%2.2s%2.2s"
		out.write(String.format(Locale.ENGLISH,"%-6s%5d %-4.4s%s%3.3s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4.4s%2.2s%2.2s\n",
			Record_ID,
			idx,
			at.getName(),
			at.getAltLoc(),
			g.getPDBName(),
			g.getParent().getName(),
			Integer.parseInt(g.getPDBCode()),
			" ",
			at.getX(),at.getY(),at.getZ(),1.0,0.0," ",getElement(at)," ")
			);
		//out.write(String.format("%-6s%5d %-4.4s %3.3s\n",Record_ID,idx,at.getName(),g.getPDBName()));
		//out.write(String.format(FORMAT_ATOM,idx,Record_ID,idx,' ',g.getPDBName(),at.getName()," ",g.getPDBName()," ",g.getParent().getName(),g.getPDBCode(),"   ",at.getX(),at.getY(),at.getZ(),0.0,0.0," ","   "));
		//out.write(at.getPDBline());
		idx++;
		//System.out.println(at.getFullName());
	    }
	    
	}
	
    }
    
    public static Vector<Group> getGroups(Structure struc, List<ResidueLocator> centerlocators)
    {
        Vector<Group> centers = new Vector<Group>();
        for(ResidueLocator rl: centerlocators)
        {
          List<Chain> chains    = struc.getChains(0);
          boolean found = false;
          
          for(Chain ch: chains)
          {
            if(!ch.getName().equals(rl.chain))
              continue;
            List<Group> agr = ch.getAtomGroups();
            for(Group g: agr)
            {
               if( g.getPDBName().equals(rl.name) && Integer.parseInt(g.getPDBCode())==rl.idx)
               {
                  centers.addElement( g );
                  found = true;
                  logger.info("found requested group "+rl);
               }
            }
          }
          if(!found)
          {
              logger.info("did not find requested group "+rl);
          }
        }
        return centers;
    }
    
    public static List<Group> getNearbyGroups(List<Group> x, List<Group> y, double cutoff)
    {
        List<Group> ret = new Vector<Group>();
        for(Group g1: x)
        {
          double minD = 10000;  
          for(Atom a: g1.getAtoms())
          {
              for(Group g2: y)
              {
                for(Atom b: g2.getAtoms())
                {
                    
                 double D = 10000;
                 try
                 {
                    D = Calc.getDistance(a,b);
                 } 
                 catch (StructureException e)
                 {
                       e.printStackTrace();
                 }
                 if(minD>D)
                    minD = D;
                }
              }
          }
          if(minD<cutoff)
              ret.add( g1 );
        }    
        return ret;
    }

    public static List<Group> getNearbyGroups2(List<Group> x, List<Vector3> y, double cutoff)
    {
        List<Group> ret = new Vector<Group>();
        for(Group g1: x)
        {
          double minD = 10000;  
          for(Atom a: g1.getAtoms())
          {
              Vector3 va = new Vector3(a.getX(),a.getY(),a.getZ());
              
              for(Vector3 v: y)
              {
                double D = va.distance(v);                 
                if(minD>D)
                    minD = D;
              }
             
          }
          if(minD<cutoff)
              ret.add( g1 );
         
        }    
        return ret;
    }
    
    public static List<SequenceAnchoredResidue> getASCResidues(Structure struc, List<ResidueLocator> centerlocators, List<Vector3> positions, double cutoff)
    {
        List<Group> centers     = getGroups(struc, centerlocators);
        
        List<Chain>    chains     = struc.getChains(0);
        Vector<SequenceAnchoredResidue> ascresidues = new Vector<SequenceAnchoredResidue>();
        
        Chain chain1       = chains.get(0);
        HashMap<Integer,Integer> seqIDX2pdbIDX = new HashMap<Integer,Integer>();
        HashMap<Integer,Integer> pdbIDX2seqIDX = new HashMap<Integer,Integer>();
        
        Sequence seqChain1 = getSequence(chain1,seqIDX2pdbIDX,pdbIDX2seqIDX);
        
        PairwiseAlignment pwa = new PairwiseAlignment();
        
        for(Chain c: chains)
        {
            Sequence seqChain = getSequence(c,null,null);
            
            double pid = 0.0;
            try
            {
              pid = pwa.align(seqChain1, seqChain);
            }
            catch( IllegalSymbolException e)
            {
                e.printStackTrace();
            }
            
            if(pid==1.0)
            {
                Vector<String> allowed = new Vector<String>();
                allowed.addElement("MSE");
                allowed.addElement("175");
                logger.info("examining chain "+c.getName()+" (which is identical to first chain)");
                
                List<Group> groups  = getGroups( c , allowed);
                List<Group> vinc    = getNearbyGroups(groups,centers, cutoff);
                List<Group> vinc2   = getNearbyGroups2(groups,positions, cutoff);
                HashSet<Group> collected = new HashSet<Group>();
                
                for(Group g: vinc)
                {
                    //AminoAcid aa = (AminoAcid) g;
                    Integer pI = Integer.parseInt(g.getPDBCode());
                    Integer pS = pdbIDX2seqIDX.get(pI);
                    SequenceAnchoredResidue ascr = new SequenceAnchoredResidue(g, pS, pI );
                    ascresidues.addElement(ascr);
                    collected.add(g);
                }
                for(Group g: vinc2)
                {
                    if(collected.contains(g))
                         continue;
                    //AminoAcid aa = (AminoAcid) g;
                    Integer pI = Integer.parseInt(g.getPDBCode());
                    Integer pS = pdbIDX2seqIDX.get(pI);
                    SequenceAnchoredResidue ascr = new SequenceAnchoredResidue(g, pS, pI );
                    ascresidues.addElement(ascr);
                }                
            }
        }
        return ascresidues;
    }
    
    public static List<Group> getASCGroups(Structure struc, List<ResidueLocator> centerlocators, List<Vector3> positions, double cutoff)
    {
        List<Group> centers     = getGroups(struc, centerlocators);
        
        List<Chain>    chains     = struc.getChains(0);
        List<Group> ascresidues = new Vector<Group>();
        
        Chain chain1       = chains.get(0);
        HashMap<Integer,Integer> seqIDX2pdbIDX = new HashMap<Integer,Integer>();
        HashMap<Integer,Integer> pdbIDX2seqIDX = new HashMap<Integer,Integer>();
        
        Sequence seqChain1 = getSequence(chain1,seqIDX2pdbIDX,pdbIDX2seqIDX);
        
        PairwiseAlignment pwa = new PairwiseAlignment();
        
        for(Chain c: chains)
        {
            Sequence seqChain = getSequence(c,null,null);
            
            double pid = 0.0;
            try
            {
              pid = pwa.align(seqChain1, seqChain);
            }
            catch( IllegalSymbolException e)
            {
                e.printStackTrace();
            }
            
            if(pid==1.0)
            {
                Vector<String> allowed = new Vector<String>();
                allowed.addElement("MSE");
                
                logger.info("examining chain "+c.getName()+" (which is identical to first chain)");
                
                List<Group> groups  = getGroups( c , allowed);
                List<Group> vinc    = getNearbyGroups(groups,centers, cutoff);
                List<Group> vinc2   = getNearbyGroups2(groups,positions, cutoff);
                HashSet<Group> collected = new HashSet<Group>();
                
                for(Group g: vinc)
                {
                    collected.add(g);
                    ascresidues.add(g);
                }
                for(Group g: vinc2)
                {
                    if(collected.contains(g))
                         continue;
                    ascresidues.add(g);
                }                
            }
        }
        return ascresidues;
    }
    
    public static Sequence getSequence(Chain cha, Map<Integer,Integer> seqIDX2pdbIDX, Map<Integer,Integer> pdbIDX2seqIDX)
    {        
        String seq = "";
   
        List<Group> agr = cha.getSeqResGroups();
        if(agr.size()==0)
            agr = cha.getAtomGroups();
        for(Group g: agr)
        {
            if(!g.has3D())
            {
                seq += "X";
                continue;
            }
            if ( g instanceof AminoAcid )
            {
                AminoAcid aa = (AminoAcid) g;
                seq += SeqTools.ThreeLetterToShort(aa.getPDBName());
                if(seqIDX2pdbIDX!=null)
                   seqIDX2pdbIDX.put(seq.length(), Integer.parseInt(aa.getPDBCode()));
                if(pdbIDX2seqIDX!=null)
                    pdbIDX2seqIDX.put(Integer.parseInt(aa.getPDBCode()),seq.length());
            }
            else
            {
              if ( g instanceof HetatomImpl )
              {
                HetatomImpl het = (HetatomImpl)g;
                if(het.getPDBName().equals("MSE"))
                  seq += "X";
                if(seqIDX2pdbIDX!=null)
                   seqIDX2pdbIDX.put(seq.length(), Integer.parseInt(het.getPDBCode()));
                if(pdbIDX2seqIDX!=null)
                    pdbIDX2seqIDX.put(Integer.parseInt(het.getPDBCode()),seq.length());
              }
            } 
            
        } 
        
        Sequence pdbseq = null;
        try
        {
            FiniteAlphabet fa = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN");
            SymbolTokenization p = null;
	    try
	    {
		p = fa.getTokenization("token");
	    } 
	    catch (BioException e)
	    {
		e.printStackTrace();
	    }
            SymbolList sl = new SimpleSymbolList(p, seq);
            // make Sequence from SymbolList
            pdbseq = new SimpleSequenceFactory().createSequence(sl,"", "pdb", new SimpleAnnotation());
        } 
        catch (IllegalSymbolException e)
        {
            e.printStackTrace();
        } 
         
        
        return pdbseq;
        
    }
    
    public static Vector<Group> getGroups(Chain c)
    {
        Vector<Group> aas = new Vector<Group>();
        List<Group>   agr = c.getAtomGroups("amino");
        for(Group g: agr)
        {
            if(!g.has3D())
                continue;
            aas.addElement( g );   
        }
        return aas;
    }
    
    public static Vector<Group> getGroups(Chain c, Vector<String> allowed)
    {
        Vector<Group> aas = new Vector<Group>();
        List<Group> agr = c.getAtomGroups();
        for(Group g: agr)
        {
            if(!g.has3D())
                continue;
            if(g.getType()=="amino")
            {   
                aas.addElement( g );
                continue;
            }
            
            if( allowed.contains(g.getPDBName()) )
            {
                aas.addElement( g );
                continue;
            }
            
        }
        return aas;
    }
    
    public static void main(String args[]) throws StructureException, FileNotFoundException
    {
	PDBFileReader pdbreader = new PDBFileReader();
        
        Structure struc = null;
        try
        {
         struc = pdbreader.getStructure("/tmp/pdbfile");
        }
        catch(Exception e)
        {
         
        }
        Group g = struc.findGroup("A","234");
        System.out.println(g);
        List<Group> gs = new ArrayList<Group>();
        gs.add(g);
        gs.add(g);
        PrintWriter raus = new PrintWriter("/tmp/raus");
        writeAtomData(raus,gs);
        raus.close();
    }
}
