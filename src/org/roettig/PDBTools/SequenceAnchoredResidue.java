package org.roettig.PDBTools;

import java.io.Serializable;

import org.biojava.bio.structure.*;

public class SequenceAnchoredResidue implements Serializable
{
    
    
   private static final long serialVersionUID = -6484969438097777417L;
   
   public SequenceAnchoredResidue(Group g, int sI, int pI)
   {
      //group = g; 
      seqIdx   = sI; pdbIdx = pI;
      pdbCode  = g.getPDBCode();
      pdbName  = g.getPDBName();
      pdbChain = g.getParent().getName();
   }
   
   public String toString()
   {
       return pdbChain+" "+pdbCode+" "+pdbName+" [pdbIdx:"+pdbIdx+",seqIdx:"+seqIdx+"]";
   }
   
   //public Group  group;
   public String pdbChain;
   public String pdbCode;
   public String pdbName;
   public int    seqIdx;
   public int    pdbIdx;
}