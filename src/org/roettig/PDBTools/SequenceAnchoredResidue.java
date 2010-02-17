package org.roettig.PDBTools;

import org.biojava.bio.structure.*;

public class SequenceAnchoredResidue
{
   public SequenceAnchoredResidue(Group g, int sI, int pI)
   {
      group = g; seqIdx = sI; pdbIdx = pI;    
   }
   
   public String toString()
   {
       return group.getPDBCode()+" "+group.getPDBName()+" ["+pdbIdx+","+seqIdx+"]";
   }
   
   public Group  group;
   public int    seqIdx;
   public int    pdbIdx;
}