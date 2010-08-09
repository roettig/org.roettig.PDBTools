package org.roettig.PDBTools;

import java.io.Serializable;

public class ResidueLocator implements Serializable
{
  private static final long serialVersionUID = 7426685724689835032L;
  
  public static ResidueLocator fromString(String s)
  {
      if(s.startsWith("[")&&s.endsWith("]"))
	  s = s.substring(1,s.length()-1);
      String toks[] = s.split("\\|");
      return new ResidueLocator(toks[0],toks[1],Integer.parseInt(toks[2]));
  }
  
  public ResidueLocator(String _chain, String _name, int _idx)
  {
     chain = _chain;
     name  = _name;
     idx   = _idx;
  }
  
  public String toString()
  {
      return "["+chain+"|"+name+"|"+idx+"]";
  }
  
  public String chain;
  public String name;
  public int idx;
  
  public static void main(String args[])
  {
     System.out.println(ResidueLocator.fromString("[A|PHE|550]").toString());
  }
  
}
