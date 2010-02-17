package org.roettig.PDBTools;

import java.io.Serializable;

public class ResidueLocator implements Serializable
{
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
}
