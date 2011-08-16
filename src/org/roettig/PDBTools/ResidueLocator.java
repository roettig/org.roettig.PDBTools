package org.roettig.PDBTools;

import java.io.Serializable;

/**
 * The ResidueLocator class stores the three identifiers (chain, name, index) needed to
 * uniquely identify a residue within a PDB structure.
 * 
 * @author roettig
 *
 */
public class ResidueLocator implements Serializable
{
	private static final long serialVersionUID = 7426685724689835032L;

	/**
	 * factory method to construct a residue locator object from a canonical string representation.
	 * 
	 * Example: [A,ASP,22]
	 * 
	 * @param s
	 * @return residue locator object
	 */
	public static ResidueLocator fromString(String s)
	{
		if(s.startsWith("[")&&s.endsWith("]"))
			s = s.substring(1,s.length()-1);
		String toks[] = s.split("\\|");
		return new ResidueLocator(toks[0],toks[1],Integer.parseInt(toks[2]));
	}

	/**
	 * construct a residue locator object.
	 * 
	 * @param _chain name of the parent chain (e.g. A)
	 * @param _name name of the residue (e.g. ASP)
	 * @param _idx index of the residue (e.g. 22)
	 */
	public ResidueLocator(String _chain, String _name, int _idx)
	{
		chain = _chain;
		name  = _name;
		idx   = _idx;
	}

	/**
	 * returns a canonical string representation of the residue locator.
	 */
	public String toString()
	{
		return "["+chain+"|"+name+"|"+idx+"]";
	}

	public String chain;
	public String name;
	public int idx;

}
