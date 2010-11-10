package org.roettig.PDBTools;

import java.io.Serializable;

public class Vector3 implements Serializable
{

	private static final long serialVersionUID = 6167184707451296274L;

	public Vector3()
	{

	}

	public Vector3(double _x, double _y, double _z)
	{
		x = _x; y= _y; z = _z;
	}

	public double distance(Vector3 v)
	{
		return Math.sqrt(Math.pow(v.x-x,2)+Math.pow(v.y-y,2)+Math.pow(v.z-z,2));
	}

	public double x;
	public double y;
	public double z;
}