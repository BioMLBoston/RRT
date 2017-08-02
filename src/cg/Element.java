package cg;

import java.util.ArrayList;
import org.jblas.DoubleMatrix;
import org.jblas.Singular;

import pdb.PDBAtom;
import pdb.PDBMolecule;
import utilities.Vector;
import cg.GordonEnergy.DihedralTypes;

/**
 * @author Dong Luo 2012-06-12
 *
 */

public class Element {
	
	public static int minimumElementSize = 3;
	
	private CGMolecule parent = null;
	private DihedralTypes type;
	private int start, end;
	private Vector vector = null, center = null;
	
	public Element(CGMolecule parent, DihedralTypes type, int start, int end) {
		this.parent = parent;
		this.type = type; this.start = start; this.end = end;
	}

	/**
	 * 
	 * @param parent is the new molecule that contains the new element
	 * @param e is the old element to be copied
	 */
	public Element(CGMolecule parent, Element e) {
		if (e == null) return;
		this.parent = parent;
		type = e.type; start = e.start; end = e.end;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	public DihedralTypes getType() {
		return type;
	}
	
	public CGMolecule getParentMolecule() {
		return parent;
	}
	
	public Vector getCenter() {
		if ((center == null) && (parent != null))
			center = parent.getCentroid(start, end);
		return center;
	}
	
	/*
	 * get representative vector for this element
	 * for helix, it is the best line fit the points
	 * for sheet, it is the line normal to the plane
	 * algorithm from http://mathforum.org/library/drmath/view/69103.html
	 */
	private Vector getVector() {
		if (vector != null) return vector;
		if (parent == null) return null;
		center = getCenter();
		if (center == null) return null;
		ArrayList<PDBAtom> l = parent.getPDBAtomList();
		if (l == null) return null;
		int size = l.size();
		if (start<0 || start>(size-1)) return null;
		if (end<0 || end>(size-1) || end<start) return null;
		DoubleMatrix M = DoubleMatrix.zeros(end-start+1, 3);
		double x0 = center.getX(), y0 = center.getY(), z0 = center.getZ();
		for (int i=0;i<=(end-start);i++) {
			PDBAtom a = l.get(start+i);
			M.put(i, 0, a.getX()-x0);
			M.put(i, 1, a.getY()-y0);
			M.put(i, 2, a.getZ()-z0);
		}
		DoubleMatrix[] USV = Singular.fullSVD(M);
		DoubleMatrix V = USV[2];
		if (type == DihedralTypes.H)
			vector = new Vector(V.get(0, 0), V.get(1, 0), V.get(2, 0));
		else
			vector = new Vector(V.get(0, 2), V.get(1, 2), V.get(2, 2));
		return vector;
	}

	public void inValidateVectors() {
		center = null; vector = null;
	}
	
	@Override
	public String toString() {
		int offset = parent.getPDBAtomList().get(0).getResidueIndex();
		return (start+offset)+"-"+(end+offset)+" "+type;
	}
	
	public static double Angle(Element e1, Element e2) {
		if (e1 == null || e2 == null) return 0;
		Vector v1 = e1.getVector(), v2 = e2.getVector();
		if (v1 == null || v2 == null) return 0;
		return v1.angle(v2);
	}
	
	public static double distance(Element e1, Element e2) {
		if (e1 == null || e2 == null) return 0;
		Vector c1 = e1.getCenter(), c2 = e2.getCenter();
		if (c1 == null || c2 == null) return 0;
		return c1.distance(c2);
	}

	/*
	 * m1 and m2 are assumed to be the same molecule in different conformation
	 */
	@SuppressWarnings("unchecked")
	public static void setSharedElements(CGMolecule m1, CGMolecule m2) {
		if (m1 == null || m2 == null) return;
		ArrayList<Element> l1 = m1.getElementList(), l2 = m2.getElementList();
		if (l1 == null || l2 == null) return;
		ArrayList<Element> nl1 = new ArrayList<Element>(), nl2 = new ArrayList<Element>();
		Element e1, e2;
		DihedralTypes t1, t2;
		int beg1, end1, beg2, end2;
		int i = 0, j= 0;
		while (i < l1.size()) {
			e1 = l1.get(i);
			t1 = e1.getType();beg1 = e1.getStart();end1 = e1.getEnd();
			while (j < l2.size()) {
				e2 = l2.get(j);
				t2 = e2.getType();beg2 = e2.getStart();end2 = e2.getEnd();
				if (beg2 >= end1) break;
				if (end2 <= beg1) {j++;continue;}
				int b = java.lang.Math.max(beg1, beg2);
				int e = java.lang.Math.min(end1, end2);
				if ((t1 == t2) && ((e-b) >= minimumElementSize)) {
					nl1.add(new Element(m1, t1, b, e));
					nl2.add(new Element(m2, t2, b, e));
				}
				if (e == end1) break;
				if (e == end2) {j++;continue;}
				j++;	// will not reach here, but keep for safety
			}
			i++;
		}
		m1.setElements(nl1);
		m2.setElements(nl2);
	}
	
}
