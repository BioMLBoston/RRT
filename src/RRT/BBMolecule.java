/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package RRT;


import java.util.ArrayList;
import java.util.ListIterator;
import pdb.PDBAtom;
import pdb.PDBMolecule;
import utilities.Vector;
import RRT.BBEnergy;
import org.netlib.lapack.Dgesvd;

/**
 *
 * @author amir
 */
public class BBMolecule extends PDBMolecule {
	private ArrayList<Integer> freeDihedralAngleList=null;	// free dihedral angle indexes, depend on elements
	
        
	
    
	private BBMolecule(ArrayList<PDBAtom> list, String chainID) {
		super(list, chainID);
	}

	
        public static BBMolecule newBBMolecule(PDBMolecule molecule) {
                
		ArrayList<PDBAtom> list = new ArrayList<PDBAtom>();
                
                  ListIterator<PDBAtom> iter;
                  iter = molecule.extractBackbone().listIterator();
              
		while (iter.hasNext()) {
			PDBAtom a = iter.next();
			if (a != null) {
				String xyzStr = a.getCoordinatesString();
				double x = Double.parseDouble(xyzStr.substring(0, 7).trim());
				double y = Double.parseDouble(xyzStr.substring(8, 15).trim());
				double z = Double.parseDouble(xyzStr.substring(16, 23).trim());
				list.add(new PDBAtom(a.getIndex(), x, y, z, a.getAtomType(), a.getAminoAcidType(), a.getChainID(), a.getResidueIndex()));
			}
		}
                BBMolecule inst = new BBMolecule(list, molecule.getChainID());
		if (molecule instanceof BBMolecule) {
			BBMolecule old = (BBMolecule)molecule;
		
                        if (old.freeDihedralAngleList != null)
				inst.freeDihedralAngleList = new ArrayList<Integer>(old.freeDihedralAngleList);
				
		}
                
		return inst;
		
        }
   
        
	
	public void setFreeDihedralAngleList(ArrayList<Integer> list) {
		freeDihedralAngleList = list;
	}

        public ArrayList<Integer> getDihedralAngleList()
        {
            return this.freeDihedralAngleList;
        }
	/**
	 * assume m1 and m2 are same molecule but in different configurations
	 * This function will set the same rigid elements for m1 and m2
	 * and their freeAngleList, freeDihedralAngleList as those differed larger than threshold value
	 * @param m1
	 * @param m2
	 * @return the summation of delta of angles
	 */
        
        public double getDihedralAngle(int index) {
		ArrayList<PDBAtom> l = getPDBAtomList();
		return BBMolecule.calculateDihedralAngle(l.get(index-1), l.get(index), l.get(index+1), l.get(index+2));
	}
	public static ArrayList<Double> setElementsAndFreeAngles(BBMolecule m1, BBMolecule m2, double threshold) {
		if (m1 == null || m2 == null) return null;
	
                ArrayList<PDBAtom> l1 = m1.getPDBAtomList();
		ArrayList<PDBAtom> l2 = m2.getPDBAtomList();
		if (l1 == null || l2 == null) return null;
		if (l1.size()<4 || l2.size()<4) return null;
		double angleDeltaThreshold = threshold*java.lang.Math.PI/180;
		ArrayList<Double>deltaAngle = new ArrayList<Double>();
		double sum = 0;
		ArrayList<Integer>  dal1 = new ArrayList<Integer>();
		ArrayList<Integer> dal2 = new ArrayList<Integer>();
		PDBAtom a1, a2, a3, a4, b1, b2, b3, b4;
	
                int nres=m1.extractBackbone().size();
		for (int i=1;i<nres-2;i++) {
			a1=l1.get(i-1);a2=l1.get(i);a3=l1.get(i+1);a4=l1.get(i+2);
			b1=l2.get(i-1);b2=l2.get(i);b3=l2.get(i+1);b4=l2.get(i+2);
			double angle1 = calculateDihedralAngle(a1, a2, a3, a4);
			double angle2 = calculateDihedralAngle(b1, b2, b3, b4);
			double d = java.lang.Math.abs(angle1-angle2);
			if (d>java.lang.Math.PI) d = java.lang.Math.PI*2 - d;
			if (d>angleDeltaThreshold) {
				dal1.add(i);
				dal2.add(i);
				sum += d;
				deltaAngle.add(sum);
			}
		}
		m1.setFreeDihedralAngleList(dal1);
		m2.setFreeDihedralAngleList(dal2);
		return deltaAngle;
	}
        
        
        	
	/**
	 * get score vector relative to given molecule
	 * @param reference
	 * @return
	 */

        
	public static double calculateDihedralAngle(PDBAtom a1, PDBAtom a2, PDBAtom a3, PDBAtom a4) {
		Vector v1v2 = a2.subtract(a1);
		Vector v2v3 = a3.subtract(a2);
		Vector v3v4 = a4.subtract(a3);

		Vector surface1Normal = v2v3.crossProduct(v3v4);
		Vector surface2Normal = v1v2.crossProduct(v2v3);

		//return surface2Normal.angle(surface1Normal);

		// The sign of the dihedral angle is determined by the angle between 
		// the normal to the plane v1v2 and the vector v3v4 (which is v3-v4). 
		// If the angle is < 90 degrees (or the dot product between the two 
		// vectors is negative), reverse the sign of the dihedral angle.
		double angle = surface2Normal.angle(surface1Normal);
		double sign = surface2Normal.dotProduct(v3v4);
		if (sign < 0)
		return (-1*angle);
		else 
		return angle;
		
	}
        	
        
        
        public double getRMSD(PDBMolecule otherProtein)
	{
		
            
		double[] x = convertToCoordinateArray(this);
		double[] y = convertToCoordinateArray(otherProtein);
		
		int xLength = x.length;
		int yLength = y.length;
		
		double[] comx = new double[3]; // center of mass of x
		double[] comy = new double[3]; // center of mass of y
		
		int i; 
		int n = x.length;
		
		// compute centers of mass
		comx[0] = comx[1] = comx[2] = comy[0] = comy[1] = comy[2] = 0.0;
	    for (i = 0; i < n; i = i+3) {
	    	comx[0] += x[i]; 
	    	comx[1] += x[i+1]; 
	    	comx[2] += x[i+2];
	    	comy[0] += y[i]; 
	    	comy[1] += y[i+1]; 
	    	comy[2] += y[i+2];
	    }
	    comx[0] *= 3./n; 
	    comx[1] *= 3./n; 
	    comx[2] *= 3./n;
	    comy[0] *= 3./n; 
	    comy[1] *= 3./n; comy[2] *= 3./n;
	    
	    // compute covariance matrix
	    double[] C = new double[9];
	    double[] v = new double[n];
	    double[] w = new double[n];
	    for (i=0; i<n; i+=3) 
	    {
	    	v[i]   = x[i] - comx[0]; 
	    	v[i+1] = x[i+1] - comx[1]; 
	    	v[i+2] = x[i+2] - comx[2];
	    	w[i]   = y[i] - comy[0]; 
	    	w[i+1] = y[i+1] - comy[1]; 
	    	w[i+2] = y[i+2] - comy[2];
	    	C[0] += v[i] * w[i];   
	    	C[1] += v[i] * w[i+1];   
	    	C[2] += v[i] * w[i+2];
	    	C[3] += v[i+1] * w[i]; 
	    	C[4] += v[i+1] * w[i+1]; 
	    	C[5] += v[i+1] * w[i+2];
	    	C[6] += v[i+2] * w[i]; 
	    	C[7] += v[i+2] * w[i+1]; 
	    	C[8] += v[i+2] * w[i+2];
	    }
	    
	    // compute SVD of C
	    double[] S  = new double[3];
	    double[] U  = new double[9];
	    double[] VT = new double[9];
	    double[] work = new double[30];
	    org.netlib.util.intW info = new org.netlib.util.intW(0);
	    //Dgesvd.dgesvd("A","A",3,3,C,3,S,U,3,VT,3,work,lwork,info);
	    Dgesvd.dgesvd("A","A",3,3,C,0,3,S,0,U,0,3,VT,0,3,work,0,work.length,info);
	    
	    // compute rotation: rot=U*VT
	    double[] rot = new double[9];
	    rot[0] = U[0]*VT[0] + U[3]*VT[1] + U[6]*VT[2];
	    rot[1] = U[1]*VT[0] + U[4]*VT[1] + U[7]*VT[2];
	    rot[2] = U[2]*VT[0] + U[5]*VT[1] + U[8]*VT[2];
	    rot[3] = U[0]*VT[3] + U[3]*VT[4] + U[6]*VT[5];
	    rot[4] = U[1]*VT[3] + U[4]*VT[4] + U[7]*VT[5];
	    rot[5] = U[2]*VT[3] + U[5]*VT[4] + U[8]*VT[5];
	    rot[6] = U[0]*VT[6] + U[3]*VT[7] + U[6]*VT[8];
	    rot[7] = U[1]*VT[6] + U[4]*VT[7] + U[7]*VT[8];
	    rot[8] = U[2]*VT[6] + U[5]*VT[7] + U[8]*VT[8];
	    
	    // make sure rot is a proper rotation, check determinant
	    if ((rot[1]*rot[5]-rot[2]*rot[4])*rot[6]
	        + (rot[2]*rot[3]-rot[0]*rot[5])*rot[7]
	        + (rot[0]*rot[4]-rot[1]*rot[3])*rot[8] < 0) {
	      rot[0] -= 2*U[6]*VT[2]; rot[1] -= 2*U[7]*VT[2]; rot[2] -= 2*U[8]*VT[2];
	      rot[3] -= 2*U[6]*VT[5]; rot[4] -= 2*U[7]*VT[5]; rot[5] -= 2*U[8]*VT[5];
	      rot[6] -= 2*U[6]*VT[8]; rot[7] -= 2*U[7]*VT[8]; rot[8] -= 2*U[8]*VT[8];
	    }
	    
	    // compute rmsd
	    double dist = 0;
	    for (i = 0; i < n; i = i+3) {
	      v[i]   -= rot[0]*w[i] + rot[1]*w[i+1] + rot[2]*w[i+2];
	      v[i+1] -= rot[3]*w[i] + rot[4]*w[i+1] + rot[5]*w[i+2];
	      v[i+2] -= rot[6]*w[i] + rot[7]*w[i+1] + rot[8]*w[i+2];
	      dist += v[i]*v[i] + v[i+1]*v[i+1] + v[i+2]*v[i+2];
	      // rsmd too large => return INF
	    }
	    
	    dist = Math.sqrt(dist*3.0/n);
	    
	    return dist;
	}
	
	private static double[] convertToCoordinateArray(PDBMolecule protein)
	{
		// create the coordinates array
		// note that the size of the coordinates array is 3 times 
		// the number of atoms in the complex because the position of each
		// atom is represented as 3 double value (x, y, z coordinates)
		double[] coordinates = new double[3 * protein.getPDBAtomList().size()];
		//ArrayList<PDBMolecule> molecules = protein.getMolecules();
		//ListIterator<PDBMolecule> moleculeIter = molecules.listIterator();
		
		int i = 0;
		
				ArrayList<PDBAtom> atoms = protein.getPDBAtomList();
				ListIterator<PDBAtom> atomIter = atoms.listIterator();
				while(atomIter.hasNext())
				{
					PDBAtom atom = atomIter.next();
					coordinates[i] = atom.getX();
					coordinates[i+1] = atom.getY();
					coordinates[i+2] = atom.getZ();
					
					i = i+3;
				}
			
		
		
		return coordinates;
	}
        @Override
        public void translate(double deltaX, double deltaY, double deltaZ)
	{
		super.translate(deltaX, deltaY, deltaZ);
		//for (Element e:getElementList()) e.inValidateVectors();
	}

	/**
	 * This function is used to correct the distance between two adjacent atoms
	 * caused by numeric error accumulated during sampling
	 * @param index is the index of first atom
	 * @param delta is the distance difference to be corrected
	 */
	public void translate(int index, double delta) {
		ArrayList<PDBAtom> l = getPDBAtomList();
		PDBAtom a1 = l.get(index), a2 = l.get(index+1);
		double ratio = delta/a1.distance(a2);
		double deltaX = (a2.getX()-a1.getX())*ratio;
		double deltaY = (a2.getY()-a1.getY())*ratio;
		double deltaZ = (a2.getZ()-a1.getZ())*ratio;
		for (int i=index+1;i<l.size();i++) l.get(i).translate(deltaX, deltaY, deltaZ);
		//for (Element e:getElementList()) e.inValidateVectors();
	}
        public  boolean isEnergyValid(PDBMolecule mol )
       {
           
          double energy=getEnergy( mol,false);
          
           if(energy<mol.numberOfResidues())
                return true;
          
           return false;
        }
        public static double getEnergy(PDBMolecule mol,boolean force)
        {
            BBEnergy EnergyObj=new BBEnergy();
            EnergyObj.setNumberOfResidues(mol);
            double potential=0;
            if(force)
            {
               
            }
             potential+=EnergyObj.calcHB(true,force,mol);
             potential+=EnergyObj.waterE(mol);
             potential+=EnergyObj.calcDeBump(force,mol);
             potential+=EnergyObj.burialE(mol);
             
            return potential;
        }
}