package pdb;

import java.util.ArrayList;
import java.util.ListIterator;
import org.netlib.lapack.Dgesvd;
import pdb.AmberPDBAtom;
import pdb.PDBAtom;
import utilities.Definitions;
import utilities.Vector;

/**
 * 
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 * This class represents a single protein chain
 *
 */
public class PDBMolecule 
{
	///////////////////////////////////////////
	//// Constructors of PDBMolecule class ////
	///////////////////////////////////////////
	public PDBMolecule(String theChainID)
	{
		this.chainID = theChainID;
		this.pdbAtomList = new ArrayList<PDBAtom>(); // create an empty pdbAtomList
	}
	
	public PDBMolecule(ArrayList<PDBAtom> theAtomList, String theChainID)
	{
		this.pdbAtomList = theAtomList;
		this.chainID = theChainID;
	}
	
	//////////////////////////////////////////////////
	//// Setters and getters of PDBMolecule class ////
	//////////////////////////////////////////////////
	
	/**
	 * @param theAtomList
	 */
	public void setPDBAtomList(ArrayList<PDBAtom> theAtomList) 
	{
		this.pdbAtomList = theAtomList;
	}
	
	/**
	 * @return the pdbAtomList
	 */
	public ArrayList<PDBAtom> getPDBAtomList() {
		return this.pdbAtomList;
	}
	
	/**
	 * 
	 * @return the chain ID of the PDB molecule
	 */
	public String getChainID()
	{
		return this.chainID;
	}
	
	//////////////////////////////////////////
	//// Operations of PDBMolecule class  ////
	//////////////////////////////////////////
	
	/**
	 * @return the backbone atoms of this molecule
	 */
	public ArrayList<PDBAtom> extractBackbone()
	{
		// create an iterator on this chain
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();
		// create a list of backbone atoms
		ArrayList<PDBAtom> backbone = new ArrayList<PDBAtom>();		
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			if(a.isCA() || a.isN() || a.isC()||a.isO())
				backbone.add(a);
		}
		return backbone;
	}
	
	/**
	 * This method returns the list of C-alpha atoms in the molecule. 
	 * The list of C-alpha atoms can be used as a coarse-grained representation of the molecule.
	 * @return the C-alpha atoms of this molecule
	 */
	public ArrayList<PDBAtom> extractCAlphas()
	{
		// create an iterator on this chain
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();
		// create a list of C-alpha atoms
		ArrayList<PDBAtom> cAlphas = new ArrayList<PDBAtom>();		
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			if(a.isCA())
				cAlphas.add(a);
		}
		return cAlphas;
	}
	
	/**
	 * This method computes the centroid of a 3D point set (i.e. a molecule)
	 * @return a Vector which holds the x, y, z coordinates of the centroid
	 */
	public Vector getCentroid()
	{
		double totalX = 0;
		double totalY = 0;
		double totalZ = 0;
		Vector centroid = new Vector();
		
		ListIterator<PDBAtom> iter = this.getPDBAtomList().listIterator(0); 
		int chainSize = this.getPDBAtomList().size();	
		while(iter.hasNext() && iter.nextIndex() <= chainSize)
		{
			PDBAtom a = iter.next();
			totalX = totalX + a.getX();
			totalY = totalY + a.getY();
			totalZ = totalZ + a.getZ();
		}
		
		centroid.setX(totalX/chainSize);
		centroid.setY(totalY/chainSize);
		centroid.setZ(totalZ/chainSize);
		
		return centroid;
	}
	
	/**
	 * This method translates this molecule (i.e. chain) by calling the 
	 * translate() method on each of its atoms with the given delta values
	 * on x-coordinate, y-coordinate and z-coordinate. 
	 */
	public void translate(double deltaX, double deltaY, double deltaZ)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.translate(deltaX, deltaY, deltaZ);
		}
	}
	
	
	/**
	 * @param bondBeginAtomIndex - the index of the atom at the beginning edge of the bond
 	 * @param bondEndAtomIndex - the index of the atom at the ending edge of the bond
	 * @param angle - the angle by which atoms will be rotated around the bond
	 */
	public void rotateAroundBond(int bondBeginAtomIndex, int bondEndAtomIndex, double angle)
	{
		// Specify the bond by the indices of two atoms at the beginning and ending edges
		// This is needed to perform the rotation on atoms to be rotated.
		PDBAtom bondBegin = this.pdbAtomList.get(bondBeginAtomIndex);
		PDBAtom bondEnd = this.pdbAtomList.get(bondEndAtomIndex);
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator(bondEndAtomIndex+1);
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(bondBegin, bondEnd, angle);
		}
	}
	
	
	/**
	 * @param bondBegin - the atom at the beginning edge of the bond
	 * @param bondEnd - the atom at the ending edge of the bond
	 * @param angle -
	 */
	public  void rotateAroundBond(PDBAtom bondBegin, PDBAtom bondEnd, double angle)
	{
		// Specify the bond by two atoms at the beginning and ending edges
		// This is needed to perform the rotation on atoms to be rotated.
		int bondEndAtomIndex = this.pdbAtomList.indexOf(bondEnd);
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator(bondEndAtomIndex+1);
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(bondBegin, bondEnd, angle);
		}
	}
	
	public void rotateRigidBodyAroundArbitraryAxis(Vector begin, Vector end, double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(begin, end, angle);
		}
	}
	
	public void rotateRigidBodyAroundArbitraryAxis(double alpha, double beta, double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotate(alpha, beta, angle);
		}
	}
	
	public void rotateAroundXAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundXAxis(angle);
		}
	}
	
	public void rotateAroundYAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundYAxis(angle);
		}
	}
	
	public void rotateAroundZAxis(double angle)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom a = iter.next();
			a.rotateAroundZAxis(angle);
		}
	}
	
	public void rotateRigidBodyAroundXAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around x axis
		this.rotateAroundXAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	public void rotateRigidBodyAroundYAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around y axis
		this.rotateAroundYAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	public void rotateRigidBodyAroundZAxis(double angle)
	{
		// First translate the centroid of the molecule to the origin
		Vector centroid = this.getCentroid();
		this.translate(-1*centroid.getX(), -1*centroid.getY(), -1*centroid.getZ());
		
		// Then rotate the molecule around z axis
		this.rotateAroundZAxis(angle);
		
		// Finally translate the centroid back to original position
		this.translate(centroid.getX(), centroid.getY(), centroid.getZ());
	}
	
	/**
	 * Checks if there are any steric clashes within the given protein conformation. 
	 * @return true if clash exists between atoms of the molecule
	 * @return false if no clash exists between atoms of the molecule
	 */
	public boolean stericClashExists()
	{
		PDBAtom a1, a2;
		boolean stericClashOccurred = false;
		double energy = 0;
		
		for(int i = 0 ; i < pdbAtomList.size() ; i++)
		{	
			a1 = pdbAtomList.get(i);
			
			// check if pdbAtom at index i clashes with any other atom in the list.
			for(int j = i+4 ; j < pdbAtomList.size() ; j++)
			{
				a2 = pdbAtomList.get(j);
				
				double dist = a1.distance(a2);
				double a1Radius = a1.getRmin2();
				double a2Radius = a2.getRmin2();
				energy = energy + a1.pairwiseLJPotential(a2);
				if( (dist < a1Radius + a2Radius) && !isClashIgnorable(a1, a2) )
				{
					System.out.println("Steric clash occurred between atom " + i + " and atom " + j);
					System.out.println("Distance between these two atoms: " + dist);
					stericClashOccurred = true;
				}
			}
		}
					
		// if the control reaches here, this means that
		// no clashes detected.
		return stericClashOccurred;
	}
	
	/**
	 * This method checks if the steric clash between a pair of atoms is small enough to be ignored.
	 * @param a1 first atom of the clashing pair
	 * @param a2 second atom of the clashing pair
	 * @return true if the distance between clashing atoms is greater than 0.7 Angstroms
	 * @return false otherwise
	 */
	private static boolean isClashIgnorable(PDBAtom a1, PDBAtom a2){
		// The threshold value of 0.7 for detecting clashes is set by experimenting 
		// with the backbone of the 1COA structure of CI2, which is provided with 
		// Lydia Kavraki's Assignment-2.
		// This distance value can be improved by more trials and observations.
		double dist = a1.distance(a2);
		if(dist > 0.7)
			return true;
		return false;
	}
	
	/**
	 * This method gets the atom with the specified atom type from the specified amino acid.
	 * @param atomTypeStr the type of the desired atom, e.g "C", "CA", "N", "O", etc.
	 * @param theResidueIndex the index of the amino acid within the molecule, which the desired atom belongs to. 
	 * @return PDBAtom with the specified name and residue index
	 * @return null if the molecule does not contain such atom
	 */
	public PDBAtom getPDBAtomInMolecule(String atomTypeStr, int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getAtomType() == Definitions.getAtomType(atomTypeStr)
					&& atom.getResidueIndex() == theResidueIndex)
				return atom;
		}
		// if the control reaches here, this means that pdbAtomList
		// does not contain the given atom so return null
		return null;
	}
        
	public ArrayList<PDBAtom> getResidueInMolecule(int theResidueIndex)
	{
                ArrayList<PDBAtom> residueObj=new ArrayList<PDBAtom>();
                int startIndex=this.pdbAtomList.get(0).getResidueIndex();
                  boolean flag=false;
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getResidueIndex() == theResidueIndex )
                        { 
                                flag=true;
				residueObj.add(atom);
                               
                        }
		}
                if(flag)
                    return residueObj;
		// if the control reaches here, this means that pdbAtomList
		// does not contain the given atom so return null
		return null;
	}
	
	/**
	 * This method gets the list of atoms with the specified atom type from the specified amino acid type.
	 * @param atomTypeStr the desired type of the list of atoms, e.g "C", "CA", "N", "O", etc.
	 * @param aminoTypeStr the desired type of the amino acids that contain the atoms 
	 * @return ArrayList<PDBAtom> an array list of atoms with the specified name and residue type
	 * @return null if the molecule does not contain such atom
	 */
	public ArrayList<PDBAtom> getPDBAtomsInMolecule(String atomTypeStr, String aminoTypeStr)
	{
		ArrayList<PDBAtom> outputList = new ArrayList<PDBAtom>();
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getAtomType() == Definitions.getAtomType(atomTypeStr)
					&& atom.getAminoAcidType() == Definitions.getAminoAcidType(aminoTypeStr))
				outputList.add(atom);
		}
		// if the control reaches here, this means that pdbAtomList
		// does not contains the given atom so return null
		return outputList;
	}
	
	/**
	 * This method gets the list of atoms with the specified AMBER atom type.
	 * @param amberTypeStr the desired AMBER type of the list of atoms.
	 * @return ArrayList<AmberPDBAtom> an array list of atoms with the specified 
	 * name and residue type (null if the molecule does not contain such atom)
	 */
	public ArrayList<AmberPDBAtom> getAmberPDBAtomsInMolecule(String amberTypeStr)
	{
		ArrayList<AmberPDBAtom> outputList = new ArrayList<AmberPDBAtom>();
		
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom instanceof AmberPDBAtom)
			{
				AmberPDBAtom amberAtom = (AmberPDBAtom) atom;
				if(amberAtom.getAmberAtomType() == Definitions.getAmberAtomType(amberTypeStr)){
					outputList.add(amberAtom);
				}
			}
		}
		// if the control reaches here, this means that amberPDBAtomList
		// does not contains the given atom so return null
		return outputList;
	}
	
	/**
	 * This method adds a new pdbAtom to the pdbAtomList of this PDBMolecule.
	 * @param pdbAtom the new atom to be added
	 */
	public void addPDBAtom(PDBAtom pdbAtom)
	{
		this.pdbAtomList.add(pdbAtom);
	}
	
	/**
	 * This method gives the index of a specific PDBAtom that belongs to this PDBMolecule.
	 * @param pdbAtom
	 * @return index of the given PDBAtom; -1 if the given pdbAtom does not exist in this molecule
	 */
	public int getAtomIndex(PDBAtom pdbAtom)
	{
		for(int i = 0; i < this.pdbAtomList.size(); i++)	
		{
			PDBAtom atom = this.pdbAtomList.get(i);
			if(atom.equals(pdbAtom))
				return i;
		}
		// otherwise pdbAtom does not exist in this molecule
		// so return -1
		return -1;
	}
	
	/**
	 * This method checks if the given pair atoms are bonded or not by looking at their indices.
	 * @param a1 first PDBAtom of the given pair	
	 * @param a2 second PDBAtom of the given pair
	 * @return true if the given pair of atoms are bonded; false otherwise.
	 */
	public boolean isBonded(PDBAtom a1, PDBAtom a2)
	{
		int indexOfA1 = this.getAtomIndex(a1);
		int indexOfA2 = this.getAtomIndex(a2);
		
		if(indexOfA1 == -1 || indexOfA2 == -1)
		{
			return false;
		}
		else
		{
			if(Math.abs(indexOfA1 - indexOfA2) < 4)
			{
				return true;	
			}
			else
			{
				return false;
			}
		}
	}
	
	/**
	 * This method checks if the given pair of atoms are bonded with a peptide bond.
	 * @param a1 first PDBAtom of the given pair
	 * @param a2 first PDBAtom of the given pair
	 * @return true if the given pair of atoms share a peptide bond; false otherwise.
	 */
	public boolean isPeptideBond(PDBAtom a1, PDBAtom a2)
	{
		if(a1.getChainID().equalsIgnoreCase(this.chainID) &&
				a2.getChainID().equalsIgnoreCase(this.chainID))
		{
			if (a1.getResidueIndex() == a2.getResidueIndex()+1)
			{
				if( a1.isN() && a2.isC() )
					return true;
				return false;
			}
			else if (a2.getResidueIndex() == a1.getResidueIndex()+1)
			{
				if( a1.isC() && a2.isN() )
					return true;
				return false;
			}
			else
				return false;
		}
		else
			return false;
	}
	
	public void setTerminalTypes()
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(this.containsResidueWithSmallerIndex(atom.getResidueIndex()) &&
					this.containsResidueWithLargerIndex(atom.getResidueIndex()))
			{
				// then the chain contains residue(s) both with an index
				// larger than the atom's residue index and smaller than
				// the atom's residue index
				atom.setTerminalType(Definitions.TerminalTypes.NONTERMINAL);
			}
			else 
			{
				if(this.containsResidueWithSmallerIndex(atom.getResidueIndex()))
				{
					// then the chain contains residue(s) with an index
					// smaller than the atom's residue index but it does not
					// contain any residues with an index larger than the 
					// atom's residue index
					atom.setTerminalType(Definitions.TerminalTypes.CTERMINAL);
				}
				else
				{
					// then the chain contains residue(s) with an index
					// larger than the atom's residue index but it does not
					// contain any residues with an index smaller than the 
					// atom's residue index
					atom.setTerminalType(Definitions.TerminalTypes.NTERMINAL);
				}
			}
		}
	}
	
	private boolean containsResidueWithSmallerIndex(int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getResidueIndex() < theResidueIndex)
				return true;
		}
		return false;
	}
	
	private boolean containsResidueWithLargerIndex(int theResidueIndex)
	{
		ListIterator<PDBAtom> iter = this.pdbAtomList.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			if(atom.getResidueIndex() > theResidueIndex)
				return true;
		}
		return false;
	}

	//get number of residues
        
	public int numberOfResidues()
	{
		ArrayList<PDBAtom> iter = this.pdbAtomList;
                
                int MoleculeSize=iter.size();
                PDBAtom firstResidue=iter.get(0);
                PDBAtom LastResidue=iter.get(MoleculeSize-1);
                
                
                
		return LastResidue.getResidueIndex()-firstResidue.getResidueIndex()+1;
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
       
	/////////////////////////////////////////////
	//// The attributes of PDBMolecule class ////
	/////////////////////////////////////////////
	
	/**
	 * The list of PDB atoms in the molecule
	 */
	private ArrayList<PDBAtom> pdbAtomList; 
	
	/**
	 * The chain ID of the PDBMolecule
	 */
	private String chainID;

}
