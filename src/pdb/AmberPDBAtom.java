package pdb;
import utilities.Definitions;
import utilities.Definitions.AmberAtomTypes;

/**
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 *
 */
public class AmberPDBAtom extends PDBAtom{
	
	////////////////////////////////////////////
	//// Constructors of AmberPDBAtom class ////
	////////////////////////////////////////////
	
	public AmberPDBAtom(PDBAtom pdbAtom, String theAmberAtomTypeStr, double theCharge){
		// int theIndex, double theX, double theY, double theZ, 
		// String theAtomTypeStr, String theAminoTypeStr, int theResidueIndex
		super(pdbAtom.getIndex(), pdbAtom.getX(), pdbAtom.getY(), 
				pdbAtom.getZ(), pdbAtom.getAtomType(), 
				pdbAtom.getAminoAcidType(), pdbAtom.getChainID(), pdbAtom.getResidueIndex());
		
		this.setAmberAtomType(Definitions.getAmberAtomType(theAmberAtomTypeStr));
		this.setCharge(theCharge);
	}
	
	////////////////////////////////////////////////
	//// Setters and getters AmberPDBAtom class ////
	////////////////////////////////////////////////
	
	/**
	 * @param amberAtomType the AMBER atom type to set
	 */
	public void setAmberAtomType(AmberAtomTypes amberAtomType) {
		this.amberAtomType = amberAtomType;
	}

	/**
	 * @return the AMBER atom type
	 */
	public AmberAtomTypes getAmberAtomType() {
		return amberAtomType;
	}
	
	/**
	 * @param charge the AMBER atom charge to set
	 */
	public void setCharge(double charge) {
		this.charge = charge;
	}
	
	/**
	 * @return the charge of AMBER atom
	 */
	public double getCharge() {
		return this.charge;
	}
	
	//////////////////////////////////////////////
	//// The operations of AmberPDBAtom class ////
	//////////////////////////////////////////////
	
	public String toString()
	{
		String pdbStr = super.toString();
		String amberPDBStr = pdbStr + " " + Definitions.getAmberAtomTypeString(this.amberAtomType);
		if(this.charge < 0)
			amberPDBStr = amberPDBStr + " " + this.charge;
		else
			amberPDBStr = amberPDBStr + "  " + this.charge;
		amberPDBStr = amberPDBStr + " " + this.getEps();
		amberPDBStr = amberPDBStr + " " + this.getRmin2();
		return amberPDBStr;
	}
	
	//////////////////////////////////////////////
	//// The attributes of AmberPDBAtom class ////
	//////////////////////////////////////////////

	/**
	 * Computes the electrostatic potential between two atoms by using Coulomb's Law,
	 * which assumes that atoms behave as point charges located at their centers. 
	 * This method is supposed to be called on a AmberPDBAtom object (say a1), with another AmberPDBAtom object (say a2) as parameter.
	 * @param a2 atom 2
	 * @return electrostaticEnergy the computed electrostatic energy betweeen this atom and the parameter atom.
	 */
	public double pairwiseElectrostaticEnergy(AmberPDBAtom a2)
	{
		// Coulomb's Law: (q_i * q_j)/(e * r_ij) 
		// where q_i and q_j are charges of atoms, 
		// e is the dielectric constant, and
		// r_ij is the distance between atoms.
		
		double qa1 = this.getCharge();
		double qa2 = a2.getCharge();
		double dc = 1; 	// dielectric constant, use vacuum constant 1 for now
						// Later, value can be experimented with, trying the values
						// from the distance b/w atoms to 80 (water's dielectric constant)
		double dist = this.distance(a2); 
		
		double electrostaticEnergy = 332 * (qa1 * qa2) / (dc * dist);
		
		return electrostaticEnergy;
	}

	/**
	 * The type of the atom
	 */
	private AmberAtomTypes amberAtomType;
	private double charge;
}
