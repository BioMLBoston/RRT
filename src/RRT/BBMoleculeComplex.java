package RRT;

/**
 * @author Dong Luo 2012-05-24
 *
 */

import RRT.*;
import io.Parser;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.ListIterator;
import pdb.PDBMolecule;
import pdb.PDBMoleculeComplex;

public class BBMoleculeComplex extends PDBMoleculeComplex {

	private BBMoleculeComplex(ArrayList<PDBMolecule> theList) {
		super(theList);
	}
	
	public static BBMoleculeComplex newCGMoleculeComplex(PDBMoleculeComplex moleculeComplex) {
		ArrayList<PDBMolecule> moleculeList=new ArrayList<PDBMolecule>();
		ListIterator<PDBMolecule> iter = moleculeComplex.getMolecules().listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null)
				moleculeList.add(BBMolecule.newBBMolecule(pdbMolecule));
		}
		return new BBMoleculeComplex(moleculeList);
	}
	
	public static BBMoleculeComplex newCGMoleculeComplex(String filePath) {
		BBMoleculeComplex inst = newCGMoleculeComplex(Parser.parsePDBFile(filePath));
		inst.secondaryStructureFromPDBFile(filePath);
		return inst;
	}
	
	/**
	 * dihedral types corrected by PDB SS section but ignore first and last residues
	 * rigid element re-assign is automatically done
	 * TURN SS is ignored
	 * PDB SS Format reference:
	 * http://www.wwpdb.org/documentation/format23/sect5.html
	 * @param filePath
	 */
	public void secondaryStructureFromPDBFile(String filePath){
		try{
			BufferedReader reader = new BufferedReader(new FileReader(filePath));
			String pdbLine;
			BBMolecule old = null, molecule = null;
				while( (pdbLine = reader.readLine()) != null){
				if (pdbLine.length() < 4) continue;
				if (pdbLine.substring(0,6).equals("HELIX ")) {
					molecule = (BBMolecule) getMolecule(pdbLine.substring(19, 20));
					if (molecule == null) continue;
					if (molecule != old) {
						if (old != null ) 
						old = molecule;
					}
					int offset = Integer.parseInt(pdbLine.substring(21,25).trim())
							-molecule.getPDBAtomList().get(0).getResidueIndex();
					int length = Integer.parseInt(pdbLine.substring(71,76).trim());
				
				} else if (pdbLine.substring(0,6).equals("SHEET ")) {
					molecule = (BBMolecule) getMolecule(pdbLine.substring(21, 22));
					if (molecule == null) continue;
					if (molecule != old) {
						if (old != null )
						old = molecule;
					}
					int offset = Integer.parseInt(pdbLine.substring(22,26).trim())-1;
					int length = Integer.parseInt(pdbLine.substring(33,37).trim())-offset;
					offset -= molecule.getPDBAtomList().get(0).getResidueIndex()-1;
				
				}
			}
			// last molecule redefine dihedral types
			if (old != null ) 
			
			reader.close();
		}
		catch(FileNotFoundException e){
			System.out.println("Cannot find file " + filePath);
		}
		catch(IOException e){
			System.out.println("Error occured while reading from file " + filePath);
		}
	}
	
	/**
	 * pair wise energy is only calculated within each chain for now
	 * ToDo: define the way to calculate inter-chain energy
	 */
	public double getEnergy() {
		double energy=0;
		ListIterator<PDBMolecule> iter = getMolecules().listIterator();
		while(iter.hasNext())
		{
			BBMolecule molecule = (BBMolecule) iter.next();
			if(molecule != null) energy+=molecule.getEnergy(molecule, true);
		}
		return energy;
	}

	/**
	 * pair wise energy is only calculated within each chain for now
	 * ToDo: define the way to calculate inter-chain energy
	 */
	

	public double getMoleculeEnergy(String chainID) {
		BBMolecule molecule=(BBMolecule) getMolecule(chainID);
		if (molecule == null) return 0;
		return molecule.getEnergy(molecule, true);
	}


	@Override
	public String toString() {
		String info = new String();
		ListIterator<PDBMolecule> iter = getMolecules().listIterator();
		while (iter.hasNext()) info += (BBMolecule) iter.next();
		return info;
	}
}
