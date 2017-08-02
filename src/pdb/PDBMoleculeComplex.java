package pdb;

import pdb.PDBMolecule;
import io.Parser;
import io.Writer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.ListIterator;

import ms.MSSurfaceEntry;

import utilities.InterfaceEnergyTuple;
import utilities.OS;
import utilities.Statistics;
import utilities.UnsupportedOSException;
import utilities.Vector;

import org.netlib.lapack.Dgesvd;
import pdb.AmberPDBAtom;
import pdb.PDBAtom;

/**
 * 
 * @author Bahar Akbal-Delibas (abakbal@cs.umb.edu)
 * This class represents a protein consisting of multiple chains.
 *
 */

public class PDBMoleculeComplex 
{
	public PDBMoleculeComplex()
	{
		// create an empty list of PDB molecules
		this.pdbMoleculeList = new ArrayList<PDBMolecule>();
	}
	
	public PDBMoleculeComplex(ArrayList<PDBMolecule> thePDBMoleculeList)
	{
		this.pdbMoleculeList = thePDBMoleculeList;
	}
	
	/**
	 * This method extracts the single chain with the given id from the molecule complex.
	 * @param theChainID the name of the desired chain, eg. "A", "B", "C", "D", etc.
	 * @return pdbMolecule - a PDBMolecule object, if a molecule exists with the given chain id; false otherwise.
	 */
	public PDBMolecule getMolecule(String theChainID)
	{
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null && pdbMolecule.getChainID().equalsIgnoreCase(theChainID))
			{
				return pdbMolecule;
			}
		}
		return null;
	}
	
	/**
	 * This method returns the molecules that belong to this molecule complex
	 * @return pdbMoleculeList - A list of PDBMolecule objects within this complex
	 */
	public ArrayList<PDBMolecule> getMolecules()
	{
		return this.pdbMoleculeList;
	}
	
	public Vector getCentroid()
	{
		double totalX = 0;
		double totalY = 0;
		double totalZ = 0;
		
		int complexSize = 0;
		
		Vector centroid = new Vector();
		
		
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule chain = iter.next();
			if(chain != null)
			{
				Vector chainCentroid = chain.getCentroid();
				int chainSize = chain.getPDBAtomList().size();
				complexSize = complexSize + chainSize;
				
				totalX = totalX + chainCentroid.getX() * chainSize;
				totalY = totalY + chainCentroid.getY() * chainSize;
				totalZ = totalZ + chainCentroid.getZ() * chainSize;
			}
		}
		
		centroid.setX(totalX/complexSize);
		centroid.setY(totalY/complexSize);
		centroid.setZ(totalZ/complexSize);
		
		return centroid;
	}
	
	public void addPDBAtom(PDBAtom pdbAtom)
	{
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null 
					&& pdbMolecule.getChainID().equalsIgnoreCase(pdbAtom.getChainID()))
			{
				pdbMolecule.addPDBAtom(pdbAtom);
				return;
			}
		}
		// If the control comes here, this means that 
		// pdbAtom is not added to any molecules yet.
		// In other words, this PDB Molecule Complex
		// does not have a molecule with the chain ID 
		// of the given pdbAtom.
		// Therefore, first create a new molecule,
		// then add pdbAtom to this newly created molecule.
		PDBMolecule pdbMolecule = new PDBMolecule(pdbAtom.getChainID());
		pdbMolecule.addPDBAtom(pdbAtom);
		this.pdbMoleculeList.add(pdbMolecule);
	}
	
	/**
	 * This method searches this molecule complex for the PDBAtom objects with the specified atom type and residue index. 
	 * @param atomTypeStr the name of the desired kind of atom, eg. "C", "CA", "N", etc.
	 * @param theResidueIndex the index of the amino acid that the desired atom belongs to
	 * @return pdbAtom - a PDBAtom object if the atom being searched for is found; null otherwise
	 */
	public PDBAtom getPDBAtomInMoleculeComplex(String atomTypeStr, int theResidueIndex)
	{
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null )
			{
				PDBAtom pdbAtom = pdbMolecule.getPDBAtomInMolecule(atomTypeStr, theResidueIndex);
				if(pdbAtom != null)
				{
					return pdbAtom;
				}
			}
		}
		return null;
	}
	
	/**
	 * This method searches this molecule complex for the list of PDBAtom objects with the specified atom type and the amino acid type. 
	 * @param atomTypeStr the name of the desired atom kind, eg. "C", "CA", "N", etc.
	 * @param aminoTypeStr the type of the amino acid that the desired atoms belong to
	 * @return outputList - a list of PDBAtom objects; null if no such atoms exist in the molecule complex.
	 */
	public ArrayList<PDBAtom> getPDBAtomsInMoleculeComplex(String atomTypeStr, String aminoTypeStr)
	{
		ArrayList<PDBAtom> outputList = new ArrayList<PDBAtom>();
		
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null )
			{
				ArrayList<PDBAtom> pdbAtoms = pdbMolecule.getPDBAtomsInMolecule(atomTypeStr, aminoTypeStr);
				ListIterator<PDBAtom> iter2 = pdbAtoms.listIterator();
				while(iter2.hasNext())
				{
					PDBAtom atom = iter2.next();
					outputList.add(atom);
				}
			}
		}
		return outputList;
	}
	/**
	 * This method computes and returns a list of PDBAtom objects 
	 * that are on the interface of the given pair of molecules 
	 * belonging to this molecule complex.
	 * @param chain1List first chain (or combined set of chains) of the given pair of molecules
	 * @param chain2List second chain (or combined set of chains) of the given pair of molecules
	 * @return interfaceAtoms a list of PDBAtom objects that are on the interface of these two chains.
	 */
	public ArrayList<PDBAtom> computeInterfaceAtoms(ArrayList<String> chain1List, ArrayList<String> chain2List)
	{
		this.interfaceAtoms = new ArrayList<PDBAtom>();
		
		ArrayList<PDBAtom> chain1ListAtoms = getPDBAtomsInMoleculeList(chain1List);
		ListIterator<PDBAtom> chain1Iterator = chain1ListAtoms.listIterator();
		
		while(chain1Iterator.hasNext())
		{
			PDBAtom chain1Atom = chain1Iterator.next();
			
			ArrayList<PDBAtom> chain2ListAtoms = getPDBAtomsInMoleculeList(chain2List);
			ListIterator<PDBAtom> chain2Iterator = chain2ListAtoms.listIterator();
			
			while(chain2Iterator.hasNext())
			{				
				PDBAtom chain2Atom = chain2Iterator.next();
				double dist = chain1Atom.distance(chain2Atom);
				if( dist < 6)
				{										
					// then add chain1Atom and chain2Atom
					// to the interface atom list if they are not
					// already in the interface atom list
					if(interfaceAtoms.contains(chain1Atom) == false)
					{
						interfaceAtoms.add(chain1Atom);
					}
					if(interfaceAtoms.contains(chain2Atom) == false)
					{
						interfaceAtoms.add(chain2Atom);
					}
				}
			}
		}
		return interfaceAtoms;
	}
	
	/**
	 * This method computes and returns the list of PDBAtom objects 
	 * contained within a given list of chains
	 * @param chainList a list of chain IDs
	 * @return atomsInMolecules a list of PDBAtom objects that are included within chainList.
	 */
	private ArrayList<PDBAtom> getPDBAtomsInMoleculeList(ArrayList<String> chainList)
	{
		ArrayList<PDBAtom> atomsInMolecules = new ArrayList<PDBAtom>();
		ListIterator<String> chainListIterator = chainList.listIterator();
		while(chainListIterator.hasNext())
		{
			String chain = chainListIterator.next();
			PDBMolecule chainMolecule = this.getMolecule(chain);
			if(chainMolecule != null)
			{
				ListIterator<PDBAtom> chainIterator = 
					chainMolecule.getPDBAtomList().listIterator();
				while(chainIterator.hasNext())
				{
					PDBAtom chainAtom = chainIterator.next();
					atomsInMolecules.add(chainAtom);
				}
			}
		}
		
		return atomsInMolecules;
	}
	
	/**
	 * This method returns the list of PDBAtom objects, which 
	 * correspond to the latest computed list of interface atoms 
	 * of an arbitrary pair of chains in this complex 
	 * @return interfaceAtoms - the latest computed atoms on the interface of arbitrarily picked two chains.
	 */
	public ArrayList<PDBAtom> getInterfaceAtoms()
	{
		return this.interfaceAtoms;
	}
	
	public static boolean isIncludedInChainList(ArrayList<String> chainList, String theChain)
	{
		ListIterator<String> chainListIterator = chainList.listIterator();
		while(chainListIterator.hasNext())
		{
			String chain = chainListIterator.next();
			if(chain.equalsIgnoreCase(theChain))
				return true;
		}
		return false;
	}
	
	public InterfaceEnergyTuple getInterfaceEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths, String msroll_path)
			throws UnsupportedOSException
	{
		PDBAtom a, a1, a2;
		double energy = 0;
		double vdWTerm = 0;
		double electroStaticTerm = 0;
		double conservationTerm = 0;
		double distanceRestraintTerm = 0;
		
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = new ArrayList<PDBAtom>();
		ArrayList<PDBAtom> interfaceAtoms_chain2 = new ArrayList<PDBAtom>();
		
		for(int i = 0 ; i < interfaceAtoms.size() ; i++)
		{
			a = interfaceAtoms.get(i);
			if(isIncludedInChainList(chainList1, a.getChainID()))
			{
				interfaceAtoms_chain1.add(a);
			}
			else if(isIncludedInChainList(chainList2, a.getChainID()))
			{
				interfaceAtoms_chain2.add(a);
			}
		}
		
		// Create two temporary PDB files for chain1 interface atoms
		// and chain2 interface atoms.
		// These temporary PDB files are needed to calculate
		// the molecular surfaces of chain1 interface atoms
		// and chain2 interface atoms.
		////Writer.createPDBFile(interfaceAtoms_chain1, "pdb_files/temp_interface_chain1.pdb");
		////Writer.createPDBFile(interfaceAtoms_chain2, "pdb_files/temp_interface_chain2.pdb");
		Writer.createPDBFile(interfaceAtoms_chain1, "temp_interface_chain1.pdb");
		Writer.createPDBFile(interfaceAtoms_chain2, "temp_interface_chain2.pdb");
		
		// Using the external MS software, we can calculate
		// the molecular surfaces of chain1 and chain2 
		// interface atoms.
		try
	    {
			String /*msroll_path,*/ command, input_pdb_file, chain1SurfaceFile, chain2SurfaceFile;
			/*if(OS.isMac())
			{
				msroll_path = "msp_exe/mac/msroll";
			}
			else if(OS.isWindows())
			{
				msroll_path = "msp_exe/win/msroll.exe";
			}
			else if(OS.isSolaris())
			{
				msroll_path = "msp_exe/sun/msroll.exe";
			}
			else if(OS.isIrix())
			{
				msroll_path = "msp_exe/sgi/msroll.exe";
			}
			else
			{
				throw new UnsupportedOSException(OS.getOS());
			}
			*/
			
			// first execute MS command to find out
			// the molecular surfaces of chain1 interface atoms
			////input_pdb_file = "pdb_files/temp_interface_chain1.pdb";
			////chain1SurfaceFile = "msp_files/temp_interface_chain1_surface.sur";
			input_pdb_file = "temp_interface_chain1.pdb";
			chain1SurfaceFile = "temp_interface_chain1_surface.sur";
			command = msroll_path + " -m " + input_pdb_file + " -j " + chain1SurfaceFile;
			
			Process process = Runtime.getRuntime().exec(command);
			process.waitFor();
			
			// then execute MS command to find out
			// the molecular surfaces of chain2 interface atoms
			////input_pdb_file = "pdb_files/temp_interface_chain2.pdb";
			////chain2SurfaceFile = "msp_files/temp_interface_chain2_surface.sur";
			input_pdb_file = "temp_interface_chain2.pdb";
			chain2SurfaceFile = "temp_interface_chain2_surface.sur";
			command = msroll_path + " -m " + input_pdb_file + " -j " + chain2SurfaceFile;
			
			process = Runtime.getRuntime().exec(command);
			process.waitFor();
			
			// Now we need to parse MS surface files for chain1 and chain2
			// interface atoms.
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain1SurfaceMap = 
				Parser.parseMSSurfaceFile(chain1SurfaceFile);
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain2SurfaceMap = 
				Parser.parseMSSurfaceFile(chain2SurfaceFile);
			
			if(interfaceAtoms.isEmpty() == false)
			{
				// Parse ET rank files both for chain1 and chain2
				// this is needed for calculating conservation energy
				HashMap<Integer, Double> residueETRankMap_chain1 = 
					Parser.parseNewETRankFile(chain1ETRankFilePaths);
				HashMap<Integer, Double> residueETRankMap_chain2 = 
					Parser.parseNewETRankFile(chain2ETRankFilePaths);
				
				for(int i = 0 ; i < interfaceAtoms_chain1.size() ; i++)
				{	
					a1 = interfaceAtoms_chain1.get(i);
					for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
					{
						a2 = interfaceAtoms_chain2.get(j);
						if(!(a1.getChainID().equalsIgnoreCase(a2.getChainID())))
						{
							// increase vdWTerm
							vdWTerm = vdWTerm + a1.pairwiseLJPotential(a2);
							
							// increase electroStaticTerm
							if(a1 instanceof AmberPDBAtom && a2 instanceof AmberPDBAtom)
							{
								AmberPDBAtom amberA1 = (AmberPDBAtom) a1;
								AmberPDBAtom amberA2 = (AmberPDBAtom) a2;
								electroStaticTerm = 
									electroStaticTerm + amberA1.pairwiseElectrostaticEnergy(amberA2);
							}
							
							// increase conservationTerm
							conservationTerm = conservationTerm + conservationForAtomPair(a1, a2, 
																			residueETRankMap_chain1, 
																			residueETRankMap_chain2,
																			i+1, j+1,
																			chain1SurfaceMap,
																			chain2SurfaceMap);
							
						}
					}
				}
				//System.out.println("vdWTerm = " + vdWTerm);
				//System.out.println("electroStaticTerm = " + electroStaticTerm);
				//System.out.println("conservationTerm = " + conservationTerm);
				// energy = vdWTerm + electroStaticTerm + conservationTerm;
				
				
				//////////////////////////////////////////////////////
				//////////////////////////////////////////////////////
				//													//
				// 	REMOVE THIS FROM HERE TO ITS OWN METHOD !!!		//
				// 	IT IS CALLED HERE FOR TESTING PURPOSES FOR NOW	//
				//													//	
				//////////////////////////////////////////////////////
				//////////////////////////////////////////////////////
				distanceRestraintTerm = this.interfaceDistanceRestraintEnergy(chainList1, chainList2, residueETRankMap_chain1, residueETRankMap_chain2);
				
				energy = vdWTerm + electroStaticTerm - 5 * (conservationTerm / distanceRestraintTerm);
			}
			
	    }catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e)
	    {
	    	e.printStackTrace();
	    }
		
		return new InterfaceEnergyTuple(energy, vdWTerm, electroStaticTerm, conservationTerm, distanceRestraintTerm);
	}
	
	/**
	 * This method returns the interface energy that is computed using our own conservation term calculation;
	 * It does not use surface compatibility information like Kanamori's conservation term
	 * @return 
	 */
	public InterfaceEnergyTuple getInterfaceEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	{
		PDBAtom a, a1, a2;
		double energy = 0;
		double vdWTerm = 0;
		double electroStaticTerm = 0;
		double conservationTerm = 0;
		double distanceRestraintTerm = 0;
		
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = new ArrayList<PDBAtom>();
		ArrayList<PDBAtom> interfaceAtoms_chain2 = new ArrayList<PDBAtom>();
		
		for(int i = 0 ; i < interfaceAtoms.size() ; i++)
		{
			a = interfaceAtoms.get(i);
			if(isIncludedInChainList(chainList1, a.getChainID()))
			{
				interfaceAtoms_chain1.add(a);
			}
			else if(isIncludedInChainList(chainList2, a.getChainID()))
			{
				interfaceAtoms_chain2.add(a);
			}
		}
		
		if(interfaceAtoms.isEmpty() == false)
		{
			// Parse ET rank files both for chain1 and chain2
			// this is needed for calculating conservation energy
			HashMap<Integer, Double> residueETRankMap_chain1 = 
				Parser.parseNewETRankFile(chain1ETRankFilePaths);
			HashMap<Integer, Double> residueETRankMap_chain2 = 
				Parser.parseNewETRankFile(chain2ETRankFilePaths);
			
			for(int i = 0 ; i < interfaceAtoms_chain1.size() ; i++)
			{	
				a1 = interfaceAtoms_chain1.get(i);
				for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
				{
					a2 = interfaceAtoms_chain2.get(j);
					if(!(a1.getChainID().equalsIgnoreCase(a2.getChainID())))
					{
						// increase vdWTerm
						vdWTerm = vdWTerm + a1.pairwiseLJPotential(a2);
						
						// increase electroStaticTerm
						if(a1 instanceof AmberPDBAtom && a2 instanceof AmberPDBAtom)
						{
							AmberPDBAtom amberA1 = (AmberPDBAtom) a1;
							AmberPDBAtom amberA2 = (AmberPDBAtom) a2;
							electroStaticTerm = 
								electroStaticTerm + amberA1.pairwiseElectrostaticEnergy(amberA2);
						}
						
						// increase conservationTerm
						conservationTerm = conservationTerm + conservationForAtomPair(a1, a2, 
																		residueETRankMap_chain1, 
																		residueETRankMap_chain2);	
					}
				}
			}
			
			distanceRestraintTerm = this.interfaceDistanceRestraintEnergy(chainList1, chainList2, residueETRankMap_chain1, residueETRankMap_chain2);
			
			// energy = vdWTerm + electroStaticTerm - 5 * (conservationTerm / distanceRestraintTerm);
			energy = vdWTerm - conservationTerm;
		}
		
		return new InterfaceEnergyTuple(energy, vdWTerm, electroStaticTerm, conservationTerm, distanceRestraintTerm);
	}
	
	
	public double interfaceVDWEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2)
	{
		PDBAtom a, a1, a2;
		double vdWEnergy = 0;
		
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = new ArrayList<PDBAtom>();
		ArrayList<PDBAtom> interfaceAtoms_chain2 = new ArrayList<PDBAtom>();
		
		for(int i = 0 ; i < interfaceAtoms.size() ; i++)
		{
			a = interfaceAtoms.get(i);
			if(isIncludedInChainList(chainList1, a.getChainID()))
			{
				interfaceAtoms_chain1.add(a);
			}
			else if(isIncludedInChainList(chainList2, a.getChainID()))
			{
				interfaceAtoms_chain2.add(a);
			}
		}
		
		for(int i = 0 ; i < interfaceAtoms_chain1.size() ; i++)
		{	
			a1 = interfaceAtoms_chain1.get(i);
			for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
			{
				a2 = interfaceAtoms_chain2.get(j);
				if(!(a1.getChainID().equalsIgnoreCase(a2.getChainID())))
				{
					// increase vdWEnergy
					vdWEnergy = vdWEnergy + a1.pairwiseLJPotential(a2);
				}
			}
		}
		return vdWEnergy;
	}
	
	public double interfaceElectroStaticEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2)
	{
		PDBAtom a, a1, a2;
		double electroStaticEnergy = 0;

		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = new ArrayList<PDBAtom>();
		ArrayList<PDBAtom> interfaceAtoms_chain2 = new ArrayList<PDBAtom>();
		
		for(int i = 0 ; i < interfaceAtoms.size() ; i++)
		{
			a = interfaceAtoms.get(i);
			if(isIncludedInChainList(chainList1, a.getChainID()))
			{
				interfaceAtoms_chain1.add(a);
			}
			else if(isIncludedInChainList(chainList2, a.getChainID()))
			{
				interfaceAtoms_chain2.add(a);
			}
		}
		
		for(int i = 0 ; i < interfaceAtoms_chain1.size() ; i++)
		{	
			a1 = interfaceAtoms_chain1.get(i);
			for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
			{
				a2 = interfaceAtoms_chain2.get(j);
				if(!(a1.getChainID().equalsIgnoreCase(a2.getChainID())))
				{
					// increase electroStaticEnergy
					if(a1 instanceof AmberPDBAtom && a2 instanceof AmberPDBAtom)
					{
						AmberPDBAtom amberA1 = (AmberPDBAtom) a1;
						AmberPDBAtom amberA2 = (AmberPDBAtom) a2;
						electroStaticEnergy = 
							electroStaticEnergy + amberA1.pairwiseElectrostaticEnergy(amberA2);
					}
				}
			}
		}
		return electroStaticEnergy;
	}
	
	/**
	 * This method computes the degree of shape complementarity of two molecular surfaces 
	 * and the total conservation of neighboring vertices on each surface. The equation 
	 * is taken from Kanamori et al paper (2007). Evolutionary Trace files are needed to 
	 * parse the Evolutionary Trace scores of each amino acid when computing the conservation scores.
	 * ET rank files are retrieved from the Evolutionary Trace Server for each protein from the following link: 
	 * http://mammoth.bcm.tmc.edu/ETserver.html 
	 * 
	 * @param chainList1  first chain (or combined set of chains) of the given pair of molecules
	 * @param chainList2  second chain (or combined set of chains) of the given pair of molecules
	 * @param chain1ETRankFilePaths file path of first chain's ET rank file 
	 * @param chain2ETRankFilePaths file path of first chain's ET rank file
	 * @return The total surface complementarity value, which also embodies the interface conservation information.
	 * @throws UnsupportedOSException
	 */
	private double interfaceConservationEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	throws UnsupportedOSException
	{
		PDBAtom a1, a2;
		double conservationEnergy = 0;
		
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = getChainAtomsInInterface(interfaceAtoms, chainList1);
		ArrayList<PDBAtom> interfaceAtoms_chain2 = getChainAtomsInInterface(interfaceAtoms, chainList2);
		
		// Create two temporary PDB files for chain1 interface atoms
		// and chain2 interface atoms.
		// These temporary PDB files are needed to calculate
		// the molecular surfaces of chain1 interface atoms
		// and chain2 interface atoms.
		Writer.createPDBFile(interfaceAtoms_chain1, "pdb_files/temp_interface_chain1.pdb");
		Writer.createPDBFile(interfaceAtoms_chain2, "pdb_files/temp_interface_chain2.pdb");
		
		// Using the external MS software, we can calculate
		// the molecular surfaces of chain1 and chain2 
		// interface atoms.
		try
	    {
			String msroll_path, command, input_pdb_file, chain1SurfaceFile, chain2SurfaceFile;
			if(OS.isMac())
			{
				msroll_path = "msp_exe/mac/msroll.exe";
			}
			else if(OS.isWindows())
			{
				msroll_path = "msp_exe/win/msroll.exe";
			}
			else if(OS.isSolaris())
			{
				msroll_path = "msp_exe/sun/msroll.exe";
			}
			else if(OS.isIrix())
			{
				msroll_path = "msp_exe/sgi/msroll.exe";
			}
			else
			{
				throw new UnsupportedOSException(OS.getOS());
			}
			
			// first execute MS command to find out
			// the molecular surfaces of chain1 interface atoms
			input_pdb_file = "pdb_files/temp_interface_chain1.pdb";
			chain1SurfaceFile = "msp_files/temp_interface_chain1_surface.sur";
			command = msroll_path + " -m " + input_pdb_file + " -j " + chain1SurfaceFile;
			
			Process process = Runtime.getRuntime().exec(command);
			process.waitFor();
			
			// then execute MS command to find out
			// the molecular surfaces of chain2 interface atoms
			input_pdb_file = "pdb_files/temp_interface_chain2.pdb";
			chain2SurfaceFile = "msp_files/temp_interface_chain2_surface.sur";
			command = msroll_path + " -m " + input_pdb_file + " -j " + chain2SurfaceFile;
			
			process = Runtime.getRuntime().exec(command);
			process.waitFor();
			
			// Now we need to parse MS surface files for chain1 and chain2
			// interface atoms.
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain1SurfaceMap = 
				Parser.parseMSSurfaceFile(chain1SurfaceFile);
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain2SurfaceMap = 
				Parser.parseMSSurfaceFile(chain2SurfaceFile);
			
			if(interfaceAtoms.isEmpty() == false)
			{
				// Parse ET rank files both for chain1 and chain2
				// this is needed for calculating conservation energy
				HashMap<Integer, Double> residueETRankMap_chain1 = 
					Parser.parseNewETRankFile(chain1ETRankFilePaths);
				HashMap<Integer, Double> residueETRankMap_chain2 = 
					Parser.parseNewETRankFile(chain2ETRankFilePaths);
				
				for(int i = 0 ; i < interfaceAtoms_chain1.size() ; i++)
				{	
					a1 = interfaceAtoms_chain1.get(i);
					for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
					{
						a2 = interfaceAtoms_chain2.get(j);
						if(!(a1.getChainID().equalsIgnoreCase(a2.getChainID())))
						{
							// increase conservationTerm
							conservationEnergy = conservationEnergy + conservationForAtomPair(a1, a2, 
																			residueETRankMap_chain1, 
																			residueETRankMap_chain2,
																			i+1, j+1,
																			chain1SurfaceMap,
																			chain2SurfaceMap);
						}
					}
				}
			}
			
	    }catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e)
	    {
	    	e.printStackTrace();
	    }
		
		return conservationEnergy;
	}
	
	/**
	 * This method computes the conservation values of all atoms in a given chain 
	 * based on the equation in Kanamori et al paper (2007). 
	 * Evolutionary Trace files are needed to parse the Evolutionary Trace scores 
	 * of each amino acid when computing the conservation scores.
	 * ET rank files are retrieved from the Evolutionary Trace Server for each 
	 * protein from the following link: 
	 * 		http://mammoth.bcm.tmc.edu/ETserver.html 
	 * The computed conservation values are written to a file.
	 * @param chainList The protein chain that conservation values are computed for.
	 * @param ETRankFilePath File path of the chain's ET rank file
	 */
	public void chainConservation(ArrayList<String> chainList, ArrayList<String> ETRankFilePath)
	{
		double conservation = 0;
		ArrayList<PDBAtom> chainAtoms = getPDBAtomsInMoleculeList(chainList);
		ListIterator<PDBAtom> chainIterator = chainAtoms.listIterator();
		ArrayList<Integer> activeResidueIndices = new ArrayList<Integer>();
		System.out.println("Computing residue conservation values:");
		//try{
				// Parse ET rank file for the chain,
				// this is needed for calculating conservation energy
				HashMap<Integer, Double> residueETRankMap = 
						Parser.parseNewETRankFile(ETRankFilePath);
				
				while(chainIterator.hasNext())
				{
					PDBAtom atom = chainIterator.next();
					conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
					if ( !activeResidueIndices.contains(atom.getResidueIndex()) )
					{
					//if(isActiveResidue(residueETRankMap, atom.getResidueIndex()) && !activeResidueIndices.contains(atom.getResidueIndex()))
					//if(conservation >= 0.5 && !activeResidueIndices.contains(atom.getResidueIndex()))
						if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
						{	
							activeResidueIndices.add(atom.getResidueIndex());
							//System.out.println(atom.getChainID() + "\t" + atom.getResidueIndex() + "\t" + atom.getAminoAcidType() + "\t" + conservation);
							System.out.println(atom.getResidueIndex() + "\t" + 1);
						}
						else 
						{
							activeResidueIndices.add(atom.getResidueIndex());
							System.out.println(atom.getResidueIndex() + "\t" + 0);
						}
					}
				}
			
		//}catch (InterruptedException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		//}
		//catch (IOException e)
	    //{
	    	//e.printStackTrace();
	   // }
	}
	
	public ArrayList<PDBAtom> getActiveAtoms(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	{
		ArrayList<PDBAtom> activeAtoms = new ArrayList<PDBAtom>();
		
		// first calculate the number of active atoms in chainList1
		for(int i = 0; i < chainList1.size(); i++)
		{
			String chainID = chainList1.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain1ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				//double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
				{	
					activeAtoms.add(atom);
				}
			}
		}
		
		// then calculate the number of active atoms in chainList2
		for(int i = 0; i < chainList2.size(); i++)
		{
			String chainID = chainList2.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain2ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(conservation <= 0.5)
				{	
					activeAtoms.add(atom);
				}
			}
		}
		
		return activeAtoms;
	}
	
	public ArrayList<PDBAtom> getAtomsWithETScore(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths, 
			int stdevFromMean)
	{
		ArrayList<PDBAtom> activeAtoms = new ArrayList<PDBAtom>();
		
		// first calculate the number of active atoms in chainList1
		for(int i = 0; i < chainList1.size(); i++)
		{
			String chainID = chainList1.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain1ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<Double> rankList = new ArrayList<Double>(residueETRankMap.values());
			double mean = Statistics.mean(rankList);
			double stdev = Statistics.standardDeviation(rankList);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				Double residueRankDouble = residueETRankMap.get(new Integer(atom.getResidueIndex()));
				if(residueRankDouble != null)
				{
					double residueRank = residueRankDouble.doubleValue();
					
					if(residueRank <= (mean + stdev * stdevFromMean))
					{	
						activeAtoms.add(atom);
					}
				}
			}
		}
		
		// then calculate the number of active atoms in chainList2
		for(int i = 0; i < chainList2.size(); i++)
		{
			String chainID = chainList2.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain2ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<Double> rankList = new ArrayList<Double>(residueETRankMap.values());
			double mean = Statistics.mean(rankList);
			double stdev = Statistics.standardDeviation(rankList);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				Double residueRankDouble = residueETRankMap.get(new Integer(atom.getResidueIndex()));
				if(residueRankDouble != null)
				{
					double residueRank = residueRankDouble.doubleValue();
					
					if(residueRank <= (mean + stdev * stdevFromMean))
					{	
						activeAtoms.add(atom);
					}
				}
			}
		}
		
		return activeAtoms;
	}
	
	//public int getNumberOfActiveAtoms(ArrayList<String> chainList, ArrayList<String> allETRankFilePaths)
	public int getNumberOfActiveAtoms(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	{
		int numberOfActiveAtoms = 0;
		
		// first calculate the number of active atoms in chainList1
		for(int i = 0; i < chainList1.size(); i++)
		{
			String chainID = chainList1.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain1ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				//double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
				{	
					numberOfActiveAtoms = numberOfActiveAtoms + 1;
				}
			}
		}
		
		// then calculate the number of active atoms in chainList2
		for(int i = 0; i < chainList2.size(); i++)
		{
			String chainID = chainList2.get(i);
			PDBMolecule chain = this.getMolecule(chainID);
			
			ArrayList<String> ETRankFilePath = new ArrayList<String>();
			ETRankFilePath.add(chain2ETRankFilePaths.get(i));
			HashMap<Integer, Double> residueETRankMap = 
					Parser.parseNewETRankFile(ETRankFilePath);
			
			ArrayList<PDBAtom> chainAtoms = chain.getPDBAtomList();
			ListIterator<PDBAtom> iter = chainAtoms.listIterator();
			while(iter.hasNext())
			{
				PDBAtom atom = iter.next();
				//double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
				{	
					numberOfActiveAtoms = numberOfActiveAtoms + 1;
				}
			}
		}
		
		return numberOfActiveAtoms;
	}
	
	public int getNumberOfActiveAtomsInInterface(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	{
		int numberOfActiveAtoms = 0;
		
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		ListIterator<PDBAtom> iter = interfaceAtoms.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			
			HashMap<Integer, Double> residueETRankMap = null;		
			if(isIncludedInChainList(chainList1, atom.getChainID()))
			{
				residueETRankMap = Parser.parseNewETRankFile(chain1ETRankFilePaths);
			}
			else if(isIncludedInChainList(chainList2, atom.getChainID()))
			{
				residueETRankMap = Parser.parseNewETRankFile(chain2ETRankFilePaths);
			}
			
			if(residueETRankMap != null)
			{
				//double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
				{	
					numberOfActiveAtoms = numberOfActiveAtoms + 1;
				}
			}
		}
		return numberOfActiveAtoms;
	}
	
	public static int getNumberOfActiveAtomsInInterface(ArrayList<PDBAtom> interfaceAtoms, ArrayList<String> chainList1, 
			ArrayList<String> chainList2, ArrayList<String> chain1ETRankFilePaths, ArrayList<String> chain2ETRankFilePaths)
	{
		int numberOfActiveAtoms = 0;
		
		ListIterator<PDBAtom> iter = interfaceAtoms.listIterator();
		while(iter.hasNext())
		{
			PDBAtom atom = iter.next();
			
			HashMap<Integer, Double> residueETRankMap = null;		
			if(isIncludedInChainList(chainList1, atom.getChainID()))
			{
				residueETRankMap = Parser.parseNewETRankFile(chain1ETRankFilePaths);
			}
			else if(isIncludedInChainList(chainList2, atom.getChainID()))
			{
				residueETRankMap = Parser.parseNewETRankFile(chain2ETRankFilePaths);
			}
			
			if(residueETRankMap != null)
			{
				//double conservation = PDBMoleculeComplex.residueConservationValue(residueETRankMap, atom.getResidueIndex()); 
				if(isActiveResidue(residueETRankMap, atom.getResidueIndex()))
				{	
					numberOfActiveAtoms = numberOfActiveAtoms + 1;
				}
			}
		}
		return numberOfActiveAtoms;
	}
	
	public int getNumberOfInterfaceAtoms(ArrayList<String> chainList1, ArrayList<String> chainList2)
	{
		ArrayList<PDBAtom> interfaceAtoms = this.computeInterfaceAtoms(chainList1, chainList2);
		return interfaceAtoms.size();
	}
	
	
	private static double conservationForAtomPair(PDBAtom a1, PDBAtom a2,
			HashMap<Integer, Double> residueETRankMap_a1, 
			HashMap<Integer, Double> residueETRankMap_a2,
			int a1_index, int a2_index,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain1SurfaceMap,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain2SurfaceMap)
	{
		// n1 : normal vector on surface 1 for a1
		Vector n1 = getA1NormalVector(a1_index, a2_index, chain1SurfaceMap, chain2SurfaceMap);
		
		// n2 : normal vector on surface 2 for a2
		Vector n2 = getA2NormalVector(a1_index, a2_index, chain1SurfaceMap, chain2SurfaceMap);
	   
		double I;
		if(n1 != null || n2 != null)
		{
			I = (-1 * n1.dotProduct(n2) + 1) / 2;
		}
		else
		{
			I = 1.0;
		}
		
		double c1 = 
			PDBMoleculeComplex.residueConservationValue(residueETRankMap_a1, a1.getResidueIndex());
		double c2 = 
			PDBMoleculeComplex.residueConservationValue(residueETRankMap_a2, a2.getResidueIndex());
		double D;
		if(a1.distance(a2) <= 0.5)
			D = 1.0;
		else
			D = 1 / (4 * a1.distance(a2) * a1.distance(a2));
		
		/*
		System.out.println("I = " + I + " for a1_index = " + a1_index 
				+ " a2_index = " + a2_index);		
		*/
		return I * I * c1 * c2 * c1 * c2 * D;
		// return c1 * c2 * c1 * c2 * D;
	}
	
	private static double conservationForAtomPair(PDBAtom a1, PDBAtom a2,
			HashMap<Integer, Double> residueETRankMap_a1, 
			HashMap<Integer, Double> residueETRankMap_a2)
	{
		
		double c1 = 
			PDBMoleculeComplex.residueConservationValue(residueETRankMap_a1, a1.getResidueIndex());
		double c2 = 
			PDBMoleculeComplex.residueConservationValue(residueETRankMap_a2, a2.getResidueIndex());

		if(c1 < 0 && c2 < 0)
		{
			return -1 * c1 * c2;
		}
		else
		{
			return c1 * c2;
		}
		
		//return c1*c2;
	}
	
	public static boolean isActiveResidue(HashMap<Integer, Double> residueETRankMap, int residueNumber)
	{
		ArrayList<Double> rankList = new ArrayList<Double>(residueETRankMap.values());
		
		double mean = Statistics.mean(rankList);
		Double residueRankDouble = residueETRankMap.get(new Integer(residueNumber));
		if(residueRankDouble != null)
		{
			double residueRank = residueRankDouble.doubleValue();
			if(residueRank < mean)
				return true;
			else
				return false;
		}
		return false;
	}
	
	
	/**
	 * Effective distance restraints must be computed for each "active" amino acid of each chain.
	 * Must find out what to do with those values; sum them, multiply them, etc?
	 * Right now, the code simply computes the sum of effective distance restraints of a chain.
	 * Also, the threshold value for an atom to be considered evolutionary conserved is not known,
	 * must find this out and use here.
	 * And how about replacing the for-loops with iterators?
	 * This method needs to be modified!
	 * */
	private double interfaceDistanceRestraintEnergy(ArrayList<String> chainList1, ArrayList<String> chainList2, 
			HashMap<Integer, Double> residueETRankMap_chain1, HashMap<Integer, Double> residueETRankMap_chain2)
	{		
		PDBAtom a1, a2;
		double effectiveDistance1 = 0;
		double effectiveDistance2 = 0;
		double distance = 0;
		//ArrayList<PDBAtom> interfaceAtoms = this.getInterfaceAtoms(chain1, chain2);
		ArrayList<PDBAtom> interfaceAtoms = this.getInterfaceAtoms();
		
		// Find chain1 and chain2 atoms in the interface
		ArrayList<PDBAtom> interfaceAtoms_chain1 = getChainAtomsInInterface(interfaceAtoms, chainList1);
		ArrayList<PDBAtom> interfaceAtoms_chain2 = getChainAtomsInInterface(interfaceAtoms, chainList2);
		
		// Parse ET rank files both for chain1 and chain2
		// this is needed for calculating conservation energy
		/*HashMap<Integer, Double> residueETRankMap_chain1 = 
			Parser.parseETRankFile(chain1ETRankFilePaths);
		HashMap<Integer, Double> residueETRankMap_chain2 = 
			Parser.parseETRankFile(chain2ETRankFilePaths);*/
		
		// Compute effective distance for active atoms on one chain
		for(int i = 0 ; i <  interfaceAtoms_chain1.size() ; i++)
		{	
			a1 = interfaceAtoms_chain1.get(i);
			//double c1 = 
				//PDBMoleculeComplex.residueConservationValue(residueETRankMap_chain1, a1.getResidueIndex());
			
			if(isActiveResidue(residueETRankMap_chain1, a1.getResidueIndex()))
			{	
				for(int j = 0 ; j < interfaceAtoms_chain2.size() ; j++)
				{	
					a2 = interfaceAtoms_chain2.get(j);
					distance = a1.distance(a2);
					
					double t = 1 / (Math.pow(distance, 6));
					effectiveDistance1 = effectiveDistance1 + t ;
				}
			}
			// Computing effective distance for the current active atom 
			// against all atoms on the opposing chain ends here
		}
		effectiveDistance1 = Math.pow(effectiveDistance1, (-1/6.));
		//System.out.println("Effective distance for chain A = " + effectiveDistance1);
		
		// Compute effective distance for active atoms on the other chain
		for(int i = 0 ; i <  interfaceAtoms_chain2.size() ; i++)
		{	
			a2 = interfaceAtoms_chain2.get(i);
			//double c1 = 
				//PDBMoleculeComplex.residueConservationValue(residueETRankMap_chain2, a2.getResidueIndex());
			
			if(isActiveResidue(residueETRankMap_chain2, a2.getResidueIndex()))
			{	
				for(int j = 0 ; j < interfaceAtoms_chain1.size() ; j++)
				{	
					a1 = interfaceAtoms_chain1.get(j);
					distance = a2.distance(a1);
					double t = 1 / (Math.pow(distance, 6));
					effectiveDistance2 = effectiveDistance2 + t ;
				}
			}
			// Computing effective distance for the current active atom 
			// against all atoms on the opposing chain ends here
		}
		effectiveDistance2 = Math.pow(effectiveDistance2, (-1/6.));
		//System.out.println("Effective distance for chain B = " + effectiveDistance2);
		
		
		//System.out.println("Effective distance for chain A and B = " + (effectiveDistance1 + effectiveDistance2));
		
		return effectiveDistance1 + effectiveDistance2;
	}
	
	private static ArrayList<PDBAtom> getChainAtomsInInterface(ArrayList<PDBAtom> interfaceAtoms, 
			ArrayList<String> chainList)
	{
		// Find atoms in the interface that belong to chains in the given chainList
		ArrayList<PDBAtom> interfaceAtoms_chainList = new ArrayList<PDBAtom>();
		
		for(int i = 0 ; i < interfaceAtoms.size() ; i++)
		{
			PDBAtom a = interfaceAtoms.get(i);
			if(isIncludedInChainList(chainList, a.getChainID()))
			{
				interfaceAtoms_chainList.add(a);
			}
		}
		
		return interfaceAtoms_chainList;
	}
	
	private static Vector getA1NormalVector(int a1_index, int a2_index,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain1SurfaceMap,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain2SurfaceMap)
	{
		ArrayList<MSSurfaceEntry> a1SurfacePointEntries = 
			chain1SurfaceMap.get(new Integer(a1_index));
		ArrayList<MSSurfaceEntry> a2SurfacePointEntries = 
			chain2SurfaceMap.get(new Integer(a2_index));
		if(a1SurfacePointEntries == null || a2SurfacePointEntries == null)
		{
			return null;
		}
		else
		{
			MSSurfaceEntry optimalA1SurfacePointEntry = a1SurfacePointEntries.get(0);
			double minimumDistance = 999999999; // a very large double value
			double distance;
			for(int i = 0; i < a1SurfacePointEntries.size(); i++)
			{
				MSSurfaceEntry a1SurfacePointEntry = a1SurfacePointEntries.get(i);
				for(int j = 0; j < a2SurfacePointEntries.size(); j++)
				{
					MSSurfaceEntry a2SurfacePointEntry = a2SurfacePointEntries.get(j);
					distance = a1SurfacePointEntry.getPoint().distance(a2SurfacePointEntry.getPoint());
					if( distance < minimumDistance)
					{
						// Then update minimumDistance and optimalA1SurfacePointEntry
						minimumDistance = distance;
						optimalA1SurfacePointEntry = a1SurfacePointEntry;
					}
				}
			}
			return optimalA1SurfacePointEntry.getNormalVector();
		}
	}
	
	private static Vector getA2NormalVector(int a1_index, int a2_index,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain1SurfaceMap,
			HashMap<Integer, ArrayList<MSSurfaceEntry>> chain2SurfaceMap)
	{
		ArrayList<MSSurfaceEntry> a1SurfacePointEntries = 
			chain1SurfaceMap.get(new Integer(a1_index));
		ArrayList<MSSurfaceEntry> a2SurfacePointEntries = 
			chain2SurfaceMap.get(new Integer(a2_index));
		
		if(a1SurfacePointEntries == null || a2SurfacePointEntries == null)
		{
			return null;
		}
		else
		{
			MSSurfaceEntry optimalA2SurfacePointEntry = a2SurfacePointEntries.get(0);
			double minimumDistance = 999999999; // a very large double value
			double distance;
			for(int i = 0; i < a1SurfacePointEntries.size(); i++)
			{
				MSSurfaceEntry a1SurfacePointEntry = a1SurfacePointEntries.get(i);
				for(int j = 0; j < a2SurfacePointEntries.size(); j++)
				{
					MSSurfaceEntry a2SurfacePointEntry = a2SurfacePointEntries.get(j);
					distance = a1SurfacePointEntry.getPoint().distance(a2SurfacePointEntry.getPoint());
					if( distance < minimumDistance)
					{
						// Then update minimumDistance and optimalA1SurfacePointEntry
						minimumDistance = distance;
						optimalA2SurfacePointEntry = a2SurfacePointEntry;
					}
				}
			}
			return optimalA2SurfacePointEntry.getNormalVector();
		}
	}
	
	private static double residueConservationValue(HashMap<Integer, Double> residueETRankMap, int residueNumber)
	{
		ArrayList<Double> rankList = new ArrayList<Double>(residueETRankMap.values());
		
		double mean = Statistics.mean(rankList);
		double stdev = Statistics.standardDeviation(rankList);
		Double residueRankDouble = residueETRankMap.get(new Integer(residueNumber));
		if(residueRankDouble != null)
		{
			double residueRank = residueRankDouble.doubleValue();
			//double conservation = (MyMath.erf((residueRank - mean) / Math.sqrt(2) * stdev) + 1) / 2;
			//double conservation = (residueRank - mean) / stdev;
			double conservation = (mean - residueRank) / stdev;
			
			return conservation;
		}
		else
			return 0;
	}
	
	public double getRMSD(PDBMoleculeComplex otherProtein)
	{
		//int thisProteinAtomSize = this.getNumberOfAtoms();
		//int otherProteinAtomSize = otherProtein.getNumberOfAtoms();
		
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
	
	private static double[] convertToCoordinateArray(PDBMoleculeComplex protein)
	{
		// create the coordinates array
		// note that the size of the coordinates array is 3 times 
		// the number of atoms in the complex because the position of each
		// atom is represented as 3 double value (x, y, z coordinates)
		double[] coordinates = new double[3 * protein.getNumberOfAtoms()];
		ArrayList<PDBMolecule> molecules = protein.getMolecules();
		ListIterator<PDBMolecule> moleculeIter = molecules.listIterator();
		
		int i = 0;
		while(moleculeIter.hasNext())
		{
			PDBMolecule molecule = moleculeIter.next();
			if(molecule != null )
			{
				ArrayList<PDBAtom> atoms = molecule.getPDBAtomList();
				ListIterator<PDBAtom> atomIter = atoms.listIterator();
				while(atomIter.hasNext())
				{
					PDBAtom atom = atomIter.next();
					coordinates[i] = atom.getX();
					coordinates[i+1] = atom.getY();
					coordinates[i+2] = atom.getZ();
					
					i = i+3;
				}
			}
		}
		
		return coordinates;
	}
	
	public int getNumberOfMolecules()
	{
		ArrayList<PDBMolecule> molecules = this.getMolecules();
		return molecules.size();
	}
	
	public int getNumberOfAtoms()
	{
		int numberOfAtoms = 0;
		
		ArrayList<PDBMolecule> molecules = this.getMolecules();
		ListIterator<PDBMolecule> moleculeIter = molecules.listIterator();
		while(moleculeIter.hasNext())
		{
			PDBMolecule molecule = moleculeIter.next();
			if(molecule != null )
			{
				numberOfAtoms = numberOfAtoms + molecule.getPDBAtomList().size();
			}
		}
		
		return numberOfAtoms;
	}
	
	public ArrayList<AmberPDBAtom> getAmberPDBAtomsInMoleculeComplex(String amberTypeStr)
	{
		ArrayList<AmberPDBAtom> amberAtomsInMoleculeComplex = new ArrayList<AmberPDBAtom>();
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			if(pdbMolecule != null )
			{
				ArrayList<AmberPDBAtom> amberAtomsInMolecule =
					pdbMolecule.getAmberPDBAtomsInMolecule(amberTypeStr);
				
				ListIterator<AmberPDBAtom> iter2 =  amberAtomsInMolecule.listIterator();
				while(iter2.hasNext())
				{
					amberAtomsInMoleculeComplex.add(iter2.next());
				}
			}
		}
		return amberAtomsInMoleculeComplex;
	}
	
	public boolean isBonded(PDBAtom a1, PDBAtom a2)
	{
		if(a1.getChainID().equalsIgnoreCase(a2.getChainID()))
		{
			PDBMolecule molecule = this.getMolecule(a1.getChainID());
			if(molecule == null)
			{
				return false;
			}
			else
			{
				return molecule.isBonded(a1, a2);
			}
		}
		else
		{
			return false;
		}
	}
	
	public boolean isPeptideBond(PDBAtom a1, PDBAtom a2)
	{
		if(a1.getChainID().equalsIgnoreCase(a2.getChainID()))
		{
			PDBMolecule molecule = this.getMolecule(a1.getChainID());
			if(molecule == null)
			{
				return false;
			}
			else
			{
				return molecule.isPeptideBond(a1, a2);
			}
		}
		else
		{
			return false;
		}
	}
	
	public void setTerminalTypes()
	{
		ListIterator<PDBMolecule> iter = this.pdbMoleculeList.listIterator();
		while(iter.hasNext())
		{
			PDBMolecule pdbMolecule = iter.next();
			pdbMolecule.setTerminalTypes();
		}
	}
	
	/////////////////////////////////////////////
	//// The attributes of PDBMoleculeComplex class ////
	/////////////////////////////////////////////
	
	/**
	 * The list of PDB atoms in the molecule complex
	 */
	private ArrayList<PDBMolecule> pdbMoleculeList;
	
	/**
	 * The list of PDB atoms on the interface of an 
	 * arbitrary pair of chains in the molecule complex
	 */
	private ArrayList<PDBAtom> interfaceAtoms;
	
}
