/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package RRT;

/**
 *
 * @author amir
 */
import java.util.*;
import pdb.PDBAtom;
import pdb.PDBMolecule;
import javax.vecmath.Vector3d;
import java.util.AbstractList;

        
        
public class BBEnergy {
     int nres=5000;

   int[] resDensities;
  
  
  
    public static enum atomtype_t{
	// This is for the backbone + CB implementation
		CA, C, N, O, OT1, OT2, OXT, CB, NONE
	}
    
    
    
  public static String atom_types[] =   { "CA", "C" , "N" ,"O","OT1","OT2","OXT","CB"};
  public static double pi=3.1415926535897932;

  public void setNumberOfResidues(PDBMolecule obj )
  {
      nres=obj.numberOfResidues();
      resDensities=new int[nres];
  }
 public static enum amino_t
{
    ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, 
    MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK
}
public static String amino_names[] = 
    {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
     "HIS", "ILE", "LEU", "LYS", 
     "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "UNK"
    };
        
public static  amino_t get_amino(String  amino_name)
{
    for(int i = 0; i < NR_AMINOACIDS; i++)
    {
	if(amino_names[i].equals(amino_name) ) 
	    return amino_t.valueOf(amino_names[i]);
	
	//pfsgen renames HIS according to its state
	if((amino_name.equals("HSD" ))  || 
	   (amino_name.equals( "HSP")) ||
	   (amino_name.equals("HSE") ) )
	    return amino_t.HIS;
    }

    return amino_t.UNK;
}

public static  atomtype_t get_atype(String type_name)
{
    for(int i = 0; i < 8; i++)
    {
      if((atom_types[i].equals(type_name)) ) 
	return (atomtype_t.valueOf(atom_types[i]));
    }

    return atomtype_t.NONE;
}

// This is according to the parameters in the paper
 public static double BURIAL_GAMMA_TABLE[][] = 
{
    {0.84, 0.88, 0.57}, //ALA
    {0.94, 0.83, 0.13}, //ARG
    {0.96, 0.79, 0.25}, //ASN
    {0.98, 0.75, 0.20}, //ASP
    {0.67, 0.94, 0.66}, //CYS
    {0.96, 0.79, 0.24}, //GLN
    {0.97, 0.78, 0.16}, //GLU
    {0.94, 0.81, 0.34}, //GLY
    {0.92, 0.85, 0.13}, //HIS
    {0.78, 0.92, 0.55}, //ILE
    {0.78, 0.94, 0.46}, //LEU
    {0.98, 0.75, 0.00}, //LYS
    {0.82, 0.92, 0.46}, //MET
    {0.81, 0.94, 0.33}, //PHE
    {0.97, 0.76, 0.25}, //PRO
    {0.94, 0.79, 0.38}, //SER
    {0.92, 0.82, 0.40}, //THR
    {0.85, 0.91, 0.34}, //TRP
    {0.83, 0.92, 0.34}, //TYR
    {0.77, 0.93, 0.55}  //VAL
};
  
    
    
public static double  VDW_PENETRATION_FACTOR =0.9;
public static double  HUGE_VALUE =10000;
public static double  NO_ATOMS_RES= 5;
//for assigning hydrogen bonds
public static double  FREE_DONOR      =                     -1;
public static double  FREE_ACCEPTOR                        =-1;
public static double  DONOR_ACCEPTOR_RES_DISTANCE           =3;     //at least three residues apart
public static double  MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE =5.0;   //maximum distance in angstroms   
public static double  OPT_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE =3.5;   //maximum distance in angstroms
public static double  MAX_OUT_OF_PLANE_DIHEDRAL_FOR_HBOND  =0.698 ;  //maximum of 40 degrees
public static double  MIN_TWOBOND_ANGLE                    =1.5708;  //minimum of 90 degrees

public static double  BB_BB_HB_MAX_ENERGY                  =-0.5;   //kcal/mol
public static double  BB_SD_HB_MAX_ENERGY                  =-1.0;   //kcal/mol
public static double  DONOR_ACCEPTOR_SHORT_RANGE_INTER     =4;
public int[] hBondAssignedNDonors=new int[nres];

//for additional wolynes burial term
public static double  BURIAL_K        =5.0;
public static double  CONTACT_DENSITY =20.25 ;//angstroms. This is the original
// square contact density
public static double  LOW_BURIAL_MIN = 0 ;  //0 <= x
public static double  LOW_BURIAL_MAX = 3 ;  //x <= 3
public static double  MED_BURIAL_MIN = 3 ;  //3 <= x
public static double  MED_BURIAL_MAX = 6  ; //x <= 6
public static double  HIGH_BURIAL_MIN= 6 ;  //6 <= x 
public static double  HIGH_BURIAL_MAX =9 ;  //x <= 9

public static double  NR_AMINOACIDS       =  21; //20 + unknown type

//for additional wolynes water term
public static double  WATER_K                        =  5.0;
public static double  FIRST_WELL_WATER_LEFT_ENDPOINT   =4.5;
public static double  FIRST_WELL_WATER_RIGHT_ENDPOINT = 6.5;
public static double  SECOND_WELL_WATER_LEFT_ENDPOINT  =6.5;
public static double  SECOND_WELL_WATER_RIGHT_ENDPOINT =9.5;
public static double  RO_THRESHOLD_WATER       =  2.6;
public static double  NO_ATOM_RES= 5 ;// Number of atoms per residue
    
public static double maxDist = 40; // Maximum sq. distance for considered interactions

/* Force field parameters */
public static double CNSpringpublic  = 490;
public static double CACBSpringpublic  = 310;
public static double CACSpringpublic  = 317;
public static double COSpringpublic  = 570;
public static double CANSpringpublic  = 307;
public static double NCOAnglepublic  = 80;
public static double CCACBAnglepublic  = 63;
public static double CNCAAnglepublic  = 50;
public static double NCACBAnglepublic  = 80;
public static double CACOAnglepublic  = 80;
public static double CACNAnglepublic  = 70;
public static double NCACAnglepublic  = 63;
public static double CNDist = 1.335;
public static double CACBDist = 1.526;
public static double CACDist = 1.522;
public static double CODist =  1.229;
public static double CANDist = 1.449;
public static double NCOAngle = 2.145; // 122.6 degrees
public static double CCACBAngle = 1.939; //111.1 degrees
public static double CNCAAngle = 2.128; // 121.9 degrees
public static double NCACBAngle = 1.915; // 109.7 degrees
public static double CACOAngle = 2.101; // 120.4 degrees
public static double CACNAngle = 2.035; // 116.6 degrees
public static double NCACAngle = 1.922; // 110.1 degrees 




// for wolynes water energy term
//gamma_ij prot
public static double WATER_GAMMA_IJ_PROT_TABLE[][] = 
{ 
   //ALA   ARG   ASN   ASP   CYS   GLN    GLU   GLY   HIS   ILE
    {0.09, 0.04, 0.00, 0.00, 0.26, -0.02, 0.02, 0.04, 0.02, 0.12, 
     0.09, 0.02, 0.15, 0.31, -0.00, -0.00, 0.05, 0.08, 0.18, 0.32}, //ALA_ALL
   //LEU   LYS   MET   PHE   PRO    SER    THR   TRP   TYR   VAL

    {0.04, -0.04, -0.05, 0.02, 0.42, -0.03, -0.03, -0.00, -0.05, -0.03, //ARG ALL
     -0.07, -0.08, -0.16, -0.13, 0.01, 0.01, -0.01, -0.20, 0.14, 0.01},
   
    {0.00, -0.05, -0.03, -0.00, 0.16, -0.02, -0.03, 0.00, 0.00, -0.22,  
     -0.13, -0.05, -0.09, -0.11, -0.00, 0.00, -0.02, 0.08, 0.13, -0.10}, //ASN ALL

    {0.00, 0.02, -0.00, 0.00, -0.23, -0.03, -0.04, -0.02, 0.00, -0.18,
     -0.19, -0.02, -0.18, -0.19, -0.02, -0.00, -0.00, -0.12, 0.04, -0.14}, //ASP ALL
 
    {0.26, 0.42, 0.16, -0.23, 0.38, 0.16, 0.15, 0.38, 0.02, 0.32,
     0.31, -0.01, 0.73, 0.88, 0.39, 0.52, 0.33, 0.58, 0.51, 0.62}, //CYS ALL

    {-0.02, -0.03, -0.02, -0.03, 0.16, 0.03, -0.03, 0.01, 0.03, -0.09, //GLN ALL
     -0.13, -0.06, -0.12, 0.04, 0.01, -0.01, -0.03, -0.05, -0.09, 0.09},

    {0.02, -0.03, -0.03, -0.04, 0.15, -0.03, -0.04, -0.01, -0.04, -0.10, //GLU ALL
     -0.25, -0.03, -0.23, -0.15, -0.02, -0.01, -0.00, 0.00, -0.04, -0.12},

    {0.04, -0.00, 0.00, -0.02, 0.38, 0.01, -0.01, 0.08, -0.02, -0.05, 
     -0.05, -0.03, 0.21, -0.08, 0.06, 0.03, 0.01, -0.03, 0.07, 0.05}, //GLY ALL

    {0.02, -0.05, 0.00, 0.00, 0.02, 0.03, -0.04, -0.02, 0.11, -0.00,
     -0.00, -0.10, 0.08, 0.38, 0.03, 0.05, 0.02, 0.48, 0.34, 0.02}, //HIS ALL

    {0.12, -0.03, -0.22, -0.18, 0.32, -0.09, -0.10, -0.05, -0.00, 1.00,
     1.00, -0.09, 0.72, 0.93, -0.18, 0.01, 0.10, 0.34, 0.22, 0.68}, //ILE ALL
    
    {0.09, -0.07, -0.13, -0.19, 0.31, -0.13, -0.25, -0.05, -0.00, 1.00, 
     1.00, -0.07, 0.74, 0.69, -0.15, -0.13, 0.12, 0.38, 0.34, 0.63}, //LEU ALL

    {0.02, -0.08, -0.05, -0.02, -0.01, -0.06, -0.03, -0.03, -0.10, -0.09, 
     -0.07, -0.05, -0.14, -0.16, -0.02, -0.03, -0.03, -0.20, -0.28, -0.13}, //LYS ALL

    {0.15, -0.16, -0.09, -0.18, 0.73, -0.12, -0.23, 0.21, 0.08, 0.72, 
     0.74, -0.14, 0.32, 0.71, 0.01, 0.07, 0.05, 0.50, 0.27, 0.40}, //MET ALL

    {0.31, -0.13, -0.11, -0.19, 0.88, 0.04, -0.15, -0.08, 0.38, 0.93, 
     0.69, -0.16, 0.71, 1.00, -0.20, 0.00, 0.12, 0.65, 0.27, 0.82}, //PHE ALL

    {-0.00, 0.01, -0.00, -0.02, 0.39, 0.01, -0.02, 0.06, 0.03, -0.18, 
     -0.15,  -0.02, 0.01, -0.20, -0.00, -0.00, -0.01, 0.47, -0.07, -0.10}, // PRO ALL

    {-0.00, 0.01, 0.00, -0.00, 0.52, -0.01, -0.01, 0.03, 0.05, 0.01,
     -0.13, -0.03, 0.07, 0.00, -0.00, 0.02, -0.01, 0.10, -0.03, -0.00}, // SER ALL

    {0.05, -0.01, -0.02, -0.00, 0.33, -0.03, -0.00, 0.01, 0.02, 0.10, 
     0.12, -0.03, 0.05, 0.12, -0.01, -0.01, -0.01, 0.12, -0.05, -0.10}, //THR ALL

    {0.08, -0.20, 0.08, -0.12, 0.58, -0.05, 0.00, -0.03, 0.48, 0.34,
     0.38, -0.20, 0.50, 0.65, 0.47, 0.10, 0.12, 0.42, 0.15, 0.43}, //TRP ALL

    {0.18, 0.14, 0.13, 0.04, 0.51, -0.09, -0.04, 0.07, 0.34, 0.22, 
     0.34, -0.28, 0.27, 0.27, -0.07, -0.03, -0.05, 0.15,  0.21, 0.59}, //TYR ALL 
    
    {0.32, 0.01, -0.10, -0.14, 0.62, 0.09, -0.12, 0.05, 0.02, 0.68, 
     0.63, -0.13, 0.40, 0.82, -0.10, -0.00, -0.10, 0.43, 0.59, 0.73} //VAL ALL
};

//gamma_ij wat
public static double WATER_GAMMA_IJ_WAT_TABLE[][] = 
{ 
   //ALA   ARG   ASN   ASP   CYS   GLN    GLU   GLY   HIS   ILE
    {0.09, 0.04, 0.00, 0.00, 0.26, -0.02, 0.02, 0.04, 0.02, 0.12,
     0.09, 0.02, 0.15, 0.31, -0.00, -0.00, 0.05, 0.08, 0.18, 0.32}, //ALA ALL
   //LEU   LYS   MET   PHE   PRO    SER    THR   TRP   TYR   VAL

    {0.04, -0.04, -0.05, 0.02, 0.42, -0.03, -0.03, -0.00, -0.05, -0.03,
     -0.07, -0.08, -0.16, -0.13, 0.01, 0.01, -0.01, -0.20, 0.14, 0.01}, //ARG ALL

    {0.00, -0.05, -0.03, -0.00, 0.16, -0.02, -0.03, 0.00, 0.00, -0.22,
     -0.13, -0.05, -0.09, -0.11, -0.00, 0.00, -0.02, 0.08, 0.13, -0.10}, //ASN ALL

    {0.00, 0.02, -0.00, 0.00, -0.23, -0.03, -0.04, -0.02, 0.00, -0.18,
     -0.19, -0.02, -0.18, -0.19, -0.02, -0.00, -0.00, -0.12, 0.04, -0.14}, //ASP ALL

    {0.26, 0.42, 0.16, -0.23, 0.38, 0.16, 0.15, 0.38, 0.02, 0.32,
     0.31, -0.01, 0.73, 0.88, 0.39, 0.52, 0.33, 0.58, 0.51, 0.62},  //CYS ALL

    {-0.02, -0.03, -0.02, -0.03, 0.16, 0.03, -0.03, 0.01, 0.03, -0.09,
     -0.13, -0.06, -0.12, 0.04, 0.01, -0.01, -0.03, -0.05, -0.09, 0.09}, //GLN ALL

    {0.02, -0.03, -0.03,  -0.04, 0.15, -0.03, -0.04, -0.01, -0.04, -0.10,
     -0.25, -0.03, -0.23, -0.15, -0.02, -0.01, -0.00, 0.00, -0.04, -0.12},//GLU ALL

    {0.04, -0.00, 0.00, -0.02, 0.38, 0.01, -0.01, 0.08, -0.02, -0.05,
     -0.05, -0.03, 0.21, -0.08, 0.06, 0.03, 0.01, -0.03, 0.07, 0.05}, //GLY ALL

    {0.02, -0.05, 0.00, 0.00, 0.02, 0.03, -0.04, -0.02, 0.11, -0.00,
     -0.00, -0.10, 0.08, 0.38, 0.03, 0.05, 0.02, 0.48, 0.34, 0.02}, //HIS ALL

    {0.12, -0.03, -0.22, -0.18, 0.32, -0.09, -0.10, -0.05, -0.00, 1.00,
     1.00, -0.09, 0.72, 0.93, -0.18, 0.01, 0.10, 0.34, 0.22, 0.68}, //ILE ALL

    {0.09, -0.07, -0.13, -0.19, 0.31, -0.13, -0.25, -0.05, -0.00, 1.00,
     1.00, -0.07, 0.74, 0.69, -0.15, -0.13, 0.12, 0.38, 0.34, 0.63}, //LEU ALL

    {0.02, -0.08, -0.05, -0.02, -0.01, -0.06, -0.03, -0.03, -0.10, -0.09,
     -0.07, -0.05, -0.14, -0.16, -0.02, -0.03, -0.03, -0.20, -0.28, -0.13}, //LYS ALL

    {0.15, -0.16, -0.09, -0.18, 0.73, -0.12, -0.23, 0.21, 0.08, 0.72, 
     0.74, -0.14, 0.32, 0.71, 0.01, 0.07, 0.05, 0.50, 0.27, 0.40}, //MET ALL
    
    {0.31, -0.13, -0.11, -0.19, 0.88, 0.04, -0.15, -0.08, 0.38, 0.93,
     0.69, -0.16, 0.71, 1.00, -0.20, 0.00, 0.12, 0.65, 0.27, 0.82}, //PHE ALL

    {-0.00, 0.01, -0.00, -0.02, 0.39, 0.01, -0.02, 0.06, 0.03, -0.18,
     -0.15, -0.02, 0.01, -0.20, -0.00, -0.00, -0.01, 0.47, -0.07, -0.10}, //PRO ALL

    {-0.00, 0.01, 0.00, -0.00, 0.52, -0.01, -0.01, 0.03, 0.05, 0.01,
     -0.13, -0.03, 0.07, 0.00, -0.00, 0.02, -0.01, 0.10, -0.03, -0.00}, //SER ALL

    {0.05, -0.01, -0.02, -0.00, 0.33, -0.03, -0.00, 0.01, 0.02, 0.10, 
     0.12, -0.03, 0.05,  0.12, -0.01, -0.01, -0.01, 0.12, -0.05, -0.10}, //THR ALL

    {0.08, -0.20, 0.08, -0.12, 0.58, -0.05, 0.00, -0.03, 0.48, 0.34,
     0.38, -0.20, 0.50, 0.65, 0.47, 0.10, 0.12, 0.42, 0.15, 0.43}, //TRP ALL

    {0.18, 0.14, 0.13, 0.04, 0.51, -0.09, -0.04, 0.07, 0.34, 0.22, 
     0.34,  -0.28, 0.27, 0.27, -0.07, -0.03, -0.05, 0.15, 0.21, 0.59}, //TYR ALL
    
    {0.32, 0.01, -0.10, -0.14, 0.62, 0.09, -0.12, 0.05, 0.02, 0.68, 
     0.63, -0.13, 0.40, 0.82, -0.10, -0.00, -0.10, 0.43, 0.59, 0.73} //VAL ALL
};



public PDBAtom GetithResidueofList(PDBMolecule molecule, int i)
{
                ArrayList<PDBAtom> iter = molecule.getPDBAtomList();
                    PDBAtom firstAtom=iter.get(0);
                   int residueNumber= firstAtom.getResidueIndex();
                    
                for(int temp=0; temp < iter.size();temp++)
		{
			PDBAtom atom = iter.get(temp);
			if(atom.getResidueIndex() ==(i+residueNumber))
				return atom;
		}
		// if the control reaches here, this means that pdbAtomList
		// does not contain the given atom so return null
		return null;
}

//water energy function.
// Originally this was a method inside class PDBMolecule. 
// Therefore, all the indices (like in the dij) should be *(this)[resIndex] etc. 
public double waterE(PDBMolecule molcule )
{    
    Double water_energy = 0.0;
    int amino_i_code=0, amino_j_code=0;
   
    Double theta_two_ij = 0.0, theta_one_ij = 0.0;
    Double sigma_ij_prot = 0.0, sigma_ij_wat = 0.0;
    Double gamma_ij_prot = 0.0, gamma_ij_wat = 0.0;

    //compute ro_i (first well) for each residue i
    ArrayList<Double> ro_i= new ArrayList<Double>();
    for(int index=0;index<nres;index++)
    {
        ro_i.add(0.0) ;
    }
    
    Double dij = 0.0;
    //get number of residues
    
     nres=molcule.numberOfResidues();
    
    for(int resi = 0; resi < nres-1; resi++)
    {	
	for(int resj = resi + 1; resj < nres; resj++)
	{
		// This can be ignored for now. It's only if we need to update the water calculation or something. 
	//  if(!tree[sseIndices[resi]].isDirty() && 
	  //   !tree[sseIndices[resj]].isDirty()) continue;

	// it should be the molecule's resIndex. 
                 ArrayList<PDBAtom> ResidueObj=molcule.getResidueInMolecule(GetithResidueofList(molcule, resi).getResidueIndex());
                 PDBAtom atomObj=ResidueObj.get(3);
                 ResidueObj=molcule.getResidueInMolecule(GetithResidueofList(molcule, resj).getResidueIndex());
                 PDBAtom atomj=ResidueObj.get(0);
                 
    
            
	  dij =    atomObj.distance(atomj);
          

	// That's the actual formula. 
	    theta_one_ij = 
	      0.25 * (1 + Math.tanh(WATER_K * (dij - FIRST_WELL_WATER_LEFT_ENDPOINT))) * 
		(1 + Math.tanh(WATER_K * (FIRST_WELL_WATER_RIGHT_ENDPOINT - dij)));
            
            ro_i.set(resi,ro_i.get(resi)+  theta_one_ij);
             ro_i.set(resj,ro_i.get(resj)+  theta_one_ij);
           
            
	    
	}
    }
        
    double energy_term;
    for(int resi = 0; resi < nres; resi++)
    {
        
        
 
        ArrayList<PDBAtom> residue=molcule.getResidueInMolecule(GetithResidueofList(molcule, resi).getResidueIndex());
        PDBAtom  atomObj=residue.get(0);
        String AminoiName = atomObj.getAminoAcidType().toString();
	
	//amino_i_code = (*this)[resIndex[resi]].code;
	for(int resj = resi + 1; resj < nres; resj++)
	{
	 //if(!tree[sseIndices[resi]].isDirty() && 
	   //  !tree[sseIndices[resj]].isDirty()) continue;
	    //amino_j_code = (*this)[resIndex[resj]].code;
	      ArrayList<PDBAtom> residueJ=molcule.getResidueInMolecule(GetithResidueofList(molcule, resj).getResidueIndex());
                 PDBAtom  atomObj2=residueJ.get(0);
                 String AminoJName= atomObj2.getAminoAcidType().toString();
                 for(int i=0 ; i < 20;i++)
                 {
                     if(amino_names[i].equals(AminoJName))
                         amino_j_code=i;
                     if(amino_names[i].equals(AminoiName))
                         amino_i_code=i;
                 }
	
	    gamma_ij_prot = WATER_GAMMA_IJ_PROT_TABLE[amino_i_code][amino_j_code];
	    gamma_ij_wat = WATER_GAMMA_IJ_WAT_TABLE[amino_i_code][amino_j_code];
	    
	    sigma_ij_wat  = 0.5 * (1 - Math.tanh(WATER_K * (ro_i.get(resi) - RO_THRESHOLD_WATER))) * 
		            0.5 * (1 - Math.tanh(WATER_K * (ro_i.get(resj) - RO_THRESHOLD_WATER)));
	    sigma_ij_prot = 1.0 - sigma_ij_wat;

	    //theta_two_ij
            
                 
    
              //dij = (*this)[resIndex[resi]+3].dist((*this)[resIndex[resj]]);;
	  
	  dij =    residue.get(3).distance(atomObj2);
          
	    theta_two_ij = 0.25 * (1 + Math.tanh(WATER_K * (dij - SECOND_WELL_WATER_LEFT_ENDPOINT))) * 
		(1 + Math.tanh(WATER_K * (SECOND_WELL_WATER_RIGHT_ENDPOINT - dij)));
	    
	    energy_term = theta_two_ij * (sigma_ij_wat * gamma_ij_wat + sigma_ij_prot * gamma_ij_prot);
	    water_energy += energy_term;
	}
    }
    water_energy *= (-0.5);
    return water_energy;
}


// Soft vdW energy
public double calcDeBump(boolean force, PDBMolecule molecule)
{
    
    
  double d, result = 0;
  int starti,endi,startj,endj;
  int i,j,k,l,sz = molecule.numberOfResidues();
  for(i=0;i<nres-2;++i) {
     
    starti = GetithResidueofList(molecule, i).getResidueIndex();
    endi = GetithResidueofList(molecule, i+1).getResidueIndex();
  
    for(j=i+2;j<nres;++j) {
        ArrayList<PDBAtom> residuei=molecule.getResidueInMolecule(GetithResidueofList(molecule, i).getResidueIndex());
      ArrayList<PDBAtom> residuej=molecule.getResidueInMolecule(GetithResidueofList(molecule, j).getResidueIndex());
      
      
      d = residuei.get(0).distance(residuej.get(0));
      
      if(d > maxDist) continue;
      startj = GetithResidueofList(molecule, j).getResidueIndex();
      if(j < nres-1) endj =GetithResidueofList(molecule, j+1).getResidueIndex();
      else endj = sz;
    
      for(k=starti;k<endi;++k)
	for(l=startj;l<endj-1;++l)
	  result += pairwiseVDW(GetithResidueofList(molecule, k),GetithResidueofList(molecule, l), force);
           
     
    }   
  }
  
  return result;
}

// Find VDW between two atoms. If force = true, calculate force too (for minimization)
public double pairwiseVDW(PDBAtom i, PDBAtom j, boolean force)
{
  double r_a = 0.0, r_b = 0.0; //vdw radii
  double d_ab = HUGE_VALUE;      //actual Euclidean distance

  // added to make it as in amber
  // Make sure you understand these parameters. 
  double eps_a = 0.0, eps_b = 0.0; //epsilon values
  double eps_ab = 0.0;
  double common_term = 0.0, attraction = 0.0, repulsion = 0.0;
  r_a = i.getRmin2();
  eps_a  = i.getEps();
  //a and b vdw radii
  r_b = j.getRmin2();
	    
  //added to make it as in amber
  eps_b  = j.getEps();
  eps_ab = (eps_a * eps_b); // removed sqrt
		    
  //actual distance between atoms a and b
  //d_ab = i.dist(j);
    d_ab=i.distance(j);
    
  // This is just an approximation. If distance between atoms > than the sum of the VDW radii
  // do nothing. Basically check for collisions. 
  if(d_ab > r_a + r_b)
    return 0;
      
  //else, soft sphere potential is applied. 
  // Calculation is done to avoid exp function.
  common_term = (VDW_PENETRATION_FACTOR * (r_a + r_b))/d_ab;
  common_term = common_term * common_term * common_term; //3rd power
  attraction = common_term * common_term;          //6th power
  repulsion = attraction * attraction;          //12th power
      
  // If not force, that's it. Just calculate LJ potential
  double result = eps_ab * (repulsion - 2 * attraction);

 // Do I need to apply in reference object
  if(force) {
   // Vector3 axis = new Vector3();
  
    // Again, why do we need to copy over i and j? 
  
      Vector3d axis = new Vector3d(i.getX()-j.getX(),i.getY()-j.getY(),i.getZ()-j.getZ());
    
    // Remember that force is the derivative of the potential.
    double norm =axis.length();
    double inv = 1.0/norm;
    double mag = inv*12*eps_ab*(repulsion-attraction);

    axis.scale(mag*inv);
  /* // force  
    Vector3d force_i = new Vector3d(i.force);
    force_i.sub( axis);

    Vector3d force_j = new Vector3d(j.force);
    force_j.sub( axis);   */
  }
  return result;
}


public double calcHB(boolean updateHBondsFlag,boolean force,PDBMolecule molecule) {
  if(updateHBondsFlag) {
    resetNoHBonds();
    assignHBondsManyAcceptorstoOneDonor(molecule); //assignHBonds();
  }
  double hb_energy      = 0.0;

  double scaling_denom  = MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE - 
    OPT_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE;
  double energy_unit    = 0.0;    
  double scaling_factor = 1.0; //portion of maximum energy
  double proportionality_constant = 1.0; 
  
  //for force calculations
  // donor and acceptor positions, their distance
  int DonorAtIdx;

  // Go over all nitrogen atoms
 /*
  for(int DonorResIdx = 1; DonorResIdx <= nres; DonorResIdx++)
    {
        
      findN(DonorResIdx,DonorAtIdx,p_d);
      
      int AcceptorAtIdx = hBondAssignedNDonors[DonorResIdx-1];
      if(AcceptorAtIdx == FREE_DONOR)
	continue;
      if(!tree[sseIndices[DonorResIdx-1]].isDirty() && 
	 !tree[sseIndices[[AcceptorAtIdx].resid-1]].isDirty()) continue;
      
      int AcceptorResIdx = [AcceptorAtIdx].resid;
      double DonorAcceptorDistance         = HUGE_VALUE;
      double AcceptorWithDonorPeptideAngle = HUGE_VALUE;
      double twoBondAngle                  = HUGE_VALUE;
      if(!hbondGeomCriteriaMet(DonorAtIdx, AcceptorAtIdx, 
				  DonorResIdx, AcceptorResIdx,
				  DonorAcceptorDistance, 
				  AcceptorWithDonorPeptideAngle,
			       twoBondAngle)) continue ; //Should be assert!!!
      
      // Scale energy according to distance criteria
      if(DonorAcceptorDistance < OPT_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
	scaling_factor = 1.0;
      // too far away - discard
      else if(DonorAcceptorDistance > MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
	scaling_factor = 0.0;
      //else, linear scaling
      else
	scaling_factor = (MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE - DonorAcceptorDistance)/
	  scaling_denom;
      
      if([AcceptorAtIdx].isO())
	energy_unit = BB_BB_HB_MAX_ENERGY;
      else
        energy_unit = BB_SD_HB_MAX_ENERGY;	
      
      //a proportionality constant to favor long-range interactions
      if(abs(AcceptorResIdx - DonorResIdx) <= DONOR_ACCEPTOR_SHORT_RANGE_INTER)
        proportionality_constant = 1.0;
      else
        proportionality_constant = 2.0;
      
      hb_energy += proportionality_constant * scaling_factor * energy_unit;
      // Calculate forces
      if(force) {
	//OOPSMPdeclare_vector3(v_da);
        Vector3 v_da = new Vector3 ();
        Vector3 donor_coords = new Vector3([DonorAtIdx].m_coords);
        Vector3 acceptor_coords = new Vector3([AcceptorAtIdx].m_coords);
        v_da = vect_sub(donor_coords, acceptor_coords);
	double tmpval = proportionality_constant/(scaling_denom * DonorAcceptorDistance);
     
	//OOPSMPmultiply_scalar_vector3(v_da,tmpval,v_da);
	v_da = scalar_mult(v_da, tmpval);
	//OOPSMPadd_vector3((*this)[DonorAtIdx].force,v_da,
			  //(*this)[DonorAtIdx].force);
        Vector3 donor_force = new Vector3([DonorAtIdx].force);
        donor_force = vect_add(donor_force, v_da);
	//OOPSMPsub_vector3((*this)[AcceptorAtIdx].force,v_da,
			  //(*this)[AcceptorAtIdx].force);
       Vector3 acceptor_force = new Vector3([DonorAtIdx].force);
       acceptor_force = vect_add(acceptor_force, v_da);
      }
    }*/
  for(int DonorResIdx = 1; DonorResIdx < nres-1; DonorResIdx++)
    {
        
      PDBAtom AtomObj=GetithResidueofList(molecule, DonorResIdx);
      ArrayList<PDBAtom> ResidueObj=molecule.getResidueInMolecule(AtomObj.getResidueIndex());
      
      
             
      PDBAtom  AtomN= findN(ResidueObj);
      
      DonorAtIdx=AtomN.getIndex();
      Vector3d p_d = new Vector3d(AtomN.getX(), AtomN.getY(), AtomN.getZ());
       
      int AcceptorAtIdx = hBondAssignedNDonors[DonorResIdx-1];
      if(AcceptorAtIdx == FREE_DONOR || AcceptorAtIdx>nres-1)
	continue;
     // if(!tree[sseIndices[DonorResIdx-1]].isDirty() && 
	// !tree[sseIndices[[AcceptorAtIdx].resid-1]].isDirty()) continue;
      
      int AcceptorResIdx=GetithResidueofList(molecule, AcceptorAtIdx).getResidueIndex();
     
      if(AcceptorResIdx<0 || AcceptorResIdx>nres -1  )
          continue;
      double DonorAcceptorDistance         = HUGE_VALUE;
      double AcceptorWithDonorPeptideAngle = HUGE_VALUE;
      double twoBondAngle                  = HUGE_VALUE;
      if(!hbondGeomCriteriaMet(DonorAtIdx, AcceptorAtIdx, 
				  DonorResIdx, AcceptorResIdx,
				  DonorAcceptorDistance, 
				  AcceptorWithDonorPeptideAngle,
			       twoBondAngle, molecule)) continue ; //Should be assert!!!
      
      // Scale energy according to distance criteria
      if(DonorAcceptorDistance < OPT_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
	scaling_factor = 1.0;
      // too far away - discard
      else if(DonorAcceptorDistance > MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
	scaling_factor = 0.0;
      //else, linear scaling
      else
	scaling_factor = (MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE - DonorAcceptorDistance)/
	  scaling_denom;
      
      if(GetithResidueofList(molecule, AcceptorAtIdx).isO())
	energy_unit = BB_BB_HB_MAX_ENERGY;
      else
        energy_unit = BB_SD_HB_MAX_ENERGY;	
      
      //a proportionality constant to favor long-range interactions
      if(Math.abs(AcceptorResIdx - DonorResIdx) <= DONOR_ACCEPTOR_SHORT_RANGE_INTER)
        proportionality_constant = 1.0;
      else
        proportionality_constant = 2.0;
      
      hb_energy += proportionality_constant * scaling_factor * energy_unit;
      // Calculate forces
      if(force) {
	//OOPSMPdeclare_vector3(v_da);
        PDBAtom AtomObj2=GetithResidueofList(molecule, DonorAtIdx);
        
        Vector3d donor_coords = new Vector3d(AtomObj2.getX(),AtomObj2.getY(),AtomObj2.getZ() );
         AtomObj2=GetithResidueofList(molecule, AcceptorAtIdx);
        Vector3d acceptor_coords = new Vector3d(AtomObj2.getX(),AtomObj2.getY(),AtomObj2.getZ() );
        double tmpval = proportionality_constant/(scaling_denom * DonorAcceptorDistance);
     
        Vector3d v_da = new Vector3d ((donor_coords.x-acceptor_coords.x)*tmpval ,(donor_coords.y-acceptor_coords.y)*tmpval,(donor_coords.z-acceptor_coords.z )*tmpval);
       
      
        
	//OOPSMPmultiply_scalar_vector3(v_da,tmpval,v_da);
	//OOPSMPadd_vector3((*this)[DonorAtIdx].force,v_da,
			  //(*this)[DonorAtIdx].force);
/*        Vector3d donor_force = new Vector3d([DonorAtIdx].force);
        donor_force = vect_add(donor_force, v_da);
	//OOPSMPsub_vector3((*this)[AcceptorAtIdx].force,v_da,
			  //(*this)[AcceptorAtIdx].force);
       Vector3d acceptor_force = new Vector3d([DonorAtIdx].force);
       acceptor_force = vect_add(acceptor_force, v_da);*/
      }
    }
    
  return hb_energy;
}
public PDBAtom findN(ArrayList<PDBAtom> residue)
{
    if(residue.isEmpty())
        return null;
    for(PDBAtom obj :residue)
    {
        if(obj.isN())
            return obj;
        System.out.println(String.valueOf(obj.getResidueIndex()));
    }
    return null;
}
public PDBAtom findO(ArrayList<PDBAtom> residue)
{
    for(PDBAtom obj :residue)
    {
        if(obj.isO())
            return obj;
    }
    return null;
}
public void resetNoHBonds()
{
    for(int i = 0; i < nres; i++)
	hBondAssignedNDonors[i] = -1;	    
}
public PDBAtom findCA(ArrayList<PDBAtom> residue)
{
    for(PDBAtom obj :residue)
    {
        if(obj.isCA())
            return obj;
    }
    return null;
}
public PDBAtom findC(ArrayList<PDBAtom> residue)
{
    for(PDBAtom obj :residue)
    {
        if(obj.isC())
            return obj;
    }
    return null;
}
public PDBAtom findCB(ArrayList<PDBAtom> residue)
{
    for(PDBAtom obj :residue)
    {
        if(obj.isCB())
            return obj;
    }
    return null;
}
public void assignHBondsManyAcceptorstoOneDonor(PDBMolecule molecule)
{
    double DonorAcceptorDistance         = HUGE_VALUE;
    double AcceptorWithDonorPeptideAngle = HUGE_VALUE;
    double twoBondAngle                  = HUGE_VALUE;
  
    double lowestDonorAcceptorDistance = HUGE_VALUE;
    int Nidx, Oidx, CBidx; 
  
   
  
    // For each residue find closest acceptor to given N donor
    for(int i = 0; i < nres-1; i++)
    {
	//this donor has been assigned already
	if(hBondAssignedNDonors[i] != FREE_DONOR)
	    continue;	
	// Find the N atom of residue i
         PDBAtom AtomObj=GetithResidueofList(molecule, i+1);
         ArrayList<PDBAtom> ResidueObj=molecule.getResidueInMolecule(AtomObj.getResidueIndex());
        PDBAtom  AtomN= findN(ResidueObj);
      
        Nidx=AtomN.getIndex();
        Vector3d NidxPos=new Vector3d(AtomN.getX(), AtomN.getY(), AtomN.getZ());
	
            
	lowestDonorAcceptorDistance = HUGE_VALUE;

	for(int j = 0; j < nres-1; j++)
	  {
	    // Donors and acceptors must be at least 3 residues apart
	    if(Math.abs(i - j) < DONOR_ACCEPTOR_RES_DISTANCE)
		continue;

	    //check criteria
               AtomObj=GetithResidueofList(molecule, j+1);
               ResidueObj=molecule.getResidueInMolecule(AtomObj.getResidueIndex());
               PDBAtom  AtomO= findO(ResidueObj);
      
                Oidx=AtomO.getIndex();  
                Vector3d OidxPos=new Vector3d(AtomO.getX(), AtomO.getY(), AtomO.getZ());
	
	    // Check that criteria met for this pair
	    if(hbondGeomCriteriaMet(Nidx, Oidx, i+1, j+1,
				    DonorAcceptorDistance,
				    AcceptorWithDonorPeptideAngle,
				    twoBondAngle, molecule))
	    {
		if(DonorAcceptorDistance < lowestDonorAcceptorDistance)
		{
		    lowestDonorAcceptorDistance = DonorAcceptorDistance;
		    hBondAssignedNDonors[i]    = Oidx;
		}
	    }
	  }

	//if N was assigned, do not bother with side chains
	if(hBondAssignedNDonors[i] != FREE_DONOR)
	    continue;

	for(int j = 0; j < nres-1; j++) //else
	{	 
	    if(Math.abs(i - j) < DONOR_ACCEPTOR_RES_DISTANCE)
		continue;
	    // Just to get an atom
            AtomObj=GetithResidueofList(molecule, j+1);
            ResidueObj=molecule.getResidueInMolecule(AtomObj.getResidueIndex());
             PDBAtom  AtomO2= findO(ResidueObj);
      
           int Oidx2=AtomO2.getIndex();
           Vector3d OidxPos2=new Vector3d(AtomO2.getX(), AtomO2.getY(), AtomO2.getZ());
	
	   //go over Ser, Thr, Asn, Asp, Gln, and Glu
	    //if(!(*this)[Oidx2].isSideChainConsideredForHBond())
	if(!(AtomO2.isSideChainConsideredForHBond()))
		           continue;

	    //check criteria
              AtomObj=GetithResidueofList(molecule, j+1);
              ResidueObj=molecule.getResidueInMolecule(AtomObj.getResidueIndex());
             PDBAtom  AtomCB= findO(ResidueObj);
      
             CBidx=AtomCB.getIndex();
             Vector3d CBidxPos=new Vector3d(AtomCB.getX(), AtomCB.getY(), AtomCB.getZ());
	
	    assert(CBidx >= 0);

	    if(hbondGeomCriteriaMet(Nidx, CBidx, i, j,
				    DonorAcceptorDistance, 
				    AcceptorWithDonorPeptideAngle,
				    twoBondAngle, molecule))
	    {
		hBondAssignedNDonors[i]     = CBidx;
	    }
	}
    }
}

// Find out whether distance/angle criteria are satisfied. 
public boolean hbondGeomCriteriaMet( int DonorAtIdx,
		int AcceptorAtIdx, int DonorResIdx, int AcceptorResIdx,
		double  DonorAcceptorDistance, double  AcceptorWithDonorPeptideAngle, double  twoBondAngle,PDBMolecule molecule)
{

// m_coords are the coordinates of the atoms. Replace it with the function that retrieves the coordinates. 
// Second thought - why do you need to copy them over? 
// You only need them for distance or angle, it's just a waste of time and space. 
    
    
   ArrayList<PDBAtom> DonorResidueObj=  molecule.getResidueInMolecule(GetithResidueofList(molecule, DonorResIdx).getResidueIndex());
   
  PDBAtom NAtom=findN(DonorResidueObj);
  
 Vector3d pDonor = new Vector3d(NAtom.getX(),NAtom.getY(),NAtom.getZ());
  
   
   ArrayList<PDBAtom> AcceptorResidueObj=  molecule.getResidueInMolecule(GetithResidueofList(molecule, AcceptorResIdx).getResidueIndex());
   
  PDBAtom OAtom=findO(AcceptorResidueObj);
  
 
 Vector3d pAcceptor = new Vector3d(OAtom.getX(),OAtom.getY(),OAtom.getZ());
 
  //can be O of bb or CB of side chain
  // get Distance between atoms
 
  DonorAcceptorDistance = NAtom.distance(OAtom);
  // Too far apart to be a proper H-bond
  if(DonorAcceptorDistance > MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
    return false; 

  // Check angle between Donor, acceptor and atom attached to acceptor. For example - C=O and N. 
  int idx;
  
  Vector3d CAofDonorRes = new Vector3d(findCA(DonorResidueObj).getX(),findCA(DonorResidueObj).getY(),findCA(DonorResidueObj).getZ());
  Vector3d CofPrevDonorRes = new Vector3d();
 
  // Find the CA atom of the donor residue, and the C of the residue or the one before. 

  PDBAtom Catom;
  
  if(DonorResIdx > 1)
  {
      ArrayList<PDBAtom> PrevDonorResidueObj=  molecule.getResidueInMolecule(GetithResidueofList(molecule, DonorResIdx-1).getResidueIndex());
         Catom=findC(PrevDonorResidueObj);
     
  }
  else
  {
        Catom=findC(DonorResidueObj);
 
       
  }
     CofPrevDonorRes.x=Catom.getX();
        CofPrevDonorRes.y =Catom.getY();
        CofPrevDonorRes.z=Catom.getZ() ;
// Please check you indeed look at the correct atoms making up the angle

  
   pAcceptor.sub(pDonor);
    CofPrevDonorRes.sub(pDonor);
  CAofDonorRes.sub(CAofDonorRes,pDonor);
  //additional test from Linus not mentioned in paper

 // Calculate angle
  twoBondAngle = pAcceptor.angle(CAofDonorRes);
  AcceptorWithDonorPeptideAngle = dihedral_vector3(pAcceptor,CofPrevDonorRes,CAofDonorRes);

  if(AcceptorWithDonorPeptideAngle > MAX_OUT_OF_PLANE_DIHEDRAL_FOR_HBOND ||
     twoBondAngle < MIN_TWOBOND_ANGLE)
    return false;
  return true;
}

// Calculate dihedral angle for a residue. 
// The phi boolean variable determines if it's phi or psi.
 public double getDihe(int res, boolean phi, PDBMolecule molecule) {
       
     ArrayList<PDBAtom> PreviouseResidue=molecule.getResidueInMolecule(res-1);
     ArrayList<PDBAtom> CurrentResidue=molecule.getResidueInMolecule(res);
     ArrayList<PDBAtom> NextResidue=molecule.getResidueInMolecule(res+1);
       
       Vector3d v1 = new Vector3d();                        
       Vector3d v2 = new Vector3d();
       Vector3d v3 = new Vector3d();
       Vector3d v4 = new Vector3d();
       int pos;
       if(phi) { // the angle around N-CA
	// Find relevant atoms
        v1.x=findC(PreviouseResidue).getX();        
        v1.y=findC(PreviouseResidue).getY();
        v1.z=findC(PreviouseResidue).getZ();
        v2.x=findN(CurrentResidue).getX();
        v2.y=findN(CurrentResidue).getY();
        v2.z=findN(CurrentResidue).getZ();
        v3.x=findCA(CurrentResidue).getX();
        v3.y=findCA(CurrentResidue).getY();
        v3.z=findCA(CurrentResidue).getZ();
        v4.x=findC(CurrentResidue).getX();
        v4.y=findC(CurrentResidue).getY();
        v4.z=findC(CurrentResidue).getZ();
        
        
      }
      else { // Psi - the angle around CA-C
         
        v1.x=findN(CurrentResidue).getX();
        v1.y=findN(CurrentResidue).getY();
        v1.z=findN(CurrentResidue).getZ();
        v2.x=findCA(CurrentResidue).getX();
        v2.y=findCA(CurrentResidue).getY();
        v2.z=findCA(CurrentResidue).getZ();
        v3.x=findC(CurrentResidue).getX();
        v3.y=findC(CurrentResidue).getY();
        v3.z=findC(CurrentResidue).getZ();
         v4.x=findN(NextResidue).getX();
        v4.y=findN(NextResidue).getY();
        v4.z=findN(NextResidue).getZ();
        
      }
       v1.sub(v2);
       v2.sub(v3);
       v3.sub(v4);
	// Please make sure it works, because it doesn't look like it's normalized!
       return dihedral_vector3(v1,v2,v3);
}

// Calculate dihedral angle between three vectors (A-B, B-C, C-D?)
// The convention is that the values go from [-pi,pi]
// You have to make sure the vectors are normalized.
 double dihedral_vector3(Vector3d v1, Vector3d v2, Vector3d v3){
				      
    Vector3d normal = new Vector3d();
    //normal = cross_prod(v2,v3);
    normal.cross(v2,v3);
    
    
    
    Vector3d normal1 = new Vector3d();
    normal1.cross(v1,v2);

    double alpha = normal1.angle(normal);
    return alpha > pi ? -alpha : alpha;
}                                                  

// This is the harmonic force between atoms. It's needed for energy minimization
 public double calcForceAngle(PDBAtom NorCBAtom,PDBAtom CAAtom, PDBAtom CBorC, double springConst, double eqangle){
    
  
     double e;
       
    //   Vector3d pos1 = new Vector3d(NorCBAtom.getX(),NorCBAtom.getY(),NorCBAtom.getZ());                        
      // Vector3d pos2 = new Vector3d(pos2.m_coords);
       //Vector3d pos3 = new Vector3d(pos2.m_coords);
       double d12 = NorCBAtom.distance(CAAtom);
       double d23 =CAAtom.distance(CBorC);
       
       Vector3d r12 = new Vector3d(NorCBAtom.getX()-CAAtom.getX(),NorCBAtom.getY()-CAAtom.getY(),NorCBAtom.getZ()-CAAtom.getZ());
      Vector3d  r23 = new Vector3d(CAAtom.getX()-CBorC.getX(),CAAtom.getY()-CBorC.getY(),CAAtom.getZ()-CBorC.getZ());
        Vector3d tmp1 =r12;
       Vector3d tmp2 = r23;
     
       double angle= r12.angle(r23);
       double d12inv = 1.0/d12, d23inv = 1.0/d23;
       double sin_theta = Math.sin(angle);
       double cos_theta = Math.cos(angle);
       double middle = angle-eqangle;
       middle *= (-2*springConst)/(sin_theta);
       double c1 = middle*d12inv;
       double c2 = middle*d23inv;
       e = springConst*(eqangle-angle)*(eqangle-angle);
  // Calculate the actual forces
       tmp1.scale(d12inv*cos_theta);
       tmp2.scale(d23inv);
       tmp2.sub(tmp1);
       Vector3d f1 = tmp2;
       f1.scale(c1);
      
      r23.scale(d23inv*cos_theta);
       r12.scale(d12inv);
       tmp2.sub(tmp1);
       Vector3d f3 = tmp2;
       f3.scale(c2);
        Vector3d f2 = f1;
                     
    
       
      f2.add(f3);
      /*
       NorCBAtom.force = vect_sub(NorCBAtom.force,f1);
       CAAtom.force = vect_sub(CAAtom.force,f2);
       CBorC.force = vect_sub(CBorC.force,f3);*/
       return e;                                                                                                          }


// burial energy function
public double burialE( PDBMolecule mol)
{
  double burial_energy = 0.0;

  amino_t amino_code;

     updateResDensities(mol);
  double low_S = -1, med_S = -1, high_S = -1;
  double ref_low_S = 0, ref_med_S = 0, ref_high_S = 0;
  for(int i = 0; i < nres; i++) {
    int ro_i = resDensities[i];
    low_S = Math.tanh(BURIAL_K * (ro_i - LOW_BURIAL_MIN)) + Math.tanh(BURIAL_K * (LOW_BURIAL_MAX - ro_i));
    med_S = Math.tanh(BURIAL_K * (ro_i - MED_BURIAL_MIN)) + Math.tanh(BURIAL_K * (MED_BURIAL_MAX - ro_i));
    high_S = Math.tanh(BURIAL_K * (ro_i - HIGH_BURIAL_MIN)) + Math.tanh(BURIAL_K * (HIGH_BURIAL_MAX - ro_i));
	
    int amino_codes=0;
        ArrayList<PDBAtom> residueJ=mol.getResidueInMolecule(GetithResidueofList(mol, i).getResidueIndex());
                 PDBAtom  atomObj2=residueJ.get(0);
                 String AminoJName= atomObj2.getAminoAcidType().toString();
                 for(int k=0 ; k < 20;k++)
                 {
                     if(amino_names[k].equals(AminoJName))
                         amino_codes=k;
         
                 }
    double low_gamma = BURIAL_GAMMA_TABLE[amino_codes][0];
    double med_gamma = BURIAL_GAMMA_TABLE[amino_codes][1];
    double high_gamma = BURIAL_GAMMA_TABLE[amino_codes][2];

	
    //ref is extended, where each residue averages density of 2 contact residues
    ref_low_S = Math.tanh(BURIAL_K * (2 - LOW_BURIAL_MIN)) + Math.tanh(BURIAL_K * (LOW_BURIAL_MAX - 2));
    ref_med_S = Math.tanh(BURIAL_K * (2 - MED_BURIAL_MIN)) + Math.tanh(BURIAL_K * (MED_BURIAL_MAX - 2));
    ref_high_S = Math.tanh(BURIAL_K * (2 - HIGH_BURIAL_MIN)) + Math.tanh(BURIAL_K * (HIGH_BURIAL_MAX - 2));		
    double ref_resContr = low_gamma * ref_low_S + med_gamma * ref_med_S + high_gamma * ref_high_S;
	
    double resContr = low_gamma * low_S + med_gamma * med_S + high_gamma * high_S - ref_resContr;

    burial_energy += resContr;
  }    
  return burial_energy;
}


/**
 ** for Wolynes burial and contact term
 **/
public void resetResDensities()
{
    for(int i = 0; i < nres; i++)
	resDensities[i] = 0;
}

//compute densities
public void updateResDensities(PDBMolecule mol)
{    
  resetResDensities();
  for(int resi = 0; resi < nres-1; resi++)
    {	
      for(int resj = resi + 1; resj < nres; resj++)
	{  
            try{
                if(areResInContact(resi, resi,mol))
                {


                  resDensities[resi]++;
                  resDensities[resj]++;
                }
            }
            catch(Exception es){
            }
	}
    }
}

public boolean areResInContact(int resi,  int resj,PDBMolecule mol)
{
    
    
  int starti =   GetithResidueofList(mol, resi).getResidueIndex();
  int startj =  GetithResidueofList(mol, resj).getResidueIndex();
  int endi,endj;
  if(resi >= nres-1)
    endi = nres-1;
  else
    endi = GetithResidueofList(mol, resi+1).getResidueIndex();
  if(resj >= nres-1)
    endj =nres-1;
  else endj =  GetithResidueofList(mol, resj+1).getResidueIndex();
  for(int i=starti;i<endi;++i)
    for(int j=startj;j<endj;++j) {
      double d =  GetithResidueofList(mol, i).getResidueIndex();
       
      if (d <= CONTACT_DENSITY) return true;
    }
  return false;
}
}

       
 
