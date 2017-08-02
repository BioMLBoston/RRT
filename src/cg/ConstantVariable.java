/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cg;

import java.util.ArrayList;
import pdb.PDBAtom;
import pdb.PDBMolecule;

/**
 *
 * @author amir
 */
public class ConstantVariable  {
     int nres=10;
    public static enum atomtype_t{
		// B is for strong attraction
		// V is for weak attraction
		// N is for weak repulsion
		// L is for strong repulsion
		CA, C, N, O, OT1, OT2, OXT, CB, NONE
	}
    
    
    

  public static String atom_types[] =   { "CA", "C" , "N" ,"O","OT1","OT2","OXT","CB"};

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



//water energy function
public double waterE()
{    
   double water_energy = 0.0;
    
    amino_t amino_i_code, amino_j_code;
   
    double theta_two_ij = 0.0, theta_one_ij = 0.0;
    double sigma_ij_prot = 0.0, sigma_ij_wat = 0.0;
    double gamma_ij_prot = 0.0, gamma_ij_wat = 0.0;

    //compute ro_i (first well) for each residue i
    ArrayList<Double> ro_i= new ArrayList<Double>();
    double dij = 0.0;
    for(int resi = 0; resi < nres; resi++)
    {	
	for(int resj = resi + 1; resj < nres; resj++)
	{
	  if(!tree[sseIndices[resi]].isDirty() && 
	     !tree[sseIndices[resj]].isDirty()) continue;
	  dij = (*this)[resIndex[resi]+3].dist((*this)[resIndex[resj]]);

	    theta_one_ij = 
	      0.25 * (1 + tanh(WATER_K * (dij - FIRST_WELL_WATER_LEFT_ENDPOINT))) * 
		(1 + tanh(WATER_K * (FIRST_WELL_WATER_RIGHT_ENDPOINT - dij)));
	    ro_i[resi] += theta_one_ij;
	    ro_i[resj] += theta_one_ij;
	}
    }
        
    double energy_term;
    for(int resi = 0; resi < nres; resi++)
    {
	amino_i_code = (*this)[resIndex[resi]].code;
	for(int resj = resi + 1; resj < nres; resj++)
	{
	  if(!tree[sseIndices[resi]].isDirty() && 
	     !tree[sseIndices[resj]].isDirty()) continue;
	    amino_j_code = (*this)[resIndex[resj]].code;
	    
	    gamma_ij_prot = WATER_GAMMA_IJ_PROT_TABLE[amino_i_code][amino_j_code];
	    gamma_ij_wat = WATER_GAMMA_IJ_WAT_TABLE[amino_i_code][amino_j_code];
	    
	    sigma_ij_wat  = 0.5 * (1 - tanh(WATER_K * (ro_i[resi] - RO_THRESHOLD_WATER))) * 
		            0.5 * (1 - tanh(WATER_K * (ro_i[resj] - RO_THRESHOLD_WATER)));
	    sigma_ij_prot = 1.0 - sigma_ij_wat;

	    //theta_two_ij
	    dij = (*this)[resIndex[resi]+3].dist((*this)[resIndex[resj]]);;
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
public double calcDeBump(boolean force)
{
    
  double d, result = 0;
  int starti,endi,startj,endj;
  int i,j,k,l,sz = size();
  for(i=0;i<nres-2;++i) {
    starti = resIndex[i];
    endi = resIndex[i+1];
    for(j=i+2;j<nres;++j) {
      if(!force && !tree[sseIndices[i]].isDirty() && 
	 !tree[sseIndices[j]].isDirty()) continue;
      // Only non bonded
      d = OOPSMPsquare_distance_vector3((*this)[resIndex[i]].m_coords,
					(*this)[resIndex[j]].m_coords);
      if(d > maxDist) continue;
      startj = resIndex[j];
      if(j != nres-1) endj = resIndex[j+1];
      else endj = sz;
      for(k=starti;k<endi;++k)
	for(l=startj;l<endj;++l)
	  result += pairwiseVDW((*this)[k],(*this)[l], force);
    }   
  }
  return result;
}


public double pairwiseVDW(PDBAtom i, PDBAtom j,boolean force)
{
  double r_a = 0.0, r_b = 0.0; //vdw radii
  double d_ab = HUGE_VALUE;      //actual Euclidean distance

  //added to make it as in amber
  double eps_a = 0.0, eps_b = 0.0; //epsilon values
  double eps_ab = 0.0;
  double common_term = 0.0, attraction = 0.0, repulsion = 0.0;
  r_a = i.rmin2;
  eps_a  = i.eps;
  //a and b vdw radii
  r_b = j.rmin2;
	    
  //added to make it as in amber
  eps_b  = j.eps;
  eps_ab = (eps_a * eps_b); // removed sqrt
		    
  //actual distance between atoms a and b
  d_ab = i.dist(j);
		    
  if(d_ab > r_a + r_b)
    return 0;
      
  //else, soft sphere potential is applied
  common_term = (VDW_PENETRATION_FACTOR * (r_a + r_b))/d_ab;
  common_term = common_term * common_term * common_term; //3rd power
  attraction = common_term * common_term;          //6th power
  repulsion = attraction * attraction;          //12th power
      
  double result = eps_ab * (repulsion - 2 * attraction);
  if(force) {
    OOPSMPdeclare_vector3(axis);
    OOPSMPsub_vector3(i.m_coords,j.m_coords,axis);
    double norm = OOPSMPnorm_vector3(axis);
    double inv = 1.0/norm;
    double mag = inv*12*eps_ab*(repulsion-attraction);
    OOPSMPmultiply_scalar_vector3(axis,mag*inv,axis);
    OOPSMPsub_vector3(i.force,axis,i.force);
    OOPSMPadd_vector3(j.force,axis,j.force);   
  }
  return result;
}
public double calcHB(boolean updateHBondsFlag,boolean force) {
  if(updateHBondsFlag) {
    resetNoHBonds();
    assignHBondsManyAcceptorstoOneDonor(); //assignHBonds();
  }
  double hb_energy      = 0.0;

  double scaling_denom  = MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE - 
    OPT_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE;
  double energy_unit    = 0.0;    
  double scaling_factor = 1.0; //portion of maximum energy
  double proportionality_constant = 1.0; 
  
  //for force calculations
  // donor and acceptor positions, their distance
  OOPSMPdeclare_vector3(p_d);
  int DonorAtIdx;

  // Go over all nitrogen atoms
  for(int DonorResIdx = 1; DonorResIdx <= nres; DonorResIdx++)
    {
      findN(DonorResIdx,DonorAtIdx,p_d);
      int AcceptorAtIdx = hBondAssignedNDonors[DonorResIdx-1];
      if(AcceptorAtIdx == FREE_DONOR)
	continue;
      if(!tree[sseIndices[DonorResIdx-1]].isDirty() && 
	 !tree[sseIndices[(*this)[AcceptorAtIdx].resid-1]].isDirty()) continue;
      
      int AcceptorResIdx = (*this)[AcceptorAtIdx].resid;
      double DonorAcceptorDistance         = HUGE_VAL;
      double AcceptorWithDonorPeptideAngle = HUGE_VAL;
      double twoBondAngle                  = HUGE_VAL;
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
      
      if((*this)[AcceptorAtIdx].isO())
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
	OOPSMPdeclare_vector3(v_da);
	OOPSMPsub_vector3((*this)[DonorAtIdx].m_coords,
			  (*this)[AcceptorAtIdx].m_coords,v_da);
	double tmpval = proportionality_constant/(scaling_denom * DonorAcceptorDistance);
	OOPSMPmultiply_scalar_vector3(v_da,tmpval,v_da);
	
	OOPSMPadd_vector3((*this)[DonorAtIdx].force,v_da,
			  (*this)[DonorAtIdx].force);
	OOPSMPsub_vector3((*this)[AcceptorAtIdx].force,v_da,
			  (*this)[AcceptorAtIdx].force);
      }
    }
  return hb_energy;
}

public void resetNoHBonds()
{
    for(int i = 0; i < nres; i++)
	hBondAssignedNDonors[i] = FREE_DONOR;	    
}


public void assignHBondsManyAcceptorstoOneDonor()
{
    double DonorAcceptorDistance         = HUGE_VAL;
    double AcceptorWithDonorPeptideAngle = HUGE_VAL;
    double twoBondAngle                  = HUGE_VAL;
  
    double lowestDonorAcceptorDistance = HUGE_VAL;
    int Nidx, Oidx, CBidx;
    OOPSMPdeclare_vector3(NidxPos); 
    OOPSMPdeclare_vector3(OidxPos); 
    OOPSMPdeclare_vector3(CBidxPos); 
    // For each residue find closest acceptor to given N donor
    for(int i = 0; i < nres; i++)
    {
	//this donor has been assigned already
	if(hBondAssignedNDonors[i] != FREE_DONOR)
	    continue;	
	// Find the N atom of residue i
	findN(i+1,Nidx,NidxPos);

	lowestDonorAcceptorDistance = HUGE_VALUE;

	for(int j = 0; j < nres; j++)
	  {
	    // Donors and acceptors must be at least 3 residues apart
	    if(Math.abs(i - j) < DONOR_ACCEPTOR_RES_DISTANCE)
		continue;

	    //check criteria
	    findO(j+1,Oidx,OidxPos);

	    // Check that criteria met for this pair
	    if(hbondGeomCriteriaMet(Nidx, Oidx, i+1, j+1,
				    DonorAcceptorDistance,
				    AcceptorWithDonorPeptideAngle,
				    twoBondAngle))
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

	for(int j = 0; j < nres; j++) //else
	{	 
	    if(Math.abs(i - j) < DONOR_ACCEPTOR_RES_DISTANCE)
		continue;
	    // Just to get an atom
	    findO(j+1,Oidx,OidxPos);
	    //go over Ser, Thr, Asn, Asp, Gln, and Glu
	    if(!(*this)[Oidx].consideredSidechainForHBond())
		continue;

	    //check criteria
	    findCB(j+1,CBidx,CBidxPos);
	    assert(CBidx >= 0);

	    if(hbondGeomCriteriaMet(Nidx, CBidx, i, j,
				    DonorAcceptorDistance, 
				    AcceptorWithDonorPeptideAngle,
				    twoBondAngle))
	    {
		hBondAssignedNDonors[i]     = CBidx;
	    }
	}
    }
}

public boolean hbondGeomCriteriaMet( int DonorAtIdx, 
				        int AcceptorAtIdx, 
				        int DonorResIdx, 
				        int AcceptorResIdx,
				       double& DonorAcceptorDistance,
				       double& AcceptorWithDonorPeptideAngle,
				       double& twoBondAngle)
{
  OOPSMPdeclare_vector3(pDonor);
  OOPSMPdeclare_vector3(pAcceptor);
  OOPSMPcopy_vector3(pDonor,(*this)[DonorAtIdx].m_coords);
  OOPSMPcopy_vector3(pAcceptor,(*this)[AcceptorAtIdx].m_coords);
  //can be O of bb or CB of side chain
  // get Distance between atoms
  DonorAcceptorDistance = OOPSMPdistance_vector3(pDonor, pAcceptor);

  // Too far apart to be a proper H-bond
  if(DonorAcceptorDistance > MAX_DONOR_ACCEPTOR_EUCLIDEAN_DISTANCE)
    return false; 

  // Check angle
    //peptide plane
  int idx;
  OOPSMPdeclare_vector3(CAofDonorRes);
  OOPSMPdeclare_vector3(CofPrevDonorRes);
  findCA(DonorResIdx,idx,CAofDonorRes);
  if(DonorResIdx > 1)
    findC(DonorResIdx-1,idx,CofPrevDonorRes);
  else
    findC(DonorResIdx,idx,CofPrevDonorRes);
  OOPSMPdeclare_vector3(v1);
  OOPSMPdeclare_vector3(v2);
  OOPSMPdeclare_vector3(v3);
  OOPSMPsub_vector3(pAcceptor,pDonor,v1);
  OOPSMPsub_vector3(CofPrevDonorRes,pDonor,v2);
  OOPSMPsub_vector3(CAofDonorRes,pDonor,v3);

  //additional test from Linus not mentioned in paper
  twoBondAngle = OOPSMPangle_vector3(v1,v3);
  AcceptorWithDonorPeptideAngle = dihedral_vector3(v1,v2,v3);
 

  if(AcceptorWithDonorPeptideAngle > MAX_OUT_OF_PLANE_DIHEDRAL_FOR_HBOND ||
     twoBondAngle < MIN_TWOBOND_ANGLE)
    return false;
  return true;
}


 

}
