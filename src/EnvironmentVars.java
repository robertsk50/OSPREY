/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University
	
	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as 
	published by the Free Software Foundation, either version 3 of 
	the License, or (at your option) any later version.
	
	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.
	
	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.
		
	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.
	
	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129 
			USA
			e-mail:   www.cs.duke.edu/brd/
	
	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
*/
///////////////////////////////////////////////////////////////////////////////////////////////
//	EnvironmentVars.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//     KER        Kyle E. Roberts           Duke University                 ker17@duke.edu
//     MAH        Mark A. Hallen            Duke University                 mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////


import java.math.BigInteger;
import java.util.HashMap;


public class EnvironmentVars {
	
	public enum FORCEFIELD {
		AMBER, CHARMM22, CHARMM19NEUTRAL, CHARMM19
	}
	
	static String dataDir = "./"; //Directory where all of the data files are stored
	static String localDir = "./"; //Local directory where we can quickly write to and delete files
	static FORCEFIELD forcefld; 

	//static BigInteger maxKSconfs;
	//static boolean useMaxKSconfs = false;
	
	// PGC 2014: If set to true, rotamers that are in the voxel of other rotamers will not be added to the rotamer library.  
	//  This is useful if you want to prevent wildtype rotamers that are in the rotamer library 
	//  from being added as repeats.  
	public static boolean CLUSTER_SIMILAR_ROTAMERS;
	
	// PGC 2014: DOPE_ROTAMERS: activate to study the effect of rotamer "doping", which means 
	// adding rotamers that are epsilon from other rotamers just to see how the K* score changes.
	public static boolean DOPE_ROTAMERS; 
	
	// PGC 2014: Store all the atom coordinates for the wildtype rotamer and not just the dihedral. Useful for when there are bond angle and bond lengths deviations.
	public static boolean STORE_FULL_WT_ROT = true;
	
	// PGC 2014: Compute Delta H and Delta S.
	public static boolean COMPUTE_ENTHALPY_AND_ENTROPY = false;
	
	// PGC 2014: we multiply RT by a scale value.
	public static double GLOBAL_CONST_RT_MULT;
	
	// Gas constant times the temperature; this is assigned in the KSParser, when  GLOBAL_CONST_RT_MULT is read.
	public static double constRT;
    // ILP interaction cutoff between residues
    public static double ILP_INTERACTION_CUTOFF = 0.2;
	
	public static boolean useEntropy = false;
	private static HashMap<String, Double> entropy = null;
	public static double entropyScaling = 0;
	
	public static String ksConfDir = "ksConfs";
	
//	public static RotamerLibrary aaRotLib;
	public static String aaRotLibFile;
	
    public static boolean autoFix = true;//Should structures being read in be autofixed?

	// PGC 2015: Use the Dunbrack Rotamer library
	public static boolean USE_DUNBRACK_ROTAMER_LIBRARY = true;
	public static String dunbrackRotamerLibraryLocation = "/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/dataFiles/ALL.bbdep.rotamers.lib";
	// Only accept rotamers tha that have a higher probability than this: 
	public static double DUNBRACK_PROBABILTY_CUTOFF = 0.001;
	// Remove the rotamer that is closest to the wildtype rotamer to avoid "double" counting.
	public static boolean REMOVE_CLOSEST_ROTAMER_WT_ROTAMER = true;
	// If we are removing the closest rotamer to the wildtype rotamer, set a cutoff of Dihedral RMSD
	public static double REMOVE_CLOSEST_ROTAMER_RMS_CUTOFF = 30;
	// Diversify dunbrack rotamers by Chi1
	public static boolean DIVERSIFY_DUNBRACK_ROTAMERS_CHI1 = false;
	// Use dunbrack rotamer probabilities as energies.
	public static boolean USE_DUNBRACK_ROTAMER_PROBABILITIES = false;
    
	//Save Human Readable Emat
	public static boolean SAVE_HUMAN_READABLE_EMAT = false;
    
    public static boolean useMPLP = false;
    public static int MPLP_iterations = 100;
    
        
        
	public static String getDataDir() {
		return dataDir;
	}

	public static void setDataDir(String dd) {
		if(!dd.endsWith("/") || !dd.endsWith("\\")){
			dd = dd.concat("/");
		}
		EnvironmentVars.dataDir = dd;
	}
	
	public static void setLocalDir(String ld) {
		if(!ld.endsWith("/") || !ld.endsWith("\\")){
			ld = ld.concat("/");
		}
		EnvironmentVars.localDir = ld;
	}
	
	
	public static void setForcefld(String frcefld) {
		forcefld = FORCEFIELD.valueOf(frcefld.toUpperCase());
	}
	
	/*public static void setMaxKSconfs(BigInteger confs){
		maxKSconfs = confs;
	}
	
	public static void setUseMaxKSConfs(boolean val){
		useMaxKSconfs = val;
	}*/
	
	public static double getEntropyTerm(String rName) {
		if(entropy == null){
			HashMap<String, Double> entr = new HashMap<String,Double>();
			
			entr = new HashMap<String,Double>();
			entr.put("ALA",0.0*entropyScaling);
			entr.put("CYS",1.14*entropyScaling);
			entr.put("ASP",0.61*entropyScaling);
			entr.put("GLU",1.65*entropyScaling);
			entr.put("PHE",0.58*entropyScaling);
			entr.put("GLY",0.0*entropyScaling);
			entr.put("HIS",0.99*entropyScaling);
			entr.put("HIP",0.99*entropyScaling);
			entr.put("HIE",0.99*entropyScaling);
			entr.put("HID",0.99*entropyScaling);
			entr.put("ILE",0.75*entropyScaling);
			entr.put("IIL",0.75*entropyScaling);
			entr.put("LYS",2.21*entropyScaling);
			entr.put("LEU",0.75*entropyScaling);
			entr.put("MET",1.53*entropyScaling);
			entr.put("ASN",0.81*entropyScaling);
			entr.put("PRO",0.0*entropyScaling);
			entr.put("GLN",2.02*entropyScaling);
			entr.put("ARG",2.13*entropyScaling);
			entr.put("SER",1.19*entropyScaling);
			entr.put("THR",1.12*entropyScaling);
			entr.put("VAL",0.50*entropyScaling);
			entr.put("TRP",0.97*entropyScaling);
			entr.put("TYR",0.99*entropyScaling);
			entr.put("DALA",0.0*entropyScaling);
			entr.put("DCYS",1.14*entropyScaling);
			entr.put("DASP",0.61*entropyScaling);
			entr.put("DGLU",1.65*entropyScaling);
			entr.put("DPHE",0.58*entropyScaling);
			entr.put("DGLY",0.0*entropyScaling);
			entr.put("DHIS",0.99*entropyScaling);
			entr.put("DHIP",0.99*entropyScaling);
			entr.put("DHIE",0.99*entropyScaling);
			entr.put("DHID",0.99*entropyScaling);
			entr.put("DILE",0.75*entropyScaling);
			entr.put("DIIL",0.75*entropyScaling);
			entr.put("DLYS",2.21*entropyScaling);
			entr.put("DLEU",0.75*entropyScaling);
			entr.put("DMET",1.53*entropyScaling);
			entr.put("DASN",0.81*entropyScaling);
			entr.put("DPRO",0.0*entropyScaling);
			entr.put("DGLN",2.02*entropyScaling);
			entr.put("DARG",2.13*entropyScaling);
			entr.put("DSER",1.19*entropyScaling);
			entr.put("DTHR",1.12*entropyScaling);
			entr.put("DVAL",0.50*entropyScaling);
			entr.put("DTRP",0.97*entropyScaling);
			entr.put("DTYR",0.99*entropyScaling);
			
			//double[] entr = {0,0.5,0.75,0.75,0.58,0.99,0.97,1.14,1.53,1.19,1.12,2.21,2.13,0.99,0.99,0.99,0.61,1.65,0.81,2.02,0,0,0.5,0.75,0.75,0.58,0.99,0.97,1.14,1.53,1.19,1.12,2.21,2.13,0.99,0.99,0.99,0.61,1.65,0.81,2.02,0};
			entropy = entr;	
			
			//for(int i=0;i<entropy.length;i++)
			//	entropy[i] *= entropyScaling;
		}
		
		Double result = entropy.get(rName);
		if(result != null){
			return result;
		}else{
			System.out.println("Could not find entropy for "+rName);
			return 0.0;
		}
		
	}
	
	
	public static void setEntropyScale(double es){
		if(es>0)
			useEntropy = true;
		entropyScaling = es;
	}
	
	public static void setAArotLibFile(String rlFile){
		aaRotLibFile = rlFile;//new RotamerLibrary(rl, true);
	}


}
