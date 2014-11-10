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
// KSParser.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 * 
 */

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;
import java.util.Vector;

import mpi.MPI;
import mpi.MPIException;

/**
 * 
 * The main class that sets up and handles the basic OSPREY computation and related functions.
 * 
 * The OSPREY functions include:
 * 		doDEE - perform DEE/A* redesign (this includes MinDEE, BD, and BRDEE);
 * 		genStructDEE - generate structures for a selected set of the top doDEE conformations;
 * 		precomputeBackrubs - precompute a list of allowed backrubs for each flexible residue position (used by BRDEE);
 * 		KSMaster - perform K* redesign;
 * 		doSinglePartFn - generate (bound or unbound) structures for the K* ensemble of a given protein-ligand complex;
 * 		doResEntropy - use SCMF to compute the residue entropy for each (non-Pro) residue in a protein.
 *
 */
public class KSParser
{

	static Metrics metrics = new Metrics();
	
	boolean printSegID = false;

	boolean hElect = true; // should hydrogens be used in electrostatic energy calculations
	boolean hVDW = true; // should hydrogens be used in vdw energy calculations
	boolean hSteric = false; // should hydrogens be used in steric checks

	final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)

	// Config file name
	String cfgName = "KStar.cfg";

	ParamSet rParams = null; //the main KStar parameters

	// For soft vdw potential
	double softvdwMultiplier = 1.0;
	// For electrostatics
	boolean distDepDielect = true;
	double dielectConst = 1.0;

	//Determine if dihedral and solvation energies should be computed
	boolean doDihedE = false;
	boolean doSolvationE = false;
	double solvScale = 1.0;
	double stericThresh = -10000.0f; // allowed overlap between the vdW radii of two atoms, steric clash if larger overlap
	double softStericThresh = -10000.0f; // soft steric overlap threshold

	HBondSettings hbonds;

	public enum DEEMETHOD {
		GOLDSTEIN, SPLITFLAGS1, SPLITFLAGS2, MBPAIRS, FULLPAIRS
	}

	public enum ASTARMETHOD{
		ORIG,PGREORDER,ASGUROBI,MIN,LPGUROBI,ASMPLP,BYSEQ,WCSP,BYSUBROT,ASWCSP,ASGUROBIREORDER,ASWCSPREORDER, BYSEQREORDER
	}

	int numThreads = 1;

	//Tags that are passed to denote MPI status 
	final static int regTag = 1; //regular tag for MPI messages
	final static int updateTag = 2; //used in DACS for updating the best energy found for the different partitions
	final static int doneTag = 3; //Used to tell the energy computation that we're done
	
	
	static int numProc = 1; //number of processors for MPI

	//Constant used when the whole protein complex is being computed
	final static int COMPLEX = -1;

	//Constant used to determine if a duplicate sequence has been found in K* calculations
	static final int DUPFOUND = -100;
	
	static boolean mpiRun = false; //determines if this is an MPI run


	//The printstream where we write our output
	PrintStream outPS = System.out;
	
	Molecule[] mols = new Molecule[5];
	RotamerLibrary[][] rotLibs = new RotamerLibrary[5][];
	//make room to cache molecules and rotLibs for up to 5 configurations of strands
	//though we may not need them all

	/** 
	 * Checks if this is an MPI run and calls the respective functions
	 * 
	 */
	public void checkMPI(String[] args) {

		if ((args.length>0)&&(args[0].equalsIgnoreCase("mpi"))) { //MPI run

			mpiRun = true;
			MPItoThread.initialize(mpiRun, 0);

			String tmp[] = new String[args.length-1]; //remove the mpi argument
			System.arraycopy(args, 1, tmp, 0, tmp.length);
			args = tmp;

			args = parseArgs(args);

			try{ handleDoMPI(args);} catch (Exception e){};
		}
		else { //Threaded run
			mpiRun = false;

			args = parseArgs(args);

			MPItoThread.initialize(mpiRun, numThreads); 
			KSParser.numProc = MPItoThread.numProc;
			//Store all the threads that are available for Kyle's "thread mpi"
			MPItoThread.threadEle.put(Thread.currentThread(), new ThreadElement(0));
			
			//KER: If it isn't an mpiRun start extra threads so that we can simulate
			//an mpiRun
			MPItoThread.startThreads(this,Thread.currentThread());

			outputProgInfo(); //output program information
			setConfigPars(); //set the parameters from the configuration file

			outPS = System.out;
			parse(args); //parse the arguments
		}
	}

	/**
	 * Parse and remove command line flags
	 * @param args command line arguments
	 * @return command line arguments with flags removed
	 */
	private String[] parseArgs(String[] args) {
		while(args.length>0 && args[0].startsWith("-")){

			if (args[0].equalsIgnoreCase("-c")){
				cfgName = args[1];
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
			else if(args[0].equalsIgnoreCase("-t")){
				numThreads = new Integer(args[1]);
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
		}

		return args;

	}

	/**
	 * The main function which handles the OSPREY commands
	 */
	public void parse(String[] args) {

		boolean commandLineScript = false;
		boolean firstCommandLine = false;
		byte bytebuff[];
		String s = new String("");  // line being parsed

		if (args.length > 0) {
			commandLineScript = true;
			firstCommandLine = true;
		}

		bytebuff = new byte[150];
		if (!commandLineScript) {
			System.out.print("> ");
			try {
				System.in.read(bytebuff);
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}
			s = new String(bytebuff);  // create a string from bytebuff
		}
		else if (commandLineScript && !firstCommandLine) {
			// If you were running a command line script and the file is over then quit
			s = new String("quit");
		}
		else if (firstCommandLine) {
			s = new String("");
			for(int i=0;i<args.length;i++)
				s = s.concat(args[i] + " ");
			firstCommandLine = false;
		}			

		s = s.trim();  // remove whitespace from beginning and end of line

		StringTokenizer st = new StringTokenizer(s," ;\t\n\r\f");
		String firstToken = new String("");
		if (st.hasMoreTokens())
			firstToken = st.nextToken();  // snag a copy of the first token

//		if (firstToken.equalsIgnoreCase("doSinglePairE"))
//			doSinglePairE(s,null);
		if (firstToken.equalsIgnoreCase("doResEntropy"))
			handleDoResEntropy(s,null);
		else if (firstToken.equalsIgnoreCase("selectResidues"))
			selectResidues(s);
		else if (firstToken.equalsIgnoreCase("compStericOverlap"))
			handleCompStericOverlap(s);
		else if (firstToken.equalsIgnoreCase("precomputeBackrubs"))
			handlePrecomputeBackrubs(s);

		else if (firstToken.equalsIgnoreCase("doDEE"))
			handleDoDEE(s);
		else if (firstToken.equalsIgnoreCase("doIMINDEE"))
			handlePartitionedDEE(s);
		else if (firstToken.equalsIgnoreCase("doExpandedDEE"))
			handleExpandedIMinDEE(s,true);
		else if (firstToken.equalsIgnoreCase("doExpandedMatrix"))
			handleExpandedIMinDEE(s,false);
		else if (firstToken.equalsIgnoreCase("genStructDEE"))
			handleMinDEEApplyRot(s);
		else if (firstToken.equalsIgnoreCase("generateRandConfs"))
			generateRandConfs(s);
		else if (firstToken.equalsIgnoreCase("fitEparams"))
			fitEparams(s);

		//KER: It is better to do a threaded debug than continually
		//update this function IMO
		/*else if (firstToken.equalsIgnoreCase("doSinglePartFn"))
			handleKSTest(s);*/
		else if (firstToken.equalsIgnoreCase("computeEnergyMol"))
			handleComputeEnergyMol(s);
		else if (firstToken.equalsIgnoreCase("KSMaster"))
			handleKSMaster(s);
		else if (firstToken.equalsIgnoreCase("computeEmats"))
			handleComputeAllPairwiseRotamerEnergies(s);

		else if (firstToken.equalsIgnoreCase("genBackbones"))
			generateBackbones(s);
		else if (firstToken.equalsIgnoreCase("identifyRots"))
			identifyRotamers(s);
		else if (firstToken.equalsIgnoreCase("makeStericShell"))
			makeStericShell(s);
		else if (firstToken.equalsIgnoreCase("fixStruct"))
			fixStruct(s);
		else{
			String output = "The function "+firstToken+" was not recognized\n"
					+ "The available functions are: \n"
					+ "doResEntropy, selectResidues, compStericOverlap, precomputeBackrubs,\n"
					+ "doDEE, genStructDEE, generateRandConfs, fitEparams, computeEnergyMol,\n"
					+ "KSMaster, computeEmats, genBackbones, identifyRots, makeStericShell,\n"
					+ "and fixStruct.\n"
					+ "Exiting...";
			System.out.println(output);	
		}
		//exit from all slave nodes
		cleanUpNodes();


	} // End parse function	

	public static void cleanUpNodes(){
		if (mpiRun || MPItoThread.exe != null){ //exit from all slave nodes
			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = null;
			for (int curProc=1; curProc<numProc; curProc++){
				try {MPItoThread.Send(cObj, 0, 1, ThreadMessage.OBJECT, curProc, regTag);} catch (Exception e){}
			}
			if(MPItoThread.exe!=null)
				MPItoThread.exe.shutdown();
		}
	}


	/**
	 * Displays the program version and citations
	 */
	public void outputProgInfo() {

		System.out.println();
		System.out.println("OSPREY Protein Redesign Software Version 2.1 beta");
		System.out.println("Copyright (C) 2001-2012 Bruce Donald Lab, Duke University");
		System.out.println("");
		System.out.println("This program is free software: you can redistribute it and/or modify");
		System.out.println("it under the terms of the GNU Lesser General Public License as");
		System.out.println("published by the Free Software Foundation, either version 3 of the"); 
		System.out.println("License, or (at your option) any later version.");
		System.out.println("");
		System.out.println("This program is distributed in the hope that it will be useful,");
		System.out.println("but WITHOUT ANY WARRANTY; without even the implied warranty of");
		System.out.println("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the");
		System.out.println("GNU Lesser General Public License for more details.");
		System.out.println("");
		System.out.println("There are additional restrictions imposed on the use and distribution");
		System.out.println("of this open-source code, including: (A) this header must be included");
		System.out.println("in any modification or extension of the code; (B) you are required to");
		System.out.println("cite our papers in any publications that use this code.  The citation");
		System.out.println("for the various different modules of our software, together with a");
		System.out.println("complete list of requirements and restrictions are found in the");
		System.out.println("document license.pdf enclosed with this distribution.");
		System.out.println("");

		if(mpiRun)
			System.out.println("OSPREY running on "+numProc+" processor(s)");
		else
			System.out.println("OSPREY running on "+numThreads+" thread(s)");
		System.out.println();

	}

	//Sets the parameters from the configuration file
	public void setConfigPars() {

		rParams = new ParamSet();
		rParams.addParamsFromFile(cfgName);

		// PGC Added default values for all variables. 
		hElect = (new Boolean((String)rParams.getValue("HELECT", "true"))).booleanValue();
		hVDW = (new Boolean((String)rParams.getValue("HVDW", "true"))).booleanValue();
		hSteric = (new Boolean((String)rParams.getValue("HSTERIC","false"))).booleanValue();
		distDepDielect = (new Boolean((String)rParams.getValue("DISTDEPDIELECT","true"))).booleanValue();
		dielectConst = (new Double((String)rParams.getValue("DIELECTCONST","6.0"))).doubleValue();
		doDihedE = (new Boolean((String)rParams.getValue("DODIHEDE","false"))).booleanValue();
		doSolvationE = (new Boolean((String)rParams.getValue("DOSOLVATIONE","true"))).booleanValue();
		solvScale = (new Double((String)rParams.getValue("SOLVSCALE","0.5"))).doubleValue();
		softvdwMultiplier = (new Double((String)rParams.getValue("VDWMULT","0.95"))).doubleValue();
		stericThresh = (new Double((String)rParams.getValue("STERICTHRESH","0.4"))).doubleValue();
		softStericThresh = (new Double((String)rParams.getValue("SOFTSTERICTHRESH","1.5"))).doubleValue();
		EnvironmentVars.setDataDir(rParams.getValue("DATADIR","./"));
		EnvironmentVars.setForcefld(rParams.getValue("FORCEFIELD","AMBER"));
		double entropyScale = (new Double((String)rParams.getValue("ENTROPYSCALE","0.0"))).doubleValue();
		EnvironmentVars.setEntropyScale(entropyScale);

		EnvironmentVars.setLocalDir(rParams.getValue("LOCALDIR","./"));

		double hbondScale = (new Double((String)rParams.getValue("HBONDSCALE","0"))).doubleValue();
		String dsspFile = rParams.getValue("DSSPFILE","");
		hbonds = new HBondSettings(hbondScale, dsspFile);

		EnvironmentVars.setAArotLibFile(EnvironmentVars.getDataDir().concat(rParams.getValue("ROTFILE","LovellRotamer.dat")));
		/*numAAallowed = rl.getNumAAallowed();
		resAllowed = rl.getAAtypesAllowed();
		rotamerIndexOffset = rl.getRotamerIndexOffset();
		totalNumRotamers = rl.getTotalNumRotamers();*/

		EnvironmentVars.STORE_FULL_WT_ROT = new Boolean((String)rParams.getValue("STORE_FULL_WT_ROT","true")).booleanValue();
		EnvironmentVars.autoFix = new Boolean((String)rParams.getValue("AUTOFIX","true")).booleanValue();

		String ramaGlyFile = (String)rParams.getValue("RAMAGLYFILE","rama500-gly-sym.data");

		if( ! ramaGlyFile.equalsIgnoreCase("none") ){
			String ramaFiles[] = { EnvironmentVars.dataDir + ramaGlyFile,
					EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPROFILE","rama500-pro.data"),
					EnvironmentVars.dataDir + (String)rParams.getValue("RAMAGENFILE","rama500-general.data"),
					EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPREPROFILE","rama500-prepro.data")
			};
			RamachandranChecker.getInstance().readInputFiles( ramaFiles );
		}
	}

	/******************************/
	static // This function returns the number of tokens in string s
	int numTokens(String s) {

		int curNum = 0;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

		while (st.hasMoreTokens()) {
			curNum++;
			st.nextToken();
		}
		return(curNum);
	}


	/******************************/
	static // This function returns the xth token in string s
	String getToken(String s, int x) {

		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
				st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	// Computes the factorial of input n
	public static BigInteger factorial(int n){

		if (n==0)
			return BigInteger.valueOf(1);
		if (n<0){
			System.out.println("ERROR: trying to find factorial of: "+n);
			System.exit(0);
		}

		return (factorial(n-1).multiply(BigInteger.valueOf(n)));
	}

	// This function generates all possible combinations of n choose m
	public void generateCombinations(int residueMutatable[][], int n, int m) {

		int curIndex[] = new int[1];
		int curComb[] = new int[n];
		curIndex[0] = 0;
		generateCombHelper(0,n,curIndex,residueMutatable,curComb,0,m);
	}
	private void generateCombHelper(int depth, int maxDepth, int curIndex[], int
			residueMutatable[][], int curComb[], int numUsed, int maxToUse){

		if (depth >= maxDepth){
			if (numUsed == maxToUse) {
				for (int i=0; i<maxDepth; i++) {
					residueMutatable[curIndex[0]][i] = curComb[i];
				}
				curIndex[0]++;
			}
			return;
		}

		curComb[depth] = 0;
		generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed,maxToUse);

		if (numUsed < maxToUse) {
			curComb[depth] = 1;
			generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed+1,maxToUse);
		}
	}
	// end combination code

	//Sets the allowables for AS residue position curPos
	private void setAllowablesHelper(ParamSet sParams, boolean addWT, Residue r){
		String tempResAllow = (String)sParams.getValue("RESALLOWED"+r.getResNumberString(), "");
		if(numTokens(tempResAllow) <= 0 && !addWT){
			System.out.println("Warning: addWT is false, but no amino acid type set. Using WT for residue "+r.fullName);
			r.setAllowable(r.defaultAA);
		}
		for(int q=0;q<numTokens(tempResAllow);q++)
			r.setAllowable(getToken(tempResAllow,q+1));
		if (addWT)
			r.setAllowable(r.defaultAA); //the default type is set last
	}

	//Sets up the molecule system and returns the number of ligand rotamers
	private Molecule setupMolSystem(Molecule m, ParamSet sParams, boolean[] strandPresent, String[][] strandLimits){
		return setupMolSystem(m,sParams,strandPresent,strandLimits,true);
	}

	//Sets up the molecule system and returns the number of ligand rotamers
	private Molecule setupMolSystem(Molecule m, ParamSet sParams, boolean[] strandPresent, String[][] strandLimits, boolean keepCofactor){

		try{
			FileInputStream is = new FileInputStream((String)sParams.getValue("PDBNAME"));
			new PDBChemModel(m, is);
		}
		catch (Exception e){
			System.out.println("WARNING: An error occurred while reading file");
			System.out.println(e);
			e.printStackTrace();
			System.exit(1);
		}

		int strNum = 0; //the current strand number; 0 is reserved for the protein strand, the ligand strand is 1 (if present)

		//get the number of strands that are present
		int numPresent = 0;
		for(int str=0;str<strandPresent.length;str++)
			if(strandPresent[str])
				numPresent++;

		Molecule newMol = new Molecule();

		//int curPresStr = 0;
		for(int i=0; i<strandLimits.length;i++){
			String pdbStart = strandLimits[i][0]; 
			String pdbEnd   = strandLimits[i][1];
			//if (pdbEnd>=0){ //with ligand in PDB
			int molStartNum = m.mapPDBresNumToMolResNum(pdbStart);
			int molEndNum   = m.mapPDBresNumToMolResNum(pdbEnd);
			if(molStartNum <0 || molEndNum < 0){
				System.out.println("Please make sure strand "+i+"'s begin and end are set properly.");
				System.exit(0);
			}
			int numInStrand = molEndNum-molStartNum+1;
			//strandLength[i] = numInStrand;
			Residue myStrand[] = new Residue[numInStrand];
			for (int j=(numInStrand-1); j>=0; j--){
				myStrand[j] = m.residue[molStartNum+j]; // pull out the ligand
				m.deleteResidue(molStartNum+j);
				myStrand[j].renumberResidue();
			}
			if (strandPresent[i]) { //ligand will be used in design

				newMol.addStrand(""+i);
				strNum = newMol.strand.length-1;
				newMol.strand[strNum].rotTrans = (new Boolean((String)sParams.getValue("STRANDROTTRANS"+i))).booleanValue();

				for(int j=0; j<numInStrand;j++)
					newMol.addResidue(strNum,myStrand[j],false);

				newMol.strand[strNum].isProtein = (new Boolean((String)sParams.getValue("STRANDAA"+i))).booleanValue();
				if (newMol.strand[strNum].isProtein) //use the AA rotamer library for the ligand
					newMol.setAARotLib(EnvironmentVars.aaRotLibFile);
				else //use the non-AA rotamer library for the ligand
					newMol.setGenRotLib(sParams.getValue("GROTFILE","GenericRotamers.dat"));

				//change ligand to the specified residue type
				//if ( m.strand[ligStrNum].isProtein && !m.strand[ligStrNum].residue[0].name.equalsIgnoreCase(ligType) ) //not the same ligand type
				//	(new StrandRotamers(grl,m.strand[ligStrNum])).changeResidueType(m,0,ligType,true);

				strNum++;
				//curPresStr++;

			}
			//}
			/*else if (strandPresent[i]){
			System.out.println("ERROR: Attempting to use a ligand, but ligand not found in system config file");
			System.exit(1);
		}*/
		}


		//Get the cofactors (if present)
		if(keepCofactor){
			int str=-1; 
			for(int i=0; i<strandPresent.length; i++){
				if(strandPresent[i]){
					str++;
					String cofMapString = sParams.getValue("COFMAP"+i, "-1");
					//int numCofactorRes = (new Integer((String)sParams.getValue("NUMCOFRES"))).intValue();

					Residue cof;
					for(int res=0;res<numTokens(cofMapString);res++){
						int cofactorRes = m.mapPDBresNumToMolResNum(getToken(cofMapString, res));
						if(cofactorRes >= 0){
							cof = m.residue[cofactorRes];
							cof.cofactor = true;
							m.deleteResidue(cofactorRes);
							newMol.addResidue(str,cof,false);
						}
					}
				}
			}
		}

		newMol.determineBonds();
		newMol.establishConnectivity(false);

		return newMol;

		//Determine the number of rotamers for the ligand (if used)
		/*int numLigRotamers = 0;
		if (useLig) {
			numLigRotamers = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
			if (numLigRotamers == 0)
				numLigRotamers = 1;
		}
		return numLigRotamers;*/
	}

	// Computes the bound or unbound partition function and
	//  can compute the energyminimized structure for a specified
	//  rotameric conformation
	/*public void handleKSTest(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutataion search parameter filename (string)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters		

		// Pull search parameters
		//int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String eMatrixNameMin = (String)sParams.getValue("MINENERGYMATRIXNAME");
		String eMatrixNameMax = (String)sParams.getValue("MAXENERGYMATRIXNAME");
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		boolean repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT"))).booleanValue();
		double maxIntScale = (new Double((String)sParams.getValue("MAXINTSCALE"))).doubleValue();
		double initEw = (new Double((String)sParams.getValue("INITEW"))).doubleValue();
		double pruningE = (new Double((String)sParams.getValue("PRUNINGE"))).doubleValue();
		double stericE = (new Double((String)sParams.getValue("STERICE"))).doubleValue();
		//boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		//String ligType = (String)sParams.getValue("LIGTYPE");		
		boolean saveConfs = (new Boolean((String)sParams.getValue("OUTPUTPDBS"))).booleanValue();
		boolean KSGMEC = (new Boolean((String)sParams.getValue("KSGMEC"))).booleanValue();
		boolean KSCONFTHRESH = (new Boolean((String)sParams.getValue("KSCONFTHRESH"))).booleanValue();
		String tmpNumKSconfs = (String)sParams.getValue("numKSconfs");
		BigInteger numKSconfs = new BigInteger(tmpNumKSconfs);

		String fName = (String)sParams.getValue("PDBPREFIX");
		String resMut = (String)sParams.getValue("RESMUT");
		//String mutNum = (String)sParams.getValue("mutNum");
		//double epsilon = (new Double((String)sParams.getValue("EPSILON"))).doubleValue();

		if (!doMinimize)
			minimizeBB = false;
		if (!minimizeBB)
			doBackrubs = false;

		if ( (!ligPresent) && ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"))).booleanValue()) ) { //ligPresent, or a different input structure is used for the unbound partition function computation
			sParams.setValue("PDBNAME",sParams.getValue("UNBOUNDPDBNAME"));
			sParams.setValue("PDBLIGNUM","-1");
			eMatrixNameMin = sParams.getValue("MINENERGYMATRIXNAMEUNBOUND");
			eMatrixNameMax = sParams.getValue("MAXENERGYMATRIXNAMEUNBOUND");
		}

		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();

		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);

		// Define the mutation amino-acid sequence
		System.out.print("Mutation Sequence:");
		String curSeq[] = new String[numInAS];
		for(int i=0;i<numInAS;i++){
			curSeq[i] = getToken(resMut,i+1);
			System.out.print(" "+curSeq[i]);
		}
		System.out.println();

		System.out.println("Beginning setAllowables");
		for(int i=0;i<numInAS;i++){
			rs.setAllowable(residueMap[i],curSeq[i]);
		}

		System.out.print("Loading precomputed min energy matrix...");
		loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin+".dat",doMinimize,eMatrixNameMax+".dat");
		System.out.println("done");

		BigDecimal q_L = BigDecimal.ZERO;
		if (ligPresent)
			q_L = getLigPartFn(m,numInAS,ligType,eMatrixNameMin+".dat"); //compute the ligand partition function

		System.out.println("Before start");		

		boolean prunedRotAtRes[] = new boolean[numInAS*totalNumRotamers+numLigRotamers];
		for (int i=0; i<prunedRotAtRes.length; i++)
			prunedRotAtRes[i] = false;

		//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
		prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers, 
				residueMap, rotamerIndexOffset, prunedRotAtRes, stericE);

		if (doMinimize) //precompute the interval terms in the MinDEE criterion
			rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, 
					rotamerIndexOffset, prunedRotAtRes, scaleInt, maxIntScale);

		prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
				rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, false, minimizeBB);

		//Prune with MinBounds
		prunedRotAtRes = rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, false);

		//Compute the Ec value and prunedIsSteric[]
		rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, true);

		BigDecimal initialBest = BigDecimal.ZERO;
		if (ligPresent)
			initialBest =  q_E.multiply(bestScore.multiply(q_L)).multiply(new BigDecimal(gamma * epsilon));

		//Do the rotamer search
		rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
				residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);

		if ((repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done

			System.out.println();
			System.out.println("Repeating search..");
			rs.repeatSearch = false; //reset the flag
			rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
					residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);
		}
	}*/


	// Finds the energy for a given input system (a molecule with specified flexible residues)
	public void handleComputeEnergyMol(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Ligand (boolean), is true if present
		// 3: Amino acid type for ligand (if ligand is absent, write none or anything)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		String runName = ((String)sParams.getValue("RUNNAME"));
		MolParameters mp = new MolParameters();
		loadStrandParams(sParams, mp, COMPLEX);

		//Check to see if the pdb name is a directory
		Object pdbName = (String)sParams.getValue("PDBNAME");
		File f = new File((String) pdbName);
		Object[] pdbFiles = null;
		if(f.isDirectory()){
			pdbFiles = getPdbFiles(f);
		}
		else{
			pdbFiles = new String[1];
			pdbFiles[0] = pdbName;
		}

		System.out.println("Starting energy computation");
		for(int q = 0; q<pdbFiles.length ; q++){
			//Change the pdbfile to be looked at;
			sParams.setValue("PDBNAME", (String) pdbFiles[q]);
			//Setup the molecule system
			Molecule m = new Molecule();
			m = setupMolSystem(m,sParams,mp.strandPresent,mp.strandLimits);
			mp.m = m;
			Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier,hbonds);
			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			//TODO: Fix this so that ligand numbers can get set in the energy function
			/*if (ligPresent)
				a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());
			 */

			/*boolean specificInt = true;
			if(specificInt){

				loadMutationParams(sParams, mp);

				String flag = "SHL-AS";int str1 = 1; int res1 = 0;int str2 = -1; int res2 = -1;
				RotamerSearch rs = new RotamerSearch(m,-1, mp.strandsPresent, hElect, hVDW, hSteric, true,
						true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, 
						doSolvationE, solvScale, softvdwMultiplier, rl, grl, forcefield);
				//rs.strandMut = strandMut;
				double energy = rs.getPairE(flag, str1, res1, str2, res2,mp.strandMut);
				System.out.println(""+pdbFiles[q]+": Energy: "+energy);
			}
			else{*/

			HashMap<String,double[]> eRef = null;
			String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
			boolean useEref = (new Boolean((String)sParams.getValue("USEEREF","true"))).booleanValue();
			if(useEref){
				eRef = Emat.loadErefMatrix(runNameEMatrixMin+"_COM.dat.eref");
				if(eRef == null){
					System.out.println("Eref not calculated...exiting.");
					System.exit(0);
					//					handleComputeOnlyEref(sParams);
					//					eRef = Emat.loadErefMatrix(runNameEMatrixMin+"_COM.dat.eref");
				}
			}

			double totEref = 0.0;
			double totEntropy = 0.0;
			if(useEref || EnvironmentVars.useEntropy){
				loadMutationParams(sParams, mp);
				FullConf conf = new FullConf(mp.strandMut.allMut.length);
				for(int i=0; i<mp.strandMut.allMut.length;i++){
					conf.pdbNums[i] = mp.m.residue[mp.strandMut.allMut[i]].getResNumberString();
					conf.AAnames[i] = mp.m.residue[mp.strandMut.allMut[i]].name;
				}
				if(useEref)
					totEref = getTotSeqEref(eRef,conf,m);
				if (EnvironmentVars.useEntropy)
					totEntropy = getTotSeqEntropy(conf.AAnames,mp.strandMut, mp.m );

			}





			double energy[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
			energy[0] -= (totEref - totEntropy);
			System.out.println("System energy: " + energy[0]+" (elect: "+energy[1]+" vdW: "+energy[2]+" solvation: "+energy[3]+") Eref: "+totEref+" Entropy: "+totEntropy);
			//}

		}
	}

	private double getTotSeqEref(HashMap<String, double[]> eRef, FullConf conf, Molecule m) {
		double totEref = 0;
		for(int i=0; i<conf.AAnames.length;i++){
			double tmpE = eRef.get(conf.pdbNums[i])[m.residue[m.mapPDBresNumToMolResNum(conf.pdbNums[i])].rl.getAAType(conf.AAnames[i]).index];
			totEref += tmpE;

		}
		return totEref;
	}

	private double getTotSeqEntropy(String[] AAnames, MutableResParams strandMut, Molecule m) {

		double totEref = 0;
		for(int i=0; i<AAnames.length;i++){
			int str = strandMut.resStrand[i];
			if(m.strand[str].isProtein){
				double tmpE = m.aaRotLib.getAAType(AAnames[i]).entropy;
				totEref += tmpE;
				//System.out.println("Entropy: "+tmpE);
			}
		}
		return totEref;
	}

	//input should be a directory
	public Object[] getPdbFiles(File f){
		Vector<String> pdbFiles = new Vector<String>();
		File allMyFolderObjects[]  = f.listFiles();
		for(int i =0; i<allMyFolderObjects.length; i++){
			String filename = allMyFolderObjects[i].getName();
			String ext = (filename.lastIndexOf(".")==-1)?"":filename.substring(filename.lastIndexOf(".")+1,filename.length());
			if(ext.equals("pdb"))
				pdbFiles.add(""+f.getPath()+"\\"+filename);
		}

		return pdbFiles.toArray();

	}

	/**
	 * Performs K* redesign; sets up the K* computation from the input model and configuration files and distributes the
	 * candidate mutants for evaluation by the set of available processors.
	 */
	public void handleKSMaster(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		Settings settings = new Settings();
		
		/******** Load all of the settings for DEE *******/
		// Pull search parameters
		String runName = Settings.getRunName(sParams);
		
		//DEE Settings
		Settings.DEE deeSettings = settings.new DEE(sParams);
		double difference = deeSettings.Ival;
		
		//Minimization Settings
		Settings.Minimization minSettings = settings.new Minimization(sParams);
		
		
		//EPICSettings
		EPICSettings es = new EPICSettings(sParams);
		if(deeSettings.Ival+deeSettings.initEw>es.EPICThresh2){
			System.out.println("EPICThresh2 must be at least Ival+Ew: raising to Ival="+(deeSettings.Ival+deeSettings.initEw));
			es.EPICThresh2 = deeSettings.Ival+deeSettings.initEw;
		}
		
		//Enumeration Settings
		Settings.Enum enumSettings = settings.new Enum(sParams);
		
		//Emat Settings
		Settings.Emat ematSettings = settings.new Emat(sParams, runName, minSettings.doPerturbations);
		
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		
		//Output Settings
		Settings.Output outputSettings = settings.new Output(sParams, runName);
		
		//KStar Settings
		Settings.KStar kstarSettings = settings.new KStar(sParams, runName);
		
		//Unclassified Settings
		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
				
		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
		String resumeFilename ="";
		if(resumeSearch){
			resumeFilename = ((String)sParams.getValue("RESUMEFILENAME"));
		}

		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Min Energy Matrix: "+ematSettings.runNameEMatrixMin);
		System.out.println("Volume Center: "+kstarSettings.targetVol);
		System.out.println("Volume Window Size: "+kstarSettings.volWindow);
		System.out.println("Num Residues Allowed to Mutate: "+kstarSettings.numMutations);

		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}

		MolParameters mp = loadMolecule(sParams, COMPLEX, graphSettings.neighborList, graphSettings.distCutoff,true);
		//KER: This is a placeholder so I don't have to change all the variables in the code
		Molecule m = mp.m;
		int numberMutable = mp.strandMut.numMutPos();
		int strandsPresent = mp.strandsPresent;
		String[][] strandLimits = mp.strandLimits;
		boolean[] strandPresent = mp.strandPresent;
		MutableResParams strandMut = mp.strandMut;
		String[][] strandDefault = mp.strandDefault;

		// Create the mutation list with estimated energies
		Set<OneMutation> mutSet = new TreeSet<OneMutation>();


		// Generate all combinations (include (n choose m), (n choose m-1), ... , (n choose 1), and (n choose 0) )
		int numCombAll = 0;
		int numMutations = Math.min(kstarSettings.numMutations, numberMutable);
		for (int i=numMutations; i>=0; i--)
			numCombAll += factorial(numberMutable).divide(factorial(numberMutable-i).multiply(factorial(i))).intValue();
		int residueMutatableAll[][] = new int[numCombAll][numberMutable];
		int curInd = 0;
		for (int i=numMutations; i>=0; i--){
			int numCombCur = factorial(numberMutable).divide(factorial(numberMutable-i).multiply(factorial(i))).intValue();
			int residueMutatableCur[][] = new int[numCombCur][numberMutable];
			generateCombinations(residueMutatableCur,numberMutable,i);
			for (int j=0; j<numCombCur; j++){
				residueMutatableAll[curInd] = residueMutatableCur[j];
				curInd++;
			}
		}

		// At this point each row of residueMutatble is a 0/1 array, 1 indicates
		//  that that residues can mutate

		if(minSettings.selectPerturbations)//Need to run the automatic perturbation selection
			//This only needs to be done once though: after that the perturbations can be read from pertFile
			selectPerturbations(mp, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, ematSettings.addWTRot, sParams);
		//We'll need to do this once for each strand!  They have different perturbations...easiest to keep separate files
		//this is the complex here

		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(mp.strandPresent, sParams);
		for(int resID:mp.strandMut.allMut){
				setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}
		

		System.out.print("Checking if precomputed energy matrix is already computed...");
		RotamerSearch rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
				kstarSettings.epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
				softvdwMultiplier, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, false, false, es,hbonds,mp.strandMut);


		rs.setupRCs( minSettings.doPerturbations);
		

		for(int i=0; i<m.numberOfStrands;i++){

			String strandPertFile = "STR"+i+"."+minSettings.pertFile;

			if ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"+i, "false"))).booleanValue()){ //a different input structure is used for the unbound partition function computation
				ParamSet ubParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
				ubParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
				ubParams.setValue("PDBNAME",sParams.getValue("UNBOUNDPDBNAME"+i));
				ubParams.setValue("NUMOFSTRANDS", "1");
				ubParams.setValue("STRAND0", sParams.getValue("STRAND"+i));
				ubParams.setValue("STRANDAA0", sParams.getValue("STRANDAA"+i));
				ubParams.setValue("STRANDROTRANS0", "FALSE");
				ubParams.setValue("STRANDMUTNUMS", getToken(sParams.getValue("STRANDMUTNUMS"),i+1));
				ubParams.setValue("STRANDMUT0", sParams.getValue("STRANDMUT"+i));
				ubParams.setValue("COFMAP0", sParams.getValue("COFMAP"+i,"-1"));
				ubParams.setValue("UNBOUNDSTRAND", String.valueOf(i) );

				ubParams.setValue("PERTURBATIONFILE", strandPertFile);

				System.out.print("Checking if precomputed energy matrix (unbound) is already computed...");
				rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
						kstarSettings.epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
						softvdwMultiplier, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, false, false, es,hbonds,mp.strandMut);

				rs.setupRCs(minSettings.doPerturbations);

				loadUnboundPairwiseEnergyMatrices(ubParams,rs,ematSettings.runNameEMatrixMin,true,i);

				if(es.useEPIC)
					loadUnboundCETMatrix(ubParams,rs,Double.POSITIVE_INFINITY,false,i);

				rs = null;
				ubParams = null;
				System.out.println("done");
				//BAD CODE //setupMolSystem(m,sParams,ligPresent,ligType); //re-initialize, since some molecule-relative variables have changed (e.g., ligStrNum)
				//VERY BAD CODE
				m = setupMolSystem(m,sParams,strandPresent,strandLimits);
			}
			else{
				
				mp = loadMolecule(sParams, i, graphSettings.neighborList, graphSettings.distCutoff,true); //Load new molecule for curStrForMatrix
				
				//Set the allowable AAs for each AS residue
				if(!addWT)
					mp.strandMut.checkWT(mp.strandPresent, sParams);
				for(int resID:mp.strandMut.allMut){
						setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
				}
				
//				if(rs==null || minSettings.selectPerturbations){//MH: this can happen if we just handled a special unbound-strand structure
					rs = new RotamerSearch(mp.m,mp.strandMut.numMutPos(), mp.strandsPresent,hElect,hVDW,hSteric,true,true,
							kstarSettings.epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
							softvdwMultiplier, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, false, false, es,hbonds, mp.strandMut);

					rs.setupRCs(minSettings.doPerturbations);
//				}

				sParams.setValue("PERTURBATIONFILE", strandPertFile);

				rs.resetMatrices();
				Emat emat = loadPairwiseEnergyMatrices(sParams,ematSettings.runNameEMatrixMin,minSettings.doMinimize,i, es,mp.m, false);


				if(es.useEPIC)
					loadCETMatrix(sParams,rs,i,Double.POSITIVE_INFINITY,false, emat);
			}
		}

		rs = null;
		System.out.println("done");

//		m = new Molecule();

		//revert this for the complex
		sParams.setValue("PERTURBATIONFILE", minSettings.pertFile);
		mp = loadMolecule(sParams, COMPLEX, graphSettings.neighborList, graphSettings.distCutoff,true);
//		m = setupMolSystem(m,sParams,strandPresent,strandLimits);
		rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
				kstarSettings.epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
				softvdwMultiplier, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, false, false, es,hbonds,mp.strandMut);
		//Compute the matrices for all strands (-1 is for the complex)
		//KER: Put this last so the correct Eref matrix is kept
		rs.resetMatrices();
		Emat emat = loadPairwiseEnergyMatrices(sParams,ematSettings.runNameEMatrixMin,minSettings.doMinimize,COMPLEX, es,m, false);


		if(es.useEPIC)
			loadCETMatrix(sParams,rs,COMPLEX,Double.POSITIVE_INFINITY,false, emat);

		//Load mutation list for distribution
		mutSet = hybridKSLoadMutList(kstarSettings.mutFileName, mp, numCombAll, residueMutatableAll,
				sParams, kstarSettings.targetVol, kstarSettings.volWindow);

		OneMutation[] mutArray = mutSet.toArray(new OneMutation[1]);
		BigDecimal bestScore = new BigDecimal("0.0"); //for the resume results		
		// If doing a resume, read the initial results into a bunch of OneMutations
		HashMap<String, OneMutation> resumeResults = new HashMap<String,OneMutation>();
		if (resumeSearch) {

			bestScore = new BigDecimal("0.0"); //high scores are better

			//OneMutation resumeResults[] = new OneMutation[mutArray.length];
			readResumeFile(resumeResults,resumeFilename,numberMutable,false,false,-1, mp.m, mp.strandMut);
			if(resumeResults != null){
				System.out.println("Read "+resumeResults.size()+" completed mutations");

				// Now filter removed mutations (already computed results
				//  are NOT written to file since you already have them)
				// We do need to maintain the best score
				int newIndex = 0;
				LinkedList<OneMutation> newVector2 = new LinkedList<OneMutation>();
				Iterator<OneMutation> iter = mutSet.iterator();
				while(iter.hasNext()){
					//for(int q=0;q<mutArray.length;q++) {
					OneMutation curMut = iter.next();
					String seq = "";
					for(int i=0; i<curMut.resTypes.length;i++ )
						seq += curMut.resTypes[i]+" ";
					//int w = findMutationIndex(resumeResults,curMut.resTypes);
					if (resumeResults.containsKey(seq)){
						OneMutation oneMut = resumeResults.get(seq);
						bestScore = bestScore.max(oneMut.score); //higher scores are better for Hybrid MinDEE-K*
					}
					else{
						newIndex++;
						newVector2.add(curMut);
					}
				}
				//				mutSet = new TreeSet<OneMutation>();
				//				Iterator<OneMutation> q = newVector2.iterator();
				//				while(q.hasNext()){
				//					//for(int q=0; q<newArray2.length; q++){
				//					mutSet.add(q.next());
				//				}
				mutArray = newVector2.toArray(new OneMutation[0]);
				//System.arraycopy(newArray2,0,mutArray,0,newIndex);
				System.out.println("Length of mutArray after removing already computed mutations: "+mutArray.length);
				if(mutArray.length == 0){
					System.out.println("Length of mutArray is 0, so there are no mutations to compute.");
					System.out.println("Quitting.....");
					return;
				}
			}
		}

		//Check to see if any of the strands sequences are repeated so we don't have to compute them more than once
		checkDuplicateMutations(mutArray, mp.m);

		MutationManager mutMan = new MutationManager(runName,mutArray,false);

//		mutMan.setMutationSearch(true);
		mutMan.setMolecule(mp.m);
		mutMan.setDEEsettings(deeSettings.deeSettings);
//		mutMan.setIMinDEE(doIMinDEE);
		mutMan.setSaveTopConfs(kstarSettings.saveTopConfs);
		mutMan.setPrintTopConfs(kstarSettings.printTopConfs);
		mutMan.setNumTopConfs(kstarSettings.numTopConfs);
		mutMan.setStrandMut(strandMut);
		mutMan.setStrandDefault(strandDefault);
		mutMan.setStrandPresent(strandPresent);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setStrandLimits(strandLimits);
		mutMan.setAddOrigRots(ematSettings.addOrigRots);
		mutMan.setNumMutations(numMutations);
		mutMan.setarpFilenameMin(ematSettings.runNameEMatrixMin+".dat");
		mutMan.setDoMinimization(minSettings.doMinimize);
		mutMan.setMinimizeBB(minSettings.minimizeBB);
		mutMan.setDoBackrubs(minSettings.doBackrubs);
		mutMan.setBackrubFile(minSettings.backrubFile);
		mutMan.setRepeatSearch(kstarSettings.repeatSearch);
		mutMan.setInitEw(deeSettings.initEw);
		mutMan.setGamma(kstarSettings.gamma);
		mutMan.setEpsilon(kstarSettings.epsilon);
		mutMan.setStericE(deeSettings.stericE);
		mutMan.setPruningE(deeSettings.pruningE);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setComputeEVEnergy(true);
		mutMan.setCalculateVolumes(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(deeSettings.scaleInt);
		mutMan.setMaxIntScale(deeSettings.maxIntScale);
		mutMan.setRotamerLibrary(m.aaRotLib);
		mutMan.setUseMaxKSconfs(kstarSettings.useMaxKSconfs);
		mutMan.setNumKSconfs(kstarSettings.maxKSconfs);
		//DEEPer
		mutMan.setDoPerturbations(minSettings.doPerturbations);
		mutMan.setPertFile(minSettings.pertFile);
		mutMan.setMinimizePerts(minSettings.minimizePerts);
		mutMan.setAddWTRot(ematSettings.addWTRot);
		mutMan.setIdealizeSC(Perturbation.idealizeSC);
		mutMan.setUseFlagsAStar(false);
		mutMan.setEnumSettings(enumSettings);
		mutMan.setPDBoutDir(outputSettings.pdbOutDir);

		mutMan.setES(es);

		if (resumeSearch)
			mutMan.setBestScore(bestScore);	// Set the current best score from the partial results
		else
			mutMan.setBestScore(new BigDecimal("0.0")); //the initial best score is 0.0

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			e.printStackTrace();
			System.exit(1);
		}

		System.out.println("DONE: K* computation");
	}

	/**
	 * Computes the partition function for the ligand using the rotamers from the (ligand) rotamer library
	 */
	/*private BigDecimal getLigPartFn(Molecule m, int numInAS, String ligType, String eMatrixNameMin){

		double minMatrix[][][][][][] = (double [][][][][][])readObject(eMatrixNameMin);

		int numRot = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
		if (numRot==0) //ALA or GLY
			numRot = 1;

		ExpFunction efunc = new ExpFunction();

		BigDecimal q_L = new BigDecimal("0.0");

		for (int i=0; i<numRot; i++)
			q_L = q_L.add(efunc.exp(-minMatrix[numInAS][grl.getAARotamerIndex(ligType)][i][numInAS][0][0]/constRT));

		System.out.println("Ligand partition function (double): "+q_L.doubleValue());

		return q_L;
	}*/

	//Loads the mutation sequence list for Hybrid MinDEE-K*; computes a list if one cannot be loaded
	//	private Set<OneMutation> handleHybridKSLoadMutList (String mutFileName, int numMutable,
	//			Molecule m, int numComb, int residueMutatable[][], ParamSet sParams,int strandMut[][], String strandDefault[][],
	//			StrandRotamers[] strandRot, double targetVol,double volWindow,int strandsPresent, boolean strandPresent[]){
	//
	//		// Look for previous mutation file
	//		System.out.println();
	//		System.out.print("Looking for mutation list file ");
	//		Set<OneMutation> mutSet = loadMutationList(mutFileName,numMutable,false);
	//
	//		if (mutSet == null) {
	//
	//			// Create the mutation list with estimated energies
	//			mutSet = new TreeSet<OneMutation>();
	//			RotamerSearch rs = new RotamerSearch(m, numMutable, strandsPresent, hElect, hVDW, hSteric, true,
	//					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst,doDihedE,
	//					doSolvationE,solvScale,softvdwMultiplier,grl,
	//					false,null,false,false,false,new EPICSettings(),hbonds);//Not going to consider perturbations, etc. at this level of approximation
	//
	//			//KER: load vol file for each rotamer library
	//			for (int strNum=0; strNum<rs.strandRot.length; strNum++){
	//				//KER: Only load volfile if residue is allowed to mutate
	//				//if(rs.strandRot[strNum].rl.getNumAAallowed() > 1)
	//				rs.strandRot[strNum].rl.loadVolFile(); //load the rotamer volume file
	//
	//			}
	//
	//			int curNumSeq = 0;
	//			boolean valid;
	//			for(int i=0; i<numComb; i++) {
	//				valid = true;
	//				// Reset each amino acid type
	//				System.out.print("Starting mutation combination " + i + " ... ");
	//				for(int str=0;str<strandMut.length;str++)
	//					rs.refreshStrand(str); // clears allowables and does some other stuff
	//
	//				boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
	//				if(!addWT)
	//					checkWT(strandDefault, strandPresent, sParams);
	//				int molStrand = 0;
	//				int mutIndex = 0;
	//				for (int strNum=0; strNum<strandPresent.length; strNum++){
	//					if(strandPresent[strNum]){
	//						for (int k=0; k<strandMut[molStrand].length; k++){ 
	//							if (residueMutatable[i][mutIndex] == 1)
	//								setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, strandMut, strandDefault);
	//							else{
	//								valid = false;
	//								String tempResAllow = (String)sParams.getValue("RESALLOWED"+strNum+"_"+k,"");
	//								if(numTokens(tempResAllow) <= 0 && !addWT){
	//									System.out.println("Error: resAllowed not set for strand");
	//									System.exit(1);
	//								}
	//								for(int q=0;q<numTokens(tempResAllow);q++){
	//									if(getToken(tempResAllow,q+1).equalsIgnoreCase(strandDefault[molStrand][k])){
	//										rs.setAllowable(strandMut[molStrand][k],strandDefault[molStrand][k],molStrand); //the default type is set last
	//										valid = true;
	//									}
	//								}
	//								if(addWT){
	//									rs.setAllowable(strandMut[molStrand][k],strandDefault[molStrand][k],molStrand); //the default type is set last
	//									valid = true;
	//								}								
	//							}
	//							mutIndex++;
	//						}
	//						molStrand++;
	//					}
	//				}
	//
	//				// Perform simple mutation search for this set of mutatable residues
	//				if(valid){
	//					curNumSeq = rs.simpleMasterMutationSearch(strandMut,numMutable,
	//							curNumSeq,mutSet,targetVol-volWindow,
	//							targetVol+volWindow);
	//				}
	//				System.out.println("finished");
	//			}
	//
	//			System.out.println("Sequences remaining after volume filter "+curNumSeq);
	//
	//			// We now have all the mutations in mutArray, collapse the mutArray
	//			//  to the actual number of mutations we have.
	//
	//			//KER: mutArray is now a set so we don't need to reallocate a smaller array
	//			//OneMutation newArray[] = new OneMutation[curNumSeq];
	//			//System.out.println("Allocated newArray");
	//			//System.out.println("Initial Length of mutArray: "+mutArray.length);
	//			//System.arraycopy(mutArray,0,newArray,0,curNumSeq);
	//			//mutArray = newArray;
	//			System.out.println("Trimmed Length of mutArray: "+mutSet.size());
	//
	//			//KER: mutArray is now a set so there should be no duplicates to begin with
	//			/*System.out.print("Removing duplicates...");
	//			mutArray = removeDuplicates(mutArray);
	//			System.out.println("done");*/
	//
	//			System.out.println(mutSet.size()+" unique mutation sequences found in volume range "+(targetVol-volWindow)+" to "+(targetVol+volWindow));
	//			BigInteger numConfs = BigInteger.ZERO;
	//			Iterator<OneMutation> i = mutSet.iterator();
	//			while(i.hasNext()){
	//				OneMutation tmp = i.next();
	//				numConfs = numConfs.add(tmp.numConfUB.add(tmp.numConfB));
	//			}
	//			System.out.println("Total number of conformations (bound and unbound) for all sequences: "+numConfs);// Save mutation list
	//			saveMutationList(mutSet,mutFileName,false);
	//		}
	//
	//
	//		// Sort the mutation list
	//		// System.out.print("Sorting mutation list ... ");
	//		//KER: mutSet is a TreeSet which is sorted so no need to sort
	//		//RyanQuickSort rqs = new RyanQuickSort();
	//		//rqs.Sort(mutArray);
	//		//rqs = null;
	//		// System.out.println("done");
	//
	//		return mutSet;
	//	}

	//Loads the mutation sequence list for Hybrid MinDEE-K*; computes a list if one cannot be loaded
	private Set<OneMutation> hybridKSLoadMutList (String mutFileName, MolParameters mp,
			int numComb, int[][] residueMutatable,ParamSet sParams,double targetVol,double volWindow ){

		// Look for previous mutation file
		System.out.println();
		System.out.print("Looking for mutation list file ");
		Set<OneMutation> mutSet = loadMutationList(mutFileName,mp.strandMut.numMutPos(),mp);


		if (mutSet == null) {

			mp.m.aaRotLib.loadVolFile(); //load the rotamer volume file

			// Create the mutation list with estimated energies
			mutSet = new TreeSet<OneMutation>();
			RotamerSearch rs = new RotamerSearch(mp.m,mp.strandMut.numMutPos(), mp.strandsPresent, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, 
					false, "", false, false, false, null,hbonds, mp.strandMut);


			boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","false"))).booleanValue();
			if(!addWT)
				mp.strandMut.checkWT(mp.strandPresent, sParams);

			int curNumSeq = 0;
			for(int i=0; i<numComb; i++) {
				// Reset each amino acid type
				System.out.print("Starting mutation combination " + i + " ... ");


				//for(int str=0;str<strandMut.length;str++)
				//	rs.refreshStrand(str); // clears allowables and does some other stuff

				for(int j=0; j< mp.strandMut.allMut.length;j++){
					Residue r = mp.m.residue[mp.strandMut.allMut[j]];
					r.clearAllowable(); //clear the allowables before we set them
					if (residueMutatable[i][j] == 1)
						setAllowablesHelper(sParams, addWT, r);
					else
						r.setAllowable(r.defaultAA); //the default type is set last
				}


				// Perform simple mutation search for this set of mutatable residues
				curNumSeq = rs.simpleMasterMutationSearch(mp.strandMut, mp.strandMut.allMut.length,curNumSeq,mutSet, targetVol-volWindow,
						targetVol+volWindow);

				System.out.println("finished");
			}

			System.out.println("Sequences remaining after volume filter "+curNumSeq);

			// We now have all the mutations in mutArray, collapse the mutArray
			//  to the actual number of mutations we have.

			//KER: mutArray is now a set so we don't need to reallocate a smaller array
			//OneMutation newArray[] = new OneMutation[curNumSeq];
			//System.out.println("Allocated newArray");
			//System.out.println("Initial Length of mutArray: "+mutArray.length);
			//System.arraycopy(mutArray,0,newArray,0,curNumSeq);
			//mutArray = newArray;
			System.out.println("Trimmed Length of mutArray: "+mutSet.size());

			//KER: mutArray is now a set so there should be no duplicates to begin with
			/*System.out.print("Removing duplicates...");
				mutArray = removeDuplicates(mutArray);
				System.out.println("done");*/

			System.out.println(mutSet.size()+" unique mutation sequences found"); //in volume range "+(targetVol-volWindow)+" to "+(targetVol+volWindow));
			BigInteger numConfs = BigInteger.ZERO;


			//Set Seq Num
			Iterator<OneMutation> i = mutSet.iterator();
			int ctr=1;
			while(i.hasNext()){
				OneMutation tmp = i.next();
				tmp.mutNum = ctr;
				ctr++;
			}

			//while(i.hasNext()){
			//	OneMutation tmp = i.next();
			//	numConfs = numConfs.add(tmp.numConfUB.add(tmp.numConfB));
			//}
			//System.out.println("Total number of conformations (bound and unbound) for all sequences: "+numConfs);// Save mutation list
			saveMutationList(mutSet,mutFileName,mp);
		}


		return mutSet;
	}



	/**
	 * Reads the results of a partially completed run into an array of CommucObj. The MutationManager then queries 
	 * this array before sending out a task.
	 */
	public HashMap<String,OneMutation> readResumeFile(HashMap<String,OneMutation> resumeResults, String resumeFilename, 
			int numMutable, boolean distrDACS, boolean PEMcomp, int initDepth, Molecule m, MutableResParams strandMut) {

		BufferedReader bufread = null;
		try {
			File file = new File(resumeFilename);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Resume File Not Found");
			return(null);
		}

		boolean done = false;
		String str = null;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).equalsIgnoreCase("Completed")) {
				//				if (PEMcomp) {//PEM computation
				//					resumeResults[resultNum].resMutNum = new int[numMutable];
				//					for(int q=0;q<numMutable;q++) 
				//						resumeResults[resultNum].resMutNum[q] = new Integer(getToken(str,8+q)).intValue();
				//					resumeResults[resultNum].flagMutType = getToken(str,8+numMutable);
				//				}
				//				else { //mutation search
				if (!distrDACS){ //Hybrid-K* or MinDEE/A* resume
					OneMutation oneMut = new OneMutation();
					oneMut.score = new BigDecimal(getToken(str,5));
					//String[] AANames =  new String[numMutable];
					oneMut.resTypes = new int[numMutable];
					String seq = "";
					for(int q=0;q<numMutable;q++) {
						String aaName = getToken(str,20+q);
						oneMut.resTypes[q] = m.residue[strandMut.allMut[q]].rl.getAARotamerIndex(aaName);
						seq += oneMut.resTypes[q]+" ";
					}
					if(resumeResults.containsKey(seq))
						System.out.println("Dupl Seq: "+seq);
					resumeResults.put(seq, oneMut);
				}	
				else {//distributed DACS resume
					System.out.println("DACS NOT IMPLEMENTED");
					OneMutation oneMut = new OneMutation();
					oneMut.mutNum = new Integer(getToken(str,3)).intValue();
					oneMut.score = new BigDecimal(getToken(str,7));
					oneMut.index = new Index3[initDepth];
					//KER: index is not a tuple
					for(int q=0;q<initDepth;q++)
						oneMut.index[q] = new Index3(getToken(str,9+q));
				}
				//}
			}
		}

		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		//OneMutation temp[] = new OneMutation[resultNum];
		//System.arraycopy(resumeResults,0,temp,0,resultNum);
		//resumeResults = temp;
		return (resumeResults);
	}


	//	// Finds the index of the mutation in resumeResults with the same
	//	//  mutation sequence as the targetMutation. If none are found, -1
	//	//  is returned.
	//	public int findMutationIndex(OneMutation resumeResults[],
	//			String targetMutation[]) {
	//
	//		for(int q=0;q<resumeResults.length;q++) {
	//			if (resumeResults[q].isSame(targetMutation))
	//				return(q);
	//		}
	//		return(-1);
	//	}


	// Attempts to read a list of mutations from file
	public Set<OneMutation> loadMutationList(String fName, int numMutable, MolParameters mp) {

		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... no mutation list file found. Computing one.");
			return(null);
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;
		Set<OneMutation> mutList = new TreeSet<OneMutation>();

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else {

				OneMutation tmp = new OneMutation();
				//mutList[resultNum] = new OneMutation();
				tmp.mutNum = new Integer(getToken(str,1));
				tmp.score = BigDecimal.ZERO;
				tmp.vol = new Double(getToken(str,2)).doubleValue();
				tmp.resTypes = new int[numMutable];
				for(int q=0;q<numMutable;q++) {
					Residue r = mp.m.residue[mp.strandMut.allMut[q]];
					String AAtype = getToken(str,3+q);
					tmp.resTypes[q] = r.rl.getAARotamerIndex(AAtype); 	
					if(tmp.resTypes[q] < 0){
						System.out.println("When loading "+fName+" did not recognize "+AAtype+" on line "+resultNum);
					}
				}

				boolean added = mutList.add(tmp);
				if(!added)
					System.out.println("Duplicate Seq: "+str);

				resultNum++;
				/*if (resultNum >= mutList.length){
					OneMutation newArray[] = new OneMutation[mutList.length+1000];
					System.arraycopy(mutList,0,newArray,0,resultNum);
					mutList = newArray;
				}*/
			}
		}

		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		/*OneMutation temp[] = new OneMutation[resultNum];
		System.arraycopy(mutList,0,temp,0,resultNum);
		mutList = temp;*/
		System.out.println(" ... read "+mutList.size()+" mutations from mutation list "+fName);
		return(mutList);
	}	
	// Saves the list of mutations so that a PEM computation/mutation search
	//  doesn't need to recompute these during a resume. Thus
	//  the resume can go more quickly.
	public void saveMutationList(Set<OneMutation> mutList, String fName, MolParameters mp) {

		if (mutList.size() == 0)
			return;

		int numInAS = 0;
		Iterator<OneMutation> i = mutList.iterator();
		numInAS = i.next().resTypes.length;

		PrintStream printStream = setupOutputFile(fName);
		i = mutList.iterator();
		while(i.hasNext()){
			//for(int q=0;q<mutList.length;q++) {
			OneMutation tmp = i.next();

			printStream.print(tmp.mutNum +" "+ tmp.vol);
			for(int w=0;w<numInAS;w++) {
				Residue r = mp.m.residue[mp.strandMut.allMut[w]];
				printStream.print(" "+r.rl.getAAName(tmp.resTypes[w]));
			}
			printStream.println();
		}
		printStream.close();
	}

	//Removes duplicate mutations (for which the mutation sequence is the same) from a given list
	public OneMutation [] removeDuplicates(OneMutation mutArray[]){

		//First, sort the list alphabetically, according to the mutation sequence
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;

		//Copy mutArray into nArray, excluding duplicate entries
		OneMutation nArray[] = new OneMutation[mutArray.length];

		//Copy the first element
		nArray[0] = mutArray[0];
		int nAIndex = 1;

		//Compare each mutation with the previous one in the list
		for (int i=1; i<mutArray.length; i++){ //for each mutation
			if (!(mutArray[i].isSame(mutArray[i-1].resTypes))){ //different sequence
				nArray[nAIndex] = mutArray[i];
				nAIndex++;
			}
		}

		mutArray = new OneMutation[nAIndex];
		System.arraycopy(nArray,0,mutArray,0,mutArray.length);

		return mutArray;//return the reduced list
	}

	// Mutation search Slave function
	public CommucObj handleKSSlave(CommucObj cObj) {
//		EnvironmentVars.aaRotLib = cObj.rl;
		if (cObj.PEMcomp){ //PEM computation
			if (!cObj.entropyComp) //PEM computation
				cObj = handleComputeAllPairwiseRotamerEnergiesSlave(cObj);
			else //entropy E matrix computation
				cObj = handleDoResEntropySlave(cObj);
		}else if(cObj.gurobiCalc || cObj.wcspCalc){
			try {
				gurobiSlave(cObj);
			} catch (MPIException | InterruptedException e) {
				System.out.println("Something went wrong with the slave bounds computation.");
				e.printStackTrace();
			}
		}
		else { //distributed mutation search
			if (cObj.distrDACS){ //running distributed DACS
				cObj = doDistrDACSSlave(cObj);
			}
			else if (cObj.distrDEE){ //running distributed DEE
				cObj = doDistrDEESlave(cObj);
			}
			else { //running Hybrid MinDEE-K*
				cObj = hybridKScompute(cObj);
			}
		}
		return cObj;
	}

	/**
	 * Handles the computation of the K* score for a single mutation sequence with the target ligand.
	 * The 'cObj' parameter contains the mutation search input distributed by the main processor.
	 * Returns the results of the computation to the main processor.
	 */
	private CommucObj hybridKScompute(CommucObj cObj){

		//System.out.println("Start of hybridKScompute");

		Molecule fullMol = cObj.m;

		//KER: We only need to find the partition function for each unbound
		//KER: strand and then the whole complex
		//int[][] allCombos = generateCombos(cObj.strandMut.length);
		int[][] allCombos = new int[fullMol.numberOfStrands+1][];
		for(int com=0; com<allCombos.length;com++){
			allCombos[com] = new int[fullMol.numberOfStrands];
			for(int i=0; i<allCombos[com].length;i++)
				if(com==allCombos.length-1)
					allCombos[com][i] = 1;
				else if(i==com)
					allCombos[com][i] = 1;
				else
					allCombos[com][i] = 0;
		}

		cObj.numComplexes = allCombos.length;
		cObj.repeatEW = new boolean[allCombos.length];
		cObj.allPruned = new boolean[allCombos.length];
		cObj.searchNumConfsTotal = new int[allCombos.length];
		cObj.searchNumConfsPrunedByE = new int[allCombos.length];
		cObj.searchNumConfsPrunedByS = new int[allCombos.length];
		cObj.searchNumConfsEvaluated = new int[allCombos.length];
		cObj.searchNumConfsLeft = new int[allCombos.length];
		cObj.searchNumPrunedMinDEE = new int[allCombos.length];
		cObj.searchBestEnergyFound = new double[allCombos.length];
		cObj.q = new BigDecimal[allCombos.length];
		cObj.bestEMin = new double[allCombos.length];
		cObj.bestE = new double[allCombos.length];
		cObj.q_Time = new int[allCombos.length];
		//cObj.effEpsilon = new double[allCombos.length];

		int[] numMutPerStrand = new int[fullMol.numberOfStrands];
		int[] offsetMutPerStrand = new int[fullMol.numberOfStrands+1];
		for(int i=0; i<offsetMutPerStrand.length;i++){offsetMutPerStrand[i] = 0;}
		for(int i=0; i<numMutPerStrand.length;i++){
			numMutPerStrand[i] = fullMol.numMutableForStrand(i);
			for(int j=i+1; j<offsetMutPerStrand.length-1;j++)
				offsetMutPerStrand[j] += numMutPerStrand[i];
		}

		System.out.print("## CurMut: "+cObj.curMut+" Starting Sequence: ");
		for(int i=0;i<cObj.currentMutation.length;i++)
			System.out.print(" "+cObj.currentMutation[i]);
		System.out.println(" &&");


		if( Boolean.valueOf(cObj.params.getValue("USECLASHFILTER","false")) ) {
			if(doClashFilter(cObj))//sequence pruned by clash filter
				return cObj;
		}
		
		//KER: Run through all of the partition function calculations
		for(int runNum = 0; runNum<fullMol.numberOfStrands+1; runNum++) {
			long startTime = System.currentTimeMillis();

			boolean lastRun = false;
			if(runNum == fullMol.numberOfStrands)
				lastRun = true;

			int curStrForMatrix = runNum;
			if(lastRun)
				curStrForMatrix = -1;

			
			//Skip this run if this is a duplicate sequence that
			//is being calculated somewhere else
			//Code will break if lastRun is skipped, Also there should
			//never be a last run (fully bound run) that is a duplicate
			if(cObj.duplicateMut != null && !lastRun){
				if(cObj.duplicateMut[runNum]>=0){
					if(cObj.computedPartFuns != null && cObj.computedPartFuns[runNum]!=null){
						cObj.setPartitionProperties(runNum, cObj.computedPartFuns[runNum]);
					}
					continue;
				}
			}
			
			String unboundStr = Integer.toString(runNum);
			int strandsPresent = 0;
			boolean strandPresent[] = new boolean[cObj.strandPresent.length];

			boolean notFullComplex = (runNum != fullMol.numberOfStrands) ;
			int numMut = 0;
			for(int str=0;str<strandPresent.length;str++){
				strandPresent[str] = (1==allCombos[runNum][str]);

				if(strandPresent[str]){

					strandsPresent++;
					numMut += numMutPerStrand[str];
				}
			}

			ParamSet params = null;
			String minEmatrixFile = null;
//			String maxEmatrixFile = null;
			String strandPertFile = null;
			//TODO:fix this hack....
			if ( notFullComplex && ((new Boolean((String)cObj.params.getValue("USEUNBOUNDSTRUCT"+unboundStr,"false"))).booleanValue()) ) { //use a different input PDB structure for the unbound case
				params = new ParamSet();
				params.setParamsValues(cObj.params.getParams(), cObj.params.getValues(), cObj.params.getCurNum());
				params.setValue("PDBNAME",params.getValue("UNBOUNDPDBNAME"+unboundStr));
				//params.setValue("PDBLIGNUM","-1");
				minEmatrixFile = cObj.arpFilenameMin;
//				maxEmatrixFile = cObj.arpFilenameMax;
				minEmatrixFile = minEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
//				maxEmatrixFile = maxEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
				strandPertFile = "STR"+unboundStr+"."+cObj.pertFile;
				//minEmatrixFile = params.getValue("MINENERGYMATRIXNAMEUNBOUND"+unboundStr)+".dat";
				//maxEmatrixFile = params.getValue("MAXENERGYMATRIXNAMEUNBOUND"+unboundStr)+".dat";
			}
			else { //a single input PDB structure is used for the bound and unbound computations
				params = cObj.params;
				minEmatrixFile = cObj.arpFilenameMin;
//				maxEmatrixFile = cObj.arpFilenameMax;
				if(notFullComplex){
					minEmatrixFile = minEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
//					maxEmatrixFile = maxEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
					strandPertFile = "STR"+unboundStr+"."+cObj.pertFile;
				}
				else{
					minEmatrixFile = minEmatrixFile.replace(".dat", "_COM.dat");
//					maxEmatrixFile = maxEmatrixFile.replace(".dat", "_COM.dat");
					strandPertFile = cObj.pertFile;
				}
			}

			if(cObj.doPerturbations)
				Perturbation.idealizeSC = cObj.idealizeSC;

			//Setup the molecule system
			MolParameters mp = new MolParameters();
			mp.numOfStrands = cObj.numberOfStrands;	
			mp.strandLimits = cObj.strandLimits;
			mp.strandPresent = strandPresent;
			mp.strandsPresent = strandsPresent;

			//Setup the molecule system	for each of the strands
			if(mols[curStrForMatrix+1] == null){
				
				// For each unbound entity, copy the full molecule and then delete the other strands.
				if(runNum != fullMol.numberOfStrands){
					//KER: Test to not have to load the molecule from scratch
					Molecule m = null;
					try {
						m = (Molecule)MPItoThread.deepCopy(fullMol);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

					int strOffset = 0;
					for(int i=0; i<fullMol.numberOfStrands;i++){
						if(i != runNum){
							m.deleteStrand(strOffset+i);
							strOffset--;
						}
					}
					mp.m = m;

				}else{
					mp.m = fullMol;
				}
				mols[curStrForMatrix+1] = mp.m;
			}
			else{
				mp.m = mols[curStrForMatrix+1];
			}

			//			loadMutationParams(cObj.params, mp);
			//			mp.numberMutable = getNumberMutable(mp.strandMut);

//			System.out.println("NumMutable "+mp.strandMut.numMutPos());

			//KER: Create local variables so I don't have to add mp to everything...
			Molecule m = mp.m;
			
			int numMutPos = 0;
			try{
			for(Residue r: mp.m.residue)
				if(r.isMutable)
					numMutPos++;
			}catch(Exception E){
				E.printStackTrace();
			}
			
			MutableResParams strandMut = new MutableResParams(numMutPos, mp.m.numberOfStrands);

			//If we shorten the molecule we need to shorten the allMut as well
			int ctr = 0;
//			if(!lastRun){
				for(Residue r:mp.m.residue){
					if(r.isMutable){
						//Get the molResNum relative to the short molecule
						strandMut.addRes(ctr, r, m.rotLibForStrand(r.strandNumber), cObj.addOrigRots);
//						allMut[ctr] = m.mapPDBresNumToMolResNum(r.getResNumberString());
						ctr++;
					}
				}
//			}else{
//				for(Residue r:fullMol.residue){
//					if(r.isMutable){
//						//Get the molResNum relative to the short molecule
//						strandMut.addRes(ctr, r, m.rotLibForStrand(r.strandNumber), cObj.addOrigRots);
//						allMut[ctr] = m.mapPDBresNumToMolResNum(r.getResNumberString());
//						ctr++;
//					}
//				}
//			}

			double minEBound = 0.0;
			BigInteger numConfsPrunedMinDEESteric = null;
			double initEw = cObj.initEw;

			
			
			RotamerSearch rs = new RotamerSearch(m,numMutPos,strandsPresent, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, 
					cObj.dielectConst,cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,
					cObj.doPerturbations,strandPertFile, cObj.minimizePerts, false, false, cObj.es,hbonds,strandMut);


			rs.useCCD = (new Boolean((String)cObj.params.getValue("USECCD","true"))).booleanValue();
			CCDMinimizer.EConvTol = (new Double((String)cObj.params.getValue("ECONVTOL","0.01"))).doubleValue();

			boolean doDih = false;
			double lowestBound = Double.POSITIVE_INFINITY;
			Emat emat = new Emat(minEmatrixFile,doDih,m);
			rs.setMinMatrix(emat);

			outPS.print("Starting sequence: ");

			ctr=0;
			for(Residue r: m.residue){
				if(r.isMutable){
					r.clearAllowable();
					//KER: Note: the offsetMutPerStrand[fullMol.numOfStrands] is 0
					try{
						String newAA = r.rl.getAAName(cObj.currentMutation[ctr+offsetMutPerStrand[runNum]]);
						outPS.print(newAA+" ");
						r.setAllowable(newAA);
						//Set all of the allowedRCs for the current residue
						ArrayList<ResidueConformation> allowedResConf;
						if(lastRun)
							allowedResConf = fullMol.strand[r.strandNumber].rcl.getRCsPosType(r.strandResidueNumber, newAA);
						else
							allowedResConf = fullMol.strand[runNum].rcl.getRCsPosType(r.strandResidueNumber, newAA);
						
						for(ResidueConformation rc: allowedResConf)
							r.setAllowable(rc);
						
						r.flexible = true;

						//KER: Mutate to the proper amino acids
						MutUtils.changeResidueType(m, r.moleculeResidueNumber, newAA, true);
						MutUtils.applyRC(m, r, r.getRCsForType(newAA).get(0));//We set the rotamer for the minimizer initialization
					}catch(Exception E){
						E.printStackTrace();
						System.out.println("DELETE ME");
					}
					ctr++;
				}else{
					r.flexible = false;
				}
			}
			outPS.println("");
			outPS.flush();

			//Prune all rotamers that don't belong to the allowed position
			emat.pruneNotAllowed(m);
			emat.removePrunedRotReducedMem(false);
			
			
			//Prune Steric
			int[] prunedStericPerPos = RotamerSearch.DoPruneStericTemplate(emat,cObj.stericE,true,outPS);

			//			if(cObj.doPerturbations)
			//				rs.setupRCs(cObj.addWTRot);

			boolean KSCONFTHRESH = cObj.numKSconfs.compareTo(BigInteger.ZERO) > 0;
			double Ival = cObj.Ival;
			boolean finished = false;
			boolean finalRun = true;
			//DEE Section
			boolean runDEE = true;
			int loopCtr = 1;
			
			int numSplits = 0;
			while(!finished || !finalRun){
				double DEEIval = Ival;
				
				//Run DEE with the current Ival
				int maxLoopNum = Integer.MAX_VALUE;

				if(runDEE)
					runDEE(cObj.useSF, cObj.doMinimization, cObj.minimizeBB, cObj.scaleInt,
							cObj.initEw,cObj.maxIntScale,false,DEEIval, emat,
							true, false,cObj.stericE, cObj.params,0,maxLoopNum,cObj.deeSettings,
							rs.strandRot,strandMut,rs.m,rs.doPerturbations);

				//Compute the Ec value and prunedIsSteric[] (last parameter is true)
				rs.DoMinBounds(cObj.pruningE,initEw, cObj.useSF, true, false);


				//TODO: Fix usingInitialBest and setting initialBest
				boolean usingInitialBest = (!notFullComplex);
				BigDecimal initialBest = (new BigDecimal("0.0"));

				//Inter-mutation pruning
				if (usingInitialBest){
					initialBest = cObj.bestScore.multiply(new BigDecimal(cObj.gamma * cObj.epsilon));
					for(int i=0;i<allCombos.length-1;i++){
						if(cObj.q[i] == null) 
							//This can happen when one of the sequences is a duplicate and the original sequences
							//hasn't been computed yet. Since we don't know q, initialBest = 0
							initialBest = BigDecimal.ZERO;
						else
							initialBest = initialBest.multiply(cObj.q[i]);
					}
				}


				if(cObj.es.useEPIC){
					//Load the CET Matrix and then make the smaller version that matches the emat
					//This should always already be calculated
					loadCETMatrix(cObj.params,rs,curStrForMatrix,0,true, emat);
				
				}
				//strandPresent from loadCETMatrix is the strand number if positive,
				//-1 for the complex
				//which corresponds to curStr-1 here
				//THIS SHOULD ALREADY BE CALCULATED, ONCE FOR ALL MUTANTS...

				SaveConfsParams saveConfsParams = new SaveConfsParams(cObj.numTopConfs, cObj.saveTopConfs, cObj.printTopConfs, false);
				AStarResults asr = rs.slaveDoRotamerSearch(runNum, cObj.computeEVEnergy,cObj.doMinimization,numMutPos,
						strandMut,usingInitialBest,initialBest,cObj,cObj.minimizeBB,cObj.doBackrubs,cObj.backrubFile,
						saveConfsParams, cObj.curMut, cObj.useMaxKSconfs, cObj.numKSconfs,prunedStericPerPos, DEEIval, cObj.enumSettings);

				if(KSCONFTHRESH && rs.numConfsEvaluated.compareTo(cObj.numKSconfs) >= 0){
					finished = true;
				}else if (asr.status == AStarResults.DONE){ //we have reached an appropriate epsilon score
					finished = true;
				}else if(rs.numConfsEvaluated.compareTo(rs.numConfsTotal.subtract(rs.numConfsPrunedByS)) == 0){
					finished = true;
				}
				else{
					//Ival += 0.5;
					if(Ival < 0.5)
						Ival = 0.5;
					else
						Ival *= 2;
					emat.unPrune();
				}


				if(finished){
					rs.printTopConfs(runNum, cObj.seqNum, cObj.saveTopConfs, cObj.printTopConfs,
							cObj.minimizeBB,cObj.doBackrubs,cObj.pdbOutDir);
				}

			}


			long stopTime = System.currentTimeMillis();
			cObj.q_Time[runNum] = Math.round((stopTime - startTime) / 1000.0f);
		} // end for(runNum)

		System.out.print("## CurMut: "+cObj.curMut+" Finished Sequence: ");
		for(int i=0;i<cObj.currentMutation.length;i++)
			System.out.print(" "+cObj.currentMutation[i]);
		System.out.println(" &&");

		cObj.m = null; //We don't need to pass back the molecule
		return cObj;
	}



	boolean doClashFilter(CommucObj cObj){
		//demand to have alternate ligand clash in all confs,
		//and desired ligand to have >0 non-clashing confs
		//if demands are not meant prune the sequences (by returning true)

		System.out.println("Running clash filter...");


		String thisCETMName = cObj.params.getValue("CETMATRIXNAME",cObj.params.getValue("RUNNAME")+"CETM") + "_COM.dat";
		String altCETMName = cObj.params.getValue("ALTCETMNAME");
		//precomputed CETM w/ alternate ligand 
		//must have same mutable residues

		CETMatrix altCETM = (CETMatrix)KSParser.readObject(altCETMName);
		CETMatrix thisCETM = (CETMatrix)KSParser.readObject(thisCETMName);//main ligand


		int curSeq[] = new int[altCETM.numRes];

		for(int res=0; res<altCETM.numRes; res++){

			curSeq[res] = -1;
			for(int index=0; index<altCETM.resAATypes[res].length; index++){
				if( altCETM.resAATypes[res][index].index == (cObj.currentMutation[res]) ){
					curSeq[res] = index;
					break;
				}
			}

			if(curSeq[res]==-1){
				System.err.println("ERROR: ClashFilter can't identify residue "+cObj.currentMutation[res]+" for res# "+res);
				System.exit(1);
			}
		}



		double cutoff = Double.valueOf(cObj.params.getValue("CLASHFILTERCUTOFF","10"));
		//total energy except for shell-shell interactions (so 10 means substantial net unfavorable contacts in active site)

		ClashFilter cfAlt = new ClashFilter(altCETM,curSeq,cutoff);
		if(cfAlt.hasNonClashing()){
			System.out.println("Sequence accommodates alt ligand");
			return true;
		}
		else {
			ClashFilter cfHere = new ClashFilter(thisCETM,curSeq,cutoff);
			if(!cfHere.hasNonClashing()){
				System.out.println("Sequence does not accommodate ligand");
				return true;
			}
		}

		System.out.println("Clash filter passed.");
		return false;
	}



	//Load the pairwise energy matrices; if not computed, compute, and the load
	//KER: strandPresent represents the strand that is present for this matrix (complex is -1)
	private Emat loadPairwiseEnergyMatrices(ParamSet sParams, String minMatrixFile, boolean doMinimize, 
			int strandPresent, EPICSettings es, Molecule m, boolean storeDih){
		//StrandNumbers are appended to PEM matrix, Complex is denoted with "_COM"
		// PGC
		String runName = sParams.getValue("RUNNAME");
		String suffix = "_" + strandPresent+".dat";
		if(strandPresent == COMPLEX )
			suffix = "_COM.dat";
		//rs.loadPairwiseEnergyMatrices(minMatrixFile + suffix);

		String ematFileLoc = minMatrixFile+suffix;
		//Load Emat
		Emat emat = new Emat(ematFileLoc,storeDih,m);

		/*if (doMinimize)
					rs.loadPairwiseEnergyMatrices(maxMatrixFile + suffix);*/

		if ( emat.resByPos==null ) { //at least one of the matrices not computed, so compute

			System.out.println("Precomputed energy matrices not available..");

			long startTime = System.currentTimeMillis();
			ParamSet newParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
			newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
			newParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
			newParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);
			newParams.setValue("ONLYSINGLESTRAND", new Integer(strandPresent).toString());
			newParams.setValue("DODIH", new Boolean(storeDih).toString());
			newParams.setValue("DOMINIMIZE", new Boolean(doMinimize).toString());
			handleComputeAllPairwiseRotamerEnergiesMaster(newParams,false,0,false,es);

			long stopTime = System.currentTimeMillis();
			System.out.println("PEM execution time: "+((stopTime-startTime)/(60.0*1000.0)));
			System.out.println("DONE: Pairwise energy matrix precomputation");

			emat = new Emat(ematFileLoc,storeDih,m);

		}

		return emat;

	}


	private void loadCETMatrix( ParamSet sParams, RotamerSearch rs, int strandPresent, double ival, boolean local, Emat emat ){
		String runName = sParams.getValue("RUNNAME");
		String suffix = "_" + strandPresent+".dat";
		if(strandPresent == COMPLEX )
			suffix = "_COM.dat";

		String CETMatrixName = (String)sParams.getValue("CETMATRIXNAME",runName+"CETM") + suffix;
		rs.loadCETMatrix(CETMatrixName);

		if( rs.cetm != null ){//Check if the continuous energy matrix is present and complete enough
			if(rs.cetm.ivalCutoff >= ival){
				if(rs.cetm.ivalCutoff > ival){
					CETMatrix reducedcetm = new CETMatrix(emat, rs.m);
					reducedcetm.copyAllTerms(rs.cetm);
					reducedcetm.DOFList = rs.cetm.DOFList;
					rs.cetm = reducedcetm;	
				}
				return;
			}
		}

		//If we get here we have to compute the matrix
		System.out.println("Precomputed continuous energy term matrix not available..");

		long startTime = System.currentTimeMillis();
		ParamSet newParams = new ParamSet();
		newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
		newParams.setValue("CETMATRIXNAME", CETMatrixName);

		newParams.setValue("ONLYSINGLESTRAND", new Integer(strandPresent).toString());
		computeCETM(newParams, emat, ival, local, rs.es, rs.cetm);

		long stopTime = System.currentTimeMillis();
		System.out.println("Continuous energy term matrix execution time: "+((stopTime-startTime)/(60.0*1000.0)));
		System.out.println("DONE: Pairwise energy matrix precomputation");

		rs.loadCETMatrix(CETMatrixName);
	}


	private void loadUnboundCETMatrix( ParamSet sParams, RotamerSearch rs, double ival, boolean local, int unboundStrand ){

		String runName = sParams.getValue("RUNNAME");
		String suffix = "_" + unboundStrand + ".dat";

		String CETMatrixName = (String)sParams.getValue("CETMATRIXNAME",runName+"CETM") + suffix;
		rs.loadCETMatrix(CETMatrixName);

		if( rs.cetm != null ){//Check if the continuous energy matrix is present and complete enough
			if(rs.cetm.ivalCutoff >= ival)
				return;
		}

		//If we get here we have to compute the matrix
		System.out.println("Precomputed unbound continuous energy term matrix not available..");

		long startTime = System.currentTimeMillis();
		ParamSet newParams = new ParamSet();
		newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
		newParams.setValue("CETMATRIXNAME", CETMatrixName);

		newParams.setValue("ONLYSINGLESTRAND", "-1");
		handleComputeAllPairwiseRotamerEnergiesMaster( newParams, true, ival, local, rs.es );

		long stopTime = System.currentTimeMillis();
		System.out.println("Continuous energy term matrix execution time: "+((stopTime-startTime)/(60.0*1000.0)));
		System.out.println("DONE: Pairwise energy matrix precomputation");

		rs.loadCETMatrix(CETMatrixName);
	}


	//Load the pairwise energy matrices; if not computed, compute, and the load
	//KER: strandPresent represents the strand that is present for this matrix (complex is -1)
	private void loadUnboundPairwiseEnergyMatrices(ParamSet sParams, RotamerSearch rs, String minMatrixFile, boolean doMinimize,
			int unboundStrand){
		//StrandNumbers are appended to PEM matrix, Complex is denoted with "_COM"
		// PGC
		String runName = sParams.getValue("RUNNAME");

		String suffix = "_" + unboundStrand+".dat";

		Emat emat = new Emat(minMatrixFile + suffix,false,rs.m);

		if ( emat.resByPos == null ) { //at least one of the matrices not computed, so compute

			System.out.println("Precomputed energy matrices not available..");

			long startTime = System.currentTimeMillis();
			ParamSet newParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
			newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
			newParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
			newParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);
			newParams.setValue("ONLYSINGLESTRAND", "-1");
			handleComputeAllPairwiseRotamerEnergiesMaster(newParams,false,0,false,rs.es);

			long stopTime = System.currentTimeMillis();
			System.out.println("PEM execution time: "+((stopTime-startTime)/(60.0*1000.0)));
			System.out.println("DONE: Pairwise energy matrix precomputation");

			emat = new Emat(minMatrixFile + suffix,false,rs.m);
		}
	}


	/////////////////////////////////////////////////////////////////////////
	// MIN and MAX pairwise energy matrices computation
	/////////////////////////////////////////////////////////////////////////

	//KER: wrapper for computing the pairwise energy matrices
	public void handleComputeAllPairwiseRotamerEnergies(String s){

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		// PGC
		String runName = sParams.getValue("RUNNAME");
		int curStrForMatrix = new Integer(sParams.getValue("ONLYSINGLESTRAND","-1"));

		String suffix = "_" + curStrForMatrix + ".dat";
		if(curStrForMatrix == COMPLEX )
			suffix = "_COM.dat";

		sParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
		sParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);

		handleComputeAllPairwiseRotamerEnergiesMaster(sParams, false, 0, false, new EPICSettings());
	}


	/**
	 * This function sets up the computation for all min and max pairwise rotamer interaction energies and stores 
	 * these energies into user-specified precomputed energy matrices. The computation is distributed to the available processors.
	 * For each rotamer, for each allowed amino acid type, for each flexible residue position, the following energies
	 * are computed: intra-rotamer, rotamer-to-template, rotamer-rotamer, and template-only energies. If a ligand is
	 * present, all energies involving the ligand are also computed.
	 */
	public void handleComputeAllPairwiseRotamerEnergiesMaster(ParamSet sParams, boolean compCETM, double ival, 
			boolean local, EPICSettings es ) {
		//If compCETM is true we calculate a continuous energy matrix, using the given pruning information (this can prunedRot=splitFlags=null and ival=arbitrary otherwise)

		int numMutations = 2; //pairwise energies are computed
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB", "false"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
		String backrubFile ="";
		if(doBackrubs){
			backrubFile = (String)sParams.getValue("BACKRUBFILE");
		}
		String runName = (String)sParams.getValue("RUNNAME");
		String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
		String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM");
		String eRefMatrix = (String)(sParams.getValue("EREFMATRIXNAME","Eref"+runName));
		RotamerSearch.MINIMIZATIONSCHEME minScheme = RotamerSearch.MINIMIZATIONSCHEME.valueOf(sParams.getValue("MINIMIZATIONSCHEME","PAIRWISE").toUpperCase()); 

		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
		boolean typeDep = (new Boolean((String)sParams.getValue("TYPEDEP","false"))).booleanValue();
		//KER: Do we keep the dihedrals the rotamers minimize to
		boolean doDih = (new Boolean((String)sParams.getValue("DODIH", "false"))).booleanValue();
		double stericE = (new Double((String)sParams.getValue("STERICE","10000.0"))).doubleValue();

		boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
		String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
		boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
		Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


		boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
		boolean addOrigRots = false, addWTRot = false;
		if( addWTRotsSomehow ){
			if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
				addWTRot = true;
			else//Otherwise just add them to the rotamer library
				addOrigRots = true;
		}

		boolean neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
		double distCutoff = 0;
		if(neighborList){
			distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
		}

		boolean useCCD = (new Boolean((String)sParams.getValue("USECCD","true"))).booleanValue();
		String CETMatrixName = (String)sParams.getValue("CETMATRIXNAME",runName+"CETM");

		CCDMinimizer.EConvTol = (new Double((String)sParams.getValue("ECONVTOL","0.01"))).doubleValue();
		//ContSCObjFunction.gradStep = (new Double((String)sParams.getValue("GRADSTEP","0.0002"))).doubleValue();

		if( ( es.useEPIC && (!useCCD) ) ) {
			System.err.println("ERROR: CCD needed for only-pairwise minimization, which in turn is needed for EPIC");
			System.exit(1);
		}

		if (!doMinimize) //no minimization
			minimizeBB = false;
		if (!minimizeBB) //not backbone minimization
			doBackrubs = false;

		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams, curStrForMatrix, neighborList, distCutoff,true);
		Molecule m = mp.m;
		int numberMutable = mp.strandMut.allMut.length;
		int strandsPresent = mp.strandsPresent;
		String[][] strandLimits = mp.strandLimits;
		boolean[] strandPresent = mp.strandPresent;
		MutableResParams strandMut = mp.strandMut;

		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);


		System.out.println("Computing _All_ Rotamer-Rotamer Energies");

		System.out.println("Starting minimum and maximum bound energy computation");

		RotamerSearch rs = new RotamerSearch(m, numberMutable, strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, 
				doSolvationE, solvScale, softvdwMultiplier,  doPerturbations, pertFile, 
				minimizePerts, false, false, es, hbonds, strandMut);

		int resMut[] = new int[numberMutable];
		for (int i=0; i<resMut.length; i++)
			resMut[i] = 1;

		//System.out.println("Beginning setAllowables");
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(strandPresent, sParams);
		for(int resID:mp.strandMut.allMut){
				setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}


		//select perturbations for strand if appropriate
		//(for K* this needs to be done here, with strand-only molecule)
		if(curStrForMatrix!=COMPLEX){
			boolean selectPerturbations = (new Boolean((String)sParams.getValue("SELECTPERTURBATIONS","false"))).booleanValue();//Should perturbations be automatically selected?
			if(selectPerturbations)//Need to run the automatic perturbation selection
				selectPerturbations(mp, doPerturbations, pertFile, minimizePerts, addWTRot, sParams);
		}

		if(!mp.loadedFromCache)
			rs.setupRCs(doPerturbations);

		//initialize the pairwise energy matrices (full initialization - for all residues in residueMap[], the ligand, and the template)
		//KER: See if the singles mat has already been computed
		Emat singlesEmat = new Emat(minEMatrixName+"_single",doDih,m);

		MutationManager mutMan;
		OneMutation mutArray[];

		if(singlesEmat.singles == null || singlesEmat.singles.E == null){

			//initialize the pairwise energy matrices (full initialization - for all residues in residueMap[], the ligand, and the template)
			//PEMHandler pemH = new PEMHandler();
			//float mutationEnergiesMin[][][][][][] = pemH.initializePairEMatrix(resMut,strandMut,true,true,true);
			//float mutationEnergiesMax[][][][][][] = null;
			singlesEmat = new Emat(m,strandMut,true,doDih);
			//Emat mutationEnergiesMin = new Emat(m,strandMut,false);
			//Emat mutationEnergiesMax = null;
			//if(doMinimize)
			//	mutationEnergiesMax = new Emat(strandMut);//mutationEnergiesMax = pemH.copyMultiDimArray(mutationEnergiesMin);


			//KER: Only do the AS-SHL entries first, then prune the matrix using sterics before moving
			//KER: on to pairs calculation
			mutArray = getMutArraySingleEcomp(numberMutable,minimizeBB);

			mutMan = new MutationManager(null,mutArray,true);
			mutMan.setDoDih(doDih);
			mutMan.setMolecule(m);
			mutMan.setMinScheme(minScheme);
			mutMan.setStrandMut(strandMut);
			mutMan.setStrandPresent(strandPresent);
			mutMan.setStrandLimits(strandLimits);
			mutMan.setStrandsPresent(strandsPresent);
			mutMan.setAddOrigRots(addOrigRots);
			mutMan.setTypeDep(typeDep);
//			mutMan.setMutableSpots(numberMutable);
			mutMan.setarpFilenameMin(minEMatrixName);
			mutMan.setPairEMatrixMin(singlesEmat);
			//mutMan.setErefMatrix(eRef);
			//mutMan.setErefMatrixName(eRefMatrix);
			mutMan.setParams(sParams);
			mutMan.setComputeEVEnergy(true);
			mutMan.setDoMinimization(doMinimize);
			mutMan.setMinimizeBB(minimizeBB);
			mutMan.setDoBackrubs(doBackrubs);
			mutMan.setBackrubFile(backrubFile);
			mutMan.setCalculateVolumes(false);
			mutMan.setStericThresh(stericThresh);
			mutMan.setSoftStericThresh(softStericThresh);
			mutMan.setDistDepDielect(distDepDielect);
			mutMan.setDielectConst(dielectConst);
			mutMan.setDoDihedE(doDihedE);
			mutMan.setDoSolvationE(doSolvationE);
			mutMan.setSolvScale(solvScale);
			mutMan.setVdwMult(softvdwMultiplier);
			mutMan.setRotamerLibrary(m.aaRotLib);
			mutMan.setcurStrForMatrix(curStrForMatrix);
			mutMan.setNeighborList(neighborList);
			mutMan.setDistCutoff(distCutoff);

			mutMan.setDoPerturbations(doPerturbations);
			mutMan.setMinimizePerts(minimizePerts);
			mutMan.setPertFile(pertFile);
			mutMan.setIdealizeSC(Perturbation.idealizeSC);
			mutMan.setAddWTRot(addWTRot);

			mutMan.setUseCCD(useCCD);

			mutMan.setES(es);

			try{
				handleDoMPIMaster(mutMan,mutArray.length);
			}
			catch (Exception e){
				System.out.println("ERROR: "+e);
				e.printStackTrace();
				System.exit(1);
			}

			//KER: prune steric rot and shrink emat accordingly
			rs.setMinMatrix(singlesEmat);

			mutMan.getMinEmatrix().save(minEMatrixName+"_single",m);
		}

		rs.setMinMatrix(singlesEmat);

		RotamerSearch.DoPruneStericTemplate(singlesEmat, stericE, false,System.out);
		singlesEmat.removePrunedRotReducedMem(true);
		Emat pairsEmat = new Emat(singlesEmat,m,doDih);

		HashMap<String,double[]> eRef = new HashMap<String,double[]>();
		for(int i=0; i<strandMut.numMutPos();i++){
			eRef.put(m.residue[strandMut.allMut[i]].getResNumberString(), new double[m.residue[strandMut.allMut[i]].rl.getNumAAallowed()]);
		}

		//KER: now do all the remaining calculations
		mutArray = getMutArrayPairEcomp(numberMutable,minimizeBB,pairsEmat);

		mutMan = new MutationManager(null,mutArray,true);
		mutMan.setMolecule(m);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStrandPresent(strandPresent);
		mutMan.setStrandLimits(strandLimits);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setAddOrigRots(addOrigRots);
		mutMan.setTypeDep(typeDep);
//		mutMan.setMutableSpots(numberMutable);
		mutMan.setarpFilenameMin(minEMatrixName);
		mutMan.setPairEMatrixMin(pairsEmat);
		mutMan.setErefMatrix(eRef);
		mutMan.setErefMatrixName(eRefMatrix);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setRotamerLibrary(m.aaRotLib);
		mutMan.setcurStrForMatrix(curStrForMatrix);

		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setMinimizePerts(minimizePerts);
		mutMan.setPertFile(pertFile);
		mutMan.setIdealizeSC(Perturbation.idealizeSC);
		mutMan.setAddWTRot(addWTRot);

		mutMan.setUseCCD(useCCD);

		mutMan.setES(es);


		if(compCETM){
			mutMan.setCompCETM(true);
			m.DOFs = DegreeOfFreedom.makeDOFArray(rs.strandRot, strandMut,m);


			CETMatrix lsbm = new CETMatrix(pairsEmat,m);
			lsbm.ivalCutoff = ival;
			mutMan.setCETM(lsbm);

//			mutMan.setPrunedRot(prunedRot);
//			mutMan.setSpFlags(splitFlags);
		}

		mutMan.setEConvTol(CCDMinimizer.EConvTol);

		if(local){
			handleMasterLocally(mutMan,mutArray.length);
		}
		else{
			try{
				handleDoMPIMaster(mutMan,mutArray.length);
			}
			catch (Exception e){
				System.out.println("ERROR: "+e);
				e.printStackTrace();
				System.exit(1);
			}
		}

		if(compCETM){
			outputObject(mutMan.cetm,CETMatrixName);

			if(es.PBTest || es.quantumTest){
				//just testing CET matrix computation...analyze and leave
				mutMan.cetm.analyzeFitTypes();
				System.exit(0);
			}
		}
		else{
			mutMan.getMinEmatrix().save(minEMatrixName,m);
		}

		System.out.println("DONE: Pairwise energy matrix precomputation..");
	}
	
	/**
	 * This function sets up the computation for the CETM and stores 
	 * these energies into a user-specified precomputed matrix. The computation is distributed to the available processors.
	 * @param cetm 
	 * 
	 */
	public void computeCETM(ParamSet sParams, Emat emat, double ival, boolean local, EPICSettings es, CETMatrix cetm ) {
		//compCETM is true we calculate a continuous energy matrix, using the given pruning information (this can prunedRot=splitFlags=null and ival=arbitrary otherwise)

		int numMutations = 2; //pairwise energies are computed
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB", "false"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
		String backrubFile ="";
		if(doBackrubs){
			backrubFile = (String)sParams.getValue("BACKRUBFILE");
		}
		String runName = (String)sParams.getValue("RUNNAME");
		String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
		RotamerSearch.MINIMIZATIONSCHEME minScheme = RotamerSearch.MINIMIZATIONSCHEME.valueOf(sParams.getValue("MINIMIZATIONSCHEME","PAIRWISE").toUpperCase());

		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
		boolean typeDep = (new Boolean((String)sParams.getValue("TYPEDEP","false"))).booleanValue();
		//KER: Do we keep the dihedrals the rotamers minimize to
		boolean doDih = (new Boolean((String)sParams.getValue("DODIH", "false"))).booleanValue();
		double stericE = (new Double((String)sParams.getValue("STERICE","10000.0"))).doubleValue();

		boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
		String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
		boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
		Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


		boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
		boolean addOrigRots = false, addWTRot = false;
		if( addWTRotsSomehow ){
			if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
				addWTRot = true;
			else//Otherwise just add them to the rotamer library
				addOrigRots = true;
		}

		boolean neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
		double distCutoff = 0;
		if(neighborList){
			distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
		}

		boolean useCCD = (new Boolean((String)sParams.getValue("USECCD","true"))).booleanValue();
		String CETMatrixName = (String)sParams.getValue("CETMATRIXNAME",runName+"CETM");

		CCDMinimizer.EConvTol = (new Double((String)sParams.getValue("ECONVTOL","0.01"))).doubleValue();
		//ContSCObjFunction.gradStep = (new Double((String)sParams.getValue("GRADSTEP","0.0002"))).doubleValue();

		if( ( es.useEPIC && (!useCCD) )) {
			System.err.println("ERROR: CCD needed for  EPIC");
			System.exit(1);
		}

		if (!doMinimize) //no minimization
			minimizeBB = false;
		if (!minimizeBB) //not backbone minimization
			doBackrubs = false;

		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams, curStrForMatrix, neighborList, distCutoff,true);
		Molecule m = mp.m;
		int numberMutable = mp.strandMut.allMut.length;
		int strandsPresent = mp.strandsPresent;
		String[][] strandLimits = mp.strandLimits;
		boolean[] strandPresent = mp.strandPresent;
		MutableResParams strandMut = mp.strandMut;

		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);


		System.out.println("Computing _All_ Rotamer-Rotamer Energies");

		System.out.println("Starting minimum and maximum bound energy computation");

		RotamerSearch rs = new RotamerSearch(m, numberMutable, strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, 
				doSolvationE, solvScale, softvdwMultiplier,  doPerturbations, pertFile, 
				minimizePerts, false, false, es, hbonds, strandMut);

		int resMut[] = new int[numberMutable];
		for (int i=0; i<resMut.length; i++)
			resMut[i] = 1;

		//System.out.println("Beginning setAllowables");
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(strandPresent, sParams);
		for(int resID:mp.strandMut.allMut){
			setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}


		//select perturbations for strand if appropriate
		//(for K* this needs to be done here, with strand-only molecule)
		if(curStrForMatrix!=COMPLEX){
			boolean selectPerturbations = (new Boolean((String)sParams.getValue("SELECTPERTURBATIONS","false"))).booleanValue();//Should perturbations be automatically selected?
			if(selectPerturbations)//Need to run the automatic perturbation selection
				selectPerturbations(mp, doPerturbations, pertFile, minimizePerts, addWTRot, sParams);
		}

		if(!mp.loadedFromCache)
			rs.setupRCs(doPerturbations);

		//KER: now do all the remaining calculations
		OneMutation[] mutArray1 = getMutArraySingleEcomp(numberMutable,minimizeBB);
		OneMutation[] mutArray2 = getMutArrayPairEcomp(numberMutable,minimizeBB,emat);
		
		OneMutation[] mutArray = new OneMutation[mutArray1.length+mutArray2.length];
		System.arraycopy(mutArray1, 0, mutArray, 0, mutArray1.length);
		System.arraycopy(mutArray2, 0, mutArray, mutArray1.length, mutArray2.length);

		MutationManager mutMan = new MutationManager(null,mutArray,true);
		mutMan.setMolecule(m);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStrandPresent(strandPresent);
		mutMan.setStrandLimits(strandLimits);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setAddOrigRots(addOrigRots);
		mutMan.setTypeDep(typeDep);
//		mutMan.setMutableSpots(numberMutable);
		mutMan.setarpFilenameMin(minEMatrixName);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setRotamerLibrary(m.aaRotLib);
		mutMan.setcurStrForMatrix(curStrForMatrix);

		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setMinimizePerts(minimizePerts);
		mutMan.setPertFile(pertFile);
		mutMan.setIdealizeSC(Perturbation.idealizeSC);
		mutMan.setAddWTRot(addWTRot);

		mutMan.setUseCCD(useCCD);

		mutMan.setES(es);

		mutMan.setCompCETM(true);
		m.DOFs = DegreeOfFreedom.makeDOFArray(rs.strandRot, strandMut,m);


		
		CETMatrix lsbm = new CETMatrix(emat,m); 
		if(cetm != null){
			//We've already computed some CET, so only compute the new ones that we need and just copy over the rest
			lsbm.copyOverTerms(cetm);
			//Setup the mutArray run params to only compute the terms we need.
			lsbm.updateMutArray(mutArray);
		}
		lsbm.ivalCutoff = ival;
		mutMan.setCETM(lsbm);

		mutMan.setEConvTol(CCDMinimizer.EConvTol);

		if(local){
			handleMasterLocally(mutMan,mutArray.length);
		}
		else{
			try{
				handleDoMPIMaster(mutMan,mutArray.length);
			}
			catch (Exception e){
				System.out.println("ERROR: "+e);
				e.printStackTrace();
				System.exit(1);
			}
		}

		
		outputObject(mutMan.cetm,CETMatrixName);

		if(es.PBTest || es.quantumTest){
			//just testing CET matrix computation...analyze and leave
			mutMan.cetm.analyzeFitTypes();
			System.exit(0);
		}
		
		System.out.println("DONE: Continuous energy matrix precomputation..");
	}

	//Generates and saves to file the mutation list for the pairwise energy matrix computation
	private OneMutation[] getMutArrayPairEcomp(int numberMutable, boolean minimizeBB,Emat emat){

		int numMutations = 2; //pairwise energy computation

		//KER: sometimes there will only be one thing flexible (for example if strand only has one residue flexible)
		if(numberMutable < numMutations)
			numMutations = numberMutable;

		// Generate all combinations
		int numComb = factorial(numberMutable).divide(factorial(numberMutable-numMutations).multiply(factorial(numMutations))).intValue();
		int residueMutatable[][] = new int[numComb][numberMutable];
		generateCombinations(residueMutatable,numberMutable,numMutations);
		// At this point each row of residueMutatble is a 0/1 array which specifies a mutation 
		//  pair combination, 1 indicates that that residue can mutate in the specified combination

		System.out.println("Number of possible mutation combinations: "+numComb);

		// Create the mutation list with estimated energies
		ArrayList<OneMutation> mutList = new ArrayList<OneMutation>();			

		//Set the AS-AS mutations
		int curMutNum = 0;
		//KER: if we only have one mutable position there can't be pairwise interactions
		if(numMutations >= 2){
			for(int i=0; i<numComb; i++) {
				OneMutation curMut = new OneMutation();

				curMut.flagMutType = "AS-AS";
				curMut.resMut = new int[numberMutable];
				curMut.runParams = new EmatCalcParams(-1,-1);
				boolean firstMut = true;
				for(int j=0; j<numberMutable; j++) {
					curMut.resMut[j] = residueMutatable[i][j];
					if(firstMut && curMut.resMut[j] == 1){
						firstMut = false;
						curMut.mutNum = j;
						curMut.runParams.pos1 = j;
					}
					else if(curMut.resMut[j] == 1){
						curMut.runParams.pos2 = j;
					}
				}

				if(emat.areNeighbors(curMut.runParams.pos1, curMut.runParams.pos2)){
					mutList.add(curMut);
					curMutNum++;
				}

				// Perform simple mutation search for this set of mutatable residues

			}
		}

		//Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
		int numOtherMut;
		/*if (ligPresent)
				numOtherMut = 2+t+2*numInAS;
			else
				numOtherMut = 1+t+numInAS;*/
		numOtherMut = 2;//+numberMutable;
		OneMutation otherMutArray[] = new OneMutation[numOtherMut];

		for (int i=0; i<numOtherMut; i++){
			otherMutArray[i] = new OneMutation();
			otherMutArray[i].resMut = new int[numberMutable];
		}

		//Set the AS-shell mutations
		/*for (int i=0; i<numberMutable; i++){
				otherMutArray[i].flagMutType = "SHL-AS";
				for (int j=0; j<numberMutable; j++){
					if (i==j)
						otherMutArray[i].resMut[j] = 1;
					else
						otherMutArray[i].resMut[j] = 0;
				}
			}*/

		//Set the intra-residue energies run
		otherMutArray[0].flagMutType = "INTRA";
		for (int j=0; j<numberMutable; j++)
			otherMutArray[0].resMut[j] = 1;

		//Set the template energy run
		otherMutArray[1].flagMutType = "TEMPL";
		for (int j=0; j<numberMutable; j++)
			otherMutArray[1].resMut[j] = 0;


		/*if (ligPresent){//if the ligand is present, set the corresponding runs

				//Set the AS-ligand mutations
				for (int i=1+t+numInAS; i<=2*numInAS+t; i++){
					otherMutArray[i].flagMutType = "LIG-AS";
					for (int j=0; j<numInAS; j++){
						if ((i-1-t-numInAS)==j)
							otherMutArray[i].resMut[j] = 1;
						else
							otherMutArray[i].resMut[j] = 0;
					}
				}

				//Set the ligand-shell run
				otherMutArray[1+t+2*numInAS].flagMutType = "LIG-SHL";
				for (int j=0; j<numInAS; j++)
					otherMutArray[1+t+2*numInAS].resMut[j] = 0;
			}*/

		// We now have all the mutations in mutArray, collapse the mutArray
		//  to the actual number of mutations we have.
		OneMutation mutArray[] = mutList.toArray(new OneMutation[0]);
		OneMutation newArray[] = new OneMutation[curMutNum+numOtherMut];
		System.arraycopy(otherMutArray,0,newArray,0,numOtherMut);//add the other mutations first
		System.arraycopy(mutArray,0,newArray,numOtherMut,curMutNum);//then add the AS-AS mutations

		mutArray = newArray;
		System.out.println("Length of mutArray: "+mutArray.length);

		return mutArray;
	}



	//Generates and saves to file the mutation list for the pairwise energy matrix computation
	private OneMutation[] getMutArraySingleEcomp(int numberMutable, boolean minimizeBB){

		//Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
		int numOtherMut;
		numOtherMut = numberMutable;
		OneMutation otherMutArray[] = new OneMutation[numOtherMut];

		for (int i=0; i<numOtherMut; i++){
			otherMutArray[i] = new OneMutation();
			otherMutArray[i].resMut = new int[numberMutable];
		}

		//Set the AS-shell mutations
		for (int i=0; i<numberMutable; i++){
			otherMutArray[i].flagMutType = "SHL-AS";
			for (int j=0; j<numberMutable; j++){
				if (i==j){
					otherMutArray[i].resMut[j] = 1;
					otherMutArray[i].mutNum = j;
					otherMutArray[i].runParams = new EmatCalcParams(j, -1);
				}
				else
					otherMutArray[i].resMut[j] = 0;
			}
		}

		return otherMutArray;
	}

	/**
	 * Computes a specific part of the pairwise energy matrices, as specified by the parameters in the 'cObj' parameter,
	 * distributed by the main processor. Returns the results of the computation to the main processor.
	 */
	public CommucObj handleComputeAllPairwiseRotamerEnergiesSlave(CommucObj cObj) {

		long startTime = System.currentTimeMillis();

		//CurStrForMatrix starts at -1 so push everything up one;
		int curComplex = cObj.curStrForMatrix+1;
		int unboundStrand = Integer.valueOf( cObj.params.getValue("UNBOUNDSTRAND","-1") );//MH: If using a special unbound structure,
		//curStrForMatrix is -1 (as needed since we need to use the whole unbound structure),
		//so we need to handle the strand number for molecule reloading and .dat file numbering separately
		if(unboundStrand!=-1)
			curComplex = unboundStrand+1;

		MolParameters mp = loadMolecule(cObj.params, cObj.curStrForMatrix, cObj.m);

		RotamerSearch rs = new RotamerSearch(mp.m,cObj.strandMut.allMut.length, cObj.strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst, 
				cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,
				cObj.doPerturbations,cObj.pertFile,cObj.minimizePerts,false,false,cObj.es,hbonds,cObj.strandMut);

		rs.useCCD = cObj.useCCD;

		//for use in K* CETM calculations
		rs.templateSt = (new Double((String)cObj.params.getValue("STERICE","30.0"))).doubleValue();
		rs.pairSt = (new Double((String)cObj.params.getValue("PAIRST", "100.0"))).doubleValue();


		CCDMinimizer.EConvTol = cObj.EConvTol;
		//ContSCObjFunction.gradStep = cObj.gradStep;

		boolean shellRun = false;
		boolean intraRun = false;
		boolean templateOnly = false;		

		if (cObj.flagMutType.compareTo("TEMPL")==0){

			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = true;/*ligPresent = true*/;intraRun = false;templateOnly = true;
		}
		else if (cObj.flagMutType.compareTo("AS-AS")==0){

			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = false;/*ligPresent = true*/;intraRun = false;templateOnly = false;
		}	
		else if (cObj.flagMutType.compareTo("SHL-AS")==0){

			// Then shell runs for the active site residues
			// Computes the active site residue rotamers to shell energies					
			shellRun = true;/*ligPresent = true*/;intraRun = false;templateOnly = false;				
		}
		else if (cObj.flagMutType.compareTo("INTRA")==0){

			// Compute all intra-residue energies					
			shellRun = false;intraRun = true;templateOnly = false;
		}				
		else if (cObj.flagMutType.compareTo("LIG-AS")==0){

			// **** Ligand present runs ****
			// This section computes the inter-residue energies between
			//  active site residues and the ligand
			shellRun = false;intraRun = false;templateOnly = false;			
		}
		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)

			// Computes ligand rotamer to shell energies
			shellRun = true; intraRun = false;templateOnly = false;
		}

		// The goal is that the total energy of a system can be bounded by the sum of 
		//  all pairwise active site residue entries plus the entry for each active site
		//  residue's shell run plus each active site residue's self intra-residue energy.
		//  If a ligand is present then one should add the ligand to shell energy, the
		//  ligand to each active site residue pairwise energy, and the ligand self intra-
		//  residue energy.

//		if(cObj.compCETM){
//			//Need to provide rs with pruning information so only unpruned RCs or pairs get energies computed
//			rs.eliminatedRotAtRes = cObj.prunedRot;
//			rs.splitFlags = cObj.splitFlags;
//		}
		
		Emat minEmatrix = cObj.emat;
		rs.setMinMatrix(minEmatrix);


		//Compute the corresponding matrix entries
		rs.simplePairwiseMutationAllRotamerSearch(cObj.strandMut,cObj.strandMut.allMut.length,cObj.doMinimization,shellRun,intraRun,
				cObj.resMut,cObj.minimizeBB,cObj.doBackrubs,
				templateOnly,cObj.backrubFile, cObj.minScheme, cObj.runParams, cObj.compCETM);


		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);


		if(cObj.compCETM)
			cObj.cetm = rs.cetm;
		else{
			cObj.compEE = new ArrayList<EMatrixEntrySlim>();

			if(templateOnly){
				cObj.compEE.add(new EMatrixEntrySlim(minEmatrix.getTemplMinE(),null));
			}
			else if(intraRun){
				Iterator<EMatrixEntryWIndex> iter = minEmatrix.singlesIterator();
				while(iter.hasNext()){
					EMatrixEntryWIndex eme = iter.next();
					cObj.compEE.add(new EMatrixEntrySlim(eme.eme, eme.index));
				}
			}
			else{

				if(shellRun){
					Iterator<EMatrixEntryWIndex> iter = minEmatrix.singlesIterator(cObj.runParams.pos1);
					while(iter.hasNext()){
						EMatrixEntryWIndex eme = iter.next();
						if(cObj.runParams.rotamers==null || cObj.runParams.rotamers.contains(new Index3(eme.index[0],eme.index[1],eme.index[2])))
							cObj.compEE.add(new EMatrixEntrySlim((RotamerEntry)eme.eme, eme.index));
					}
				}
				else{
					Iterator<EMatrixEntryWIndex> iter = null;
					if(cObj.runParams.AAs1 == null || cObj.runParams.AAs2 == null)
						iter = minEmatrix.pairsIterator(cObj.runParams.pos1,cObj.runParams.pos2);
					else
						iter = minEmatrix.pairsIterator(cObj.runParams.pos1,cObj.runParams.pos2,cObj.runParams.AAs1,cObj.runParams.AAs2);
					while(iter.hasNext()){
						EMatrixEntryWIndex eme = iter.next();
						if(cObj.runParams.rotamers==null || cObj.runParams.rotamers.contains(new Index3(eme.index[0],eme.index[1],eme.index[2])))
							cObj.compEE.add(new EMatrixEntrySlim((RotamerPairEntry)eme.eme, eme.index));
					}
				}
			}
		}

		//KER: reset evalAtoms
		mp.m.setAllEval(); 

		cObj.emat = null;
		cObj.m = null;
		outPS.print("Finished Mutation\n");
		
		return cObj;
	}

	public MolParameters loadMolecule(ParamSet sParams, int curStrForMatrix, Molecule m){	

		MolParameters mp = new MolParameters();

		loadStrandParams(sParams, mp, curStrForMatrix);

		//Setup the molecule system
		if(m != null){
			mols[curStrForMatrix+1] = m;
			mp.m = mols[curStrForMatrix+1];
		}
		else{
			mp.m = mols[curStrForMatrix+1];
		}
		
		mp.loadedFromCache = true;
		return mp;

	}

	//Generates and saves to file the mutation list for the pairwise energy matrix computation
	private OneMutation[] getMutArrayPairEcomp(int numberMutable, boolean minimizeBB){

		int numMutations = 2; //pairwise energy computation

		//KER: sometimes there will only be one thing flexible (for example if strand only has one residue flexible)
		if(numberMutable < numMutations)
			numMutations = numberMutable;

		// Generate all combinations
		int numComb = factorial(numberMutable).divide(factorial(numberMutable-numMutations).multiply(factorial(numMutations))).intValue();
		int residueMutatable[][] = new int[numComb][numberMutable];
		generateCombinations(residueMutatable,numberMutable,numMutations);
		// At this point each row of residueMutatble is a 0/1 array which specifies a mutation 
		//  pair combination, 1 indicates that that residue can mutate in the specified combination

		System.out.println("Number of possible mutation combinations: "+numComb);

		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[numComb];			

		//Set the AS-AS mutations
		int curMutNum = 0;
		//KER: if we only have one mutable position there can't be pairwise interactions
		if(numMutations >= 2){
			for(int i=0; i<numComb; i++) {

				mutArray[i] = new OneMutation();
				mutArray[i].flagMutType = "AS-AS";
				mutArray[i].resMut = new int[numberMutable];
				for(int j=0; j<numberMutable; j++) {
					mutArray[i].resMut[j] = residueMutatable[i][j];
				}

				// Perform simple mutation search for this set of mutatable residues
				curMutNum++;
			}
		}

		//Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
		int numOtherMut;
		/*if (ligPresent)
			numOtherMut = 2+t+2*numInAS;
		else
			numOtherMut = 1+t+numInAS;*/
		numOtherMut = 2+numberMutable;
		OneMutation otherMutArray[] = new OneMutation[numOtherMut];

		for (int i=0; i<numOtherMut; i++){
			otherMutArray[i] = new OneMutation();
			otherMutArray[i].resMut = new int[numberMutable];
		}

		//Set the AS-shell mutations
		for (int i=0; i<numberMutable; i++){
			otherMutArray[i].flagMutType = "SHL-AS";
			for (int j=0; j<numberMutable; j++){
				if (i==j)
					otherMutArray[i].resMut[j] = 1;
				else
					otherMutArray[i].resMut[j] = 0;
			}
		}

		//Set the intra-residue energies run
		otherMutArray[numberMutable].flagMutType = "INTRA";
		for (int j=0; j<numberMutable; j++)
			otherMutArray[numberMutable].resMut[j] = 1;

		//Set the template energy run
		otherMutArray[numberMutable+1].flagMutType = "TEMPL";
		for (int j=0; j<numberMutable; j++)
			otherMutArray[numberMutable+1].resMut[j] = 0;


		/*if (ligPresent){//if the ligand is present, set the corresponding runs

			//Set the AS-ligand mutations
			for (int i=1+t+numInAS; i<=2*numInAS+t; i++){
				otherMutArray[i].flagMutType = "LIG-AS";
				for (int j=0; j<numInAS; j++){
					if ((i-1-t-numInAS)==j)
						otherMutArray[i].resMut[j] = 1;
					else
						otherMutArray[i].resMut[j] = 0;
				}
			}

			//Set the ligand-shell run
			otherMutArray[1+t+2*numInAS].flagMutType = "LIG-SHL";
			for (int j=0; j<numInAS; j++)
				otherMutArray[1+t+2*numInAS].resMut[j] = 0;
		}*/

		// We now have all the mutations in mutArray, collapse the mutArray
		//  to the actual number of mutations we have.
		OneMutation newArray[] = new OneMutation[curMutNum+numOtherMut];
		System.arraycopy(otherMutArray,0,newArray,0,numOtherMut);//add the other mutations first
		System.arraycopy(mutArray,0,newArray,numOtherMut,curMutNum);//then add the AS-AS mutations

		mutArray = newArray;
		System.out.println("Length of mutArray: "+mutArray.length);

		return mutArray;
	}

	// Finds the index of the mutation in resumeResults with the same
	//  mutation sequence as resMut. If none are found, -1 is returned.
	public int sampFindMutationIndex(OneMutation resumeResults[], String flMutType, int mutResidues[]) {

		for(int q=0;q<resumeResults.length;q++) 
			if ((resumeResults[q].flagMutType.compareTo(flMutType)==0) && (sameSeq(resumeResults[q].resMut,mutResidues)))
				return(q);

		return(-1);
	}

	//Determines if the residues that can mutate are the same for two mutation sequences
	private boolean sameSeq (int computedResMut[], int allResMut[]){

		boolean found = true;
		for (int i=0; i<computedResMut.length; i++){
			if (computedResMut[i]!=allResMut[i])
				found = false;
		}
		return found;
	}
	///////////////////////////////////////////////////////////////////////////
	//	End of MIN and MAX Pairwise Energy Precomputation
	///////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////
	//	 Compute minimized-GMEC section
	////////////////////////////////////////////////////////////////	
	/** 
	 * Computes the (energy-minimized) structure for the rotameric conformations specified by the input file.
	 * This function is used to generate structures for a set of output conformations from a DEE/A* search.
	 */
	public void handleMinDEEApplyRot(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		// Pull search parameters
		String runName = (String)sParams.getValue("RUNNAME");
		String confResFile = (String)sParams.getValue("CONFRESFILE", "c_"+runName);
		int numResults = (new Integer((String)sParams.getValue("NUMRESULTS"))).intValue();
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB","false"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
		String backrubFile ="";
		if(doBackrubs){
			backrubFile = (String)sParams.getValue("BACKRUBFILE");
		}
		boolean outputPDB = (new Boolean((String)sParams.getValue("OUTPUTPDBS", "true"))).booleanValue();

		boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
		String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
		boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
		Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();
		boolean useCCD = (new Boolean((String)sParams.getValue("USECCD","true"))).booleanValue();

		boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
		boolean addOrigRots = false, addWTRot = false;
		if( addWTRotsSomehow ){
			if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
				addWTRot = true;
			else//Otherwise just add them to the rotamer library
				addOrigRots = true;
		}

		boolean neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
		double distCutoff = 0;
		if(neighborList){
			distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
		}


		MolParameters mp = loadMolecule(sParams, COMPLEX, neighborList, distCutoff,true);
		Molecule m = mp.m;
		int numberMutable = mp.strandMut.numMutPos();
		MutableResParams strandMut = mp.strandMut;
		int numOfStrands = mp.m.numberOfStrands;

		StrandRotamers[] strandRot = new StrandRotamers[mp.m.numberOfStrands];
		for(int i=0; i<mp.m.numberOfStrands;i++){
			if(doPerturbations)
				strandRot[i] = new StrandRCs(m.rotLibForStrand(i),m.strand[i]);
			else
				strandRot[i] = new StrandRotamers(m.rotLibForStrand(i),m.strand[i]);
		}
		
		//Load rotamer and residue conformation libraries
		String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
		
		//TODO: Right now I just check to see if an expanded matrix exists and then use it if it does
		//It might be better to have the user denote the expanded matrix somehow instead of just blindly using it
		File expandedRCfile = new File(runNameEMatrixMin+"_expanded_COM.dat.rcl_0");
		if(expandedRCfile.exists())
			runNameEMatrixMin = runNameEMatrixMin+"_expanded";
		
		m.aaRotLib.loadGlobalRots(runNameEMatrixMin+"_COM.dat.aaRots");
		if(m.genRotLib != null)
			m.genRotLib.loadGlobalRots(runNameEMatrixMin+"_COM.dat.genRots");
		if(doPerturbations)
			PertFileHandler.readPertFile(pertFile, m, strandRot,true);
		for(Strand strand : m.strand){
			if(strand.isProtein)
				strand.rcl.loadGlobalRCs(runNameEMatrixMin+"_COM.dat.rcl_"+strand.number, m.aaRotLib);
			else
				strand.rcl.loadGlobalRCs(runNameEMatrixMin+"_COM.dat.rcl_"+strand.number, m.genRotLib);
		}
		
		int numResidues = numberMutable;

		if(outputPDB){
			//KER: make pdbs directory if it doesn't exist
			File pdbDir = new File("pdbs");
			if(!pdbDir.exists())
				pdbDir.mkdir();
		}
				
		//Read the results file into the AA and rot matrices
		ArrayList<FullConf> confs = new ArrayList<FullConf>();

		readRotResFile(confResFile,confs,numberMutable);
		
		int numSaved = 0;
		for (int curResult=0; curResult<numResults; curResult++){


			if( doPerturbations && curResult>0 ){//Reload the molecule
//				mp.m.origMol();
				mp = loadMolecule(sParams, COMPLEX, neighborList, distCutoff,false);
//				m = mp.m;
			}


			System.out.print("Starting minimization of result "+(curResult+1)+"..");

			Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier,hbonds);

			

			//the starting index for the current result within AAtypes[] and rotNums[]
			int startInd = numResidues*curResult;
			a96ff.calculateTypesWithTemplates();

			ArrayList<Integer> mutRes = new ArrayList<Integer>();
			
			for(int i=0; i<confs.get(curResult).pdbNums.length;i++){
				int resID = mp.m.mapPDBresNumToMolResNum(confs.get(curResult).pdbNums[i]);
				Residue r = mp.m.residue[resID];
				mutRes.add(resID);
				mp.m.residue[resID].flexible = true;
				mp.m.residue[resID].curRC = null; 
				String AAname = confs.get(curResult).AAnames[i];
				
				r.setAllowable(AAname);
				if(r.canMutate)
					MutUtils.changeResidueType(mp.m, resID, AAname, true);
				ResidueConformation rc = mp.m.strand[r.strandNumber].rcl.getRC(confs.get(curResult).rcNums[i]);
				
				MutUtils.applyRC(mp.m, r, rc);
			}
			System.out.println("");


//			if(doPerturbations){//This has to be done after setting the allowables
//				PertFileHandler.readPertFile(pertFile, m, strandRot);
//				for(int str=0;str<numOfStrands;str++){
//					((StrandRCs)strandRot[str]).addUnperturbedRCs(addWTRot,m);
//					((StrandRCs)strandRot[str]).countRCs();
//				}
//			}

			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			//TODO: Fix the energy function so we can set each strand
			//a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());

//			int curAA[] = new int[m.numberOfResidues];
//			for(int str=0;str<numOfStrands;str++){
//				for(int j=0;j<m.strand[str].numberOfResidues;j++){
//					int molResNum = m.strand[str].residue[j].moleculeResidueNumber;
//					curAA[molResNum] = strandRot[str].getIndexOfNthAllowable(j,0);
//				}
//			}

			SimpleMinimizer simpMin = null;
			BBMinimizer bbMin = null;
			BackrubMinimizer brMin = null;
			CCDMinimizer ccdMin = null;

			EnergyFunction ef = new ForceFieldEnergy(m, a96ff);
			ContSCObjFunction of = null;


			if ( doMinimize && (!minimizeBB) ){

				if(useCCD){
					m.DOFs = DegreeOfFreedom.makeDOFArray(strandRot, strandMut, m);
					of = new ContSCObjFunction(m,numOfStrands,ef,strandRot,doDihedE,null);
					ef = of.efunc;//Dihedral energies now included if needed
					ccdMin = new CCDMinimizer(of,false);
				}
				else{

					if(doPerturbations)
						simpMin = new PMinimizer(minimizePerts);
					else
						simpMin = new SimpleMinimizer();
					simpMin.initialize(m,numOfStrands,a96ff,strandRot,doDihedE);
				}
			}
			else if (minimizeBB) {
				if (!doBackrubs){
					bbMin = new BBMinimizer();
					/*if (ligPresent)
						bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
					else*/
					bbMin.initialize(m, a96ff, strandMut, numOfStrands);
				}
				else {
					brMin = new BackrubMinimizer();
					/*if (ligPresent)
						brMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum, backrubFile, hSteric, stericThresh);
					else*/
					brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, stericThresh,numOfStrands, true);
				}
			}

			//In DEEPer we will need to reload the molecule for each result instead of backing up
			//(the backup restoration would mess up the perturbations)
			if(!doPerturbations)
				m.backupAtomCoord();	

			double unMinE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy before minimization
			if (outputPDB){ //save molecule
				String fName = runName;
				String filename = String.format("pdbs/%1$s_%2$03d_unmin.pdb",fName,numSaved);
				m.saveMolecule(filename,(double)unMinE[0]);
			}
			if ( doMinimize && (!minimizeBB) ){
				long minimizeStartTime = System.currentTimeMillis();
				if(useCCD){
					of.updateIdealDihedrals();//Dihedrals are now in the ideal rotameric states
					of.setDOFs( ccdMin.minimize() );
				}
				else
					simpMin.minimize(35);
				System.out.println( "Minimization time (ms): "+(System.currentTimeMillis()-minimizeStartTime) );
			}
			else if (minimizeBB) {
				if (!doBackrubs)
					bbMin.minimizeFull(false);
				else{
					brMin.minimizeFull();
				}
			}

			/*double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy after minimization
			if ((doMinimize)&&(!minimizeBB)&&(doDihedE)) //add dihedral energies
				minE[0] += simpMin.computeDihedEnergy();
			 */
			double minE[] = { ef.getEnergy() };

			if (outputPDB){ //save molecule
				String fName = runName;
				String filename = String.format("pdbs/%1$s_%2$03d_min.pdb",fName,numSaved);
				m.saveMolecule(filename,(double)minE[0]);
				numSaved++;
			}

			if(!doPerturbations){
				m.restoreAtomCoord();
				m.updateCoordinates();
			}

			if (minE[0]>unMinE[0])
				minE = unMinE;

			System.out.println("done");
		}
		System.out.println("done");
	}

	private void readRotResFile (String resFile, ArrayList<FullConf> confs, int numResidues){

		BufferedReader bufread = null;
		try {
			File file = new File(resFile);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 


			boolean done = false;
			String str = null;
			int resultNum = 0;

			while (!done) {
				try {
					str = bufread.readLine();
				}
				catch ( Exception e ){
					System.out.println("ERROR: An error occurred while reading input");
					System.exit(0);
				}

				if (str == null) // stop if we've reached EOF
					done = true;
				else {
					FullConf conf = new FullConf(numResidues);
					int offset = 4;
					for (int i=0; i<numResidues; i++){
						conf.pdbNums[i] = getToken(str,2+i*offset);
						conf.AAnames[i] = getToken(str,2+i*offset+1);
						conf.rcNums[i]  = new Integer((String)getToken(str,2+i*offset+2)).intValue();
						conf.rotNums[i] = new Integer((String)getToken(str,2+i*offset+3)).intValue();
					}
					conf.minE = new Double((String)getToken(str,numResidues*offset+5)).doubleValue();
					resultNum++;
					confs.add(conf);
				}
			}

			// We're done reading them in

			bufread.close();

		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Results File Not Found");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		//Collections.sort(confs);
	}

	private void readRotResFile (String resFile, String AAtypes[], int rotNums[], int numResults, int numResidues){

		BufferedReader bufread = null;
		try {
			File file = new File(resFile);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Results File Not Found");
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if(resultNum >= numResults)
				done = true;
			else {
				for (int i=0; i<numResidues; i++){
					AAtypes[resultNum*numResidues+i] = getToken(str,2+i);
					rotNums[resultNum*numResidues+i] = new Integer((String)getToken(str,2+numResidues+i)).intValue();
				}
				resultNum++;
			}
		}

		if (numResults!=resultNum){
			System.out.println("Error: Not all results available for reading");
			System.exit(0);
		}

		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
	}



	//Identify the current rotamers in the flexible residues
	//This function is sort of the inverse of handleMinDEEApplyRot
	public void identifyRotamers(String s){

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters


		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams, COMPLEX, false,0.0,true);
		int numMutable = mp.strandMut.numMutPos();

		//Data to output
		String AAtypes[] = new String[numMutable];
		int curRot[] = new int[numMutable];
		double maxDev[] = new double[numMutable];//Maximum deviation of a dihedral from the rigid-rotamer value


		//Find the rotamer at each position with the minimax deviation from the real structure
		for(int j=0; j<numMutable; j++){

			maxDev[j] = Double.POSITIVE_INFINITY;

			//Get information on the residue
			int str = mp.strandMut.resStrand[j];
			int strResNum = mp.strandMut.resStrandNum[j];
			Residue res = mp.m.residue[mp.strandMut.allMut[j]];
			AAtypes[j] = res.name;

			//Choose the correct rotamer library for this residue
			RotamerLibrary curRL = mp.m.rotLibForStrand(str);

			AARotamerType AAType = curRL.getAAType(AAtypes[j]);
			int AATypeInd = AAType.index;
			int numDihedrals = AAType.numDihedrals();
			int symmCount = ResSymmetry.getSymmetry( AAtypes[j] );//Get the number of symmetry states of the residue

			double curDih[][] = new double[numDihedrals][symmCount];//Current values for the dihedrals
			//in each of the symmetry states

			for(int a=0; a<numDihedrals; a++){

				for(int symm=0; symm<symmCount; symm++){

					int molAtNum[] = new int[4];//Molecule atoms numbers for atoms involved in the dihedral

					for( int b=0; b<4; b++ ){
						String atName = AAType.dihedralAtomNames[a][b];
						atName = ResSymmetry.getPermutedAtomName(atName, AAtypes[j], symm);
						molAtNum[b] = res.getAtomNameToMolnum(atName);
					}

					curDih[a][symm] = mp.m.getTorsion( molAtNum[0], molAtNum[1], molAtNum[2], molAtNum[3] );
				}
			}



			int numRot = curRL.getNumRotamers( AAtypes[j] );

			for(int rot=0; rot<numRot; rot++){//Check all the possible rotamers

				for(int symm=0; symm<symmCount; symm++){//Check each symmetry state against them

					double rotMaxDev = Double.NEGATIVE_INFINITY;
					//Maximum dihedral deviation between the structure (in the current symmetry)
					//and the current rotamer

					for(int k=0; k<numDihedrals; k++){

						double idealDih = AAType.rotamers.get(rot).values[k];

						double dev = (double)Math.abs( idealDih - curDih[k][symm] );//the angles' absolute values can't exceed 180 so this is at most 360
						dev = Math.min(dev, 360-dev);//dev is the absolute difference between the angles, which are mod 360

						if( dev > rotMaxDev )
							rotMaxDev = dev;
					}

					if( rotMaxDev < maxDev[j] ){//This rotamer is the best so far
						maxDev[j] = rotMaxDev;
						curRot[j] = rot;
					}
				}
			}

			if(numRot == 0){
				curRot[j] = 0;
				maxDev[j] = 0;
			}

		}


		//Now write out the answer
		//First all the amino-acid types
		System.out.print( "AA types: " );
		for(int j=0; j<numMutable; j++)
			System.out.print( AAtypes[j] + " " );
		System.out.println();

		System.out.print( "Rotamers: " );
		for(int j=0; j<numMutable; j++)
			System.out.print( curRot[j] + " " );
		System.out.println();

		System.out.print( "Max dihedral deviations: " );
		for(int j=0; j<numMutable; j++)
			System.out.print( maxDev[j] + " " );
		System.out.println();
	}







	//Make a reduced PDB file containing just the designed residues
	//specified with regular configuration files
	//plus any residues within shellThickness
	//(extra parameters required are shellThickness and shellPDBFile)
	public void makeStericShell(String s){

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters


		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams, COMPLEX, false, 0.0, true);
		int numMutable = mp.strandMut.numMutPos();
		String shellPDBFile = (String)sParams.getValue("SHELLPDBFILE");
		double shellThickness = (new Double((String)sParams.getValue("SHELLTHICKNESS"))).doubleValue();

		Molecule newMolec = new Molecule();//only designed residues + those within shellThickness of them
		newMolec.addStrand();

		HashSet<Integer> des = new HashSet<Integer>();//molecule residue numbers of designed residues

		//Add the designed residues to the shell
		for(int j=0; j<numMutable; j++){

			//Get information on the residue
			int str = mp.strandMut.resStrand[j];
			int strResNum = mp.strandMut.resStrandNum[j];
			Residue res = mp.m.residue[mp.strandMut.allMut[j]];
			des.add(res.moleculeResidueNumber);
		}

		ArrayList<Integer> newMolecRes = new ArrayList<Integer>();
		//molecule residue numbers of residues in the new molecule

		for(Residue res : mp.m.residue){
			int molResNum = res.moleculeResidueNumber;
			if(des.contains(molResNum))
				newMolecRes.add(molResNum);
			else{
				for(int molResNum2 : des){//get distances to shell residues
					if(res.getDist(mp.m.residue[molResNum2], true) <= shellThickness){
						newMolecRes.add(molResNum);
						break;
					}
				}
			}
		}

		for(int molResNum : newMolecRes)
			newMolec.addResidue(0, mp.m.residue[molResNum]);

		newMolec.saveMolecule(shellPDBFile,0);
	}

	////////////////////////////////////////////////////////////////
	//	 End of Compute minimized-GMEC section
	////////////////////////////////////////////////////////////////




	///////////////////////////////////////////////////////////////////////////
	//	DEE section
	///////////////////////////////////////////////////////////////////////////
	/**
	 * Performs a DEE (Traditional, MinDEE, BD, BRDEE, or DEEPer) pruning with A* search, with or without DACS;
	 * the only parameter 's' (String) includes the command-line arguments specifying the filenames of the two input configuration files.
	 * If distributed DACS is performed, the computation is distributed to the available processors for evaluation.
	 */
	public void handleDoDEE(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		handleDoDEE(sParams);
	}
	
	public void handleDoDEE(ParamSet sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)

		metrics.setStartTime();

		System.out.println("Performing DEE");

		Settings settings = new Settings();
		
		/******** Load all of the settings for DEE *******/
		// Pull search parameters
		String runName = Settings.getRunName(sParams);
		
		//DEE Settings
		Settings.DEE deeSettings = settings.new DEE(sParams);
		double difference = deeSettings.Ival;
		
		//Minimization Settings
		Settings.Minimization minSettings = settings.new Minimization(sParams);
		
		
		//EPICSettings
		EPICSettings es = new EPICSettings(sParams);
		if(deeSettings.Ival+deeSettings.initEw>es.EPICThresh2){
			System.out.println("EPICThresh2 must be at least Ival+Ew: raising to Ival="+(deeSettings.Ival+deeSettings.initEw));
			es.EPICThresh2 = deeSettings.Ival+deeSettings.initEw;
		}
		
		//Enumeration Settings
		Settings.Enum enumSettings = settings.new Enum(sParams);
		
		//Emat Settings
		Settings.Emat ematSettings = settings.new Emat(sParams, runName, minSettings.doPerturbations);
		
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		
		//Output Settings
		Settings.Output outputSettings = settings.new Output(sParams, runName);
				
		//Unclassified Settings
		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
				
		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
		String resumeFilename ="";
		if(resumeSearch){
			resumeFilename = ((String)sParams.getValue("RESUMEFILENAME"));
		}
		
		/***************** Done Loading Settings ************/
				
		//Check that settings are valid
		if( !deeSettings.useFlags && (enumSettings.useFlagsAStar || deeSettings.useTriples || deeSettings.algOption >=4) ){//These all rely heavily on the split flags
			System.err.println("ERROR: Options requiring split flags (flags in A*, triples pruning, or algOption>=4) are set but split flags are turned off");
			System.exit(1);
		}
		
		if( minSettings.doMinimize && (!deeSettings.useMinDEEPruningEw) && ( deeSettings.useTriples || deeSettings.algOption >=4 ) ){
			System.err.println("ERROR: Options requiring iMinDEE (triples pruning and/or algOption >=4) are set but iMinDEE is turned off");
			System.exit(1);
		}

		if ((!mpiRun)&&((deeSettings.distrDACS)||deeSettings.distrDEE)){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}

		if (!minSettings.doMinimize) //no minimization
			minSettings.minimizeBB = false;
		if (!minSettings.minimizeBB) //not backbone minimization
			minSettings.doBackrubs = false;

		if (graphSettings.genInteractionGraph) //DACS is not performed when generating the interaction graph
			deeSettings.doDACS = false;

		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams,curStrForMatrix, graphSettings.neighborList, graphSettings.distCutoff,true);

		
		boolean reload=false;

		if(minSettings.selectPerturbations){//Need to run the automatic perturbation selection
			//This only needs to be done once though: after that the perturbations can be read from pertFile
			selectPerturbations(mp, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, ematSettings.addWTRot, sParams);
			reload = true;
		}
		
		


		// 2010: If useMinDEEPruningEw is set to false, this cycle never repeats itself.
		//  If it is set to true, it can repeat at most once: if none of the rotamer vectors
		//  between the conformation of lowest energy (i.e. lowestBound) 
		//  and lowestBound+InitEw can minimize to a lower energy than lowestBound+InitEw, 
		//   then let minimumEnergy be the minimum nergy found among the enumerated conformations,
		//   we set a new Ew equal to minimumEnergy - lowestBount and repeat this cycle.  We
		//   only have to do it at most twice.   
		double interval = 0;
		double actualLowestBound = Double.POSITIVE_INFINITY;
		double bestScore = Double.POSITIVE_INFINITY;
		do{

			if( reload ){
//				reloadMolecule(mp, sParams, curStrForMatrix, neighborList, distCutoff, pertFile);
				mp = loadMolecule(sParams,curStrForMatrix,graphSettings.neighborList,graphSettings.distCutoff,false);
//				mp.m.origMol();
			}
			else if(minSettings.doPerturbations)
				reload = true;


			RotamerSearch rs = new RotamerSearch(mp.m,mp.strandMut.numMutPos(), mp.strandsPresent, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, 
					minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, deeSettings.useTriples, enumSettings.useFlagsAStar, es,hbonds, mp.strandMut);

			rs.useCCD = minSettings.useCCD;

			/////////////////////////////////////////////////////////////
			// DEE section

			long startTime = System.currentTimeMillis();

			//If molecule was reloaded we need to update the mutable information
			//Set the allowable AAs for each AS residue
			boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
			if(!addWT)
				mp.strandMut.checkWT(mp.strandPresent, sParams);
			int molStrand = 0;
			for(int resID:mp.strandMut.allMut){
				setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
			}
			
				
			if((!minSettings.selectPerturbations || reload) && !mp.loadedFromCache)
					rs.setupRCs(minSettings.doPerturbations);

			long startEmat = System.currentTimeMillis();
			System.out.print("Loading precomputed energy matrix...");
			Emat emat = loadPairwiseEnergyMatrices(sParams,ematSettings.runNameEMatrixMin,minSettings.doMinimize,curStrForMatrix,es,mp.m, false);
			rs.setMinMatrix(emat);
			long endEmat = System.currentTimeMillis();
			metrics.EmatTime += (endEmat-startEmat);
			metrics.setLoopStart();

			//Prune all rotamers that don't belong to the allowed position
			emat.pruneNotAllowed(mp.m);
			emat.removePrunedRotReducedMem(false);

			addErefAndEntropy(ematSettings.useEref, rs, emat);

			int[] numRotForRes = compNumRotForRes(emat);
			int totalNumRot = 0;
			for(int q:numRotForRes){totalNumRot+=q;	System.out.print(q+" ");};System.out.println("");

			//first prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
			System.out.println("Pruning all rotamers incompatible with the template..");			
			RotamerSearch.DoPruneStericTemplate(emat, deeSettings.stericE,true,outPS);
			System.out.println();

			emat.removePrunedRotReducedMem(false);
			int[] unprunedRotForRes = emat.remainingRot();
			for(int q:unprunedRotForRes){System.out.print(q+" ");}System.out.println();

			//preprocess pairs of rotamers (mark pairwise energies greater than the cutoff as steric clashes)
			//This is only useful if these pairs will not be pruned by rs.prunedRidiculousPairs (i.e. if !useFlags)
			if (deeSettings.preprocPairs && (!deeSettings.useFlags) ){
				System.out.println("Preprocessing pairs of rotamers, cutoff of "+deeSettings.pairSt);
				emat.preprocessPairs(deeSettings.pairSt, deeSettings.stericE);
				System.out.println();
			}

			if(deeSettings.useFlags){
				double cutoff = deeSettings.stericE;
				if ( deeSettings.preprocPairs && (deeSettings.pairSt < deeSettings.stericE) )
					cutoff = deeSettings.pairSt;
				emat.pruneRidiculousPairs(cutoff);
			}


			boolean localUseMinDEEpruningEw = deeSettings.useMinDEEPruningEw;
			if(deeSettings.doDACS){ //If we are doing dacs we don't want to prune stuff too early so turn off
				localUseMinDEEpruningEw = false;   //iMinDEE until the last depth is reached

				//Without iMinDEE we can't do triples or indirect pruning though
				if(deeSettings.useTriples){
					System.out.println("Warning: Can't prune triples with DACS.  Turning off triples pruning.");
					deeSettings.useTriples = false;
					rs.useTriples = false;
				}
				if(deeSettings.algOption >= 4){
					System.out.println("Warning: Can't do indirect pruning (deeSettings.algOption 4) with DACS.  Reverting to deeSettings.algOption 3");
					deeSettings.algOption = 3;
				}
			}

			if(es.useEPIC){
				String suffix = "_" + curStrForMatrix+".dat";
				if(curStrForMatrix == COMPLEX )
					suffix = "_COM.dat";

				String CETMatrixName = (String)sParams.getValue("CETMATRIXNAME",runName+"CETM") + suffix;
				rs.loadCETMatrix(CETMatrixName);

				if( rs.cetm != null ){//Check if the continuous energy matrix is present and complete enough
					if(rs.cetm.ivalCutoff >= deeSettings.Ival)
						deeSettings.Ival = rs.cetm.ivalCutoff;
				}
				rs.cetm = null;
			}

//			wcspUpdateEmat(mp, bestScore, emat, numRotForRes);
			
			runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
					deeSettings.initEw, deeSettings.maxIntScale, deeSettings.typeDep, deeSettings.Ival,
					emat, localUseMinDEEpruningEw, true,deeSettings.stericE, sParams,deeSettings.maxFullPairs,
					deeSettings.maxDEELoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);

			if(ematSettings.saveLatestPruned)
				emat.save(ematSettings.runNameEMatrixMin+"_latestPruned_COM.dat",mp.m);
			


			long pruneTime = System.currentTimeMillis();

			if(minSettings.pertScreen){//No A*, so we're done
				PertFileHandler.writePertFile(minSettings.screenOutFile, rs.m, emat, rs.strandRot, mp.strandMut, true);
				System.out.println("Screening time: "+((pruneTime-startTime)/(60.0*1000.0)));
				System.out.println("Screening done");
				return;
			}


			if(es.useEPIC){
				loadCETMatrix(sParams,rs,curStrForMatrix,deeSettings.Ival,false, emat);
			}

			long startAStarTime = System.currentTimeMillis();

			if (!deeSettings.doDACS){ //DACS will not be performed

				if (graphSettings.genInteractionGraph) //generate interaction graph
					genInteractionGraph(mp.strandMut.numMutPos(), rs, emat, runName, mp.strandMut, graphSettings.eInteractionCutoff, graphSettings.distCutoff, mp.m, 
							deeSettings.preprocPairs, deeSettings.pairSt);

				else { //perform A* search to enumerate conformations
					
					// 2010: A* now returns a new value for Ew.  Note that right now useMinDEEdeeSettings.pruningEw
					//   is incompatible with the traditional usage of initEw.  This can be easily 
					//   fixed by adding a different Ew for this method.  If A* returns an initEw of 
					//    -1 that means that an error occured somewhere.  If it returns 0 it means 
					//    that the GMEC or minGMEC was found.  If it returns 'Ew'>0 then it means that
					//    useMinDEEdeeSettings.pruningEw = true and that the energy window must be "enlarged" 
					//    to 'Ew'
					//  useTopKHeuristic: Only the top Kvalue conformations are enumerated.  If 
					//     not enough, a new initEw is returned. Note that we may have to do this 
					//     several times.

					if(es.useEPIC){
						System.out.println("Starting lowestBound calculation for EPIC");
						rs.es.gettingLowestBound = true;
						//run A* without polynomial fits to get the lowestBound
						//this is to ensure that our deeSettings.Ival is going to be sufficient
						rs.doAStarGMEC(outputSettings.outputConfInfo,true,minSettings.doMinimize,
								mp.strandMut.numMutPos(),mp.strandMut,deeSettings.initEw,
								bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
								ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
								localUseMinDEEpruningEw, deeSettings.Ival,enumSettings);

						rs.es.gettingLowestBound = false;
					}

					AStarResults asr = rs.doAStarGMEC(outputSettings.outputConfInfo,true,minSettings.doMinimize,
							mp.strandMut.allMut.length,mp.strandMut,deeSettings.initEw,
							bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
							ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
							localUseMinDEEpruningEw, deeSettings.Ival,enumSettings);

					
					interval = asr.bestE - asr.lowestBound;
					bestScore = Math.min(bestScore, asr.bestE);
					actualLowestBound = Math.min(asr.lowestBound, actualLowestBound);
					if(es.useEPIC){
						//we want to raise the Ival to make sure that
						//we have all conformations up to minE, where minE may vary a little bit from run to run
						//depending on the EpicFitter samples (hence we need a tolerance).  
						//But the variation should be (generally much) less than thermal energy, so we use that as tolerance
						interval += RotamerSearch.constRT;
					}
					difference = interval - deeSettings.Ival;
					deeSettings.Ival = interval;

					if(difference > 0.001 && deeSettings.useMinDEEPruningEw && minSettings.doMinimize){//we're going to need a repeat run
						System.out.println("Raising deeSettings.Ival to "+deeSettings.Ival);

						if(es.useEPIC){
							if(deeSettings.Ival+deeSettings.initEw>es.EPICThresh2){//need to raise EPICThresh2 for repeat run
								//this way it will always stay at least deeSettings.Ival+initEw, so we can trust enumerations
								//by doAStarGMECHelper
								System.out.println("Raising EPICThresh2 to "+(deeSettings.Ival+deeSettings.initEw));
								es.EPICThresh2 = deeSettings.Ival+deeSettings.initEw;
							}
						}
					}
				}
			}
			else { //DACS		

				if(es.useEPIC)
					System.out.println("Warning: EPIC not supported for DACS.  Running w/o EPIC");

				numRotForRes = compNumRotForRes(emat);
				BigInteger numInitUnprunedConfs = compNumUnprunedConfs(emat);			

				int msp[] = new int[mp.strandMut.numMutPos()]; //split positions (for DACS)
				for (int i=0; i<msp.length; i++)
					msp[i] = -1;

				if (deeSettings.distrDACS) { //distributed DACS (only for level 0)

					if (deeSettings.initDepth<=0){
						System.out.println("ERROR: distributed DACS called with 'deeSettings.initDepth="+deeSettings.initDepth+"' partitioning positions; use a positive value");
						System.exit(1);
					}
					if (deeSettings.subDepth<0)
						deeSettings.subDepth = 0;

					deeSettings.distrDEE = false; //do not perform both distributed DACS and distributed DEE

					//choose the major splitting positions
					for (int i=0; i<deeSettings.initDepth; i++)
						msp[i] = chooseSplitPos(mp.strandMut.numMutPos(),emat,numRotForRes, msp, i, deeSettings.minRatioDiff);

					int maxNumPartitions = 1;
					for (int i=0; i<deeSettings.initDepth; i++)
						maxNumPartitions *= numRotForRes[msp[i]];

					OneMutation resumeResults[] = null;		
					//TODO: Fix the resume search for DACS
//					if (resumeSearch){ //read resume results
//						System.out.println("Reading resume results..");
//						resumeResults = new OneMutation[maxNumPartitions];
//						for(int q=0;q<resumeResults.length;q++)
//							resumeResults[q] = new OneMutation();
//						resumeResults = readResumeFile(resumeResults,resumeFilename,mp.numberMutable,true,false,deeSettings.initDepth);
//						System.out.println("Read "+resumeResults.length+" completed partitions.");
//						System.out.println();
//					}

					doDistrDACSMaster(runName, mp.strandMut.numMutPos(), rs, mp.strandMut, emat,
							deeSettings.algOption, deeSettings.useFlags, deeSettings.initEw, deeSettings.pruningE, deeSettings.initDepth, msp,
							numInitUnprunedConfs, deeSettings.diffFact, outputSettings.outputPruneInfo, outputSettings.outputConfInfo, 
							deeSettings.minRatioDiff, minSettings.doMinimize,
							ematSettings.runNameEMatrixMin, sParams, enumSettings.approxMinGMEC, 
							enumSettings.lambda, numRotForRes, resumeResults, resumeFilename, minSettings.minimizeBB, 
							enumSettings.numMaxMut, deeSettings.scaleInt, deeSettings.maxIntScale, 
							ematSettings.useEref, minSettings.doBackrubs, minSettings.backrubFile, deeSettings.subDepth,mp.strandPresent,
							mp.strandLimits,mp.strandsPresent,ematSettings.addWTRot, enumSettings);
				}
				else { //single-processor DACS

					deeSettings.initDepth = 0; //only used for distributed DACS
					if (deeSettings.subDepth<=0){
						System.out.println("ERROR: single-processor DACS called with 'deeSettings.subDepth="+deeSettings.subDepth+"' partitioning positions; use a positive value");
						System.exit(1);
					}

					PrintStream logPS = setupOutputFile(outputSettings.outputPruneInfo);

					PrunedRotamers<Boolean> prunedRotamers = new PrunedRotamers<Boolean>(emat.singles.pruned,false);
					doDACS(mp.strandMut.numMutPos(), rs, mp.strandMut,
							emat, prunedRotamers, deeSettings.algOption, deeSettings.useFlags, deeSettings.initEw, deeSettings.pruningE, deeSettings.initDepth, 0, logPS, msp,
							numInitUnprunedConfs, deeSettings.diffFact, outputSettings.outputConfInfo, deeSettings.minRatioDiff, minSettings.doMinimize, 
							ematSettings.runNameEMatrixMin,	deeSettings.distrDEE, sParams, enumSettings.approxMinGMEC, enumSettings.lambda, null, 
							null, minSettings.minimizeBB, enumSettings.numMaxMut, deeSettings.scaleInt, deeSettings.maxIntScale, ematSettings.useEref, 
							minSettings.doBackrubs, minSettings.backrubFile, deeSettings.subDepth,
							deeSettings.typeDep, ematSettings.addWTRot, deeSettings.Ival, deeSettings.deeSettings,enumSettings);
				}
			}

			long stopTime = System.currentTimeMillis();

			System.out.println("Pruning time: "+((pruneTime-startTime)/(60.0*1000.0)));
			//the difference between pruneTime and startAStarTime is CET matrix loading is in between (if applicable)
			if (graphSettings.genInteractionGraph)
				System.out.println("Graph generation time: "+((stopTime-startAStarTime)/(60.0*1000.0)));
			else
				System.out.println("Enumeration/DACS time: "+((stopTime-startAStarTime)/(60.0*1000.0)));
			System.out.println("DEE execution time: "+((stopTime-startTime)/(60.0*1000.0)));
			System.out.println("DEE done");
			System.out.println("Total execution time: "+((stopTime-metrics.startTime)/(60.0*1000.0)));
			//end of DEE section
			/////////////////////////////////////////////////////////////
		}
		while(difference > 0.001 && deeSettings.useMinDEEPruningEw && minSettings.doMinimize); // 2010: if I1-I0 >0 we must repeat the cycle with the new energy 
		// window.  This can only happen if useMinDEEdeeSettings.pruningEw is true
		//   and not topK
		
		metrics.setEndTime();

		metrics.trueIval = bestScore - actualLowestBound; 
		metrics.print();

	}

	private double wcspUpdateEmat(MolParameters mp, double bestScore,
			Emat emat, int[] numRotForRes) {
		boolean WCSPupdate = true;
		if(WCSPupdate){
			
//			MCSearch mcs = new MCSearch(mp.m,emat);
//			double bestMCE = mcs.doCalculation(1000000,10);
//			bestMCE -= emat.templ_E;
//			
			int[] conf = new int[emat.resByPos.size()];
			for(int i=0; i<conf.length;i++){conf[i] = -1;}
			PGQueueNode dummy = new PGQueueNode (emat.resByPos.size(), conf, Double.NEGATIVE_INFINITY,0,-1);
			WCSPOptimization wcspOpt = new WCSPOptimization(dummy, emat, null, null, emat.numRotPerPos(),Double.POSITIVE_INFINITY);
			String postprocOutFile = wcspOpt.outDir+File.separator+"postproc_"+wcspOpt.filename;
			String[] additionalCommands = {"-z=2", "-filetodump="+postprocOutFile, "-maxarity=2","-wcsponly"};
			wcspOpt.getBound(additionalCommands);
//				tbOpt.optimize(additionalCommands);
			
//				//Print the best conformation
//				PrintStream logPS = null; //the output file for conf info
//				try {			
//					FileOutputStream fileOutputStream = new FileOutputStream(outputConfInfo,true); //append file if more than 1 partition
//					BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
//					logPS = new PrintStream( bufferedOutputStream );
//				}
//				catch (Exception ex) {
//					System.out.println("ERROR: An exception occured while opening log file");
//				}
//				
//				logPS.print(1+" ");
//				for (int i=0; i<tbOpt.bestConf.conf.length; i++){
//					
//					logPS.print(tbOpt.bestConf.conf[i].eme.printRes(mp.m,emat.resByPos));	
//				}
//
//				logPS.print("unMinE: "+(tbOpt.bestConf.E+emat.templ_E)+" ");
//				logPS.print("minE: "+(tbOpt.bestConf.E+emat.templ_E)+" ");
//				logPS.println();
//				logPS.close();
//				//Done Printing conformation
			
//			double pruneE = initEw + bestMCE;
//			tbOpt.pruneEmat(postprocOutFile, pruneE);
//			bestScore = Math.min(bestScore, bestMCE );
////				bestScore = Math.min(bestScore, tbOpt.bestConf.E);
			wcspOpt.updateEmat(postprocOutFile);
			File f = new File(postprocOutFile);
			f.delete();
			wcspOpt.cleanUp();
		}
		return bestScore;
	}

	/**
	 * Given a MolParameters, reload the molecule and copy over all of the information about mutating
	 * the molecule.
	 */
	private MolParameters reloadMolecule(MolParameters mp, ParamSet sParams,
			int curStrForMatrix, boolean neighborList, double distCutoff, String pertFile) {
		
		//Reload the molecule
//		MolParameters retMP = mp.copy();
//
//		//Setup the molecule system
//		retMP.m = new Molecule();
//		retMP.m = setupMolSystem(retMP.m,sParams,mp.strandPresent,mp.strandLimits);
//		mols[curStrForMatrix+1] = mp.m;
//		
//		//Copy over the mutable information
//		retMP.m.copyMutableInfo(mp.m);
//		PertFileHandler.readPertFile( pertFile, retMP.m ,null,false);
		
		mp.m.origMol();
		return mp;
//		return retMP;
	}

	private void runDEE(boolean useFlags, boolean doMinimize,
			boolean minimizeBB, boolean scaleInt, double initEw,
			double maxIntScale, boolean typeDep, double Ival, Emat emat,
			boolean localUseMinDEEPruningEw, boolean removeRot, double stericE,
			ParamSet sParams, int timesRunFullPairs, int maxLoopNum,
			DEEsettings deeSettings, StrandRotamers[] strandRot, 
			MutableResParams strandMut, Molecule m, boolean doPerturbations) {

		long deeStart = System.currentTimeMillis();
		boolean deeDone;
		int tmpUnprunedBefore,tmpUnprunedAfter;
		int tmpUnprunedPairsBefore,tmpUnprunedPairsAfter;
		int fullPairsCtr = 0;
		deeDone = false;
		int loopNum = 0;
		int rank = -1;

		boolean doGold = deeSettings.gold;
		boolean doSplit1 = deeSettings.split1;
		boolean doSplit2 = deeSettings.split2;
		boolean doMagicBullets = deeSettings.mb;

		try{
			rank = MPItoThread.Rank();
		}catch(Exception e){
			System.out.println("Can't determine node rank.");
		}

		RotamerSearch.DoPruneStericTemplate(emat, stericE,true,outPS);
		while (!deeDone && loopNum < maxLoopNum){ //repeat the pruning cycle until no more rotamers are pruned	

			int unprunedBefore = countUnprunedRot(emat);
			int unprunedPairsBefore = countUnprunedPairs(emat);

			tmpUnprunedBefore = unprunedBefore;
			tmpUnprunedPairsBefore = unprunedPairsBefore;

			outPS.println("Starting DEE cycle run: "+loopNum);

			//Depending on the chosen algorithm option, apply the corresponding pruning criteria;			
			outPS.println("Starting Goldstein (singles): ");
			outPS.println("UnprunedRot: "+ unprunedBefore + " UnprunedPairs: "+unprunedPairsBefore);
			//Goldstein
			//KER: Always do as much goldstein as possible
			if(doGold){
				boolean goldDone = false;
				long goldStart = System.currentTimeMillis();
				while(!goldDone){

					//unprunedBefore = countUnprunedRot(emat);

					//Depending on the chosen algorithm option, apply the corresponding pruning criteria;
					if(rank == 0 && numProc > 2 && emat.numMutPos() > 2 )
						doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
								minimizeBB, scaleInt, maxIntScale,Ival,typeDep,false,Settings.DEEMETHOD.GOLDSTEIN);
					else
						RotamerSearch.DoDEEGoldstein(emat,initEw, useFlags, typeDep,
								localUseMinDEEPruningEw,Ival,false,null,null,removeRot);

					if(removeRot)
						emat.removePrunedRotReducedMem(false);

					//check how many rotamers/pairs are pruned this run
					tmpUnprunedAfter = countUnprunedRot(emat);

					outPS.println("Goldstein (singles) Pruned: "+(tmpUnprunedBefore - tmpUnprunedAfter)+" rotamers");
					outPS.println("UnprunedRot: "+tmpUnprunedAfter );

					if(tmpUnprunedBefore == tmpUnprunedAfter)
						goldDone = true;

					tmpUnprunedBefore = tmpUnprunedAfter;

				}
				long goldEnd = System.currentTimeMillis();
				outPS.println("Goldstein Pruning took: "+((goldEnd-goldStart)/1000)+" seconds");
			}

			//Split1f
			if(doSplit1){
				outPS.println("Starting SplitFlags1: ");
				long startSplit1 = System.currentTimeMillis();

				if(rank == 0 && numProc > 2 && emat.numMutPos() > 2)
					doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
							minimizeBB, scaleInt, maxIntScale,Ival,typeDep,false,Settings.DEEMETHOD.SPLITFLAGS1);
				else
					RotamerSearch.DoDEEConfSplitting(emat,initEw, null, doMinimize, useFlags, 1, 
							false, minimizeBB, typeDep,localUseMinDEEPruningEw,Ival,null);

				if(removeRot)
					emat.removePrunedRotReducedMem(false);	

				tmpUnprunedAfter = countUnprunedRot(emat);
				tmpUnprunedPairsAfter = countUnprunedPairs(emat);
				outPS.println("Split1 pruned: "+(tmpUnprunedBefore - tmpUnprunedAfter)+" rotamers "+(tmpUnprunedPairsBefore-tmpUnprunedPairsAfter)+" pairs");
				outPS.println("UnprunedRot: "+tmpUnprunedAfter + " UnprunedPairs: "+tmpUnprunedPairsAfter);

				tmpUnprunedBefore = tmpUnprunedAfter;
				tmpUnprunedPairsBefore = tmpUnprunedPairsAfter;

				long endSplit1 = System.currentTimeMillis();
				outPS.println("Split1 Pruning took: "+((endSplit1-startSplit1)/1000)+" seconds");
			}

			//Split2f
			if(doSplit2){
				outPS.println("Starting SplitFlags2: ");
				long startSplit2 = System.currentTimeMillis();
				if(rank == 0 && numProc > 2 && emat.numMutPos() > 2)
					doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
							minimizeBB, scaleInt, maxIntScale,Ival,typeDep,false,Settings.DEEMETHOD.SPLITFLAGS2);
				else
					RotamerSearch.DoDEEConfSplitting(emat,initEw, null, doMinimize, useFlags, 2, 
							false, minimizeBB, typeDep,localUseMinDEEPruningEw,Ival,null);

				if(removeRot)
					emat.removePrunedRotReducedMem(false);		

				tmpUnprunedAfter = countUnprunedRot(emat);
				tmpUnprunedPairsAfter = countUnprunedPairs(emat);
				outPS.println("Split2 pruned: "+(tmpUnprunedBefore - tmpUnprunedAfter)+" rotamers "+(tmpUnprunedPairsBefore-tmpUnprunedPairsAfter)+" pairs");
				outPS.println("UnprunedRot: "+tmpUnprunedAfter + " UnprunedPairs: "+tmpUnprunedPairsAfter);

				long endSplit2 = System.currentTimeMillis();
				outPS.println("Split2 Pruning took: "+((endSplit2-startSplit2)/1000)+" seconds");
			}

			//Magic Bullet Pairs
			if(doMagicBullets){
				outPS.println("Starting Magic Bullet Pairs: ");
				long startMB = System.currentTimeMillis();
				if(rank == 0 && numProc > 2 && emat.numMutPos() > 2)
					doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
							minimizeBB, scaleInt, maxIntScale,Ival,typeDep,true,Settings.DEEMETHOD.MBPAIRS);
				else
					RotamerSearch.DoDEEPairs(emat,initEw, null, doMinimize, useFlags, true, 
							false, minimizeBB, scaleInt, maxIntScale,typeDep,localUseMinDEEPruningEw,Ival,null);

				if(removeRot)
					emat.removePrunedRotReducedMem(false);				

				tmpUnprunedAfter = countUnprunedRot(emat);
				tmpUnprunedPairsAfter = countUnprunedPairs(emat);
				outPS.println("MB Pairs pruned: "+(tmpUnprunedBefore - tmpUnprunedAfter)+" rotamers "+(tmpUnprunedPairsBefore-tmpUnprunedPairsAfter)+" pairs");
				outPS.println("UnprunedRot: "+tmpUnprunedAfter + " UnprunedPairs: "+tmpUnprunedPairsAfter);

				long endMB = System.currentTimeMillis();
				outPS.println("MB Pairs Pruning took: "+((endMB-startMB)/1000)+" seconds");
			}

			//check how many rotamers/pairs are pruned this run
			int unprunedAfter = countUnprunedRot(emat);
			int unprunedPairsAfter = countUnprunedPairs(emat);

			if(unprunedBefore == unprunedAfter && unprunedPairsBefore == unprunedPairsAfter)
				deeDone = true;

			if(removeRot)
				emat.removePrunedRotReducedMem(false);

			//			if ((useFlags)||(algOption>=3)){
			//				System.out.println("Starting pruning with Bounding Flags");			
			//				rs.DoBoundFlags(mp.numberMutable, mp.strandMut,
			//						pruningE, prunedRotAtRes, initEw, useFlags);
			//				System.out.println();
			//			}
			//
			//			System.out.println("Starting pruning with Bounds");			
			//			prunedRotAtRes = rs.DoMinBounds(mp.numberMutable, mp.strandMut,
			//					pruningE, prunedRotAtRes, initEw, useFlags, false);
			//			System.out.println();


			//Goldstein Full Pairs
			if(deeDone && (timesRunFullPairs < 0 || fullPairsCtr < timesRunFullPairs) && emat.numMutPos() > 2){
				unprunedBefore = countUnprunedRot(emat);
				unprunedPairsBefore = countUnprunedPairs(emat);



				fullPairsCtr++;

				outPS.println("Starting pruning with DEE (full pairs)");			

				
				if(rank == 0 && numProc > 2 && emat.numMutPos() > 2)
					doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
							minimizeBB, scaleInt, maxIntScale,Ival,typeDep,true,Settings.DEEMETHOD.FULLPAIRS);
				else
					RotamerSearch.DoDEEPairs(emat,initEw, null, doMinimize, useFlags, false, 
						false, minimizeBB, scaleInt, maxIntScale,typeDep,localUseMinDEEPruningEw,Ival,null);
				
				outPS.println();




				//				if( useTriples ){
				//					//Prune Goldstein triples
				//
				//					System.out.println("Starting triples pruning");
				//					rs.DoDEETriples(mp.numberMutable, mp.strandMut, initEw, prunedRotAtRes, null,
				//							doMinimize, magicBulletTriples, magicBulletNumTriples, false,
				//							minimizeBB, typeDep, localUseMinDEEPruningEw, Ival);
				//
				//					System.out.println();
				//				}
				//
//				if(algOption >= 4){
					//Indirect pruning goes at the end because it benefits strongly from prior pruned pairs and triples
//					System.out.println("Starting indirect pruning");
//
//					RotamerSearch.DoDEEIndirect(emat, strandMut, initEw, 
//							null, doMinimize, false, false,
//							minimizeBB, typeDep, localUseMinDEEPruningEw, Ival,doPerturbations,m, strandRot);
//
//					System.out.println();
					//This prunes full rather than magic-bullet pairs
//				}



				//check how many rotamers/pairs are pruned this run
				unprunedAfter = countUnprunedRot(emat);
				unprunedPairsAfter = countUnprunedPairs(emat);

				outPS.println("DEE (full pairs) pruned: "+(unprunedBefore - unprunedAfter)+" rotamers "+(unprunedPairsBefore-unprunedPairsAfter)+" pairs");
				outPS.println("UnprunedRot: "+unprunedAfter + " UnprunedPairs: "+unprunedPairsAfter);

				if(!(unprunedPairsBefore == unprunedPairsAfter) ||
						!(unprunedAfter == unprunedBefore)){
					deeDone = false;
				}


			}
			loopNum++;
		}
		//KER: This is merely informative but takes a while to compute so we might want to take it out
		emat.checkIfAllPruned(null);

		long deeEnd = System.currentTimeMillis();
		metrics.DEEtime += (deeEnd - deeStart);

		

		//checkBest(sParams,rs, emat);
	}

	private int countUnprunedRot(Emat emat){
		int countPruned = 0;
		int totalRot = 0;
		Iterator<EMatrixEntryWIndex> iter = emat.singlesIterator();
		while(iter.hasNext()){
			EMatrixEntry re = iter.next().eme;
			if (re.isPruned())
				countPruned++;
			totalRot++;
		}
		int retArray = totalRot - countPruned; 
		return retArray;
	}

	private int countUnprunedPairs(Emat emat){
		int countPruned = 0;
		int totalPairs = 0;
		Iterator<EMatrixEntryWIndex> iter = emat.pairsIterator();
		while(iter.hasNext()){
			EMatrixEntry re = iter.next().eme;
			if (re.isPruned())
				countPruned++;
			totalPairs++;
		}
		int retArray = totalPairs - countPruned;
		return retArray;
	}

	private void addErefAndEntropy(boolean useEref, RotamerSearch rs, Emat emat) {
		if(useEref){
			if(emat.eRef == null){
				System.out.println("Reference energies turned on but not calculated");
				System.out.println("Exiting....");
				System.exit(0);
			}
			if( !emat.hasEref)
				rs.addEref(emat.eRef);
		}
		if(EnvironmentVars.useEntropy && !emat.hasEntropy)
			rs.addEntropyTerm();
	}

	//KER: If addWT isn't true and there is only one AA type for a mutable
	//position we treat that position as WT
//	private void checkWT(String[][] strandDefault,boolean[] strandPresent, ParamSet sParams) {
//		int ctr = 0;
//		for(int str=0; str<strandPresent.length; str++){
//			if(strandPresent[str]){
//				for(int i=0; i<strandDefault[ctr].length;i++){
//					String tempResAllow = (String)sParams.getValue("RESALLOWED"+str+"_"+i);
//					if(numTokens(tempResAllow)==1)
//						strandDefault[ctr][i] = getToken(tempResAllow, 1);
//				}
//				ctr++;
//			}
//		}	
//	}


	class MolParameters{
		String[][] strandLimits;
		int numOfStrands;
		boolean[] strandPresent;
		int strandsPresent;
		MutableResParams strandMut;
		String[][] strandDefault;
		Molecule m;
//		int numberMutable;
		boolean loadedFromCache = false;
		
		MolParameters(){
			
		}
		
		public MolParameters copy(){
			MolParameters retMP = new MolParameters();
			retMP.strandLimits = strandLimits;
			retMP.numOfStrands = numOfStrands;
			retMP.strandPresent = strandPresent;
			retMP.strandsPresent = strandsPresent;
			retMP.strandMut = strandMut;
			retMP.strandDefault = strandDefault;
			retMP.m = m;
			retMP.loadedFromCache = loadedFromCache;
			return retMP;
		}
		
	}

	//	public MolParameters loadMolecule(ParamSet sParams, int curStrForMatrix, Molecule m){	
	//		
	//		MolParameters mp = new MolParameters();
	//		
	//		loadStrandParams(sParams, mp, curStrForMatrix);
	//		
	//		//Setup the molecule system
	//		if(m != null){
	//			mols[curStrForMatrix+1] = m;
	//			mp.m = mols[curStrForMatrix+1];
	//		}
	//		else{
	//			mp.m = mols[curStrForMatrix+1];
	//		}
	//		
	//		return mp;
	//		
	//	}

	public MolParameters loadMolecule(ParamSet sParams, int curStrForMatrix, 
			boolean neighborList, double distCutoff, boolean useCache){	

		MolParameters mp = new MolParameters();

		loadStrandParams(sParams, mp, curStrForMatrix);

		//Setup the molecule system

		//Setup the molecule system
		if(mols[curStrForMatrix+1] == null || !useCache){
			mp.m = new Molecule();
			mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);
			mols[curStrForMatrix+1] = mp.m;
		}else{
			mp.m = mols[curStrForMatrix+1];
			mp.loadedFromCache = true;
		}

		loadMutationParams(sParams, mp);

		if(neighborList){
			mp.m.genDistNeighborList(distCutoff);
			mp.m.dumpNeighborList();
		}

		return mp;

	}


	private int getNumberMutable(int[][] strandMut){
		int numberMutable = 0;
		for (int i=0; i<strandMut.length;i++)
			numberMutable += strandMut[i].length;

		return numberMutable;
	}

	private void loadStrandParams(ParamSet sParams, MolParameters mp, int curStrForMatrix){
		mp.numOfStrands = (new Integer((String)sParams.getValue("NUMOFSTRANDS"))).intValue();	
		mp.strandLimits = new String[mp.numOfStrands][2];
		mp.strandPresent = new boolean[mp.numOfStrands];
		mp.strandsPresent = 0;
		for (int i=0; i<mp.numOfStrands; i++){
			String strandLimit = (String)sParams.getValue("STRAND"+i);
			String limit1 = getToken(strandLimit,1);
			String limit2 = getToken(strandLimit,2);
			mp.strandLimits[i][0] = limit1;
			mp.strandLimits[i][1] = limit2;
			if(curStrForMatrix == COMPLEX)
				mp.strandPresent[i] = true;
			else if(curStrForMatrix == i)
				mp.strandPresent[i] = true;
			else
				mp.strandPresent[i] = false;
			//strandPresent[i] = (new Boolean((String)sParams.getValue("STRANDPRESENT" + i))).booleanValue();
			if(mp.strandPresent[i]){
				mp.strandsPresent++;
			}
		}
	}



	//KER: mp must already have the following terms set:
	//mp.numOfStrands	
	//mp.strandLimits
	//mp.strandPresent
	//mp.strandsPresent
	private void loadMutationParams(ParamSet sParams, MolParameters mp) {

		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","true"))).booleanValue();
		/**********Get the regions of each strand that are mutable****/
		String strandMutNums = (String)sParams.getValue("STRANDMUTNUMS");
		boolean addOrigRots = (new Boolean((String)sParams.getValue("ADDWTROTS", "false"))).booleanValue();
		int totalNumMut = 0;
		for(int i=0; i<mp.strandPresent.length; i++){
			if(mp.strandPresent[i])
				totalNumMut += (new Integer(KSParser.getToken(strandMutNums,i+1))).intValue();
		}
		mp.strandMut = new MutableResParams(totalNumMut,mp.m.numberOfStrands);  //taking the place of resMap and ligMap

		int flatCtr = 0;
		int strCtr = 0;
		try{
			for(int i=0; i<mp.strandPresent.length; i++){
				if(mp.strandPresent[i]){
					int numberOfMutables = (new Integer(getToken(strandMutNums,i+1))).intValue();
					//mp.strandMut[ctr] = new int[numberOfMutables];
					String strandMutResNum = (String)sParams.getValue("STRANDMUT"+i);
					for(int j=0; j<numberOfMutables; j++){
						String strandMutRes = getToken(strandMutResNum,j+1);
						Residue r = mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)];
						mp.strandMut.addRes(flatCtr,r,mp.m.rotLibForStrand(strCtr),addOrigRots);
						if(mp.m.strand[r.strandNumber].isProtein)
							mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)].canMutate = true;
						flatCtr++;
					}
					strCtr++;
				}
			}
			strCtr = 0;
			//Say rotamers have been set
			for(int i=0; i<mp.strandPresent.length;i++){
				if(mp.strandPresent[i]){
					mp.m.rotLibForStrand(strCtr).setAddedRotamers(true);
					strCtr++;
				}
			}
			if(!addWT)
				mp.strandMut.checkWT(mp.strandPresent, sParams);
		}

		catch(Exception E){
			System.out.println("PROBLEM Loading the strand mut numbers (Check System.cfg) ");
			E.printStackTrace();
			return;
		}
	}

	//Implements threads for DACS: allows the current best score among all partitions to be distributed to every partition;
	//This thread performs the DACS computation, while the main thread monitors for updates of the best energy from the other partitions;
	//The communication is performed via the common RotamerSearch object (synchronized access to the bestEMin variable)


	/**
	 * Performs the DACS partition-specific computation. The parameters 'rs' and 'rs.sysLR' must be valid.
	 */
	private void doDACS(int numMutable, RotamerSearch rs, MutableResParams strandMut,
			Emat emat, PrunedRotamers<Boolean> prunedRotAtRes, int algOption, boolean useFlags, double initEw, double pruningE,
			int initDepth, int curDepth, PrintStream logPS, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			boolean distrDEE, ParamSet sParams, 
			boolean approxMinGMEC, double lambda, Index3 partIndex[], CommucObj cObj, 
			boolean minimizeBB, int numMaxMut, boolean scaleInt, double maxIntScale, boolean useEref,
			boolean doBackrubs, String backrubFile, int subDepth,boolean typeDep, boolean addWTRot, double Ival,
			DEEsettings deeSettings, Settings.Enum enumSettings){

		if (curDepth>=(initDepth+subDepth))
			return;

		System.out.println("Starting pruning with DEE (DACS).");

		//prunedRotAtRes[] should not be modified here, in order to be able to distinguish
		//	newly pruned rotamers and rotamers pruned by Bounds or DEE

		//the num rotamers for each AS residue and the ligand (if present);
		//	sysLR in rs must be valid (with all the possible AA's for each residue position)
		int numRotForRes[] = compNumRotForRes(emat);

		if (curDepth>=initDepth) {//sub-partition, so majorSplitPos[curDepth] is unknown; compute it
			majorSplitPos[curDepth] = chooseSplitPos(numMutable,emat,numRotForRes,
					majorSplitPos, curDepth, minRatioDiff);//the splitting position
		}

		int numPartitions = numRotForRes[majorSplitPos[curDepth]]; //the number of partitions
		int numPrunedPartitions = 0; //the number of partitions with no unpruned confs

		System.out.println("Current depth: "+curDepth);
		System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
		System.out.println();

		//map rotamer index to rot num for the splitting residue
		Index3[] indexMap = getIndexMap(numPartitions,emat,majorSplitPos[curDepth],strandMut);

		//count the total num confs and the num unpruned confs
		BigInteger numConfsTotalForPartition[] = new BigInteger[numPartitions];
		BigInteger numConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		//BigInteger numInitConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		BigInteger numTotalConfs = new BigInteger("0");
		BigInteger numUnprunedConfs = new BigInteger("0");
		BigInteger numEvaluatedConfs = new BigInteger("0");

		//update the best energy found
		pruningE = (double)Math.min(pruningE, rs.getBestE());
		double bestScore = pruningE;

		boolean[][][][][][] savedSplitFlags = emat.copyPairPruned();

		//determine the prunings for each of the sub-solutions (the partitions)
		PrunedRotamers<Boolean> prunedForPartition[] = new PrunedRotamers[numPartitions];//[prunedRotAtRes.length];
		for(int i=0; i<numPartitions;i++)
			prunedForPartition[i] = new PrunedRotamers<Boolean>(emat.singles.pruned,false);
		for (int i=0; i<prunedForPartition.length; i++){

			if ((curDepth>=initDepth)||(indexMap[i]==partIndex[curDepth])) { //sub-partitions or current partition is the partition distributed for computation

				//copy the prunings from before conf splitting (from Bounds and simple traditional DEE)
				//System.arraycopy(prunedRotAtRes,0,prunedForPartition[i],0,prunedRotAtRes.length);
				Iterator<EMatrixEntryWIndex> iter1 = emat.singlesIterator();
				while(iter1.hasNext()){
					EMatrixEntryWIndex emeWI = iter1.next();
					prunedForPartition[i].set(emeWI.index, emeWI.eme.isPruned());
				}

				Index3 curPartIndex = indexMap[i]; //the rotamer index of the partitioning rotamer

				//artificially set all rotamers at the splitting position, other than the rotamer 
				//	for the current partition, to pruned, so that there will be only one rotamer at
				//	that residue position when the conf splitting criterion is applied;
				//	when done, subtract the artifcially pruned rotamers from the total number of pruned
				boolean indToUnprune[] = new boolean[numPartitions];
				for (int j=0; j<indToUnprune.length; j++)
					indToUnprune[j] = false;

				//check the partition only if the current partitioning rotamer is not already pruned
				if (!emat.getSinglePruned(curPartIndex)){

					for (int j=0; j<numPartitions; j++){
						if (j!=i) {//not the rotamer for the current partition
							Index3 curInd = indexMap[j]; //the index of the current rotamer
							if (!prunedForPartition[i].get(curInd)){ //not pruned by the other DEE methods
								prunedForPartition[i].set(curInd, true);
								indToUnprune[j] = true;
							}
						}
					}

					emat.setPairPruned(savedSplitFlags);

					int numPrunedRot = countPrunedRot(prunedForPartition[i]);
					int numPrunedPairs = 0;
					if ((useFlags)||(algOption>=3))
						numPrunedPairs = countPrunedPairs(rs.getMinMatrix().pairs.pruned);
					int numPrunedRotThisRun = 0;
					int numPrunedPairsThisRun = 0;
					boolean done = false;
					int numRuns = 1;

					runDEE(useFlags, doMinimize, minimizeBB, scaleInt,
							initEw, maxIntScale, typeDep, Ival,
							emat, true, true,Double.POSITIVE_INFINITY, sParams,1,
							3,deeSettings,rs.strandRot,strandMut,rs.m,rs.doPerturbations);

				}

				//count the number of pruned rotamers for each residue (except for the splitting residue,
				//	which always has only 1 available rotamer for the current partition)
				int numPrunedRotForRes[] = new int[numMutable]; //after pruning for this partition
				//int numInitPrunedRotForRes[] = new int[numInAS]; //initial prunings, before pruning for this partition
				for (int j=0; j<numMutable; j++){
					numPrunedRotForRes[j] = 0;
					if (j!=majorSplitPos[curDepth]){
						Iterator<RotInfo<Boolean>> iter = prunedForPartition[i].iterator(j);
						RotInfo<Boolean> ri = iter.next();
						while(ri.curPos == j && iter.hasNext()){
							//for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
							if (prunedForPartition[i].get(ri))
								numPrunedRotForRes[j]++;
							//if (prunedRotAtRes[j*totalNumRotamers + k])
							//numInitPrunedRotForRes[j]++;
							ri = iter.next();
						}
					}
					else {// j==majorSplitPos
						numPrunedRotForRes[j] = 0;
						//numInitPrunedRotForRes[j] = 0;
					}
				}
				/*int numPrunedLigRot = 0;
				//int numInitPrunedLigRot = 0;
				for (int k=0; k<numLigRotamers; k++){
					if (prunedForPartition[i][numInAS*totalNumRotamers + k])
						numPrunedLigRot++;
					//if (prunedRotAtRes[numInAS*totalNumRotamers + k])
					//	numInitPrunedLigRot++;
				}*/

				//count the total num confs and the num unpruned confs
				numConfsTotalForPartition[i] = new BigInteger("1");
				numConfsUnprunedForPartition[i] = new BigInteger("1");
				//numInitConfsUnprunedForPartition[i] = new BigInteger("1");
				if (emat.getSinglePruned(curPartIndex)){ //current partitioning rotamer already pruned, so no unpruned confs for this partition
					numConfsUnprunedForPartition[i] = new BigInteger("0");
					numPrunedPartitions++;
				}

				for (int j=0; j<numMutable; j++){
					if (!(isSplitRes(j,majorSplitPos,curDepth))){ //the split residues contribute only 1 rotamer
						numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]));
						numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
						//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numInitPrunedRotForRes[j]));
					}
				}
				/*if(ligPresent){
					numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numLigRotamers));
					numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
					//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numInitPrunedLigRot));
				}*/			

				numTotalConfs = numTotalConfs.add(numConfsTotalForPartition[i]);
				numUnprunedConfs = numUnprunedConfs.add(numConfsUnprunedForPartition[i]);

				BigInteger pruneDiffFact = BigInteger.valueOf(10).pow(diffFact);

				System.out.println("Num unpruned confs: "+numConfsUnprunedForPartition[i]+" diffFact: "+pruneDiffFact);
				System.out.println();
				System.out.println();

				//output pruning info to file
				logPS.print("curDepth: "+curDepth+" curPartition: "+i+" majorSplitPos: "+majorSplitPos[curDepth]+" ");
				logPS.print(numConfsTotalForPartition[i]+" "+numInitUnprunedConfs+" "+numConfsUnprunedForPartition[i]);
				logPS.println();
				logPS.println();logPS.flush();

				//if ((curDepth+1<maxDepth)&&(numConfsUnprunedForPartition[i].compareTo(numInitUnprunedConfs.divide(pruneDiffFact))==1)){ //not enough pruned, so partition at new depth
				if ((curDepth+1<(initDepth+subDepth))&&(numConfsUnprunedForPartition[i].compareTo(pruneDiffFact)==1)){ //not enough pruned, so partition at new depth
					doDACS(numMutable, rs, strandMut, emat,
							prunedForPartition[i], algOption, useFlags, initEw, pruningE, initDepth, curDepth+1, 
							logPS, majorSplitPos, numConfsUnprunedForPartition[i], diffFact, outputConfInfo, minRatioDiff,
							doMinimize, minPEM, distrDEE, sParams, approxMinGMEC, 
							lambda, partIndex, cObj, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, doBackrubs, 
							backrubFile, subDepth,typeDep,addWTRot,Ival,deeSettings,enumSettings);
				}
				else if (!prunedRotAtRes.get(curPartIndex)){ //if enough pruned or maxDepth partitioning reached, do the rotamer search

					bestScore = Math.min(bestScore,rs.getBestE());//best E for the partitions so far

					//Do the rotamer search
					rs.doAStarGMEC(outputConfInfo,true,doMinimize,numMutable,strandMut,
							initEw,bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,doBackrubs,
							backrubFile, false, 0.0,enumSettings);

					numEvaluatedConfs = numEvaluatedConfs.add(rs.numConfsEvaluated); //add the evaluated confs for this partition
					pruningE = Math.min(pruningE,rs.getBestE());//update cutoff energy for MinBounds
				}

				//unprune the artificially pruned indices
				for (int j=0; j<indToUnprune.length; j++){
					if (indToUnprune[j]){
						prunedForPartition[i].set(indexMap[j], false);
					}
				}
			}
		}

		if (cObj==null){ //not distributed DACS

			System.out.println("numTotalConfs: "+numTotalConfs+"; numUnprunedConfs: "+numUnprunedConfs+"; numEvaluatedConfs: "+numEvaluatedConfs);
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsTotalForPartition[i]+" ");System.out.println();
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsUnprunedForPartition[i]+" ");
			System.out.println();
			System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
			System.out.println("Number of partitions: "+numPartitions);
			System.out.println("Number of non-zero partitions: "+(numPartitions-numPrunedPartitions));
			System.out.println("Additional pruned rotamers: ");

			//count the number of partitions for which a rotamer is pruned (counting only the rotamers
			//	not pruned by Bounds or simple traditional DEE)
			PrunedRotamers<Integer> countNumPartitionsPrunedRot = new PrunedRotamers<Integer>(prunedRotAtRes,0);
			Iterator<RotInfo<Boolean>> iter = prunedRotAtRes.iterator();
			while(iter.hasNext()){
				RotInfo<Boolean> ri = iter.next();
				//for (int i=0; i<prunedRotAtRes.length; i++){//for each rotamer
				//countNumPartitionsPrunedRot[i] = 0;
				if (!ri.state){ //only if not pruned by the other two methods
					for (int j=0; j<numPartitions; j++){ //check for each partition
						if (prunedForPartition[j].get(ri))
							countNumPartitionsPrunedRot.set(ri, countNumPartitionsPrunedRot.get(ri)+1);
					}
				}

				//output information
				if (countNumPartitionsPrunedRot.get(ri)>0)
					System.out.println("index: "+ri.printCoord()+"; num partitions in which pruned: "+countNumPartitionsPrunedRot.get(ri));
			}
		}
		else { //distributed DACS
			cObj.bestScore = cObj.bestScore.min(BigDecimal.valueOf(rs.getBestE()));
		}
	}

	/**
	 * Implements threads for DACS: allows the current best score among all partitions to be distributed to every partition;
	 * This thread performs the DACS computation, while the main thread monitors for updates of the best energy from the other partitions;
	 * The communication is performed via the common RotamerSearch object (synchronized access to the bestEMin variable)
	 */
	private class DACSthread implements Runnable {
		private RotamerSearch rs = null;
		private CommucObj cObj = null;
		private PrintStream logPS = null;
		String outputConfInfo = null;
		DACSthread(RotamerSearch rsP, CommucObj cObjP, PrintStream lP, String ociP){
			rs = rsP;
			cObj = cObjP;
			logPS = lP;
			outputConfInfo = ociP;
		}
		public void run(){
			//Perform DACS
			doDACS(cObj.strandMut.allMut.length, rs, cObj.strandMut,
					cObj.emat, cObj.prunedRot, cObj.algOption, cObj.useSF, cObj.initEw, cObj.pruningE,
					cObj.initDepth, 0, logPS, cObj.msp, cObj.numInitUnprunedConf,
					cObj.diffFact, outputConfInfo, cObj.minRatioDiff, cObj.doMinimization, null, 
					false, cObj.params,  
					cObj.approxMinGMEC, cObj.lambda, cObj.partIndex, cObj, cObj.minimizeBB, cObj.numMutations,
					cObj.scaleInt, cObj.maxIntScale, cObj.useEref, cObj.doBackrubs, cObj.backrubFile, 
					cObj.subDepth,cObj.typeDep,cObj.addWTRot,
					cObj.Ival, cObj.deeSettings, cObj.enumSettings);
		}
	}

	//Compute the number of unpruned conformations 
	//Compute the number of unpruned conformations 
	private BigInteger compNumUnprunedConfs(Emat emat) {


		int[] unprunedRot = emat.remainingRot();
		BigInteger numConfsUnpruned = new BigInteger("1");
		for (int j=0; j<unprunedRot.length; j++){
			numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(unprunedRot[j]));
		}

		return numConfsUnpruned;

	}

	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position)
	private int chooseSplitPosRandom(int numMutable){		
		Random randNum = new Random();
		return randNum.nextInt(numMutable);
	}

	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position);
	//	Choose the AS residue position with the smallest fraction of pruned rotamers (from MinBounds and simple MinDEE)
	private int chooseSplitPos(int numMutable, Emat emat, int numRotForRes[],
			int majorSplitPos[], int curDepth, double minRatioDiff){		

		final int minPartitions = 5; //the min number of rotamers that a splitting residue can have

		double pruneRatio[] = new double[numMutable];
		int minPos = -1;
		double minRatio = (double)Math.pow(10,38);
		for (int curRes=0; curRes<numMutable; curRes++){

			if (!(isSplitRes(curRes,majorSplitPos,curDepth))){
				if (numRotForRes[curRes]>=minPartitions){ //do not split at residues with very small number of rotamers
					int curPruned = 0;
					Iterator<EMatrixEntryWIndex> iter = emat.singlesIterator(curRes);
					while(iter.hasNext()){
						EMatrixEntryWIndex emeWI= iter.next();
						//for (int curRot=0; curRot<numTotalRotamers; curRot++){
						if (emeWI.eme.isPruned()){//prunedRotAtRes[curRes*numTotalRotamers+curRot]){ //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
							curPruned++;
						}
					}
					pruneRatio[curRes] = (double)curPruned/numRotForRes[curRes];
					if (minRatio>=pruneRatio[curRes]){
						if ((minPos==-1)||(curRes<minPos)||(minRatio>=pruneRatio[curRes]+minRatioDiff)) {//preference to split at lower-numbered residues
							minRatio = pruneRatio[curRes];
							minPos = curRes;
						}
					}
				}
			}
		}

		if (minPos!=-1){
			//System.out.println("minPos: "+minPos);
			//for (int i=0;i<numInAS;i++)System.out.print(pruneRatio[i]+" ");System.out.println();
			return minPos;
		}
		else //if split position not chosen, choose randomly
			return chooseSplitPosRandom(numMutable);
	}

	//Check if the residue curRes is one of the splitRes
	private boolean isSplitRes(int curRes, int majorSplitPos[], int curDepth){
		for (int i=0; i<=curDepth; i++){
			if (curRes==majorSplitPos[i])
				return true;
		}
		return false;
	}

	//Compute the number of rotamers for each residue position (assign to numRotForRes[])
	private int [] compNumRotForRes(Emat emat){

		int numberMutable = emat.singles.E.length;
		int numRotForRes[] = new int[numberMutable];
		//boolean ligPresent = (numLigRot==0); //ligand present
		int treeLevels = numberMutable;
		/*if (ligPresent)
			treeLevels++;*/

		numRotForRes = new int[treeLevels];

		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			for(int curAA=0; curAA<emat.singles.E[curLevel].length;curAA++)
				numRotForRes[curLevel] += emat.singles.E[curLevel][curAA].length;
		}
		return numRotForRes;
	}


	//Get the mapping between rotamer indices (into the pruning matrix) and the number of the
	//	current rotamer for the giveen residue; assumes sysLR in rs is valid (all allowables for the AS residues)
	private Index3 [] getIndexMap(int numPartitions, Emat emat, int curRes, MutableResParams strandMut){

//		Residue r = m.residue[strandMut.allMut[curRes]];
//		int str = r.strandNumber;
//		int strResNum = r.strandResidueNumber;

		Index3 indexMap[] = new Index3[numPartitions];

		SinglesIterator iter = emat.singlesIterator(curRes);
		int indNum = 0;
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			indexMap[indNum] = new Index3(emeWI.index);
			indNum++;
		}
		
//		int indNum = 0;
//		for (AARotamerType aaType: r.AATypesAllowed()){ //for each AA for the given AS residue
////			int curAA = rs.strandRot[str].getIndexOfNthAllowable(strResNum,AA);
//			int numRotForAA = r.getRCsForType(aaType);//rs.getNumRot(str, strResNum, curAA);
//
//			for (int curRot=0; curRot<numRotForAA; curRot++){ //for each rot for the given AA
//				indexMap[indNum] = new Index3(curRes,curAA,curRot);//curRes*numTotalRot + rotamerIndexOffset[curAA] + curRot;
//				indNum++;
//			}
//		}
		return indexMap;
	}

	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
	}	

	///////////////////////////////////////////////////////////////////////////
	//	End of DEE section
	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//	Distributed DACS section
	///////////////////////////////////////////////////////////////////////////
	/**
	 * Handles the distribution of the DACS computation to the set of available processors
	 */
	private void doDistrDACSMaster(String runName, int numMutable, RotamerSearch rs, MutableResParams strandMut,
			Emat emat, int algOption, boolean useFlags, double initEw, double pruningE,
			int initDepth, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputPruneInfo, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			ParamSet sParams, 
			boolean approxMinGMEC, double lambda, int numRotForRes[], OneMutation resumeResults[], String resumeFileName, boolean minimizeBB, int numMaxMut,
			boolean scaleInt, double maxIntScale, boolean useEref, boolean doBackrubs, String backrubFile, int subDepth,
			boolean [] strandPresent, String[][] strandLimits, int strandsPresent, boolean addWTRot, Settings.Enum enumSettings){

		System.out.println("Starting DACS (distributed)");		
		System.out.println("Forming DACS partitions..");

		Index3 indexMap[][] = new Index3[initDepth][];
		int numPartitions[] = new int[indexMap.length];
		int maxNumPartitions = 1;

		//map rotamer index to rot num for the splitting residues
		for (int i=0; i<indexMap.length; i++){
			numPartitions[i] = numRotForRes[majorSplitPos[i]]; //the number of partitions
			indexMap[i] = getIndexMap(numPartitions[i],emat,majorSplitPos[i],strandMut);
			maxNumPartitions *= numPartitions[i];
		}

		//get the partitions
		OneMutation mutArray[] = formDACSpartitions(maxNumPartitions, initDepth, indexMap, numPartitions, emat, majorSplitPos);		

		if (resumeResults!=null){ //remove completed partitions
			int curMut = 0;
			OneMutation tmpArray2[] = new OneMutation[mutArray.length];
			for (int i=0; i<mutArray.length; i++){
				boolean partFound = true;
				for (int j=0; j<resumeResults.length; j++){
					partFound = true;
					for (int k=0; k<initDepth; k++){
						if (mutArray[i].resMut[k]!=resumeResults[j].resMut[k]){
							partFound = false;
							break;
						}
					}
					if (partFound) //partition already computed
						break;
				}
				if (!partFound){
					tmpArray2[curMut] = new OneMutation();
					tmpArray2[curMut].mutNum = mutArray[i].mutNum;
					tmpArray2[curMut].resMut = mutArray[i].resMut;
					curMut++;
				}
			}
			OneMutation tmpArray[] = new OneMutation[curMut]; //trim the size of the partition array
			System.arraycopy(tmpArray2, 0, tmpArray, 0, curMut);
			mutArray = tmpArray;

			System.out.println("Number non-zero partitions after removing completed results: "+mutArray.length);

			//update the best energy so far
			for (int j=0; j<resumeResults.length; j++){
				pruningE = (double)Math.min(pruningE, resumeResults[j].score.doubleValue());
			}
		}

		//output the rs object
		//outputObject(rs,rotFile);

		//sort the partitions
		sortDACSpartitions(mutArray,initDepth,majorSplitPos,numPartitions,emat,indexMap,numMutable,
				strandMut,pruningE,initEw,rs);

		System.out.println();

		MutationManager mutMan = new MutationManager(runName,mutArray,false);
		mutMan.setStrandMut(strandMut);
		mutMan.setStrandPresent(strandPresent);
		mutMan.setStrandLimits(strandLimits);
		mutMan.setStrandsPresent(strandsPresent);
//		mutMan.setMutableSpots(numMutable);
		mutMan.setarpFilenameMin(minPEM);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
//		mutMan.setMutableSpots(numMutable);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setNumMutations(numMaxMut);
		mutMan.setInitEw(initEw);
		mutMan.setPruningE(pruningE);
		mutMan.setUseSF(useFlags);
		mutMan.setDistrDACS(true);
		mutMan.setDistrDEE(false);
		mutMan.setBestScore(new BigDecimal(pruningE)); //the best E initially is the pruningE read from the parameter file
		mutMan.setAlgOption(algOption);
		mutMan.setInitDepth(initDepth);
		mutMan.setSubDepth(subDepth);
		mutMan.setDiffFact(diffFact);
		mutMan.setMinRatioDiff(minRatioDiff);
		mutMan.setNumInitUnprunedConf(numInitUnprunedConfs);
		mutMan.setOutputPruneInfo(outputPruneInfo);
		mutMan.setOutputConfInfo(outputConfInfo);
		mutMan.setMSP(majorSplitPos);
		mutMan.setApproxMinGMEC(approxMinGMEC);
		mutMan.setLambda(lambda);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setUseEref(useEref);
		mutMan.setRotamerLibrary(rs.m.aaRotLib);

		mutMan.setIdealizeSC(Perturbation.idealizeSC);
		mutMan.setAddWTRot(addWTRot);
		mutMan.setEnumSettings(enumSettings);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}

	}

	// Distributed DACS Slave function
	private CommucObj doDistrDACSSlave(CommucObj cObj) {

		long startTime = System.currentTimeMillis();

		//boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits);

		RotamerSearch rs = (RotamerSearch)readObject(cObj.rotFileIn); //load the saved rs from master

		Perturbation.idealizeSC = cObj.idealizeSC;

		String fn = "";
		for (int i=0; i<cObj.partIndex.length; i++)
			fn += ("_"+cObj.partIndex[i]);
		String outputConfInfo = ("./conf_info/"+cObj.outputConfInfo+fn);
		PrintStream logPS = setupOutputFile("./conf_info/"+cObj.outputPruneInfo+fn);

		if(logPS == null){
			System.out.println("ERROR: Please make folder conf_info!!!");
			return null;
		}

		double otherBestE = cObj.pruningE; //the best energy from other partitions

		Thread t = new Thread(new DACSthread(rs,cObj,logPS,outputConfInfo)); //create a new thread for performing the DACS search
		t.start();
		long waitTime = 300000; //five minutes
		while (t.isAlive()){ //DACS search thread is still running

			try{ t.join(waitTime);} catch (Exception e){} //wait for waitTime before the next interruption of the DACS search		

			double rsBestE = rs.getBestE(); //the best energy from the current partition

			//System.out.println("partition "+cObj.partIndex+": curBestE "+bestE+" rsBestE "+rsBestE);

			if (rsBestE<otherBestE){ //new best energy for the current partition; update

				CommucObj c[] = new CommucObj[1];
				c[0] = new CommucObj();
				c[0].pruningE = rsBestE;

				//System.out.println("partition "+cObj.partIndex+": sending update to main node..");

				try { MPItoThread.Send(c, 0, 1, ThreadMessage.OBJECT, 0, updateTag);} catch (Exception e){}; //send back updated best energy

				otherBestE = rsBestE;
			}

			//check if there are updates for the best energy from the other partitions
			try {				
				//System.out.println("partition "+cObj.partIndex+": checking for update from main node..");

				double c[] = new double[1];
				while (MPItoThread.Iprobe(0, updateTag)!=null) { //new update message received

					//System.out.println("partition "+cObj.partIndex+": update from main node received..");

					MPItoThread.Recv(c, 0, 1, ThreadMessage.FLOAT, 0, updateTag);

					//System.out.println("partition "+cObj.partIndex+": updateE "+c[0]+" curBestE: "+bestE);

					rs.updateBestE(c[0]);
					otherBestE = Math.min(otherBestE,c[0]);				
				}
			}
			catch (Exception e){};
		}		

		logPS.flush();
		logPS.close();		

		rs = null;
		cObj.prunedRot = null;//smaller object, for sending back		

		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);

		//System.out.println("Partition "+cObj.partIndex+" done, time: "+cObj.elapsedTime/60.0f+" minutes; sending results to main node..");

		return cObj;
	}	

	//Forms the DACS partitions based on the splitting residue positions
	private OneMutation [] formDACSpartitions(int maxNumPartitions, int initDepth, Index3 indexMap[][], int numPartitions[],
			Emat emat, int majorSplitPos[]){

		OneMutation mutArray[] = new OneMutation[maxNumPartitions];
		Index3 curInd[] = new Index3[initDepth];
		int curMut[] = new int[1];
		curMut[0] = 0;

		formDACSpartitionsHelper(initDepth, indexMap, numPartitions, emat, mutArray, curMut, 0, curInd);

		OneMutation tmpArray[] = new OneMutation[curMut[0]]; //trim the size of the partition array
		System.arraycopy(mutArray, 0, tmpArray, 0, curMut[0]);
		mutArray = tmpArray;

		System.out.print("Partitioning residues: ");
		for (int i=0; i<initDepth; i++)
			System.out.print(majorSplitPos[i]+" ");
		System.out.println();
		System.out.println("Number of non-zero partitions: "+curMut[0]);

		return mutArray;
	}

	//Determines all non-pruned partitions for DACS deistribution
	//Called by formDACSpartitions(.)
	private void formDACSpartitionsHelper(int initDepth, Index3 indexMap[][], int numPartitions[],
			Emat emat, OneMutation mutArray[], int curMut[], int curDepth, Index3 curInd[]){

		if (curDepth>=initDepth){ //new partition
			mutArray[curMut[0]] = new OneMutation();
			mutArray[curMut[0]].mutNum = curMut[0];
			mutArray[curMut[0]].index = new Index3[initDepth];
			for (int i=0; i<initDepth; i++)
				mutArray[curMut[0]].index[i] = curInd[i];
			curMut[0]++;
		}
		else {
			for (int i=0; i<numPartitions[curDepth]; i++){
				if (!emat.getSinglePruned(indexMap[curDepth][i])){ //only consider non-pruned partitions for distribution
					curInd[curDepth] = indexMap[curDepth][i];
					formDACSpartitionsHelper(initDepth, indexMap, numPartitions, emat, mutArray, curMut, curDepth+1, curInd);
				}
			}
		}
	}

	//Sorts the DACS partitions by lower energy bounds;
	//The sorted array is returned in mutArray[]
	private void sortDACSpartitions(OneMutation mutArray[], int initDepth, int majorSplitPos[], int numPartitions[], 
			Emat emat, Index3 indexMap[][], int numMutable, MutableResParams strandMut,
			double pruningE, double initEw, RotamerSearch rsP) {

		RotamerSearch rs = rsP; //no changes should be made to the original RotamerSearch object
		rsP = null;

		System.out.print("Computing a lower bound on the conformational energy for each partition..");
		for (int m=0; m<mutArray.length; m++){ //for each partition

			if(m%100 == 0)
				System.out.println("Starting Partition.. "+m+" out of "+mutArray.length+" ");

			PrunedRotamers<Boolean> prunedForPartition = new PrunedRotamers<Boolean>(emat.singles.pruned,false); //no changes should be made to prunedRotAtRes[]
			Iterator<EMatrixEntryWIndex> iter1 = emat.singlesIterator();
			while(iter1.hasNext()){
				EMatrixEntryWIndex emeWI = iter1.next();
				prunedForPartition.set(emeWI.index, emeWI.eme.isPruned());
			}
			//System.arraycopy(prunedRotAtRes, 0, prunedForPartition, 0, prunedRotAtRes.length);

			mutArray[m].setSortScores(true); //sort by lower bounds

			//artificially set to pruned all rotamers at the splitting position, other than the current partitioning rotamer for the given partitioning position
			for (int i=0; i<initDepth; i++){ //first, set all other rotamers for the partitioning positions to pruned
				Index3 curPart = mutArray[m].index[i];
				for (int j=0; j<numPartitions[i]; j++){
					Index3 curInd = indexMap[i][j];					
					if (curInd!=curPart) { //not the rotamer for the current partition
						if (!prunedForPartition.get(curInd)) //rotamer not already pruned
							prunedForPartition.set(curInd, true);
					}
				}
			}

			//compute a lower bound on the conformational energies for this partition
			rs.DoMinBounds(pruningE, initEw, false, false, true);
			mutArray[m].score = new BigDecimal(rs.getBoundForPartition());
		}
		System.out.println("done");

		//sort the partitions
		System.out.print("Sorting partitions by their lower energy bounds..");
		//RyanQuickSort rqs = new RyanQuickSort();
		//rqs.Sort(mutArray);
		//rqs = null;
		Arrays.sort(mutArray);
		System.out.println("done");
		System.out.println("MinBoundOfMinPartition: "+mutArray[0].score);
	}

	///////////////////////////////////////////////////////////////////////////
	//	End of Distributed DACS section
	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////
	//	Distributed DEE section
	// WARNING: the Distributed DEE section is outdated and must be carefully checked before called
	///////////////////////////////////////////////////////////////////////////////////////////////
	//Handles the distribution of the DEE computation to the slave nodes
	//Handles the distribution of the DEE computation to the slave nodes
	private void doDistrDEEMaster(Emat emat, ParamSet sParams, double initEw, boolean doMinimize,  
			boolean minimizeBB, boolean scaleInt, double maxIntScale, double Ival, boolean typeDep, 
			boolean pruningPairs, Settings.DEEMETHOD deeMethod){

		//the total number of residues (active site + ligand, if present)
		int totalNumRes = emat.numMutPos();

		// Generate all combinations
		int numMutRes = 2; //pairs criterion
		int residueMutatable[][];
		if(pruningPairs){
			int numComb = factorial(totalNumRes).divide(factorial(totalNumRes-numMutRes).multiply(factorial(numMutRes))).intValue();
			residueMutatable = new int[numComb][totalNumRes];
			generateCombinations(residueMutatable,totalNumRes,numMutRes);

		}
		else{ //Pruning singles
			residueMutatable = new int[totalNumRes][totalNumRes];
			for(int i=0; i<totalNumRes;i++)
				residueMutatable[i][i] = 1;
		}


		System.out.println("Using this many procs: "+(numProc-1));
		//KER: figure out the number of pairs
		long numTotalPairs = 0;
		long numTotalSingles = 0;
		int[] rotPerPos = emat.numRotPerPos();
		for (int curMut=0; curMut<residueMutatable.length; curMut++){
			int mutPos1 = -1;
			int mutPos2 = -1;

			for (int curRes=0; curRes<totalNumRes; curRes++){
				if(residueMutatable[curMut][curRes] == 1){
					if(mutPos1 == -1)
						mutPos1 = curRes;
					else
						mutPos2 = curRes;
				}
			}
			numTotalSingles += rotPerPos[mutPos1];
			if(pruningPairs)
				numTotalPairs += rotPerPos[mutPos1]*rotPerPos[mutPos2];
		}

		int numToDistribute;
		if(pruningPairs)
			numToDistribute = (int)(numTotalPairs)/((numProc-1)*2);
		else
			numToDistribute = (int)(numTotalSingles)/((numProc-1)*2);

		if(numToDistribute <= 0)
			numToDistribute = 1;

		ArrayList<OneMutation> mutList = new ArrayList<OneMutation>();
		for (int curMut=0; curMut<residueMutatable.length; curMut++){
			int mutPos1 = -1;
			int mutPos2 = -1;

			for (int curRes=0; curRes<totalNumRes; curRes++){
				if(residueMutatable[curMut][curRes] == 1){
					if(mutPos1 == -1)
						mutPos1 = curRes;
					else
						mutPos2 = curRes;
				}
			}


			int numRuns;
			if(pruningPairs){
				int numPairs = rotPerPos[mutPos1]*rotPerPos[mutPos2];
				numRuns = numPairs/numToDistribute;
			}else{//Singles
				numRuns = rotPerPos[mutPos1]/numToDistribute;
			}


			for(int i=0; i<numRuns+1;i++){
				OneMutation oneMut = new OneMutation();
				oneMut.mutNum = curMut;
				oneMut.resMut = new int[totalNumRes];

				for (int curRes=0; curRes<totalNumRes; curRes++){
					oneMut.resMut[curRes] = residueMutatable[curMut][curRes];
				}
				oneMut.pairStartEnd = new int[2];
				oneMut.pairStartEnd[0] = i*numToDistribute;
				oneMut.pairStartEnd[1]   = i*numToDistribute+numToDistribute;

				mutList.add(oneMut);
			}



		}

		OneMutation[] mutArray = mutList.toArray(new OneMutation[0]);

		//Save Emat
		String runName = sParams.getValue("RUNNAME");
		String runNameEMatrixMin = sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");



		MutationManager mutMan = new MutationManager(null,mutArray,false);

		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
//		mutMan.setMutableSpots(totalNumRes);
		//mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.setarpFilenameMin(runNameEMatrixMin);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setInitEw(initEw);
		mutMan.setTypeDep(typeDep);
		//mutMan.setLigPresent(ligPresent);
		mutMan.setDistrDACS(false);
		mutMan.setDistrDEE(true);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setIval(Ival);
		mutMan.setDEEMethod(deeMethod);
		//mutMan.setUseEref(useEref);
		//mutMan.setEref(eRef);
		//mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}


	// Distributed DEE Slave function
	private CommucObj doDistrDEESlave(CommucObj cObj) {

		long startTime = System.currentTimeMillis();

		//Load Emat
		String runName = cObj.params.getValue("RUNNAME");
		String runNameEMatrixMin = cObj.params.getValue("MINENERGYMATRIXNAME",runName+"minM");


		Emat emat = cObj.emat;


		int totalNumRes = emat.numMutPos();


		//determine the two residues in the pair (for pairs DE) or the one residue (split-DEE)
		int[] mutRes = new int[2];
		int firstPos = -1;
		boolean resInMut[] = new boolean[totalNumRes];
		for (int i=0; i<totalNumRes; i++){
			if (cObj.resMut[i]==1){
				resInMut[i] = true;
				if(firstPos == -1){
					firstPos = i;
					mutRes[0] = i;
				}
				else{
					mutRes[1] = i;
				}
			}
			else
				resInMut[i] = false;
		}

		Emat smallEmat = null;
		//if (cObj.typeDEE==optPairs){ //simple Goldstein pairs	
		boolean useFlags = true;
		boolean localUseMinDEEPruningEw = true;
		switch(cObj.deeMethod){
		case GOLDSTEIN:
			RotamerSearch.DoDEEGoldstein(emat,cObj.initEw, useFlags, cObj.typeDep,localUseMinDEEPruningEw,
					cObj.Ival,true,resInMut,cObj.pairStartEnd,false);
			smallEmat = new Emat(emat, true,new int[0]);
			break;
		case SPLITFLAGS1:
			RotamerSearch.DoDEEConfSplitting(emat,cObj.initEw, resInMut, cObj.doMinimization, useFlags, 1, 
					true, cObj.minimizeBB, cObj.typeDep,localUseMinDEEPruningEw,cObj.Ival,cObj.pairStartEnd);
			smallEmat = new Emat(emat, true,new int[0]);
			break;
		case SPLITFLAGS2:
			RotamerSearch.DoDEEConfSplitting(emat,cObj.initEw, resInMut, cObj.doMinimization, useFlags, 2, 
					true, cObj.minimizeBB, cObj.typeDep,localUseMinDEEPruningEw,cObj.Ival,cObj.pairStartEnd);
			smallEmat = new Emat(emat, true,new int[0]);
			break;
		case MBPAIRS:
			RotamerSearch.DoDEEPairs(emat,cObj.initEw, resInMut, cObj.doMinimization, useFlags, true, 
					true, cObj.minimizeBB, cObj.scaleInt, cObj.maxIntScale,cObj.typeDep,localUseMinDEEPruningEw,cObj.Ival,cObj.pairStartEnd);
			smallEmat = new Emat(emat, false,mutRes);
			break;
		case FULLPAIRS:
			RotamerSearch.DoDEEPairs(emat,cObj.initEw, resInMut, cObj.doMinimization,
					cObj.useSF, false, true, cObj.minimizeBB, cObj.scaleInt, cObj.maxIntScale,cObj.typeDep, true, cObj.Ival,cObj.pairStartEnd);
			smallEmat = new Emat(emat, false,mutRes);
			break;
		default:
			System.out.println("I don't recognize the DEE Method");
			break;
		}


		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);

		cObj.emat = smallEmat;

		return cObj;
	}	


	///////////////////////////////////////////////////////////////////////////
	//	End of Distributed DEE section
	///////////////////////////////////////////////////////////////////////////

	static Object readObject(String inFile){
		return readObject(inFile,true);
	}

	static Object readObject(String inFile, boolean repeat){
		Object inObj = null;
		boolean done = false;
		while (!done){
			try{
				ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
				inObj = in.readObject();
				in.close();
				done = true;
			}
			catch (Exception e){
				//System.out.println(e.toString());
				//System.out.println("ERROR: An exception occurred while reading from object file");
				if (repeat)
					done = false;
				else
					done = true;
			}
		}
		return inObj;
	}

	static void outputObject(Object outObj, String outFile){
		try{
			FileOutputStream fout = new FileOutputStream(outFile);
			ObjectOutputStream out = new ObjectOutputStream(fout);
			out.writeObject(outObj);
			out.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing object file");
			System.exit(0);
		}
	}

	private int countPrunedRot(PrunedRotamers<Boolean> prunedRot){
		int countPruned = 0;
		Iterator<RotInfo<Boolean>> iter = prunedRot.iterator();
		while(iter.hasNext()){
			//for (int i=0; i<prunedRot.length; i++){
			RotInfo<Boolean> ri = iter.next();
			if (ri.state)
				countPruned++;
		}
		return countPruned;
	}

	private int countPrunedPairs(boolean prunedPairs[][][][][][]){
		int countPruned = 0;

		boolean[][][][][][] fromMatrix = prunedPairs;
		for (int p1=0; p1<fromMatrix.length; p1++){
			if (fromMatrix[p1]!=null){
				for (int a1=0; a1<fromMatrix[p1].length; a1++){
					if (fromMatrix[p1][a1]!=null){
						for (int r1=0; r1<fromMatrix[p1][a1].length; r1++){
							if (fromMatrix[p1][a1][r1]!=null){
								for (int p2=0; p2<fromMatrix[p1][a1][r1].length; p2++){
									if (fromMatrix[p1][a1][r1][p2]!=null){
										for (int a2=0; a2<fromMatrix[p1][a1][r1][p2].length; a2++){
											if (fromMatrix[p1][a1][r1][p2][a2]!=null){
												for (int r2=0; r2<fromMatrix[p1][a1][r1][p2][a2].length; r2++){
													if(fromMatrix[p1][a1][r1][p2][a2][r2])
														countPruned++;
												}
												//System.arraycopy(fromMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, fromMatrix[p1][a1][r1][p2][a2].length);
											}
										}
									}
								}
							}
						}
					}
				}
			}				
		}
		/*for (int i=0; i<prunedPairs.length; i++){
			for (int j=i+1; j<prunedPairs.length; j++){
				if (prunedPairs[i][j])
					countPruned++;
			}
		}*/
		return countPruned;
	}

	//Function taken from: http://www.javaworld.com/javaworld/javatips/jw-javatip76.html?page=2
	//Java Tip 76: An alternative to the deep copy technique
	//Author: Dave Miller
	static public Object deepCopy(Object oldObj) throws Exception {
		ObjectOutputStream oos = null;
		ObjectInputStream ois = null;
		try
		{
			ByteArrayOutputStream bos = 
					new ByteArrayOutputStream(); // A
			oos = new ObjectOutputStream(bos); // B
			// serialize and pass the object
			oos.writeObject(oldObj);   // C
			oos.flush();               // D
			ByteArrayInputStream bin = 
					new ByteArrayInputStream(bos.toByteArray()); // E
			ois = new ObjectInputStream(bin);                  // F
			// return the new object
			return ois.readObject(); // G
		}
		catch(Exception e)
		{
			System.out.println("Exception in ObjectCloner = " + e);
			throw(e);
		}
		finally
		{
			oos.close();
			ois.close();
		}
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MPI section
	//		The tutorial "One-step Tutorial: MPI: It's easy to get started" 
	//			(http://www.lam-mpi.org/tutorials/one-step/ezstart.php ; accessed Oct 23, 2006) was used as MPI code reference
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup and do MPI
	//	The parameter is only used by the master node
	public void handleDoMPI(String args[]) throws MPIException, InterruptedException {

		MPI.Init(args);

		int procRank = MPItoThread.Rank();
		numProc = MPItoThread.Size();	

		//System.out.println("Node rank: "+procRank+" of "+numProc);
		MPI.COMM_WORLD.Barrier();

		setConfigPars();
		MPI.COMM_WORLD.Barrier();

		if (procRank==0){ //master node
			outputProgInfo();
			parse(args);

		}
		else {//slave node
			for(int i=0; i<mols.length;i++)
				mols[i] = null;

			handleDoMPISlave();
		}

		MPI.Finalize();
	}

	//Do MPI for the master node
	public void handleDoMPIMaster(MutationManager mutMan, int size) throws MPIException, InterruptedException {



		CommucObj cObjArray[] = new CommucObj[size];
		int numFinished = 0;

		int curMut = 0;
		for (int curProc=1; curProc<numProc; curProc++){ //distribute a single mutation per processor, for all processors

			if (curMut<cObjArray.length){ //more mutations to distribute

				System.out.println("Retrieving "+curMut+" of "+(cObjArray.length));
				cObjArray[curMut] = mutMan.getNextComObj(curMut);

				MPItoThread.Send(cObjArray, curMut, 1, ThreadMessage.OBJECT, curProc, regTag);
				curMut++;

				System.out.println("Sent to proc "+curProc);
				System.out.println();
			}
			else
				break;
		}

		boolean distrDACS = mutMan.getDistrDACS(); //distributed DACS computation

		while (numFinished<cObjArray.length){ //distribute and receive all remaining mutations

			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = new CommucObj();

			//System.out.println("Receiving message on main node..");

			Object s = MPItoThread.Recv(cObj, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);

			//System.out.println("Received message on main node: tag "+s.tag+" source "+s.source);

			if (distrDACS){ //DACS computation, so check if the new energy is better than the best energy so far

				double curBestE = mutMan.getPruningE();

				if (cObj[0].pruningE<curBestE){ //the new energy is better

					//System.out.println("Updating best energy on main node from "+curBestE+" to "+cObj[0].pruningE+", source (partition): "+cObj[0].partIndex);

					mutMan.setBestScore(new BigDecimal(cObj[0].pruningE));
					mutMan.setPruningE(cObj[0].pruningE);

					double c[] = new double[1];
					c[0] = cObj[0].pruningE;
					for (int curProc=1; curProc<numProc; curProc++){ //update the best energy in each partition
						MPItoThread.Isend(c, 0, 1, ThreadMessage.FLOAT, curProc, updateTag);
					}
				}
			}

			if (MPItoThread.getStatusTag(s)==regTag){ //completed job
				mutMan.processFinishedMutation(cObj[0]);
				numFinished++;

				System.out.println("Finished: "+cObj[0].mutationNumber+", Time: "+(cObj[0].elapsedTime/60.0));

				if (curMut<cObjArray.length){

					System.out.print("Retrieving "+curMut+" of "+(cObjArray.length));
					cObjArray[curMut] = mutMan.getNextComObj(curMut);

					MPItoThread.Send(cObjArray, curMut, 1, ThreadMessage.OBJECT, MPItoThread.getStatusSource(s), regTag);
					curMut++;

					System.out.println(", Sent to proc "+MPItoThread.getStatusSource(s));
					System.out.println();
				}
			}
		}
	}



	//Do the work of handleDoMPIMaster on the master node
	//the MPI is not set up to let a slave have its own slaves,
	//so this is for if a slave needs to run handleDoMPIMaster (e.g. has to compute an energy matrix)
	public void handleMasterLocally(MutationManager mutMan, int size) {

		if(mutMan.getDistrDACS()){
			System.err.println("ERROR: distrDACS not supported for handleMasterLocally");
			new Exception().printStackTrace();
			System.exit(1);
		}

		System.out.println("Performing slave calculations locally...");
		for(int curMut=0; curMut<size; curMut++){
			CommucObj cObj = mutMan.getNextComObj(curMut);
			cObj = handleKSSlave(cObj);
			mutMan.processFinishedMutation(cObj);
			System.out.println("Finished: "+cObj.mutationNumber+", Time: "+(cObj.elapsedTime/60.0));
		}
	}



	//Do MPI for a slave node
	public void handleDoMPISlave() throws MPIException, InterruptedException {

		int rank = MPItoThread.Rank();

		while (true){

			Object s = MPItoThread.Probe(0, ThreadMessage.ANY_TAG);
			if (s!=null) {
				if (MPItoThread.getStatusTag(s)==regTag){ //new computation			
					CommucObj cObj[] = new CommucObj[1];
					cObj[0] = new CommucObj();

					//System.out.println("node "+rank+" receiving message from main node..");

					MPItoThread.Recv(cObj, 0, 1, ThreadMessage.OBJECT, 0, regTag);

					//System.out.println("node "+rank+" received message from main node..");

					if (cObj[0]==null) //computation is done
						return;

					cObj[0] = handleKSSlave(cObj[0]); //perform computation
					MPItoThread.Send(cObj, 0, 1, ThreadMessage.OBJECT, 0, regTag); //send back result
				}
				else { //(s.tag==updateTag), so discard
					double c[] = new double[1];
					MPItoThread.Recv(c, 0, 1, ThreadMessage.FLOAT, 0, updateTag);
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	 End of MPI section
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	//	Backbone Flexibility Section
	///////////////////////////////////////////////////////////////////////////
	public void generateBackbones(String s){
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Backbone config filename (string)

		System.out.println("Performing Backbone Generation");

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		// Pull search parameters	
		String runName = ((String)sParams.getValue("RUNNAME"));
		boolean sysSampling = (new Boolean((String)sParams.getValue("SYSSAMPLING"))).booleanValue();
		double theta = (new Double((String)sParams.getValue("THETA"))).doubleValue();
		double alpha = (new Double((String)sParams.getValue("ALPHA"))).doubleValue();
		int numSamples = (new Integer((String)sParams.getValue("NUMSAMPLES"))).intValue();
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		boolean neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
		double distCutoff = 0;
		if(neighborList){
			distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
		}

		if (theta%alpha!=0){
			System.out.println("ERROR: Choose theta = k*alpha, for k - an integer.");
			System.exit(1);
		}


		MolParameters mp = loadMolecule(sParams, COMPLEX, neighborList, distCutoff,true);
		//Setup the molecule system
		mp.m = new Molecule();
		mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);

		int numberMutable = mp.strandMut.numMutPos();
		int strandsPresent = mp.strandsPresent;
		MutableResParams strandMut = mp.strandMut;

		RotamerSearch rs = new RotamerSearch(mp.m,numberMutable,strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, 
				doSolvationE, solvScale, softvdwMultiplier, false, null, false, false, false, 
				new EPICSettings(),hbonds, strandMut);

		rs.doGenBackbones(runName, numberMutable, strandMut, theta, alpha, numSamples, sysSampling);
	}
	///////////////////////////////////////////////////////////////////////////
	//	End of Backbone Flexibility Section
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	//	Generate Random Conformations Section
	///////////////////////////////////////////////////////////////////////////

	//Generates a random set of mutations/conformations for a given system
	public void generateRandConfs(String s){
		System.out.println("RANDOM CONFS IS NOT IMPLEMENTED RIGHT NOW");
	}


	//Generates a random set of mutations/conformations for a given system
	/*public void generateRandConfs1(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Run name (for output files)
		// 3: Ligand is present (boolean)
		// 4: Ligand type (if present)
		// 5: Number of conformations to be generated

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters

		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String runName = getToken(s,3);
		boolean ligPresent = (new Boolean(getToken(s,4))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,5);
		int num = (new Integer(getToken(s,6))).intValue();

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();		

		PrintStream logPS = setupOutputFile(runName);

		Random r = new Random();
		for (int i=0; i<num; i++){
			String AAs[] = new String[numInAS];
			int rot[] = new int[numInAS+1];
			for (int a=0; a<numInAS; a++){
				AAs[a] = rl.getAAName(r.nextInt(numAAallowed));
				int n = rl.getNumRotamers(AAs[a]);
				if (n<=0)
					n = 1;
				rot[a] = r.nextInt(n);
			}
			if (ligPresent){
				int n = grl.getNumRotamers(ligType);
				if (n<=0)
					n = 1;
				rot[numInAS] = r.nextInt(n);
			}

			logPS.print(i+" ");
			for (int a=0; a<numInAS; a++)
				logPS.print(AAs[a]+" ");
			if (ligPresent)
				logPS.print(ligType+" ");

			for (int a=0; a<numInAS; a++)
				logPS.print(rot[a]+" ");
			if (ligPresent)
				logPS.print(rot[numInAS]+" ");

			logPS.println();
			logPS.flush();
		}
		logPS.close();
	}*/
	///////////////////////////////////////////////////////////////////////////
	//	End of Generate Random Conformations Section
	///////////////////////////////////////////////////////////////////////////

	//KER: It is better to debug a threaded version of doDEE
	//than to continue to keep this function up to date.
//	public void doSinglePairE(String s, ParamSet sParams) {
//
//		// Takes the following parameters
//		// 1: System parameter filename (string)
//		// 2: Mutation search parameter filename (string)		
//
//		//Only read system and mutation files if sParams is null
//
//		if (sParams==null){ //parameter files not read yet, so read them in
//			sParams = new ParamSet();
//			sParams.addParamsFromFile(getToken(s,2)); //read system parameters
//			sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
//		}
//
//		//int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
//		int numMutations = 2; //pairwise energies are computed
//		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
//		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
//		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
//		String runName = (String)sParams.getValue("RUNNAME");
//		String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
//		String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM");
//		boolean templateAlwaysOn = (new Boolean((String)sParams.getValue("TEMPLATEALWAYSON"))).booleanValue();
//		double distCutoff=0;
//		boolean neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
//		if(neighborList){
//			distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
//		}
//
//		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
//		boolean inputSysWithLig = ligPresent;
//		String ligType = null;
//		if (ligPresent)
//			ligType = (String)sParams.getValue("LIGTYPE");
//
//		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH"))).booleanValue();
//		String resumeFilename = (String)sParams.getValue("RESUMEFILENAME");
//
//
//		if( (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue() ){
//			System.err.println("ERROR: DEEPer not supported in doSinglePairE");
//			System.exit(1);
//		}
//
//		System.out.println("Run Name: "+runName);
//		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
//		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
//		System.out.println("Ligand Type: "+ligType);
//		System.out.println("Num Residues Allowed to Mutate: "+numMutations);
//
//
//		System.out.println("Computing _All_ Rotamer-Rotamer Energies");
//
//		System.out.println("Starting minimum and maximum bound energy computation");
//
//		if(resumeSearch) {
//			System.out.println("** Resuming Search **");
//			System.out.println("     resuming from file: "+resumeFilename);
//		}
//
//		MolParameters mp = loadMolecule(sParams, COMPLEX, neighborList, distCutoff,true);
//		Molecule m = mp.m;
//		int numberMutable = mp.numberMutable;
//		int strandsPresent = mp.strandsPresent;
//		String[][] strandLimits = mp.strandLimits;
//		boolean[] strandPresent = mp.strandPresent;
//		MutableResParams strandMut = mp.strandMut;
//		String[][] strandDefault = mp.strandDefault;
//
//		RotamerSearch rs = new RotamerSearch(m,numberMutable,strandsPresent, hElect, hVDW, hSteric, true,
//				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst,
//				doDihedE, doSolvationE, solvScale, softvdwMultiplier, grl,
//				false, null, false, false, false, new EPICSettings(),hbonds, strandMut);
//
//
//
//		int resMut[] = new int[numberMutable];
//		for (int i=0; i<resMut.length; i++)
//			resMut[i] = 0;
//
//		String flagMutType = "AS-AS";
//		resMut[0] = 1; resMut[1] = 1;
//
//		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","true"))).booleanValue();
//		if(!addWT)
//			mp.strandMut.checkWT(mp.strandPresent, sParams);
//		int molStrand = 0;
//		for(int resID:mp.strandMut.allMut){
//			if(mp.strandPresent[mp.m.residue[resID].strandNumber]){
//				setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
//			}
//		}
//
//
//		boolean shellRun = false;
//		boolean intraRun = false;
//		boolean templateOnly = false;		
//
//		if (flagMutType.compareTo("TEMPL")==0){
//
//			// **** Normal Active Site residue runs ****
//			// Computes active site residue to active site residue pair energies					
//			shellRun = true;ligPresent = true;intraRun = false;templateOnly = true;
//		}
//		else if (flagMutType.compareTo("AS-AS")==0){
//
//			// **** Normal Active Site residue runs ****
//			// Computes active site residue to active site residue pair energies					
//			shellRun = false;ligPresent = true;intraRun = false;templateOnly = false;
//		}	
//		else if (flagMutType.compareTo("SHL-AS")==0){
//
//			// Then shell runs for the active site residues
//			// Computes the active site residue rotamers to shell energies					
//			shellRun = true;ligPresent = true;intraRun = false;templateOnly = false;				
//		}
//		else if (flagMutType.compareTo("INTRA")==0){
//
//			// Compute all intra-residue energies					
//			shellRun = false;intraRun = true;templateOnly = false;
//		}				
//		else if (flagMutType.compareTo("LIG-AS")==0){
//
//			// **** Ligand present runs ****
//			// This section computes the inter-residue energies between
//			//  active site residues and the ligand
//			shellRun = false;intraRun = false;templateOnly = false;			
//		}
//		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)
//
//			// Computes ligand rotamer to shell energies
//			shellRun = true; intraRun = false;templateOnly = false;
//		}
//
//		// The goal is that the total energy of a system can be bounded by the sum of 
//		//  all pairwise active site residue entries plus the entry for each active site
//		//  residue's shell run plus each active site residue's self intra-residue energy.
//		//  If a ligand is present then one should add the ligand to shell energy, the
//		//  ligand to each active site residue pairwise energy, and the ligand self intra-
//		//  residue energy.
//
//		//initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
//		PairwiseEnergyMatrix minEmatrix = new PairwiseEnergyMatrix(numberMutable,resMut,strandMut,
//				rs,shellRun,intraRun,false);
//		//		PairwiseEnergyMatrix maxEmatrix = minEmatrix.copy();
//
//		//Compute the corresponding matrix entries
//		rs.simplePairwiseMutationAllRotamerSearch(strandMut,numberMutable,true,shellRun,intraRun,
//				resMut,minEmatrix,minimizeBB,doBackrubs,templateOnly,backrubFile, templateAlwaysOn,strandDefault, false);
//
//		return;
//	}

	//Computes conformation energies for different combinations of the energy function parameters
	private void fitEparams(String s){

		String firstParam = getToken(s,1);

		String sysFile = getToken(s,2);

		// Pull search parameters
		String confResFile = getToken(s,3);
		String runName = getToken(s,4);
		boolean ligPresent = (new Boolean(getToken(s,5))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,6);
		int numResults = (new Integer(getToken(s,7))).intValue();
		boolean minimizeBB = (new Boolean(getToken(s,8))).booleanValue();

		int numSteps = 10;
		double maxVdwMult = 1.05;
		double maxSolvScale = 0.3;

		solvScale = 0.0;
		double initVdwMult = 0.63;

		double vdwDelta = (maxVdwMult-initVdwMult)/numSteps;
		double solvDelta = (maxSolvScale-solvScale)/numSteps;

		for (int i1=0; i1<numSteps; i1++){
			solvScale += solvDelta;
			for (int i2=0; i2<2; i2++){
				if (i2==0)
					distDepDielect = true;
				else
					distDepDielect = false;

				for (int i3=0; i3<=numSteps; i3++){
					if (i3==0)
						dielectConst = 1.0;
					else
						dielectConst = 4*i3;

					softvdwMultiplier = initVdwMult-vdwDelta;
					for (int i4=0; i4<=numSteps; i4++){
						softvdwMultiplier += vdwDelta;

						String runNameParams = (runName+"_"+solvScale+"_"+distDepDielect+"_"+dielectConst+"_"+softvdwMultiplier);

						String s1 = (firstParam+" "+sysFile+" "+confResFile+" "+runNameParams+" false none "+numResults+" "+minimizeBB);
						handleMinDEEApplyRot(s1);

						if (ligPresent){
							runNameParams = (runNameParams+"_lig");
							s1 = (firstParam+" "+sysFile+" "+(confResFile+"_lig")+" "+runNameParams+" true "+ligType+" "+numResults+" "+minimizeBB);
							handleMinDEEApplyRot(s1);
						}
					}
				}
			}
		}
	}


	////////////////////////////////////////////////////////////////
	// Compute Residue Entropy Section
	////////////////////////////////////////////////////////////////
	/**
	 * Handles the SCMF residue entropy computation.
	 */
	public void handleDoResEntropy(String s, ParamSet sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		

		//Only read system and mutation files if sParams is null

		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new ParamSet();
			sParams.addParamsFromFile(getToken(s,2)); //read system parameters
			sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		}

		String runName = (String)sParams.getValue("RUNNAME");

		String matrixName = (String)sParams.getValue("MATRIXNAME");
		boolean useEref = (new Boolean((String)sParams.getValue("USEEREF", "true"))).booleanValue();
		double dist = (new Double((String)sParams.getValue("DIST", "8.0"))).doubleValue();
		String rotProbFile = (String)sParams.getValue("ROTPROBFILE", runName+"rotProb");
		double stericE = (new Double((String)sParams.getValue("STERICE", "30"))).doubleValue();
		double maxPairE = (new Double((String)sParams.getValue("MAXPAIRE","1000.0"))).doubleValue();

		MolParameters mp = new MolParameters();
		loadStrandParams(sParams, mp, COMPLEX);

		//Set nonprotein strands as not present
		for(int i=0; i<mp.numOfStrands;i++){
			boolean isProtein = (new Boolean((String)sParams.getValue("STRANDAA"+i))).booleanValue();
			if(!isProtein){
				mp.strandPresent[i] = false;
			}
		}

		//Setup the molecule system
		mp.m = new Molecule();
		mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits,false);
		Molecule m = mp.m;



		/*MolParameters mp = loadMolecule(sParams, COMPLEX);
		Molecule m = mp.m;
		int numberMutable = mp.numberMutable;
		int strandsPresent = mp.strandsPresent;
		String[][] strandLimits = mp.strandLimits;
		boolean[] strandPresent = mp.strandPresent;
		int[][] strandMut = mp.strandMut;
		String[][] strandDefault = mp.strandDefault;
		int[] mutRes2Strand = mp.mutRes2Strand;
		int [] mutRes2StrandMutIndex = mp.mutRes2StrandMutIndex;*/
		String[][] strandDefault = new String [mp.m.strand.length][];
		for(int str=0; str<strandDefault.length;str++){
			strandDefault[str] = new String[m.strand[str].numberOfResidues];
			for(int i=0; i<m.strand[str].numberOfResidues;i++){
				strandDefault[str][i] = m.strand[str].residue[i].name;
			}
		}

		int numRes = m.numberOfResidues;

		//KER: this function assumes we're only looking at one strand
		//KER: for now I'm assuming that's a protein strand
		RotamerLibrary rotLib = m.aaRotLib;


		rotProbFile = (rotProbFile+".dat");
		double rotProb[][] = (double [][])readObject(rotProbFile,false);
		if (rotProb==null){ //Perform SCMF to compute the rotamer probabilities		

			//read in or compute all of the energy matrices;
			//NOTE: backbone-to-backbone and rotamer-to-backbone energies are included in the pairwise rotamer energies
			double asasE[][][][] = getResEntropyEmatricesPair(matrixName, sParams, numRes, strandDefault, mp, runName, dist);
			double intraEnergies[][] = getResEntropyEmatricesIntra(matrixName, sParams, numRes, mp, runName);
			double eRef[] = getResEntropyEmatricesEref(useEref, null, null, null, intraEnergies, numRes,null,null,rotLib);

			rotProb = compRotProbSCMF(numRes, intraEnergies, asasE, eRef, rotProbFile, strandDefault, stericE, maxPairE,m,rotLib);
		}

		int numProx[] = new int[numRes]; //get the number of proximate residues for each residue position
		for (int i=0; i<numProx.length; i++)
			numProx[i] = 0;
		String asDistFile = matrixName+"_dist.dat";
		boolean as[][] = (boolean [][])readObject(asDistFile,false);
		for (int i=0; i<numRes; i++){
			for (int j=i+1; j<numRes; j++){
				if (as[i][j]){
					numProx[i]++;
					numProx[j]++;
				}
			}
		}

		//m = null;		

		PrintStream logPS = setupOutputFile(runName);

		logPS.print("resNum pdbResNum resDefault entropy"+" ");
		for (int j=0; j<rotLib.getAAtypesAllowed().length; j++)
			logPS.print(rotLib.getAAtypesAllowed()[j]+" ");
		logPS.println("numProx");
		logPS.flush();


		//Compute the AA probabilities for each residue position (as a function of the rotamer probabilities);
		//Compute the entropy at each position as a function of the amino acid probabilities for that position
		final double kB = 1.0;
		for (int i=0; i<numRes; i++){
			int str=m.residue[i].strandNumber;
			int strResNum=m.residue[i].strandResidueNumber;
			if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ){

				double aaProbBBE[] = new double[rotLib.getAAtypesAllowed().length];
				int curInd = 0;
				for (int j=0; j<aaProbBBE.length; j++){

					int numCurRot = rotLib.getAAtypesAllowed()[j].numRotamers();
					if (numCurRot==0)
						numCurRot = 1;

					aaProbBBE[j] = 0.0f;
					for (int r=0; r<numCurRot; r++){
						aaProbBBE[j] += rotProb[i][curInd];
						curInd++;
					}
				}

				//Compute the unnormalized AA probabilities as a weighted function of energies and PDB stats (if included)
				double aaProbUnNorm[] = new double[rotLib.getAAtypesAllowed().length];
				double aaNorm = 0.0;
				for (int j=0; j<rotLib.getAAtypesAllowed().length; j++){

					aaProbUnNorm[j] = 1.0;
					aaProbUnNorm[j] *= aaProbBBE[j];

					aaNorm += aaProbUnNorm[j];
				}

				//Normalize the probabilities
				double aaProb[] = new double[rotLib.getAAtypesAllowed().length];
				for (int j=0; j<rotLib.getAAtypesAllowed().length; j++){
					if (aaNorm!=0.0)
						aaProb[j] = aaProbUnNorm[j]/aaNorm;
					else
						aaProb[j] = 0.0;
				}

				//Compute the entropy for the current residue position
				double sumAA = 0.0;
				for (int j=0; j<aaProb.length; j++){
					if (aaProb[j]>0.0)
						sumAA += aaProb[j]*Math.log(aaProb[j]);
				}

				double entropy = -kB * sumAA;


				logPS.print(i+" "/*+defResNum[str][strResNum]*/+" "+strandDefault[str][strResNum]+" "+entropy+" ");
				for (int j=0; j<aaProb.length; j++)
					logPS.print(aaProb[j]+" ");

				logPS.println(numProx[i]);
				logPS.flush();
			}
			else {
				logPS.println(i+" "/*+defResNum[str][strResNum]*/+" "+strandDefault[str][strResNum]+" "+0.0); //only for residue positions with wildtype Pro
			}
		}
		logPS.close();
	}

	//Computes the rotamer probabilities for all rotamers at all residue positions using SCMF
	private double[][] compRotProbSCMF(int numRes, double intraEnergies[][],
			double asasE[][][][], double eRef[], String rotProbFile, String strandDefault[][], double stericE, double maxPairE,Molecule m,
			RotamerLibrary rotLib){

		final double constR = (double)(1.9891/1000.0);//the gas constant
		double T = 50000; //initial temperature
		final double endT = 298.15f; //the minimum temperature for annealing
		double tStepSize = 100.0f; //the temperature step size for annealing
		final double eps = 0.0001f; //the convergence threshold
		final double lambda = 0.5f; //scaling factor for updating the rotamer probabilities

		int totalNumRotamers = rotLib.getTotalNumRotamers();

		for (int i=0; i<asasE.length; i++){//Set the max energy for any element in asasE[][][][] to maxPairE
			if (asasE[i]!=null){
				for (int j=0; j<asasE[i].length; j++){
					if (asasE[i][j]!=null){
						for (int k=0; k<asasE[i][j].length; k++){
							if (asasE[i][j][k]!=null){
								for (int l=0; l<asasE[i][j][k].length; l++){
									if (asasE[i][j][k][l]>maxPairE)
										asasE[i][j][k][l] = maxPairE;
								}
							}
						}
					}
				}
			}
		}

		int numPrunedRot = 0;
		boolean prunedRot[][] = new boolean[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++){
				if ( (intraEnergies[1+i*totalNumRotamers+j][0]) > stericE){
					prunedRot[i][j] = true;
					numPrunedRot++;
				}
				else
					prunedRot[i][j] = false;
			}
		}
		System.out.println("Num rotamers pruned due to incompatibility with the template: "+numPrunedRot);


		//For each residue, compute the probability of each rotamer for that residue
		double Emf[][] = new double[numRes][totalNumRotamers];
		double rotProb[][] = new double[numRes][totalNumRotamers];
		double oldProb[][] = new double[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++) {
				if (!prunedRot[i][j])
					rotProb[i][j] = 1.0f/totalNumRotamers;
				else
					rotProb[i][j] = 0.0f;

				oldProb[i][j] = rotProb[i][j];
			}
		}

		while (T>=endT){ //perform annealing

			System.out.println("Starting run at T = "+T);

			boolean done = false;
			while (!done){

				//Compute the new mean-field energy for each rotamer
				for (int i=0; i<numRes; i++){

					if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ){

						for (int j=0; j<totalNumRotamers; j++){

							if (!prunedRot[i][j]){

								Emf[i][j] = intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j,rotLib)];

								if (asasE[i]!=null){

									for (int k=0; k<numRes; k++){
										if ( (k!=i) && (!strandDefault[m.residue[k].strandNumber][m.residue[k].strandResidueNumber].equalsIgnoreCase("PRO")) ){ //for all residues with which i_j has contact
											if ( (i<k) && (asasE[i][k]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[i][k][j][l]*rotProb[k][l];
												}
											}
											else if ( (i>k) && (asasE[k][i]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[k][i][l][j]*rotProb[k][l];
												}
											}
										}
									}
								}
							}
						}
					}
				}

				//Update the rotamer probabilities			
				for (int i=0; i<numRes; i++){

					if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ){

						double normFactor = 0.0f;
						for (int j=0; j<totalNumRotamers; j++){

							if (!prunedRot[i][j])
								normFactor += (double)Math.exp( -Emf[i][j] / (constR*T));	
						}

						for (int j=0; j<totalNumRotamers; j++){

							if (!prunedRot[i][j]){

								oldProb[i][j] = rotProb[i][j]; //the probability before the update

								if (normFactor!=0.0f)
									rotProb[i][j] = lambda*((double)Math.exp( -Emf[i][j] / (constR*T)) / normFactor) + (1-lambda)*oldProb[i][j];
								else
									rotProb[i][j] = 0.0f;
							}
						}
					}
				}

				double rms = checkRotProbConvEntropy(rotProb,oldProb,strandDefault,prunedRot,m);

				if (rms>eps)
					done = false;
				else
					done = true;
			}

			T -= tStepSize;
		}

		outputObject(rotProb,rotProbFile);

		return rotProb;
	}

	//Checks if the rotamer probabilities for the entropy computation have converged
	private double checkRotProbConvEntropy(double rotProb[][], double oldProb[][], String strandDefault[][], boolean prunedRot[][],Molecule m){

		double sum = 0.0f;
		for (int i=0; i<rotProb.length; i++){
			if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ){
				for (int j=0; j<rotProb[i].length; j++){
					if (!prunedRot[i][j])
						sum += (double)Math.pow( (rotProb[i][j]-oldProb[i][j]) , 2.0);
				}
			}
		}
		double rms = (double)Math.sqrt(sum);

		System.out.println("RMS: "+rms);

		return rms;
	}

	//Reads in (if computed) or computes the energy matrices for rot-to-rot pairwise energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private double [][][][] getResEntropyEmatricesPair(String matrixName, ParamSet sParams, int numRes, String strandDefault[][], MolParameters mp,
			String runName, double dist){

		String origPDB = (String)sParams.getValue("PDBNAME");

		//Check for the pairwise rot-to-rot file;
		String pairName = matrixName+"_pair.dat";

		double asasE[][][][] = readPairMatrixEntropy(pairName,numRes);
		if (asasE==null){ //compute the rot-to-template energy matrix

			//For each residue position i, get all residue positions that are within dist
			String asDistFile = matrixName+"_dist.dat";
			boolean as[][] = (boolean [][])readObject(asDistFile,false);
			if (as==null){
				as = new boolean[numRes][numRes];
				computeEntropyEmatrixMaster(numRes,runName+".log",asDistFile,sParams,false,false,false,null,
						null,null,true,dist,as,false, mp); //the distances are returned in as[][]
			}

			int numPairs = 0;
			asasE = new double[numRes][numRes][][];

			for (int i=0; i<numRes; i++){

				if (!strandDefault[mp.m.residue[i].strandNumber][mp.m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO")){
					for (int j=i+1; j<numRes; j++){
						if (!strandDefault[mp.m.residue[j].strandNumber][mp.m.residue[j].strandResidueNumber].equalsIgnoreCase("PRO")){
							if (as[i][j]){
								asasE[i][j] = new double[1][];
								numPairs++;
							}
						}
					}
				}
			}

			computeEntropyEmatrixMaster(numRes,runName+".log",pairName,sParams,false,false,false,null,
					null,asasE,false,0.0f,null,false,mp);

			sParams.setValue("PDBNAME", origPDB);



			asasE = readPairMatrixEntropy(pairName,numRes);
		}		

		return asasE;
	}

	//Reads in (if computed) or computes the energy matrices for intra-rot energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private double [][] getResEntropyEmatricesIntra(String matrixName, ParamSet sParams, int numRes, MolParameters mp, String runName){

		String origPDB = (String)sParams.getValue("PDBNAME");


		//Check for the intra energies file
		String intraName = matrixName+"_intra.dat";

		double intraEnergies[][] = (double [][])readObject(intraName,false); //check if already computed			

		if (intraEnergies==null) //compute the intra-rotamer energy matrix
			computeEntropyEmatrixMaster(numRes,runName+".log",intraName,sParams,false,false,false,null,
					intraEnergies,null,false,0.0f,null,true,mp); 

		sParams.setValue("PDBNAME", origPDB);


		/*if (ligStrNum>=0) //the ligand is not used here
			m.deleteStrand(ligStrNum);*/

		intraEnergies = (double [][])readObject(intraName,false);

		return intraEnergies;
	}

	//Reads in the amino acid reference energies (if used);
	private double [] getResEntropyEmatricesEref(boolean useEref, double intraEnergies[][][][][][], StrandRotamers strandRot[], int strandMut[][], double intraEnergiesEntropy[][], int numRes, 
			int mutRes2Strand[],int mutRes2StrandMutIndex[], RotamerLibrary rotLib){

		if ( (intraEnergies==null && intraEnergiesEntropy==null) || (intraEnergies!=null && intraEnergiesEntropy!=null) ){ //exactly one of the two matrices should be non-null
			System.out.println("ERROR: exactly one matrix can be used for the reference energy computation.");
			System.exit(1);
		}

		double eRef[] = new double[rotLib.getAAtypesAllowed().length];

		if (useEref){ //use AA reference energies		
			//TODO: Re-implement this function because we don't make a seperate intraEnergy matrix anymore
			System.out.println("Reference energies aren't implemented for this function");
			System.exit(0);
//			eRef = compEref(intraEnergies,strandRot,strandMut,intraEnergiesEntropy,numRes,mutRes2Strand, mutRes2StrandMutIndex,rotLib);
		}
		else {
			for(int j=0; j<eRef.length;j++)
				//for (int i=0; i<eRef[j].length; i++)
				eRef[j] = 0.0f;
		}

		return eRef;
	}

	//Reads in the pairwise energy matrix for the entropy computation
	private double [][][][] readPairMatrixEntropy(String fName, int numRes){

		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			return null;
		}

		double asasE[][][][] = new double[numRes][][][];

		boolean done = false;
		String str = null;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {				
				int i = new Integer(getToken(str,1)).intValue();
				String name = getToken(str,2);
				asasE[i] = (double [][][])readObject(name,false);
				if (asasE[i]==null){
					System.out.println("ERROR: Could not read data from file "+name);
					System.exit(1);
				}
			}
		}

		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		return asasE;
	}

	//Determines how many residue positions in the molecule numbered from (pos+1...) (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(Molecule m, int pos, double dist, boolean as[]){

		int strandMap[] = new int[2];
		strandMap[0] = pos;

		for (int i=pos+1; i<m.numberOfResidues; i++){

			if (i!=pos){

				boolean done = false;

				strandMap[1] = i;

				Molecule m1 = getASASMolEntropy(m, strandMap);

				StrandRotamers strandRot = new StrandRotamers(m.aaRotLib,m1.strand[0]);

				
				
				for (int j=0; j<m1.numberOfResidues; j++){
					for (int r=0; r<strandRot.rl.getAAtypesAllowed().length; r++){
						m1.residue[j].setAllowable(strandRot.rl.getAAtypesAllowed()[r].name);
						m1.residue[j].flexible = true;
					}
				}

				for(AARotamerType aaType1 : m.residue[0].AATypesAllowed()) {

//					int AAindex1 = strandRot.getIndexOfNthAllowable(0,q1);
					strandRot.changeResidueType(m1,0,aaType1.name,true,true);
//					AARotamerType aaType1 = strandRot.rl.getAAType(AAindex1);
					
					for(AARotamerType aaType2 : m.residue[1].AATypesAllowed()) {

//						int AAindex2 = strandRot.getIndexOfNthAllowable(1,q2);
						strandRot.changeResidueType(m1,1,aaType2.name,true,true);
//						AARotamerType aaType2 = strandRot.rl.getAAType(AAindex2);
						
						int numRot1 = aaType1.numRotamers();
						
						int w1 = 0;
						if (numRot1<=0)
							w1 = -1;

						for (Rotamer rot1: aaType1.rotamers){

							if (w1!=-1)
								strandRot.applyRotamer(m1, 0, rot1);

							int numRot2 = aaType2.numRotamers();

							int w2 = 0;
							if (numRot2<=0)
								w2= -1;

							for (Rotamer rot2: aaType2.rotamers){

								if (w2!=-1)
									strandRot.applyRotamer(m1, 1, rot2);

								Residue r1 = m1.residue[0];
								Residue r2 = m1.residue[1];

								if (r1.getDist(r2,true)<=dist){
									as[i] = true;
									done = true;
								}

								w2++;
								if(done)
									break;
							}

							w1++;
							if(done)
								break;
						}
						if (done)
							break;
					}
					if (done)
						break;
				}
				if (!done)
					as[i] = false;
			}
		}

		return as;
	}

	//Returns the AA index into rotamerIndexOffset to which rotNum belongs
	private int getAAindFromRotNum(int rotNum, RotamerLibrary rotLib){
		Rotamer r = rotLib.getRot(rotNum);
		return r.aaType.index;
		
//		int[] rotamerIndexOffset = rotLib.getRotamerIndexOffset();
//		for (int i=0; i<rotamerIndexOffset.length-1; i++){
//			if ( (rotNum>=rotamerIndexOffset[i]) && (rotNum<rotamerIndexOffset[i+1]) )
//				return i;
//		}
//		if (!(rotNum>=rotLib.getTotalNumRotamers()))
//			return (rotamerIndexOffset.length-1);
//		else
//			return -1;
	}

	//Distributes the different types of energy computation for the entropy calculation
	private void computeEntropyEmatrixMaster(int numRes, String runName, String matrixName, ParamSet sParams, boolean doMinimize, boolean minimizeBB, boolean doBackrubs, String backrubFile,
			double bbEnergies[][], double asasE[][][][], boolean compASASdist, double dist, boolean asDist[][], boolean intraRun, MolParameters mp){

		//MolParameters mp = loadMolecule(sParams, COMPLEX);
		//KER: This function assumes using the protein rotamer library
		RotamerLibrary rotLib = mp.m.aaRotLib;

		int mutEnerMatrixSize = 0;
		int residueMap[] = null;
		OneMutation mutArray[] = null;
		int numMutable = numRes;

		int numMut = 0;

		if (compASASdist){ //compute the min distance between any pair of rotamers for each pair of residue positions
			mutArray = new OneMutation[numMutable];
			for (int i=0; i<asDist.length; i++)			
				mutArray[i] = new OneMutation();

			numMut = numRes;
		}
		else {
			if (intraRun) { //computing intra energies

				System.out.println("Starting intra-rot energy computation..");

				mutEnerMatrixSize = 1 + rotLib.getTotalNumRotamers()*numMutable;

				bbEnergies = new double[mutEnerMatrixSize][1];
				for(int i=0; i<mutEnerMatrixSize; i++) {
					for(int j=0; j<1; j++){
						bbEnergies[i][j] = 0.0f;
					}
				}

				numMutable = 1;
				residueMap = new int[numMutable];

				mutArray = new OneMutation[numRes];
				for (int i=0; i<mutArray.length; i++){
					mutArray[i] = new OneMutation();
					mutArray[i].flagMutType = "INTRA";
				}

				numMut = numRes;
			}

			else { //AS-AS energies

				System.out.println("Starting rot-to-rot energy computation..");

				numMutable = 2;

				int numPairs = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null)
							numPairs++;
					}
				}
				mutArray = new OneMutation[numPairs];

				int curPair = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null){				
							mutArray[curPair] = new OneMutation();
							mutArray[curPair].flagMutType = "AS-AS";
							mutArray[curPair].resMut = new int[2];
							mutArray[curPair].resMut[0] = i; //molecule-relative numbering
							mutArray[curPair].resMut[1] = j;
							curPair++;

							asasE[i][j] = new double[rotLib.getTotalNumRotamers()][rotLib.getTotalNumRotamers()];
						}
					}
				}

				numMut = numPairs;
			}
		}


		MutationManager mutMan = new MutationManager(runName,mutArray,true);
		//mutMan.setStrandMut(strandMut);
		mutMan.setStrandLimits(mp.strandLimits);
		mutMan.setStrandPresent(mp.strandPresent);
		mutMan.setStrandsPresent(mp.strandsPresent);
		//mutMan.setLigType(null);
		mutMan.setarpFilenameMin(matrixName);
		mutMan.setIntraEntropyMatrixMin(bbEnergies);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setRotamerLibrary(rotLib);
//		mutMan.setMutableSpots(numMutable);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setLigPresent(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setEntropyComp(true);
		mutMan.setPairEntropyMatrix(asasE);
		mutMan.setASdistMatrix(asDist);
		mutMan.setASdist(dist);
		mutMan.setCompASdist(compASASdist);
		mutMan.setRotamerLibrary(mp.m.aaRotLib);


		try{
			handleDoMPIMaster(mutMan,numMut);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}

		if (compASASdist){
			asDist = mutMan.getASdistMatrix();
			outputObject(asDist,matrixName);
		}
		else {
			if (intraRun){
				bbEnergies = mutMan.getMinEmatrixEntropy();

				int numCompEntries = 0;				
				for(int i1=0; i1<mutEnerMatrixSize; i1++){
					if ((bbEnergies[i1][0]!=0.0f))
						numCompEntries++;
				}
				System.out.println("Num computed entries: "+numCompEntries);

				outputObject(bbEnergies,matrixName);
			}
			else {
				asasE = mutMan.getPairEntropyEmatrix();		

				PrintStream logPS = setupOutputFile(matrixName);
				for (int i=0; i<asasE.length; i++){
					if (asasE[i]!=null){
						String fn = ("peme/pem_entr_"+i);
						logPS.println(i+" "+fn);
						outputObject(asasE[i],fn);
						logPS.flush();
					}
				}
				logPS.close();
			}
		}
	}

	private CommucObj handleDoResEntropySlave(CommucObj cObj){
		//TODO: understand what this function does and fix it for multiple strands
		long startTime = System.currentTimeMillis();

		//Set nonprotein strands as not present
		for(int i=0; i<cObj.strandPresent.length;i++){
			boolean isProtein = (new Boolean((String)cObj.params.getValue("STRANDAA"+i))).booleanValue();
			if(!isProtein){
				cObj.strandPresent[i] = false;
			}
		}
		//Setup the molecule system
		Molecule m = new Molecule();
		//BAD CODE setupMolSystem(m,cObj.params,false,null); //the ligand is not used here
		m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits, false);

		if (cObj.compASdist){ //AS-AS distance computation
			cObj.asDist = getProxAS(m,cObj.mutationNumber,cObj.dist,cObj.asDist);
			cObj.compEE = null;//new SamplingEEntries[0];
		}
		else { //AS-AS or INTRA energy computation

			boolean shellRun = false; boolean ligPresent = false; boolean intraRun = false; boolean templateOnly = false;
			MutableResParams strandMut = null;
			int resMut[] = null;
			if (cObj.flagMutType.compareTo("AS-AS")==0){ //AS-AS run

				int numMutable = 2;
				resMut = new int[numMutable];
				strandMut = new MutableResParams(numMutable, 1);
					
				for (int i=0; i<numMutable; i++)
					resMut[i] = 1;

				m = getASASMolEntropy(m,cObj.strandMut.allMut);

				strandMut.allMut[0] = 0;
				strandMut.allMut[1] = 1;

			}

			else if (cObj.flagMutType.compareTo("INTRA")==0){

				int numMutable = 1;
				resMut = new int[numMutable];
				strandMut = new MutableResParams(numMutable, 1);
				
				intraRun = true;

				m = getASASMolEntropy(m,cObj.strandMut.allMut);

				strandMut.allMut = new int[1];
				strandMut.allMut[0] = 0;

				resMut = new int[1];
				resMut[0] = 1;
			}

			else {
				System.out.println("ERROR: only AS-AS and INTRA runs allowed for the pairwise entropy matrix precomputation.");
				System.exit(1);
			}
			
			//KER: there will only be one strand with the remade molecule
			int strandsPresent = 1;

			RotamerSearch rs = new RotamerSearch(m, strandMut.allMut.length, strandsPresent, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, 
					cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult, 
					false, null, false, false, false, new EPICSettings(),hbonds,strandMut);

			for(int j=0; j<strandMut.allMut.length; j++) {
				for(int q=0;q<rs.strandRot[0].rl.getAAtypesAllowed().length;q++)
					rs.setAllowable(strandMut.allMut[j],rs.strandRot[0].rl.getAAtypesAllowed()[q].name,0);
			}

			//initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
			Emat minEmatrix = new Emat(m,strandMut,false,false);
			rs.setMinMatrix(minEmatrix);
//			PairwiseEnergyMatrix minEmatrix = new PairwiseEnergyMatrix(numMutable,resMut,strandMut,rs,shellRun,intraRun,false);
//			PairwiseEnergyMatrix maxEmatrix = minEmatrix.copy();

			rs.simplePairwiseMutationAllRotamerSearch(strandMut,strandMut.allMut.length,cObj.doMinimization,shellRun,intraRun,
					resMut,cObj.minimizeBB,cObj.doBackrubs,templateOnly,cObj.backrubFile, RotamerSearch.MINIMIZATIONSCHEME.PAIRWISE, cObj.runParams, false);

			long stopTime = System.currentTimeMillis();
			cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);

			//Store the information in less space to allow the master node to buffer several cObj at once
//			cObj.compEE = minEmatrix.generateCompEE(maxEmatrix);
		}

		return cObj;
	}

	//Returns a molecule that contains only the residues in the system strand (sysStrNum) of molecule m that are specified by residueMap[];
	//	This function is used for the pairwise energy matrix computation in the residue entropy calculations
	private Molecule getASASMolEntropy (Molecule m, int strandMap[]){

		Molecule m1 = new Molecule();

		for (int i=0; i<strandMap.length; i++){

			Residue oldResidue = m.residue[strandMap[i]];

			Residue newResidue = new Residue();
			newResidue.name = oldResidue.name;
			newResidue.fullName = oldResidue.fullName;

			for (int j=0; j<oldResidue.numberOfAtoms; j++){

				Atom oldAtom = oldResidue.atom[j];

				Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
				newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
				newAtom.strandNumber = oldAtom.strandNumber;
				newAtom.elementType = oldAtom.elementType;
				newResidue.addAtom(newAtom);
			}

			m1.addResidue(0,newResidue);
		}
		//Determine the bonds between the atoms in the molecule
		m1.determineBonds();

		// Assign the molecule relative atom numbers
		m1.updateMoleculeAtomNumbers();

		m1.strand[0].isProtein = true;

		return m1;
	}

	//Computes the amino acid reference energies using the intra-rotamer energies from intraEnergies or intraEnergiesEntropy;
	//For each amino acid type, takes the min energy among all rotamers for that amino acid type, for all numRes residues
//	private double [] compEref(double intraEnergies[][][][][][], StrandRotamers strandRot[], int strandMut[][], double intraEnergiesEntropy[][], int numRes,
//			int mutRes2Strand[],int mutRes2StrandMutIndex[], RotamerLibrary rotLib){
//
//		int numAAallowed = rotLib.getNumAAallowed();
//
//		double bigE = (double)Math.pow(10,38);
//		double eRef[] = new double[numAAallowed];
//		for(int j=0; j<numAAallowed;j++)
//			//for (int i=0; i<eRef[j].length; i++)
//			eRef[j] = bigE;
//
//		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
//
//
//
//		for (int i=0; i<numRes; i++){
//			int str = -1;
//			int strResNum = -1;
//			if(intraEnergies!=null){
//				str = mutRes2Strand[i];
//				strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
//			}
//			int numAA = numAAallowed;
//			if (intraEnergies!=null) { //the six-dimensional, so the energies for only a subset of the amino acid types are available 
//				numAA = strandRot[str].getNumAllowable(strResNum);
//			}
//			for (int j=0; j<numAA; j++){
//				int aaInd = j;
//				if (intraEnergies!=null)
//					aaInd = strandRot[str].getIndexOfNthAllowable(strResNum,j);
//				int numRot = rotLib.getNumRotForAAtype(aaInd);
//				if (numRot==0) //ALA or GLY
//					numRot = 1;
//				double curMin = bigE;
//				for (int k=0; k<numRot; k++){
//					if (intraEnergies!=null)
//						curMin = Math.min(curMin,intraEnergies[i][aaInd][k][i][0][0]);
//					else
//						curMin = Math.min(curMin,intraEnergiesEntropy[ind][0]);
//					ind++;
//				}
//				eRef[aaInd] = Math.min(eRef[aaInd],curMin);
//			}			
//		}
//
//		//for(int j=0; j<numRes;j++)
//		for (int i=0; i<eRef.length; i++)
//			if (eRef[i]==bigE)
//				eRef[i] = 0.0f;
//
//		return eRef;
//	}
	//////////////////////////////////////////////////////
	// End Compute Residue Entropy Section
	//////////////////////////////////////////////////////

	private void selectResidues(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Residue search filename (string)		

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters	

		String runName = (String)sParams.getValue("RUNNAME");
		int numRes = (new Integer((String)sParams.getValue("NUMRES"))).intValue();
		// PGC 2013: support for kabat numbering by using strings
		String pdbRes[] = new String[numRes];
		double dist[] = new double[numRes];
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.getValue("LIGTYPE"));

		String resString = (String)sParams.getValue("RESIDUES");
		String distString = ((String)sParams.getValue("DIST"));
		for (int i=0; i<numRes; i++){
			pdbRes[i] = (String)getToken(resString,i+1);
			dist[i] = new Double((String)getToken(distString,i+1)).doubleValue();
		}

		MolParameters mp = new MolParameters();
		loadStrandParams(sParams, mp, COMPLEX);

		Molecule m = new Molecule();
		m = setupMolSystem(m,sParams,mp.strandPresent,mp.strandLimits);

		//Map from the pdb residue numbers to the residue index in m.residue[]
		int residues[] = new int[numRes];
		int curRes = 0;
		for (int i=0; i<m.numberOfResidues; i++){
			for (int j=0; j<numRes; j++){
				if (m.residue[i].getResNumber().equals(pdbRes[j])){
					residues[curRes] = i;
					curRes++;
					break;
				}
			}
		}


		boolean asProx[] = new boolean[m.numberOfResidues];
		for (int i=0; i<asProx.length; i++)
			asProx[i] = false;

		for (int res=0; res<numRes; res++){
			Residue r1 = m.residue[residues[res]];

			for (int i=0; i<asProx.length; i++){

				if (i!=residues[res]){

					Residue r2 = m.residue[i];

					if (r1.getDist(r2,true)<=dist[res])
						asProx[i] = true;
				}
				else
					asProx[i] = true;
			}
		}

		Molecule m1 = new Molecule(); //add all proximate residues; the connectivity/bonds will not be valid
		for (int i=0; i<asProx.length; i++){
			if (asProx[i]){
				m1.addResidue(0, m.residue[i]);
			}
		}
		m1.saveMolecule(runName+".pdb",0.0f);
	}

	//Computes the information necessary to generate a residue interaction graph for the given system;
	//		computes the minimum distance and minimum energy (absolute value) for each residue pair in residueMap[], 
	//		considering all possible unpruned rotamers for the given residues;
	//		ligand interactions are also computed if a ligand is present
	private void genInteractionGraph(int numMutable, RotamerSearch rs, Emat emat, String runName, 
			MutableResParams strandMut, double eInteractionCutoff, double distCutoff, Molecule m, 
			boolean usePairSt, double pairSt) {

		if (eInteractionCutoff<0.0f) //the cutoff should be non-negative, since we are comparing absolute values of energies against it
			eInteractionCutoff = 0.0f;

		double dist[][] = new double[numMutable][numMutable];
		double eInteraction[][] = new double[numMutable][numMutable];		
		for (int i=0; i<numMutable; i++){
			for (int j=0; j<numMutable; j++){
				dist[i][j] = (double)Math.pow(10, 38);
				eInteraction[i][j] = 0.0f;
			}
		}
		
		PairsIterator iter = emat.pairsIterator();
		
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			EMatrixEntry re = emeWI.eme;
			//ApplyMutation
			boolean neededMut = re.applyMutation(m,emat.resByPos,true,true);
			re.applyRC(emat.resByPos,m);
			
			Residue resi = m.residue[emat.resByPos.get(emeWI.pos1()).get(0)]; //0 Index assumes not super-rotamers
			Residue resj = m.residue[emat.resByPos.get(emeWI.pos2()).get(0)]; //0 Index assumes not super-rotamers
			
			int i = emeWI.pos1();
			int j = emeWI.pos2();
			double d = resi.getDist(resj,true);
			dist[emeWI.pos1()][emeWI.pos2()] = Math.min(dist[emeWI.pos1()][emeWI.pos2()],d);
			dist[emeWI.pos2()][emeWI.pos1()] = Math.min(dist[emeWI.pos2()][emeWI.pos1()],d);

			if ( (!usePairSt) || (re.minE()<=pairSt) ) {
				eInteraction[i][j] = Math.max(eInteraction[i][j],Math.abs(re.minE()));
				eInteraction[j][i] = Math.max(eInteraction[j][i],Math.abs(re.minE()));
			}
		}
		
		

//		for (int i=0; i<numMutable; i++){
//			Residue resi = m.residue[strandMut.allMut[i]];
//			int stri = resi.strandNumber;
//			int strResNumi = resi.strandResidueNumber;
//			for(AARotamerType aaType : resi.AATypesAllowed()) {
//
////				int AAindex1 = rs.strandRot[stri].getIndexOfNthAllowable(strResNumi,q1);
//				rs.strandRot[stri].changeResidueType(m,0,aaType.name,true,true);
//
//				int numRot1 = resi.getRCsForType(aaType).size();
//
//				for (int r1=0; r1<numRot1; r1++){
//
//					if (!emat.getSinglePruned(i,AAindex1,r1)){ //rotamer not pruned
//						
//						int globalRCid = emat.singles.getRot(i,AAindex1,r1)[0];//0 Index assumes no super-rotamers
//						
////						if(rs.doPerturbations){
//							ResidueConformation rc = ((StrandRCs)rs.strandRot[stri]).rcl.getRC(globalRCid);
//							((StrandRCs)rs.strandRot[stri]).applyRC(m, strResNumi, rc);
////						}
////						else{
////							rs.strandRot[stri].applyRotamer(m, strResNumi, r1);
////						}
//
//						for (int j=i+1; j<numMutable; j++){
//							Residue resj = m.residue[strandMut.allMut[j]];
//							int strj = resj.strandNumber;
//							int strResNumj = resj.strandResidueNumber;
//							for(int q2=0;q2<rs.strandRot[strj].getNumAllowable(strResNumj);q2++) {
//
//								int AAindex2 = rs.strandRot[strj].getIndexOfNthAllowable(strResNumj,q2);
//								rs.strandRot[strj].changeResidueType(m,strResNumj,rs.strandRot[strj].rl.getAAName(AAindex2),true,true);
//
//								int numRot2 = rs.getNumRot( strj, strResNumj, AAindex2 );
//
//								for (int r2=0; r2<numRot2; r2++){
//
//									if (!emat.getSinglePruned(j,AAindex2,r2)){
//
//										double pairE = emat.getPairwiseE( i, AAindex1, r1, j, AAindex2, r2 );
//
//										int globalRCid2 = emat.singles.getRot(j,AAindex2,r2)[0];//0 Index assumes no super-rotamers
//										
////										if(rs.doPerturbations){
//											ResidueConformation rc2 = ((StrandRCs)rs.strandRot[strj]).rcl.getRC(globalRCid2);
//											((StrandRCs)rs.strandRot[strj]).applyRC( m, strResNumj, rc2 );
////										else
////											rs.strandRot[strj].applyRotamer(m, strResNumj, r2);
//
//										double d = m.strand[stri].residue[strResNumi].getDist(m.strand[strj].residue[strResNumj],true);
//										dist[i][j] = Math.min(dist[i][j],d);
//										dist[j][i] = Math.min(dist[j][i],d);
//
//										if ( (!usePairSt) || (pairE<=pairSt) ) {
//											eInteraction[i][j] = Math.max(eInteraction[i][j],Math.abs(pairE));
//											eInteraction[j][i] = Math.max(eInteraction[j][i],Math.abs(pairE));
//										}
//									}
//								}							
//							}
//						}
//					}
//				}	
//			}
//		}

		PrintStream logPS = setupOutputFile(runName+".log");
		PrintStream logPS2 = setupOutputFile(runName);

		logPS2.println("PIG:0 "+runName); //output in Pigale-compatible ASCII format

		//Output data
		for (int i=0; i<numMutable; i++){
			Residue resi = m.residue[strandMut.allMut[i]];
			int stri = resi.strandNumber;
			int strResNumi = resi.strandResidueNumber;

			String pdbResNum1 = m.strand[stri].residue[strResNumi].getResNumber();
			for (int j=i+1; j<numMutable; j++){
				Residue resj = m.residue[strandMut.allMut[j]];
				int strj = resj.strandNumber;
				int strResNumj = resj.strandResidueNumber;
				String pdbResNum2 = resj.getResNumber();

				logPS.println(pdbResNum1+" "+pdbResNum2+" "+dist[i][j]+" "+eInteraction[i][j]);
				if ( (dist[i][j]<=distCutoff) && (eInteraction[i][j]>eInteractionCutoff) ) //these two residues interact
					logPS2.println(pdbResNum1+" "+pdbResNum2);
			}
		}
		logPS2.println("0 0");

		logPS.flush();logPS.close();
		logPS2.flush();logPS2.close();

		//outputObject(prunedRotAtRes,runName+"_pruneInfo.obj");
	}

	//Returns a molecule m1 that contains only the residues in molecule m that are specified by residueMap[] (molecul-relative residue indexing);
	private Molecule getMolRes (Molecule m, int strandMap[][]){

		Molecule m1 = new Molecule();

		for (int i=0; i<m.numberOfStrands; i++){ //create the same number of strands
			m1.addStrand(m.strand[i].name);
			m1.strand[i].isProtein = m.strand[i].isProtein;
		}

		for(int str=0; str<strandMap.length;str++){
			for (int i=0; i<strandMap[str].length; i++){

				Residue oldResidue = m.residue[strandMap[str][i]];

				Residue newResidue = new Residue();
				newResidue.name = oldResidue.name;
				newResidue.fullName = oldResidue.fullName;

				for (int j=0; j<oldResidue.numberOfAtoms; j++){

					Atom oldAtom = oldResidue.atom[j];

					Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
					newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
					newAtom.strandNumber = oldAtom.strandNumber;
					newAtom.elementType = oldAtom.elementType;
					newResidue.addAtom(newAtom);
				}

				m1.addResidue(oldResidue.strandNumber,newResidue);
			}
		}

		//Determine the bonds between the atoms in the molecule
		m1.determineBonds();

		// Assign the molecule relative atom numbers
		m1.updateMoleculeAtomNumbers();

		return m1;
	}

	//////////////////////////////////////////////////////
	// Begin Steric Overlap Check Section
	//////////////////////////////////////////////////////
	//Compute the amount of overlap between a set of structures and a reference structure
	public void handleCompStericOverlap (String s) {

		// Takes the following parameters
		// 1: System parameter filename for the reference structure (string)
		// 2: System parameter filename for the set of structures to compare (string)
		// 2: Mutation search parameter filename (string)		

		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		Molecule mRef = new Molecule();

		MolParameters mp = new MolParameters();
		loadStrandParams(sParams, mp, COMPLEX);

		mRef = setupMolSystem(mRef,sParams,mp.strandPresent,mp.strandLimits);
		int numInASref = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();		
		int posMapRef[] = new int[numInASref]; //molecule-relative numbering
		String posMapString = (String)sParams.getValue("RESIDUEMAP");
		for(int i=0;i<numInASref;i++)
			posMapRef[i] = (new Integer(getToken(posMapString,i+1))).intValue();



		// Read System parameters for the set of structures to compare
		sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,3)); //read system parameters
		sParams.addParamsFromFile(getToken(s,4)); //read system parameters

		String runName = (String)sParams.getValue("RUNNAME");

		String protPDBname = (String)sParams.getValue("PROTPDBNAME");
		int numPDBfiles = (new Integer((String)sParams.getValue("NUMPDBFILES"))).intValue();
		int numInAS2 = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();

		int posMap2[] = new int[numInAS2]; //molecule-relative numbering
		posMapString = (String)sParams.getValue("RESIDUEMAP");
		for(int i=0;i<numInAS2;i++)
			posMap2[i] = (new Integer(getToken(posMapString,i+1))).intValue();

		String pdbFiles[] = getPDBfiles(protPDBname,numPDBfiles); //get the PDB filenames
		numPDBfiles = pdbFiles.length;

		double minMaxOverlap[] = new double[numInAS2];
		int minMaxOverlapStruct[] = new int[numInAS2];
		for (int i=0; i<minMaxOverlap.length; i++){
			minMaxOverlap[i] = (double)Math.pow(10, 10);
			minMaxOverlapStruct[i] = -1;
		}

		// The idea is: given a structure from the set, for each of the residues in posMap2[], determine the largest
		//		steric overlap with an atom in the posMapRef[] residues from the reference structure;
		//		then, for each residue in posMap2[], find the minimum such largest overlap among all structures in the set
		for (int i=0; i<numPDBfiles; i++){

			System.out.println("Starting structure "+pdbFiles[i]+" ("+i+")");

			sParams.setValue("PDBNAME",pdbFiles[i]);

			//Setup the molecule system
			Molecule m2 = new Molecule();
			m2 = setupMolSystem(m2,sParams,mp.strandPresent,mp.strandLimits);

			for (int res2=0; res2<numInAS2; res2++){ //for each included residue in the given structure
				double maxOverlap = 0.0;
				for (int at2=0; at2<m2.residue[posMap2[res2]].numberOfAtoms; at2++){ //for each atom in that residue
					Atom a2 = m2.residue[posMap2[res2]].atom[at2];
					if ( hSteric || (!a2.elementType.equalsIgnoreCase("H"))){
						for (int resRef=0; resRef<numInASref; resRef++){ //for each included residue in the reference structure
							for (int atRef=0; atRef<mRef.residue[posMapRef[resRef]].numberOfAtoms; atRef++){ //for each atom
								Atom aRef = mRef.residue[posMapRef[resRef]].atom[atRef];
								if ( hSteric || (!aRef.elementType.equalsIgnoreCase("H"))){
									double overlap = ((a2.radius + aRef.radius)/100.0) - a2.distance(aRef);
									if (overlap<0.0)
										overlap = 0.0;
									maxOverlap = Math.max(maxOverlap, overlap);
								}
							}
						}
					}
				}
				if (minMaxOverlap[res2]>maxOverlap){
					minMaxOverlap[res2] = maxOverlap;
					minMaxOverlapStruct[res2] = i;
				}
			}
		}


		//Output the computed distances
		PrintStream logPS = setupOutputFile(runName);		
		for (int i=0; i<numInAS2; i++){
			logPS.println(posMap2[i]+" "+minMaxOverlap[i]+" "+pdbFiles[minMaxOverlapStruct[i]]);
		}		
		logPS.flush();
		logPS.close();
	}

	//Reads the pdb filenames
	private String [] getPDBfiles(String fName, int numFiles){

		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... pdbs config file not found");
			System.exit(1);
		}

		String pdbFiles[] = new String[numFiles];

		boolean done = false;
		String str = null;
		int curFile = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {
				if (curFile>=numFiles)
					break;
				else {
					pdbFiles[curFile] = getToken(str,1);
					curFile++;
				}
			}
		}

		if (curFile<numFiles){
			String tmp[] = new String[curFile];
			System.arraycopy(pdbFiles, 0, tmp, 0, tmp.length);
			pdbFiles = tmp;
		}

		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		return pdbFiles;
	}

	//////////////////////////////////////////////////////
	// End Steric Overlap Check Section
	//////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	// Begin Backrub Precomputation Section
	//////////////////////////////////////////////////////
	public void handlePrecomputeBackrubs (String s) {

		// Takes the following parameters
		// 1: System parameter filename
		// 2: Number of backrub samples in each direction
		// 3: Backrub step size
		// 4: Output file name

		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		int numBackrubSamples = new Integer(getToken(s,3)).intValue();
		double backrubStepSize = new Double(getToken(s,4)).doubleValue();
		String backrubFile = getToken(s,5);


		MolParameters mp = new MolParameters();
		loadMolecule(sParams, -1, false, 0.0, true);
		

		Amber96ext a96ff = new Amber96ext(mp.m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier,hbonds);
		a96ff.calculateTypesWithTemplates();

		BackrubMinimizer brMin = new BackrubMinimizer();
		brMin.initialize(mp.m, a96ff, mp.strandMut, backrubFile, hSteric, stericThresh, mp.numOfStrands, false);
		brMin.precomputeBackrubs(numBackrubSamples, backrubStepSize);
		System.out.println("DONE: Backrub angle precomputation..");
	}

	//Compute the amount of overlap between a set of structures and a reference structure
	/*public void handlePrecomputeBackrubs1 (String s) {

		// Takes the following parameters
		// 1: System parameter filename
		// 2: Number of backrub samples in each direction
		// 3: Backrub step size
		// 4: Output file name

		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		int numBackrubSamples = new Integer(getToken(s,3)).intValue();
		double backrubStepSize = new Double(getToken(s,4)).doubleValue();
		String backrubFile = getToken(s,5);

		Molecule m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();		
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();

		Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
		a96ff.calculateTypesWithTemplates();

		BackrubMinimizer brMin = new BackrubMinimizer();
		brMin.initialize(m, a96ff, residueMap, sysStrNum, backrubFile, hSteric, stericThresh);
		brMin.precomputeBackrubs(numBackrubSamples, backrubStepSize);

		System.out.println("DONE: Backrub angle precomputation..");
	}*/

	//////////////////////////////////////////////////////
	// End Backrub Precomputation Section
	//////////////////////////////////////////////////////

	int[] getCurrentMutOffset(int[][] origStrandMut, int[][] currStrandMut, boolean[] strandPresent, boolean[] oldStrandPresent){
		int offsetMat[] = new int[origStrandMut.length];
		int strCtr=0;
		int offset=0;
		for(int str=0; str<strandPresent.length;str++){
			if(strandPresent[str] && oldStrandPresent[str]){
				offsetMat[strCtr] = offset;
				strCtr++;
			}
			else if(oldStrandPresent[str]){
				offsetMat[strCtr] = offset;
				strCtr++;
				offset += origStrandMut[str].length;
			}
		}
		return offsetMat;
	}

	String[] shrinkCurrentMutation(String currentMut[], int origStrandMut[][] ,int newStrandMut[][],boolean origStrandPresent[], boolean newStrandPresent[]){
		//Find the length of the new mutable array
		int newSize=0;
		for(int i=0;i<newStrandMut.length;i++)
			newSize+=newStrandMut[i].length;
		String [] newCurMut = new String[newSize];

		int NewCtr=0;
		int OldCtr=0;
		for(int str=0;str<origStrandMut.length;str++){
			for(int i=0;i<origStrandMut[str].length;i++){
				if(newStrandPresent[str]&&origStrandPresent[str]){
					newCurMut[NewCtr] = currentMut[OldCtr];
					NewCtr++;
					OldCtr++;
				}
				else if(origStrandPresent[str]){
					OldCtr++;
				}
			}
		}
		return newCurMut;
	}

	private int[] rotamersRemaining(int numRotForRes[], PrunedRotamers<Boolean> prunedRotAtRes){


		//initialize numPrunedRotForRes[]
		int numNotPrunedForRes[] = new int[numRotForRes.length]; //after pruning for this partition
		for(int i=0; i<numNotPrunedForRes.length;i++){
			numNotPrunedForRes[i] = numRotForRes[i];
		}



		Iterator<RotInfo<Boolean>> i = prunedRotAtRes.iterator(); 

		while(i.hasNext()){
			RotInfo<Boolean> ri = i.next();
			if(ri.state){
				numNotPrunedForRes[ri.curPos]--;
				if(numNotPrunedForRes[ri.curPos] == 0){
					System.err.println("ALL OF THE ROTAMERS HAVE BEEN PRUNED at site: "+ri.curPos);
					Exception e = new Exception();
					e.printStackTrace();
					System.exit(0);
				}
			}
		}

		return numNotPrunedForRes;

	}


	private void selectPerturbations(MolParameters mp, boolean doPerturbations, String pertFile, boolean minimizePerts,
			boolean addWTRot, ParamSet sParams ){


		//Parameters for perturbation selection
		double min_rmsd = (new Double((String)sParams.getValue("MINRMSD","0"))).doubleValue();//Minimum backbone heavy-atom RMSD for perturbation states: used to select perturbations
		Shear.setParams( ((String)sParams.getValue("SHEARPARAMS","none")) );//Default shear parameters (for perturbation selection): State 1 min, state 1 max, state 2 min, state 2 max,...
		Backrub.setParams( ((String)sParams.getValue("BACKRUBPARAMS","none")) );//Default backrub parameters (same format)
		RamachandranChecker.denCutoff = (new Double((String)sParams.getValue("RAMACHANDRANCUTOFF","0.02"))).doubleValue();//Density cutoff for Ramachandran filter


		RotamerSearch rs = new RotamerSearch(mp.m,mp.strandMut.numMutPos(), mp.strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE,
				solvScale, softvdwMultiplier, doPerturbations, pertFile, minimizePerts, false, false,new EPICSettings(),hbonds,mp.strandMut);


		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(mp.strandPresent, sParams);
		int molStrand = 0;
		for (int resID: mp.strandMut.allMut){
				setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}

		//Read the starting perturbation file, if any
		//The perturbations in this file will be included along with any the selector generates
		//and their perturbation states will be selected as appropriate
		String startingPertFile = (String)sParams.getValue("STARTINGPERTURBATIONFILE","none");//Input file giving perturbation information
		boolean onlyStarting = (new Boolean((String)sParams.getValue("ONLYSTARTINGPERTURBATIONS", "false"))).booleanValue();//Use only the perturbations specified in the startingPerturbationFile
		//(select RCs for them without selecting any more perturbations)

		if( startingPertFile.equalsIgnoreCase("none") && onlyStarting ){
			System.err.println("ERROR: Perturbation selector can't use only starting perturbations if startingPerturbationFile is set to 'none'");
		}

		PerturbationSelector ps = new PerturbationSelector(mp.strandMut.numMutPos(), mp.strandMut, 
				mp.m,	rs.strandRot, min_rmsd, startingPertFile,onlyStarting);

		ps.selectPerturbations();
		rs.removeImpossibleRCs(mp.strandMut.numMutPos(), mp.strandMut);
		PertFileHandler.writePertFile(pertFile, mp.m, null, rs.strandRot, null, false);
	}



	public void fixStruct(String s){
		//This simple command just reads in a structure with the autofix feature on
		//and then outputs the fixed structure
		//the first argument is the PDB file to read in and the second is the filename to write to
		//(note that this is different from the usual use of configuration files)
		EnvironmentVars.autoFix = true;
		String fileIn = getToken(s,2);
		String fileOut = getToken(s,3);

		//Stop if configuration files are provided to avoid overwriting them or something
		if( fileIn.endsWith(".cfg") || fileOut.endsWith(".cfg") ){
			System.err.println("ERROR: Just provide PDB file names for fixStruct");
			System.exit(1);
		}

		Molecule m = new Molecule();

		try{
			FileInputStream is = new FileInputStream(fileIn);
			new PDBChemModel(m, is);
		}
		catch (Exception e){
			System.out.println("WARNING: An error occurred while reading file");
			System.out.println(e);
			e.printStackTrace();
			System.exit(1);
		}

		m.saveMolecule(fileOut, Double.NaN);
	}

	static public Object loadObject(String fileName) {
		Object o;

		try{
			//			ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(fileName)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
			o = in.readObject();
			in.close();
		}
		catch (Exception e){
			//e.printStackTrace();
			System.out.println("Couldn't find file "+fileName);
			return null;
		}		

		return o;
		//arpMatrix.reconnect(m);

	}

	class FullConf implements Comparable<FullConf>{
		String[] AAnames;
		int[] rotNums;
		int[] rcNums;
		String[] pdbNums;
		double unMinE;
		double minE;

		public FullConf(int numRes){
			AAnames = new String[numRes];
			rotNums = new int[numRes];
			rcNums = new int[numRes];
			pdbNums = new String[numRes];
		}

		@Override
		public int compareTo(FullConf conf2) {
			if(minE < conf2.minE)
				return -1;
			else if(minE > conf2.minE)
				return 1;
			else
				return 0;
		}
	}
	
	public static void gurobiSlave(CommucObj cObj) throws MPIException, InterruptedException {
		
		GurobiCalcInfo gci = cObj.gurobiCalcInfo;


		int ctr=0;
		Object s;
		while(true){

			s = MPItoThread.Probe(0, ThreadMessage.ANY_TAG);
			if (s!=null) {
				if (MPItoThread.getStatusTag(s)==regTag){ //new computation		
					PGQueueNode[] nodes = new PGQueueNode[1];
					MPItoThread.Recv(nodes, 0, 1, ThreadMessage.OBJECT, 0, regTag);
					PGQueueNode node = nodes[0];

					if(cObj.gurobiCalc){
						System.out.println("Doing Nothing Right Now");
//						gurobiSlaveHelper(node, gci);
					}
					else{
						wcspSlaveHelper(node,gci);
					}

					MPItoThread.Send(nodes, 0, 1, ThreadMessage.OBJECT, 0, updateTag); //send back result

				}
				else if (MPItoThread.getStatusTag(s)==doneTag){
					int[] nothing = new int[1];
					MPItoThread.Recv(nothing, 0, 1, ThreadMessage.INT, 0, doneTag);
					return;
				}
			}

		}
	}
	
	private static void wcspSlaveHelper(PGQueueNode node, GurobiCalcInfo gci) {
		
		WCSPOptimization optimizer;
		
		if(gci.bySubRot)
			optimizer = new WCSPOptimization(node, gci.emat, gci.numParentRotPerLvl, gci.parentRotIndexPerLvl, 
					gci.numSubRotPerParentRot, gci.subRotsPerLvlPerParent, gci.numNodesForLevel, gci.upperE);
		else
			optimizer = new WCSPOptimization(node,gci.emat,gci.seqIndicesPerlevel,gci.numRotRemainingBySeq,
				gci.numNodesForLevel, gci.upperE);
		
		
		System.out.print(".");
		//double E =optimizer.optimize(null);
		double E = optimizer.getBound(null);
		optimizer.cleanUp();
		node.fScore = E;
	}
	
	/**
	 * Given the array of K* sequences this function checks whether
	 * any of the sequences are duplicates and marks them as duplicates
	 * so their scores only need to be computed once by once by K*
	 * @param mutArray Array that holds all of the K* slave sequences
	 * @param m Molecule m
	 */
	void checkDuplicateMutations(OneMutation[] mutArray, Molecule m){
		
		//Hash table that is used to store all of the sequences for each strand
		Hashtable<String, Integer>[] sequences = new Hashtable[m.numberOfStrands];
		//Set up the hashtable
		for(int i=0; i<sequences.length;i++){
			sequences[i] = new Hashtable<String, Integer>();
		}

		//check for duplications
		int[] numMutPerStrand = new int[m.numberOfStrands];
		for(int i=0; i<numMutPerStrand.length;i++)
			numMutPerStrand[i] = m.numMutableForStrand(i);

		for(int i=0; i<mutArray.length; i++){
			OneMutation mut = mutArray[i];
			mut.duplicateMut = new int[m.numberOfStrands];
			for(int q=0; q<mut.duplicateMut.length; q++)
				mut.duplicateMut[q] = -1;
			
			//Initialize sequences for this OneMutation entry
			String[] seqs = new String[m.numberOfStrands];
			for(int j=0; j<seqs.length; j++) 
				seqs[j] = "";

			//Build the sequences
			int strOffset=0; 
			int curStr=0;
			for(int j=0; j<mut.resTypes.length;j++){
				if(j-strOffset >= numMutPerStrand[curStr]){
					strOffset += numMutPerStrand[curStr];
					curStr++;
				}
				
				seqs[curStr] += mut.resTypes[j]+" ";
				
			}

			for(int j=0; j<seqs.length; j++){
				if(sequences[j].containsKey(seqs[j])){ //If sequence is a duplicate
					//Set the current OneMutation as a duplicate so it's q calculation can be skipped
					mut.duplicateMut[j] = sequences[j].get(seqs[j]);
					//Tell the independent mutArray that is has things depending on it 
					mutArray[sequences[j].get(seqs[j])].setDupMut(DUPFOUND, j);
				}
				else{
					//Add the sequence to our list of sequences being computed so far
					sequences[j].put(seqs[j], i);
				}
			}

		}

	}
	
	/**
	 * Performs a DEE (Traditional, MinDEE, BD, or BRDEE) pruning with A* search, with or without DACS;
	 * the only parameter 's' (String) includes the command-line arguments specifying the filenames of the two input configuration files.
	 * If distributed DACS is performed, the computation is distributed to the available processors for evaluation.
	 */
	public void handleExpandedIMinDEE(String s, boolean doDEE) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)

//		metrics.setStartTime();

		System.out.println("Performing DEE");

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		Settings settings = new Settings();
		
		/******** Load all of the settings for DEE *******/
		// Pull search parameters
		String runName = Settings.getRunName(sParams);
		
		//DEE Settings
		Settings.DEE deeSettings = settings.new DEE(sParams);
		double difference = deeSettings.Ival;
		
		//Minimization Settings
		Settings.Minimization minSettings = settings.new Minimization(sParams);
		
		
		//EPICSettings
		EPICSettings es = new EPICSettings(sParams);
		if(deeSettings.Ival+deeSettings.initEw>es.EPICThresh2){
			System.out.println("EPICThresh2 must be at least Ival+Ew: raising to Ival="+(deeSettings.Ival+deeSettings.initEw));
			es.EPICThresh2 = deeSettings.Ival+deeSettings.initEw;
		}
		
		//Enumeration Settings
		Settings.Enum enumSettings = settings.new Enum(sParams);
		
		//Emat Settings
		Settings.Emat ematSettings = settings.new Emat(sParams, runName, minSettings.doPerturbations);
		
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		
		//Output Settings
		Settings.Output outputSettings = settings.new Output(sParams, runName);
				
		//Unclassified Settings
		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
				
		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
		String resumeFilename ="";
		if(resumeSearch){
			resumeFilename = ((String)sParams.getValue("RESUMEFILENAME"));
		}
		
		/***************** Done Loading Settings ************/
		
		minSettings.doMinimize = true; //Minimization is turned on for the expanded rotamer library calculation
		if (!minSettings.doMinimize) //no minimization
			minSettings.minimizeBB = false;
		if (!minSettings.minimizeBB) //not backbone minimization
			minSettings.doBackrubs = false;

		if (graphSettings.genInteractionGraph) //DACS is not performed when generating the interaction graph
			deeSettings.doDACS = false;

		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams,curStrForMatrix,graphSettings.neighborList,graphSettings.distCutoff,true);

		boolean reload=false;
		if(minSettings.selectPerturbations){//Need to run the automatic perturbation selection
			//This only needs to be done once though: after that the perturbations can be read from pertFile
			selectPerturbations(mp, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, ematSettings.addWTRot, sParams);
			reload = true;
		}
		
		// 2010: If useMinDEEPruningEw is set to false, this cycle never repeats itself.
		//  If it is set to true, it can repeat at most once: if none of the rotamer vectors
		//  between the conformation of lowest energy (i.e. lowestBound) 
		//  and lowestBound+InitEw can minimize to a lower energy than lowestBound+InitEw, 
		//   then let minimumEnergy be the minimum nergy found among the enumerated conformations,
		//   we set a new Ew equal to minimumEnergy - lowestBount and repeat this cycle.  We
		//   only have to do it at most twice.  
		//		long start;
		//		long totalEtime = 0;
		//		start = System.currentTimeMillis();	
		long totalConfsEvaluated = 0;

		RotamerSearch rs = new RotamerSearch(mp.m, mp.strandMut.numMutPos(),mp.strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, 
				doSolvationE, solvScale, softvdwMultiplier, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, 
				deeSettings.useTriples, enumSettings.useFlagsAStar, es,hbonds, mp.strandMut);

		rs.useCCD = minSettings.useCCD;
//		rs.superRotamers = superRotamers;
//		rs.tuples = doTuples || readTuples;
//		rs.subRotamers = subRotamers;

		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(mp.strandPresent, sParams);
		int molStrand = 0;
		for(int resID:mp.strandMut.allMut){
			setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}
		
		if(!minSettings.selectPerturbations || reload)
			rs.setupRCs(minSettings.doPerturbations);
		
		int numSplits = 0;

		//Load Energy Matrix
//		long startEmat = System.currentTimeMillis();
		System.out.print("Loading precomputed energy matrix...");
		Emat emat = loadPairwiseEnergyMatrices(sParams,ematSettings.runNameEMatrixMin,minSettings.doMinimize,curStrForMatrix,es,mp.m,true);
		rs.setMinMatrix(emat);
		System.out.println("done");
//		long endEmat = System.currentTimeMillis();
//		metrics.EmatTime += (endEmat-startEmat);

		
		
		
		//Recalculate energy of the matrix
		sParams.setValue("doMinimize", "false");
		File ematFile = new File(ematSettings.runNameEMatrixMin+"_expanded_COM.dat.templE");

		if(!ematFile.exists()){
			
			/********* Get rid of any minimized rotamers that we know won't be part of the GMEC ******/ 
			//Prune down the matrix to find the conformation with lowest bound
			double Ival = 0.0; //Prune as much as possible to find the lowest bound conformation
			double initEw = 0.0;
			boolean localUseMinDEEPruningEw = true;
			boolean removeRot = false;
			int timesRunFullPairs = 0;
			int maxLoopNum = 1;
			runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
					initEw, deeSettings.maxIntScale, deeSettings.typeDep, Ival,
					emat, localUseMinDEEPruningEw, removeRot,deeSettings.stericE, sParams,timesRunFullPairs,
					maxLoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);
			
			double bestScore = Double.POSITIVE_INFINITY;
			double ENUMIval = 0.0;
			int numToEnumerate = 1;
			
			//Create temporary enumSettings to do a quick GMEC search
			Settings.Enum tmpEnumSettings = enumSettings.copy();
			tmpEnumSettings.asMethod = Settings.ASTARMETHOD.WCSP;
			
			AStarResults asr = rs.doAStarGMEC("preliminaryConfs.conf",true,minSettings.doMinimize, 
					mp.strandMut.numMutPos(), mp.strandMut,0.0,
					bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
					ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
					localUseMinDEEPruningEw, ENUMIval,tmpEnumSettings);
			Ival = Math.min(enumSettings.bestE-asr.lowestBound, asr.bestE-asr.lowestBound);
			emat.unPrune();
			removeRot = true;
			runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
					deeSettings.initEw, deeSettings.maxIntScale, deeSettings.typeDep, Ival,
					emat, localUseMinDEEPruningEw, removeRot, 
					deeSettings.stericE,  sParams,timesRunFullPairs,
					maxLoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);
			
			/************* The Rotamers remaining we know are good rotamers **************/
			
			ArrayList<RCToAdd> rigidRotsToAdd = identifyRotamerDifferences(emat,mp.m);
			emat.convertToRigid(); 
			emat.pairs = null;
			//Add all of the new rigid rots
			ArrayList<ArrayList<ArrayList<ArrayList<ResidueConformation>>>> rcsToAddByPosAA = emat.getEmptyRCsToAdd();

			for(RCToAdd rcTA: rigidRotsToAdd){
				if(rcTA.rc.rot.aaType.numDihedrals()!=0){ //Don't add an ALA or GLY rotamer (also update this not to add rotamers that are too close together) 
					//					ArrayList<ArrayList<Rotamer>> rotsForPos = new ArrayList<ArrayList<Rotamer>>();
					int resid = emat.resByPos.get(rcTA.pos).get(0);
					Residue res = mp.m.residue[resid];
					String pdbNum = res.getResNumberString();
					for(double[] dihedVals: rcTA.confs){
						double[] minWidths = new double[dihedVals.length];
						for(int i=0; i<dihedVals.length;i++){//The dihedVals are just offsets so we need to convert them to actual dihedrals
							minWidths[i] = 9.0 - Math.abs(dihedVals[i]);
							dihedVals[i] += rcTA.rc.rot.values[i];
						}
						ArrayList<ResidueConformation> rcToAdd = new ArrayList<ResidueConformation>();
						//Add the new rotamer and the new Residue Conformation						
						Rotamer rot = mp.m.residue[resid].rl.addRotamer(rcTA.rc.rot.aaType.name, pdbNum,dihedVals,minWidths,false);
						ResidueConformation rc = mp.m.strand[res.strandNumber].rcl.addResidueConformation(rot, rcTA.rc.pertState, res.strandResidueNumber); 
						mp.m.residue[resid].setAllowable(rc);
						rcToAdd.add(rc);

						int a1 = emat.getRCAALoc(mp.m,rcTA.pos,rcToAdd);

						rcsToAddByPosAA.get(rcTA.pos).get(a1).add(rcToAdd);
						//						rotsForPos.add(rotToAdd);
					}

				}
				else{
					System.out.println("Skipped Rot");
				}
			}

			emat.addAllRCs(mp.m,rcsToAddByPosAA,true);

//			startEmat = System.currentTimeMillis();
			emat = recalculateAllEnergies(emat, -1,sParams,mp.m,deeSettings.stericE,minSettings.minimizeBB,minSettings.doBackrubs,
					minSettings.minScheme, mp.strandMut, mp.strandsPresent, es,
					minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts);
//			endEmat = System.currentTimeMillis();
//			metrics.EmatTime += (endEmat - startEmat);
			emat.save(ematSettings.runNameEMatrixMin+"_expanded_COM.dat",mp.m);
		}
		//else{
		emat = null;
		System.gc();
		//emat = new Emat(runNameEMatrixMin+"_expanded_COM.dat",false);
		//}

		if(doDEE){
			sParams.setValue("minEnergyMatrixName",ematSettings.runNameEMatrixMin+"_expanded" );
			sParams.setValue("doMinimize","false");
			handleDoDEE(sParams);
		}

//		metrics.setEndTime();

//		metrics.print();



	}
	
	private ArrayList<RCToAdd> identifyRotamerDifferences(Emat emat, Molecule m) {

		ArrayList<RCToAdd> rcsToAdd = new ArrayList<RCToAdd>();

		SinglesIterator iter = emat.singlesIterator();
		while(iter.hasNext()){
			EMatrixEntryWIndex rot = iter.next();

//			Rotamer r = emat.getRotamers(m, rot.rot1index3()).get(0);
			ResidueConformation rc = emat.getRCs(m,rot.rot1index3()).get(0); //0 index assumes no superRotamers
			
			//First we have to figure out the clusters for the rotamer
			ArrayList<double[]> clusters = new ArrayList<double[]>();

			clusters.add(emat.singles.rotDih[rot.pos1()][rot.aa1()][rot.rot1()]);

			for(int p2 = 0; p2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()].length;p2++){
				if(p2 != rot.pos1() && emat.areNeighbors(p2, rot.pos1())){
					for(int a2 = 0; a2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2].length;a2++){
						for(int r2 = 0; r2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2].length;r2++){
							if(emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2] < Double.POSITIVE_INFINITY){
								//For every cluster see if rotDih1 matches that cluster or we need to add a new cluster
								boolean inCluster = false;
								for(double[] cluster: clusters){
									inCluster = inCluster || clustersMatch(emat.pairs.rotDih1[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2], cluster, emat);
								}
								if(!inCluster){
									clusters.add(emat.pairs.rotDih1[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2]);
								}
							}

						}
					}
				}
			}


			//For every cluster figure out the lowest score
			ArrayList<Double> values = new ArrayList<Double>();
			for(double[] cluster: clusters){
				double minBound = 0;

				if(clustersMatch(cluster,emat.singles.rotDih[rot.pos1()][rot.aa1()][rot.rot1()],emat))
					minBound += emat.singles.E[rot.pos1()][rot.aa1()][rot.rot1()];
				else
					minBound += emat.singles.maxE[rot.pos1()][rot.aa1()][rot.rot1()];

				for(int p2 = 0; p2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()].length;p2++){
					double minEforPos = Double.POSITIVE_INFINITY;
					if(p2!= rot.pos1() && emat.areNeighbors(p2, rot.pos1())){
						for(int a2 = 0; a2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2].length;a2++){
							for(int r2 = 0; r2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2].length;r2++){
								if(emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2] < Double.POSITIVE_INFINITY)
								if(clustersMatch(cluster,emat.pairs.rotDih1[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2],emat))
									minEforPos = Math.min(minEforPos, emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2]);
								else
									minEforPos = Math.min(minEforPos, emat.pairs.maxE[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2]);
							}
						}

						minBound += minEforPos;
					}
				}
				values.add(minBound);

			}

			//calculate the most optimistic bound
			double bestBound = 0;
			bestBound += emat.singles.E[rot.pos1()][rot.aa1()][rot.rot1()];
			for(int p2 = 0; p2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()].length;p2++){
				double minEforPos = Double.POSITIVE_INFINITY;
				if(p2!= rot.pos1() && emat.areNeighbors(p2, rot.pos1())){
					for(int a2 = 0; a2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2].length;a2++){
						for(int r2 = 0; r2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2].length;r2++){
							minEforPos = Math.min(minEforPos, emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2]);
						}
					}

					bestBound += minEforPos;
				}
			}

			//calculate the most optimistic bound
			double rigidBound = 0;
			rigidBound += emat.singles.maxE[rot.pos1()][rot.aa1()][rot.rot1()];
			for(int p2 = 0; p2<emat.pairs.maxE[rot.pos1()][rot.aa1()][rot.rot1()].length;p2++){
				double minEforPos = Double.POSITIVE_INFINITY;
				if(p2!= rot.pos1() && emat.areNeighbors(p2, rot.pos1())){
					for(int a2 = 0; a2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2].length;a2++){
						for(int r2 = 0; r2<emat.pairs.E[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2].length;r2++){
							minEforPos = Math.min(minEforPos, emat.pairs.maxE[rot.pos1()][rot.aa1()][rot.rot1()][p2][a2][r2]);
						}
					}

					rigidBound += minEforPos;
				}
			}


			//Check the difference between the bestBound and the more realistic bounds for each cluster
			ArrayList<double[]> clustersToAdd = new ArrayList<double[]>();
			double minDiff = Double.POSITIVE_INFINITY;
			int ctr= 0;
			int clusterCtr=0;
			for(double bound: values){
				double diff = bound - bestBound;
				double rigidDiff = bound - rigidBound;
				if(/*diff+bestBound < 0 &&*/ rigidDiff < -0.4){
					ctr++;
					clustersToAdd.add(clusters.get(clusterCtr));
				}
				minDiff = Math.min(minDiff, diff);
				clusterCtr++;
			}

			rcsToAdd.add(new RCToAdd(clustersToAdd, rc, rot.pos1()));
			System.out.println(rot.toString()+" "+minDiff+" "+clusters.size()+" "+bestBound+" "+ctr);

		}

		return rcsToAdd;

	}
	
	public boolean clustersMatch(double[] cluster1, double[] cluster2, Emat emat){ 

		if(cluster1.length == 0 || cluster2.length == 0) //clusters can be length 0 if the rotamer or pair was pruned
			return true;
		
		for(int i=0; i<cluster1.length;i++){
			if(Math.abs(cluster1[i]-cluster2[i]) > 2)
				return false;
		}
		

		return true;
	}
	
	/**
	 * 
	 * 
	 * 
	 */
	private Emat recalculateAllEnergies(Emat emat, int curStrForMatrix,ParamSet params, Molecule m, double stericE,
			boolean minimizeBB, boolean doMinimize, RotamerSearch.MINIMIZATIONSCHEME minScheme, MutableResParams strandMut,
			int strandsPresent, EPICSettings es, boolean doPerturbations, String pertFile, boolean minimizePerts) {
		//Emat ematToSend = new Emat(emat, mutPos);
		//KER: split this into several small runs for each AA vs. Pos


		OneMutation singleMutArray[] = getMutArraySingleEcomp(emat.numMutPos(),minimizeBB);

		MutationManager mutMan = new MutationManager(null, singleMutArray, true);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setcurStrForMatrix(curStrForMatrix);
		mutMan.setParams(params);
		mutMan.setMolecule(m);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setES(es);
		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setPertFile(pertFile);
		mutMan.setMinimizePerts(minimizePerts);

		try {
			handleDoMPIMaster(mutMan, singleMutArray.length);
		} catch (Exception e) {
			e.printStackTrace();
		}

		int[] numRotPerPos = emat.numRotPerPos();
		System.out.print("RotsPerPos: ");
		for(int i:numRotPerPos)
			System.out.println(i+" ");

		
		RotamerSearch.DoPruneStericTemplate(emat,stericE, false,outPS);
		emat.removePrunedRotReducedMem(true);

		numRotPerPos = emat.numRotPerPos();
		System.out.print("RotsPerPos: ");
		for(int i:numRotPerPos)
			System.out.println(i+" ");

		emat = new Emat(emat,m,false);

		OneMutation pairMutArray[] = getMutArrayPairEcomp(emat.numMutPos(),minimizeBB,emat);

		//		OneMutation mutArray[] = new OneMutation[singleMutArray.length+pairMutArray.length];

		//		System.arraycopy(singleMutArray, 0, mutArray, 0, singleMutArray.length);
		//		System.arraycopy(pairMutArray, 0, mutArray, singleMutArray.length, pairMutArray.length);


		HashMap<String,double[]> eRef = new HashMap<String,double[]>();
		for(int res: emat.allMutRes()){
			eRef.put(m.residue[res].getResNumberString(), new double[m.residue[res].rl.getNumAAallowed()]);
		}


		mutMan = new MutationManager(null, pairMutArray, true);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setcurStrForMatrix(curStrForMatrix);
		mutMan.setParams(params);
		mutMan.setMolecule(m);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setErefMatrix(eRef);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setES(es);
		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setPertFile(pertFile);
		mutMan.setMinimizePerts(minimizePerts);
		
		try {
			handleDoMPIMaster(mutMan, pairMutArray.length);
		} catch (Exception e) {
			e.printStackTrace();
		}

		//KER: do we have to set these if we pass emat to the mutMan?
		//emat.singles = mutMan.pairEMatrixMin.singlesintraE[mutPos];
		//emat.pairs = mutMan.pairEMatrixMin.pairE[mutPos];

		return emat;

	}
	
	
	/**
	 * Performs a partitioned rotamer DEE with enumeration,
	 * the only parameter 's' (String) includes the command-line arguments specifying the filenames of the two input configuration files.
	 * 
	 */
	public void handlePartitionedDEE(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)

		System.out.println("Performing DEE");

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		Settings settings = new Settings();
		
		// Pull search parameters
		String runName = Settings.getRunName(sParams);
		
		//DEE Settings
		Settings.DEE deeSettings = settings.new DEE(sParams);
		
		//Minimization Settings
		Settings.Minimization minSettings = settings.new Minimization(sParams);
		
		
		//EPICSettings
		EPICSettings es = new EPICSettings(sParams);
		if(deeSettings.Ival+deeSettings.initEw>es.EPICThresh2){
			System.out.println("EPICThresh2 must be at least Ival+Ew: raising to Ival="+(deeSettings.Ival+deeSettings.initEw));
			es.EPICThresh2 = deeSettings.Ival+deeSettings.initEw;
		}
		
		//Enumeration Settings
		Settings.Enum enumSettings = settings.new Enum(sParams);
		
		//Emat Settings
		Settings.Emat ematSettings = settings.new Emat(sParams, runName, minSettings.doPerturbations);
		
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		
		//Output Settings
		Settings.Output outputSettings = settings.new Output(sParams, runName);
				
		//Unclassified Settings
		int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
		
		Settings.ImprovedBounds ibSettings = settings.new ImprovedBounds(sParams);
		
		
		// 2010: Use energy window MinDEE method.  If this is set to true,
		//   MinDEE will use traditional DEE with an energy window (initEw) 
		//   for pruning.  Max terms will be ignored and only the min terms for pruning and 
		//boolean useMinDEEPruningEw = (new Boolean((String)sParams.getValue("imindee", "false"))).booleanValue();
		double initIval = 0.0;
		double interval = 0;
		if(minSettings.doMinimize && !minSettings.minimizeBB)
			initIval = deeSettings.Ival;
		double difference = 0;

		if ((!mpiRun)&&((deeSettings.distrDACS)||deeSettings.distrDEE)){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}

		if (!minSettings.doMinimize) //no minimization
			minSettings.minimizeBB = false;
		if (!minSettings.minimizeBB) //not backbone minimization
			minSettings.doBackrubs = false;

		if (graphSettings.genInteractionGraph) //DACS is not performed when generating the interaction graph
			deeSettings.doDACS = false;

		//Setup the molecule system
		MolParameters mp = loadMolecule(sParams,curStrForMatrix,graphSettings.neighborList,graphSettings.distCutoff,true);

		boolean reload = false;
		if(minSettings.selectPerturbations){//Need to run the automatic perturbation selection
			//This only needs to be done once though: after that the perturbations can be read from pertFile
			selectPerturbations(mp, minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, ematSettings.addWTRot, sParams);
			reload = true;
		}


		// 2010: If useMinDEEPruningEw is set to false, this cycle never repeats itself.
		//  If it is set to true, it can repeat at most once: if none of the rotamer vectors
		//  between the conformation of lowest energy (i.e. lowestBound) 
		//  and lowestBound+InitEw can minimize to a lower energy than lowestBound+InitEw, 
		//   then let minimumEnergy be the minimum nergy found among the enumerated conformations,
		//   we set a new Ew equal to minimumEnergy - lowestBount and repeat this cycle.  We
		//   only have to do it at most twice.  
		long start;
		long totalEtime = 0;
		start = System.currentTimeMillis();	
		long totalConfsEvaluated = 0;

		RotamerSearch rs = new RotamerSearch(mp.m,mp.strandMut.numMutPos(), mp.strandsPresent, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, 
				minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts, deeSettings.useTriples, enumSettings.useFlagsAStar, es,hbonds, mp.strandMut);

		rs.useCCD = minSettings.useCCD;
		rs.partitionedRotamers = ibSettings.subRotamers;
//		rs.superRotamers = superRotamers;
//		rs.tuples = doTuples || readTuples;
//		rs.subRotamers = subRotamers;

		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(mp.strandPresent, sParams);
		int molStrand = 0;
		for(int resID:mp.strandMut.allMut){
			setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}
		
		if((!minSettings.selectPerturbations || reload) && !mp.loadedFromCache)
			rs.setupRCs(minSettings.doPerturbations);
		
		int numSplits = 0;

		//Load Energy Matrix
		long startEmat = System.currentTimeMillis();
		System.out.print("Loading precomputed energy matrix...");
		Emat emat = loadPairwiseEnergyMatrices(sParams,ematSettings.runNameEMatrixMin,minSettings.doMinimize,curStrForMatrix,es,mp.m,false);
		rs.setMinMatrix(emat);
		System.out.println("done");
		long endEmat = System.currentTimeMillis();
//		metrics.EmatTime += (endEmat-startEmat);

		
		int ctr=0;

		
		HashMap<ArrayList<Index3>,EnergyTuple> energyTuples=null;
		if(ibSettings.doTuples){
			Object o = loadObject("energyTuples.dat");
			if(o != null)
				energyTuples = (HashMap<ArrayList<Index3>,EnergyTuple>)o;
			else
				energyTuples = new HashMap<ArrayList<Index3>,EnergyTuple>();
			rs.energyTuples = energyTuples;
		}
		
		double bestE = enumSettings.bestE;
		int loopNum = 0;
		AStarResults asr;

		//KER: load lowest energy conformation found so far
		if(ibSettings.doTuples){
			ArrayList<FullConf> confs = new ArrayList<FullConf>();
			readRotResFile(outputSettings.outputConfInfo, confs, emat.allMutRes().size());

			for(FullConf conf : confs){
				if(conf.minE < bestE)
					bestE = conf.minE;
			}
			rs.bestEMin = bestE;
		}


		double lowestPairBound = Double.NEGATIVE_INFINITY;
		double actualLowestBound = Double.POSITIVE_INFINITY;
		//double Ival = initIval;
		double DEEIval = initIval;
		double ENUMIval = initIval;
		int numToEnumerate = 5;
		boolean runDEE = true;
		if(ibSettings.doTuples)
			runDEE = false;
		if(ibSettings.subRotamers)
			numToEnumerate=1;
		ArrayList<Index3> permanentlyPrune = new ArrayList<Index3>();
//		GurobiOptimization gurobiOpt = null; 



//		metrics.setLoopStart();
		while(true){
			
			loopNum++;

			if(ibSettings.superRotamers)
				DEEIval = initIval;
			

			int[] numRotForRes = compNumRotForRes(emat);
			int totalNumRot = 0;
			for(int q:numRotForRes){
				totalNumRot+=q;
				System.out.print(q+" ");
			};System.out.println("");

			if(ibSettings.subRotamers && loopNum >= 2){
				emat.unPrune();
			}

			int[] unprunedRotForRes = emat.remainingRot();
			for(int q:unprunedRotForRes){System.out.print(q+" ");}System.out.println();

			boolean localUseMinDEEPruningEw = minSettings.doMinimize && !minSettings.minimizeBB;
			if(deeSettings.doDACS) //If we are doing dacs we don't want to prune stuff too early so turn off 
				localUseMinDEEPruningEw = false;   //iMinDEE until the last depth is reached


			addErefAndEntropy(ematSettings.useEref, rs, emat);

			boolean deeDone = false;
			boolean supRotDone = false;
			double bestScore = Double.POSITIVE_INFINITY;

			int[] tmpPruned;
			int[] tmpPrunedPairs;

			int runNum = 0;


			/***runDEE*****/

			int timesRunFullPairs = -1; //-1 means run as much as possible
			int maxLoopNum = Integer.MAX_VALUE; //number of times to loop DEE
			boolean removeRot = false;
			if(ibSettings.superRotamers) //Don't spend that much time because we have to unprune anyway
				timesRunFullPairs = 0;
			else if(ibSettings.doTuples){ //Don't spend that much time because we have to unprune anyway
				timesRunFullPairs = 0;
				if(loopNum==1)
					runDEE = true;
			}
			else if(ibSettings.subRotamers){ //This is heuristic right now. You don't need to prune everytime since it
											 //just takes extra time. So only prune so often
				timesRunFullPairs = 0;
				maxLoopNum = 2;
				removeRot = true;
				//Every so often actually remove rotamers based on the current real Ival
				if(loopNum % 5 == 0 || loopNum == 2)
					runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
							deeSettings.initEw, deeSettings.maxIntScale, deeSettings.typeDep, DEEIval, 
							emat, localUseMinDEEPruningEw, removeRot,deeSettings.stericE, sParams,
							timesRunFullPairs,maxLoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);

				
				runDEE = true;
				removeRot = false;
				if(DEEIval < 1.0 && loopNum > 1){
					enumSettings.asMethod = Settings.ASTARMETHOD.BYSUBROT;
					numToEnumerate = Integer.MAX_VALUE;
					ENUMIval = Double.POSITIVE_INFINITY;
				}
				else
					DEEIval = 0;
				maxLoopNum = 1;
				timesRunFullPairs = 0;
			}



			if(runDEE)
				runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
						deeSettings.initEw, deeSettings.maxIntScale, deeSettings.typeDep, DEEIval, 
						emat, localUseMinDEEPruningEw, removeRot,deeSettings.stericE, sParams,
						timesRunFullPairs,maxLoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);


			//KER: If we are doing tuples we need to know what the lowest bound is without using tuples
			if(ibSettings.doTuples && loopNum == 1){
				rs.energyTuples = new HashMap<ArrayList<Index3>,EnergyTuple>();
				asr = rs.doAStarGMEC(outputSettings.outputConfInfo,true,minSettings.doMinimize,
						mp.strandMut.numMutPos(),mp.strandMut,0,
						bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
						ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
						localUseMinDEEPruningEw, 0,enumSettings); //Ival and initEw set to 0 to only find the lowest conformation
				totalConfsEvaluated += asr.numConfsEvaluated;
				lowestPairBound = asr.lowestBound;
				actualLowestBound = Math.min(asr.lowestBound, actualLowestBound);
				rs.energyTuples = energyTuples; 
			}





//			if(asMethod == ASTARMETHOD.LPGUROBI && gurobiOpt == null){
//				gurobiOpt = new GurobiOptimization(emat,true);
//
//				//Need to add the tuples in the right order so the parents will already
//				//be around when the children are added
//				//KER: This is super inefficient, but I just want something that works right now
//				int length = 2;
//				int numRotAtLength=0;
//				while(length==2 || numRotAtLength > 0){
//					numRotAtLength = 0;
//					for(EnergyTuple tuple: energyTuples.values()){
//						if(tuple.rots.length == length){
//							//KER: also need to make sure that rotamer in tuple isn't pruned 
//							//KER: this will happen if we are loading tuples from a previous run
//							boolean validTuple = true;
//							for(Index3  rot: tuple.rots){
//								if(emat.getSinglePruned(rot))
//									validTuple = false;
//							}
//							if(validTuple){
//								gurobiOpt.addTuple(tuple,emat);
//								numRotAtLength++;
//							}
//						}
//					}
//					length++;
//				}
//			}

			//Run Enumeration
			asr = rs.doAStarGMEC(outputSettings.outputConfInfo,true,minSettings.doMinimize,
					mp.strandMut.numMutPos(),mp.strandMut,deeSettings.initEw,
					bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
					ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
					localUseMinDEEPruningEw, ENUMIval,enumSettings);

			actualLowestBound = Math.min(asr.lowestBound, actualLowestBound);
			totalConfsEvaluated += asr.numConfsEvaluated;
			bestE = Math.min(bestE,rs.getBestE());

			if(ibSettings.superRotamers){
				interval = bestE - asr.lowestBound;
				difference = interval - DEEIval;
			}
			else if(ibSettings.doTuples){
				interval = asr.lastBound - lowestPairBound;
				difference = DEEIval-interval;
			}
			else if(ibSettings.subRotamers){
				interval = bestE - asr.lowestBound;
				difference = interval;
			}
			else{ //Normal run
				System.out.println("DON'T KNOW WHAT TO DO!!");
				System.out.println("Please run doDEE instead!!");
				System.exit(0);
			}

			if(ibSettings.superRotamers){
				DEEIval = interval;
				ENUMIval = interval;
			}
			if(ibSettings.subRotamers){
				DEEIval = interval;
			}

			System.out.println("Ival difference: "+difference);

			if(difference < 0.2 || bestE - asr.lastBound < 0.001){
				if(ibSettings.doTuples && DEEIval < bestE-lowestPairBound){
					double diff = bestE-lowestPairBound; 
					DEEIval = Math.min(diff, DEEIval*2);
					if(DEEIval == 0)
						DEEIval = 0.5;
					System.out.println("New DEEIval: "+DEEIval);
					System.out.println("CHANGING DEEIVAL");
					emat.unPrune();
					rs.MSAStarSearch = null;
//					gurobiOpt = null;
					runDEE = true;

				}
				else
					break;
			}else{
				runDEE = false; //We haven't changed the DEEIval so don't prune any more 
			}



			/*****    Post processing and update for next DEE/Enumuration ****/
			//If we're using superRotamers
			if(ibSettings.superRotamers){
				//Resetting the energy matrix
				emat.unPrune();
				//runDEE with Ival = I1

				runDEE(deeSettings.useFlags, minSettings.doMinimize, minSettings.minimizeBB, deeSettings.scaleInt,
						deeSettings.initEw, deeSettings.maxIntScale, deeSettings.typeDep, DEEIval, 
						emat, localUseMinDEEPruningEw, true,deeSettings.stericE, sParams, 
						-1,maxLoopNum,deeSettings.deeSettings,rs.strandRot,mp.strandMut,rs.m,rs.doPerturbations);

				emat.removePrunedRotReducedMem(false);

				//Save a temporary superrotamer emat 
				emat.save(ematSettings.runNameEMatrixMin+"_COM.dat.tmp",mp.m);
			

				long energyStartTime = System.currentTimeMillis();
				if(difference > 0.2){
					//Contract Rotamers
					unprunedRotForRes = emat.remainingRot();
					for(int q:unprunedRotForRes){System.out.print(q+" ");}System.out.println();
					//if(runNum>=1 && runNum <=7){
					for(int i=0; i<1;i++){
						int[] pos = new int[2];
						switch(ibSettings.contractMethod){
						case LEASTPAIRS:
							pos = emat.findPositionsToContractLeastPairs();
							break;
						case PERCENTPRUNED:
							pos = emat.findPositionsToContract();
							break;
//						case CLOSESTRES:
//							pos = emat.findPositionsToContractClosestRes(rs.contractPos,mp.m);
//							break;
						case PERCENTLEAST: 
							pos = emat.findPositionsToContractPercentLeast();
							if((pos[1] == -1 || pos[0]== -1) && emat.resByPos.size() > 1){
								pos[0] = 0;
								pos[1] = 1;
							}
							break;
//						case LARGESTDIFF:
//							if(rs.contractPos < rs.contractPos2){
//								pos[0] = rs.contractPos;
//								pos[1] = rs.contractPos2;
//							}else{
//								pos[1] = rs.contractPos;
//								pos[0] = rs.contractPos2;
//							}
//							break;
						default:
							pos[0] = -1;
							pos[1] = -1;
							break;
						}
	
	
						if(pos[0] != -1 && pos[1] != -1){
							System.out.println("Combining Pos: "+pos[0]+", "+pos[1]);
							supRotDone = false;
							emat.combinePos(pos[0], pos[1]);
							long startE = System.currentTimeMillis();
							if(minSettings.doMinimize){
								recalculateEnergies(emat, pos[0],curStrForMatrix,sParams,mp.m,
										minSettings.doMinimize,minSettings.minScheme,mp.strandMut,mp.strandsPresent,es,
										minSettings.doPerturbations, minSettings.pertFile, minSettings.minimizePerts);
								if(ematSettings.useEref)
									rs.addEref(emat.eRef, pos[0]);
								if(EnvironmentVars.useEntropy)
									rs.addEntropyTerm(pos[0]);
	
								emat.save(ematSettings.runNameEMatrixMin+"_COM.dat.tmp",mp.m);
							}
	
						}
						else{
							System.out.println("Couldn't find any positions to combine.");
						}
					}
				}
			}
			else if(ibSettings.doTuples &&  rs.worstRots != null){
				System.out.println("Tuples not implemented yet.");
				System.exit(0);;
//				addTuples(sParams, curStrForMatrix, templateAlwaysOn, mp, rs,
//						emat, useEref, energyTuples, bestE, asr, DEEIval,
//						gurobiOpt,addTuplesByDistance);
				
			}//}
			else if(ibSettings.subRotamers){ //Split only the worst rotamers
				//KER: Now we can only rely on RotamerSearch storing the conformations it found
				//Here we have to figure out the worst rotamers

				numSplits = partitionRotamers(sParams, ematSettings.runNameEMatrixMin,
						ematSettings.useEref, mp, rs, numSplits, emat,
						bestE,curStrForMatrix,minSettings,mp.strandMut,es);
			}
			else if(ibSettings.readTuples){
				System.out.println("Reading tuples is not currently implemented.");
				System.exit(0);
//				if(loopNum == 1){
//					ArrayList<ArrayList<Index3>> tuples = readTuples(ibSettings.tupleFile, mp.m, emat);
//
//					//For each tuple calculate the EnergyTuple for it
//					for(ArrayList<Index3> tuple: tuples){
//
//						//Make sure that there is a breadcrumb to the tuple
//						Iterator<Index3> iter = tuple.iterator();
//						ArrayList<Index3> curTuple = new ArrayList<Index3>();
//						ArrayList<Index3> parent = new ArrayList<Index3>();
//						Index3 i3 = iter.next();
//						curTuple.add(i3);
//						parent.add(i3);
//						while(iter.hasNext()){
//							i3 = iter.next();
//							curTuple.add(i3);
//							EnergyTuple child;
//							if(!energyTuples.containsKey(curTuple)){
//								ArrayList<Index3> computeTuple = new ArrayList<Index3>();
//								for(Index3 i2:curTuple)
//									computeTuple.add(i2);
//								child = recalculateEnergyTuple(emat, computeTuple, curStrForMatrix, sParams,energyTuples,mp.m,useEref,templateAlwaysOn);
//							}
//							else{
//								child = energyTuples.get(curTuple);
//							}
//							if(parent.size() >=2){
//								EnergyTuple parentTuple = energyTuples.get(parent);
//								parentTuple.children.add(child);
//							}
//
//							parent.add(i3);
//						}
//
//
//					}
//				}
//
//				//Save energyTuples
//				outputObject(energyTuples,"energyTuples.dat");

			}
			else{
				//Run A*
				asr = rs.doAStarGMEC(outputSettings.outputConfInfo,true,minSettings.doMinimize,
						mp.strandMut.numMutPos(),mp.strandMut,deeSettings.initEw,
						bestScore,null,enumSettings.approxMinGMEC,enumSettings.lambda,minSettings.minimizeBB,
						ematSettings.useEref,minSettings.doBackrubs,minSettings.backrubFile,
						localUseMinDEEPruningEw, ENUMIval,enumSettings);
				totalConfsEvaluated += asr.numConfsEvaluated;

				break;
			}
			long endEnergyTime = System.currentTimeMillis();
//			metrics.Etime += (endEnergyTime - energyStartTime);

			ctr++;
		}


//		metrics.setEndTime();

		if(ibSettings.doTuples){
			System.out.println("NUMTUPLES: "+ energyTuples.size());
			for(ArrayList<Index3> tuples:rs.energyTuples.keySet()){
				for(Index3 index: tuples){
					System.out.print("("+index.pos+","+index.aa+","+index.rot+"), ");
				}
				System.out.println("");
			}
		}
		
//		metrics.trueIval = bestE - actualLowestBound;
//		metrics.print();

		//		System.out.println("EmatTime: "+EmatTime);
		//		System.out.println("DEETime: "+DEETime);
		//		System.out.println("ASTime: "+ASTime);
		//		System.out.println("ENERGYtime: "+EnergyTime);
		//		System.out.println("TotalTime: "+time);


		

		//		System.out.println("TOTALNUMCONFS: "+totalConfsEvaluated);

	}

	private void recalculateEnergies(Emat emat, int mutPos,int curStrForMatrix,ParamSet params,
			Molecule m, boolean doMinimize, RotamerSearch.MINIMIZATIONSCHEME minScheme, MutableResParams strandMut,
			int strandsPresent, EPICSettings es, boolean doPerturbations, String pertFile, boolean minimizePerts) {
		//Emat ematToSend = new Emat(emat, mutPos);
		//KER: split this into several small runs for each AA vs. Pos


		//KER: plus 1 is for the SHL-AS
		OneMutation[] mutArray = new OneMutation[(emat.singles.E[mutPos].length)*(emat.singles.E.length-1) + 1];
		int ctr=0;
		//Loop through mutPosAAs
		for(int i=0; i<emat.singles.E.length;i++){
			if(i==mutPos){
				mutArray[ctr] = new OneMutation();
				mutArray[ctr].mutNum = mutPos;
				mutArray[ctr].resMut = new int[emat.singles.E.length];
				for(int j=0; j<mutArray[ctr].resMut.length;j++)
					mutArray[ctr].resMut[j] = 0;

				mutArray[ctr].resMut[mutPos] = 1;
				mutArray[ctr].runParams = new EmatCalcParams(mutPos, -1);

				mutArray[ctr].flagMutType = "SHL-AS";
				ctr++;
			}
			else{
				TreeSet<Integer> tmpAA2 = new TreeSet<Integer>();
				for(int aa=0; aa<emat.singles.E[i].length;aa++){
					tmpAA2.add(aa);
				}
				for(int aa=0; aa<emat.singles.E[mutPos].length;aa++){
					mutArray[ctr] = new OneMutation();
					mutArray[ctr].mutNum = mutPos;
					mutArray[ctr].resMut = new int[emat.singles.E.length];
					for(int j=0; j<mutArray[ctr].resMut.length;j++)
						mutArray[ctr].resMut[j] = 0;

					mutArray[ctr].resMut[mutPos] = 1;

					mutArray[ctr].flagMutType = "AS-AS";
					mutArray[ctr].resMut[i] = 1;

					TreeSet<Integer> tmpAA1 = new TreeSet<Integer>();
					tmpAA1.add(aa);

					mutArray[ctr].runParams = new EmatCalcParams(mutPos, i, tmpAA1, tmpAA2);

					ctr++;
				}
			}

		}

		MutationManager mutMan = new MutationManager(null, mutArray, true);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setcurStrForMatrix(curStrForMatrix);
		mutMan.setParams(params);
		mutMan.setMolecule(m);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setES(es);
		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setPertFile(pertFile);
		mutMan.setMinimizePerts(minimizePerts);

		try {
			handleDoMPIMaster(mutMan, mutArray.length);
		} catch (Exception e) {
			e.printStackTrace();
		}


	}
	

	private int partitionRotamers(ParamSet sParams, String runNameEMatrixMin,
			boolean useEref, MolParameters mp,
			RotamerSearch rs, int numSplits, Emat emat, double bestE, int curStrForMatrix,
			Settings.Minimization minSettings,MutableResParams strandMut, EPICSettings es) {
		
		ArrayList<Index3> rotWasSplit = new ArrayList<Index3>();
		ArrayList<ArrayList<Index3>> rotWasSplitTo = new ArrayList<ArrayList<Index3>>();
		//For each conformation we found check rotamers to split	
		ArrayList<Index3> rotsAddedThisRound = new ArrayList<Index3>();
		for(GlobalRCConf conformation: rs.generatedConfs){
			ArrayList<Index3> rotToRemove = new ArrayList<Index3>();

			//KER: For this conformation setup the differences between the bound and actual energy
			Index3[] rotamers = new Index3[conformation.resNum.length];
			ArrayList<ArrayList<Index3>> rotsPerPos = new ArrayList<ArrayList<Index3>>();
			
			//KER: First we must find the corresponding rotamer for each conformation rotamer
			//KER: Note that if we split a rotamer we will no longer be able to find it cause it'll be gone
			for(int i=0; i<rotamers.length;i++){
				rotamers[i] = emat.getGlobalRot(conformation.resNum[i],conformation.globalRCs[i]);
				ArrayList<Index3> tmpRot = new ArrayList<Index3>();
				tmpRot.add(rotamers[i]);
				rotsPerPos.add(tmpRot);
			}

			//KER: Now we must find the bound and then the differences
			double[] bound = new double[conformation.EforRes.length]; 
			double diff[] = new double[conformation.EforRes.length];
			
			for(int i=0; i<bound.length;i++){
				bound[i] += emat.getSingleMinE(rotamers[i]);
				for(int j=i+1; j<bound.length;j++){
					if(emat.areNeighbors(i, j)){
						try{
						double energy = emat.getPairMinE(rotamers[i], rotamers[j]);
						bound[i] += energy/2;
						bound[j] += energy/2;
						}catch(Exception E){
							E.printStackTrace();
						}
					}
				}
			}
			//KER: calculate the differences
			double oldDiff = 0;
			double totalE = emat.templ_E;
			for(int i=0; i<diff.length;i++){
				diff[i] = conformation.EforRes[i] - bound[i];
				oldDiff += diff[i];
				totalE += conformation.EforRes[i];
			}

			//KER: sort the max indices that we are going to split
			ArrayIndexComparator comparator = new ArrayIndexComparator(diff);
			Integer[] maxIndexes = comparator.createIndexArray();
			Arrays.sort(maxIndexes,comparator);
			
			double cumulativeDiff = 0;
			for(int pos=0; pos<rotamers.length;pos++){
				cumulativeDiff += diff[maxIndexes[pos]];
				Index3 rot = rotamers[maxIndexes[pos]];
				if(!rotsAddedThisRound.contains(rot)){
					rotsAddedThisRound.add(rot);

					ArrayList<ArrayList<ResidueConformation>> rotsToAdd = new ArrayList<ArrayList<ResidueConformation>>();
					ArrayList<ArrayList<ResidueConformation>> rotsToUpdate = new ArrayList<ArrayList<ResidueConformation>>();
					//KER: Assuming only one rotamer per position for now
					ArrayList<ResidueConformation> curRCs = emat.getRCs(mp.m, rot);
					//KER: For each rotamer make several sub-rotamers in Chi1
					int resID = emat.resByPos.get(rot.pos).get(0);
					String pdbNum = mp.m.residue[resID].getResNumberString();
					ResidueConformation curRC = curRCs.get(0);
					if(curRC.rot.values == null || curRC.rot.values.length == 0)
						continue;

					//Just split the rotamers in half for now
					//Could actually do more sophisticated splits in the future
					int dihedDepth = 0;
					for(int i=0; i<curRC.rot.minimizationWidth.length-1;i++)
						if(curRC.rot.minimizationWidth[i+1] > curRC.rot.minimizationWidth[i]){
							dihedDepth = i+1;
						}

					double rotInterval = curRC.rot.minimizationWidth[dihedDepth]/2;
					for(int i=0;i<2;i++){
						double[] dihedVals = new double[curRC.rot.values.length];
						double[] minimizationWidth = new double[curRC.rot.minimizationWidth.length];
						for(int q=0; q<dihedVals.length;q++){
							dihedVals[q] = curRC.rot.values[q];
							minimizationWidth[q] = curRC.rot.minimizationWidth[q];
						}
						int dir = i*2-1;
						dihedVals[dihedDepth] += dir*rotInterval;
						minimizationWidth[dihedDepth] = rotInterval;
						Rotamer newSubRot = null;
						ResidueConformation newSubResConf = null;
						ArrayList<ResidueConformation> rotToAdd = new ArrayList<ResidueConformation>();
						ArrayList<ResidueConformation> rotToUpdate = new ArrayList<ResidueConformation>();
						if(i != 0){ //Update (don't add) orig rot 
							newSubRot = mp.m.residue[resID].rl.addSubRotamer(curRC.rot.aaType.name, pdbNum,dihedVals,minimizationWidth,curRC.rot);
							newSubResConf = mp.m.strand[mp.m.residue[resID].strandNumber].rcl.addSubResidueConformation(newSubRot, curRC.pertState, curRC.res, curRC);
							rotToAdd.add(newSubResConf);
							rotsToAdd.add(rotToAdd);
						}
						else{
							//KER: need to create a new rotamer because allowed minimization will be different
							newSubRot = mp.m.residue[resID].rl.addSubRotamer(curRC.rot.aaType.name, pdbNum,dihedVals,minimizationWidth,curRC.rot);
							newSubResConf = mp.m.strand[mp.m.residue[resID].strandNumber].rcl.addSubResidueConformation(newSubRot, curRC.pertState, curRC.res, curRC);
							rotToUpdate.add(newSubResConf);
							emat.updateRotamer(mp.m,rot,rotToUpdate,false);
						}

					}


					ArrayList<Index3> rotamerIndices = emat.addRotamers(mp.m,rot.pos,rotsToAdd,false);
					rotamerIndices.add(rot);
					calculateRotamerEnergies(emat, sParams, rotamerIndices, rot.pos, mp.m,useEref,curStrForMatrix,
							minSettings.minimizeBB, minSettings.doMinimize,minSettings.minScheme,strandMut,
							mp.strandsPresent, es, minSettings.doPerturbations,minSettings.pertFile,minSettings.minimizePerts);
					numSplits++;
					
					rotWasSplit.add(rot);
					rotWasSplitTo.add(rotamerIndices);
					
					
					//Check the new energy for the new structures
					rotsPerPos.set(rot.pos, rotamerIndices);
					
					
					double bestNewE = getBestNewE(-1, new Index3[rotsPerPos.size()],rotsPerPos,Double.POSITIVE_INFINITY,rs);
					
					
					double newDiff = totalE - bestNewE;
					
					if((cumulativeDiff/oldDiff) > .75 || ((oldDiff-newDiff)/oldDiff)>0.3 || bestNewE > bestE )
						break;
					
				}
				
			}

		}

		return numSplits;
	}


	private void calculateRotamerEnergies(Emat emat, ParamSet sParams, ArrayList<Index3> rotamers, int mutPos,
			Molecule m, boolean useEref, int curStrForMatrix,
			boolean minimizeBB, boolean doMinimize, RotamerSearch.MINIMIZATIONSCHEME minScheme, MutableResParams strandMut,
			int strandsPresent, EPICSettings es, boolean doPerturbations, String pertFile, boolean minimizePerts) {

		//KER: plus 1 is for the SHL-AS
		OneMutation[] mutArray = new OneMutation[emat.singles.E.length];
		int ctr=0;
		//Loop through mutPosAAs
		for(int i=0; i<emat.singles.E.length;i++){
			if(i==mutPos){
				mutArray[ctr] = new OneMutation();
				mutArray[ctr].mutNum = mutPos;
				mutArray[ctr].resMut = new int[emat.singles.E.length];
				for(int j=0; j<mutArray[ctr].resMut.length;j++)
					mutArray[ctr].resMut[j] = 0;

				mutArray[ctr].resMut[mutPos] = 1;
				mutArray[ctr].runParams = new EmatCalcParams(mutPos, -1,rotamers);

				mutArray[ctr].flagMutType = "SHL-AS";
				ctr++;
			}
			else{
				mutArray[ctr] = new OneMutation();
				mutArray[ctr].mutNum = mutPos;
				mutArray[ctr].resMut = new int[emat.singles.E.length];
				for(int j=0; j<mutArray[ctr].resMut.length;j++)
					mutArray[ctr].resMut[j] = 0;

				mutArray[ctr].resMut[mutPos] = 1;

				mutArray[ctr].flagMutType = "AS-AS";
				mutArray[ctr].resMut[i] = 1;

				mutArray[ctr].runParams = new EmatCalcParams(mutPos, i,rotamers);

				ctr++;
			}

		}
		
		MutationManager mutMan = new MutationManager(null, mutArray, true);
		mutMan.setPairEMatrixMin(emat);
		mutMan.setcurStrForMatrix(curStrForMatrix);
		mutMan.setParams(sParams);
		mutMan.setMolecule(m);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinScheme(minScheme);
		mutMan.setStrandMut(strandMut);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setStrandsPresent(strandsPresent);
		mutMan.setES(es);
		mutMan.setDoPerturbations(doPerturbations);
		mutMan.setPertFile(pertFile);
		mutMan.setMinimizePerts(minimizePerts);
		
		
		
		
		try{
			if(MPItoThread.Rank() == 0){
			
				try {
					handleDoMPIMaster(mutMan, mutArray.length);
				} catch (Exception e) {
					e.printStackTrace();
				}
			
			}else{ //If this is a slave, the slave has to do the work
				for(int i=0; i<mutArray.length;i++){
					CommucObj cObj = handleComputeAllPairwiseRotamerEnergiesSlave(mutMan.getNextComObj(i));
					mutMan.processFinishedMutation(cObj);
				}
			}
		}catch (Exception e){
			System.out.println("Couldn't calculate partitioned rotamer energies");
			e.printStackTrace();
		}

		//KER: make sure the eRef is set correctly
		if(useEref){
			RotamerSearch.addEref(emat, m, emat.eRef, mutPos,rotamers);
		}
		if(EnvironmentVars.useEntropy)
			RotamerSearch.addEntropyTerm(emat, m, mutPos,rotamers);


	}
	
	/**
	 * Used for partitioned rotamer calculation
	 * @param i
	 * @param conf
	 * @param rotsPerPos
	 * @param bestE
	 * @param rs
	 * @return
	 */
	private double getBestNewE(int i, Index3[] conf, ArrayList<ArrayList<Index3>> rotsPerPos,
			double bestE,RotamerSearch rs) {
		i = i+1;
		if(i == conf.length){
			double curE = rs.computeBestRotEnergyBound(conf);
				
			if(curE < bestE)
				bestE = curE;
			return bestE;
		}
		
		
		for(Index3 rot:rotsPerPos.get(i)){
			conf[i] = rot;
			bestE = getBestNewE(i,conf,rotsPerPos,bestE,rs);
		}
		
		return bestE;
	}
	
	
	/**
	 * Compute new tuples and add them to the ones calculated so far
	 * @param sParams
	 * @param curStrForMatrix
	 * @param templateAlwaysOn
	 * @param mp
	 * @param rs
	 * @param emat
	 * @param useEref
	 * @param energyTuples
	 * @param bestE
	 * @param asr
	 * @param DEEIval
	 * @param gurobiOpt
	 * @param useDistance
	 * @return
	 */
//	private boolean addTuples(ParamSet sParams, int curStrForMatrix,
//			boolean templateAlwaysOn, MolParameters mp, RotamerSearch rs,
//			Emat emat, boolean useEref,
//			HashMap<ArrayList<Index3>, EnergyTuple> energyTuples, double bestE,
//			AStarResults asr, double DEEIval, GurobiOptimization gurobiOpt,boolean useDistance) {
//		boolean addedTuple = false;
//		boolean supRotDone;
//		HashMap<ArrayList<Index3>,Boolean> tmpHash = new HashMap<ArrayList<Index3>,Boolean>();
//		for(ArrayList<Index3wVal> rotList: rs.worstRots){
//			LinkedList<Index3> tuple = new LinkedList<Index3>();
//			Iterator<Index3wVal> iter = rotList.iterator();
//
//			int[] AAnums = new int[emat.allMutRes().size()];
//			int[] ROTnums = new int[AAnums.length];
//			//double[] boundPerPos = new double[mp.strandMut.allMut.length];
//			for(Index3wVal i3wV: rotList){
//				for(Index3 i3: i3wV.i3s){
//					AAnums[i3.pos] = i3.aa;
//					ROTnums[i3.pos] = i3.rot;
//				}
//			}
//
//			ArrayList<EnergyTuple> parents = new ArrayList<EnergyTuple>();
//
//			
//			Index3wVal nextTuple = iter.next();
//			if(rs.energyTuples.containsKey(nextTuple.i3s))
//				parents.add(rs.energyTuples.get(nextTuple.i3s));
//			tuple.addAll(nextTuple.i3s);
//
//			//KER: There's a problem with finding tuples when we combine tuples
//			//KER: If we combine singles it's fine cause there is a breadcrumb to find the tuple
//			//KER: My current fix is to just always break up the tuple we are combining and add each
//			//KER: individually. 
//			boolean done = false;
//			boolean addNextRot = false;
//			while((!done || addNextRot) && iter.hasNext()){
//
//				
//				
//				if(tuple.size() >= 2 && useDistance){
//					//Mutate the residues and rotamers and then find which one is the closest
//					//Would only need to do this once, but here for testing
//					Index3 i3s[] =  new Index3[AAnums.length]; 
//					for(Index3wVal i3wV: rotList){
//						for(Index3 i3: i3wV.i3s){
//							i3s[i3.pos] = i3;
//							RotamerEntry re =  emat.singles.getTerm(i3);
//							re.applyMutation(mp.m, emat.resByPos, true, true);
//							re.applyRotamer(emat.resByPos, mp.m);
//						}
//					}
//					
//					double distance[] =  new double[AAnums.length];
//					ArrayList<Integer> tupleRes = new ArrayList<Integer>();
//					for(Index3 i3: tuple){
//						distance[i3.pos] = Double.POSITIVE_INFINITY;
//						tupleRes.addAll(emat.resByPos.get(i3.pos));
//					}
//					for(int i=0; i<distance.length;i++){
//						if(distance[i] < Double.POSITIVE_INFINITY)
//							distance[i] = mp.m.minAvgDist(tupleRes, emat.resByPos.get(i));
//					}
//					double minDist = Double.POSITIVE_INFINITY;
//					int minIndex = -1;
////					System.out.println("Distances: ");
//					for(int i=0; i<distance.length;i++){
////						System.out.print(distance[i]+" ");
//						if(distance[i] < minDist){
//							minDist = distance[i];
//							minIndex = i;
//						}
//					}
//				
////					if(minIndex == -1){
////						System.out.print("Tuple includes: ");
////						for(Index3 i3: tuple){
////							System.out.print(i3.pos+" ");
////						}
////						System.out.println("");
////						mp.m.saveMolecule("brokenDistance.pdb", 0.0f);
////						for(int i : emat.allMutRes()){
////							System.out.print(mp.m.residue[i].name+" ");
////						}
////						System.out.println();
////					}
//					
//					ArrayList<Index3> nextI3 = new ArrayList<Index3>();
//					
//					nextI3.add(i3s[minIndex]);
//					
//					nextTuple = new Index3wVal(nextI3,0.0);
//					
//				}else{
//					nextTuple = iter.next();
//				}
//				double E = 0;
//				if(rs.energyTuples.containsKey(nextTuple.i3s))
//					parents.add(rs.energyTuples.get(nextTuple.i3s));
//
//				tuple.addAll(nextTuple.i3s);
//				Collections.sort(tuple); //Sorting based on position
//				addNextRot = false;
//
//
//				ArrayList<LinkedList<EnergyTuple>> tupleOptions = findTupleOptions(energyTuples,AAnums,ROTnums);
//				E = Double.NEGATIVE_INFINITY;
//				if(tupleOptions.size() <= 0)
//					E = rs.computeBestRotEnergyBoundWTuples(AAnums, ROTnums,null,null);
//				else{
//					for(LinkedList<EnergyTuple> tuples: tupleOptions){
//						double tmpE = rs.computeBestRotEnergyBoundWTuples(AAnums, ROTnums,tuples,null);
//						if(tmpE > E){
//							E = tmpE;
//						}
//					}
//				}
//
//				outPS.println("Old bound for conf: "+E);
//
//
//				if(/*!tmpHash.containsKey(tuple) && */rs.energyTuples.containsKey(tuple)){
//					//KER: If the child isn't the parents' child then set it
//					EnergyTuple newChild = rs.energyTuples.get(tuple);
//
//					for(EnergyTuple parent : parents){
//						boolean hasChild = false;
//						for(EnergyTuple child: parent.children){
//							if(child.equals(newChild)){
//								hasChild = true;
//								break;
//							}
//						}
//						if(!hasChild)
//							parent.children.add(newChild);
//
//					}
//
//					if(!tmpHash.containsKey(tuple) && tuple.size() < AAnums.length)
//						addNextRot = true;
//					
//
//					parents = new ArrayList<EnergyTuple>();
//					parents.add(newChild);
//				}
//				else{
//
//					/*while(!tmpHash.containsKey(tuple) && rs.energyTuples.containsKey(tuple)){
//								tuple.addAll(iter.next().i3s);
//								Collections.sort(tuple);
//								}*/
//
//					supRotDone = false;
//
//
//
//					//if(rs.worstRot1.pos < rs.worstRot2.pos)
//					//	recalculateEnergyTuple(emat, rs.worstRot1, rs.worstRot2, curStrForMatrix, sParams,energyTuples);
//					//else
//
//
//
//					if(!tmpHash.containsKey(tuple)){
//
//
//						ArrayList<Index3> computeTuple = new ArrayList<Index3>();
//						for(Index3 i2:tuple)
//							computeTuple.add(i2);
//
//
//						EnergyTuple child = recalculateEnergyTuple(emat, computeTuple, curStrForMatrix, sParams,energyTuples,mp.m,useEref,templateAlwaysOn);
//						for(EnergyTuple parent : parents)
//							parent.children.add(child);
//
//
//						//KER: For Gurobi Optimization add constraints and variables for new tuple
//						if(gurobiOpt != null){
//							gurobiOpt.addTuple(child,emat);
//						}
//
//						//Reset parents for next run
//						parents = new ArrayList<EnergyTuple>();
//						parents.add(child);
//
//						tmpHash.put(computeTuple, true);
//						outPS.print("Adding tuple: ");
//						for(Index3 i3:tuple){
//							outPS.print("("+i3.pos+", "+i3.aa+", "+i3.rot+"), ");
//						}
//						outPS.println("");
//						addedTuple = true;
//
//						tupleOptions = findTupleOptions(energyTuples,AAnums,ROTnums);
//						E = Double.NEGATIVE_INFINITY;
//						if(tupleOptions.size() <= 0)
//							E = rs.computeBestRotEnergyBoundWTuples(AAnums, ROTnums,null,null);
//						else{
//							for(LinkedList<EnergyTuple> tuples: tupleOptions){
//								double tmpE = rs.computeBestRotEnergyBoundWTuples(AAnums, ROTnums,tuples,null);
//								if(tmpE > E){
//									E = tmpE;
//								}
//							}
//						}
//						outPS.println("New bound for conf: "+E);
//
//
//
//
//
//
//					}
//				}
//
//				if((E < (asr.lastBound+DEEIval) && E<=bestE) && iter.hasNext() && tuple.size() < AAnums.length){
//					//Add another residue if the bound is too low still
//					//tuple.addAll(iter.next().i3s);
//					//Collections.sort(tuple);
//					done = false;
//				}else if(!addNextRot){
//					done = true;
//				}
//
//			}
//
//		}
//
//		//Remake the A* tree instead of just getting rid of it
//		//rs.MSAStarSearch.remakeTree(bestE);
//
//		//Save energyTuples
//		try{
//			if(MPItoThread.Rank() == 0)
//				outputObject(energyTuples,"energyTuples.dat");	
//		}catch(Exception E){
//			System.out.println("Could not output energyTuples.dat");
//			E.printStackTrace();
//		}
//		
//		return addedTuple;
//	}
//	

} // end of KSParser class
