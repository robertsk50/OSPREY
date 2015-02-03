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
// RotamerSearch.java
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
 * This class provides a variety of tools and search algorithms for
 *  doing rotamer-based searching over molecular conformations
 *
 * The system consists of one molecule containing: a protein strand, 
 *  a ligand strand (optional), and a cofactor strand (optional).
 * The protein strand does not have to contain sequential residues
 *  but it must be made of standard amino acids
 * The ligand strand can only be one 'thing'
 *  -if this 'thing' is an AA then the Penultimate Rotamer library is used
 *  -if this 'thing' is not an AA then a generic rotamer library is used
 *
 */

import cern.colt.matrix.DoubleFactory1D;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.*;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

/**
 * 
 * This class provides a variety of tools and search algorithms for
 *  doing rotamer-based searching over molecular conformations.
 * Contains functions for computing the pairwise energy matrices and
 *  for performing DEE/A* and K* (with different types of minimization)
 *  mutation searches.
 * The functions in this class are typically called after the necessary
 *  setup in KSParser.java is performed.
 *
 */
public class RotamerSearch implements Serializable
{

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = true;
	/*final*/ static double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)

	Molecule m;		// the molecule
	Amber96ext a96ff;	// the forcefield and energy function to use for energy evaluation etc...
	Amber96ext a96ffmin; // the forcefield that is used for the minimizer (might be the same as a96ff unless we want template on during minimization)

	//for the matrix, only for the minimization
	CCDMinimizer ccdMin;//CCD-based minimizer: supports sidechains only at the moment but will do other stuff soon
	boolean useCCD = true;//indicates use of ccdMin
	EnergyFunction efunc;//Energy function for ccdMin

	SimpleMinimizer simpMin;	// the simple energy minimizer (side-chains)
	BBMinimizer bbMin = null;	//the backbone minimizer
	BackrubMinimizer brMin = null;	//the backrub minimizer
	//KER: eliminatedRotAtRes is now a 3D matrix for easier indexing (I don't think this will save space like the 6D array)
	//	PrunedRotamers<Boolean> eliminatedRotAtRes = null;	// rotamers pruned by MinDEE
	//KER: splitFlags will now be a 6D array to save space
	//	boolean splitFlags[][][][][][] = null;

	boolean tripleFlags[][][][][][][][][] = null;//Pruned triples
	boolean useTriples = false;             //Indicates whether triples should be used
	boolean useFlagsAStar = false;          //Should split and potentially triple flags be used in A*?

	PrunedRotamers<Boolean> prunedIsSteric = null;	// pruned rotamers due to unallowed sterics
	double indIntMinDEE[] = null;	//the single-residue interval term in the MinDEE criterion
	double pairIntMinDEE[] = null;	//the pairwise interval term in the MinDEE criterion
	boolean repeatSearch = false;	// determines if the search must be repeated to achieve the desired accuracy
	//	int curConf[] = null;		// the current conformation returned by A*
	boolean allPruned = false;	// determines if all the rotamers for a given residue have been pruned by MinDEE;
	// sends this information to the master node in the mutation search

	final double stericE = (double)Math.pow(10,38);	// the energy stored for an unallowed steric
	private double Ec_const = stericE;	// the minimum lower energy bound for a pruned conformation
	private double boundForPartition = stericE; //a lower bound on the conformational energy for a given partition (used by DACS)

	final int samplesBB = 1000;//number samples for intra-rotamer energy matrix computation with backbone flexibility

	boolean distDepDielect = true; //distance-dependent dielectric
	double dielectConst = 1.0; //the dielectric constant

	boolean doDihedE = false; 		// if true dihedral energies are computed and used
	//  during energy minimization. Note that if
	//  energy minimization is NOT used then dihedral
	//  energies are not explicitly computed. In
	//  reality the total energy values are the same
	//  because although we're using AMBER dihedral
	//  energy terms we assume that each dihedral of
	//  each rotamer is at the bottom of an energy
	//  well. Thus without minimization the total
	//  dihedral energy is zero.
	boolean doSolvationE = false; //determines if solvation energies should be computed

	double solvScale = 1.0; //the solvation energies scaling factor

	PriorityQueue<ConfPair> topConfs = null; //Stores the top conformations of the run
	//Right now this is used to print out the conformations at the end

	//RotamerLibrary rl = null; //the standard rotamer library for the protein (the AA library)
	//RotamerLibrary[] grl = null; //the rotamer library for the ligand (could be the AA or non-AA library)
	//StrandRotamers sysLR = null;// the rotamers object for the system strand
	//StrandRotamers ligROT = null;	// the rotamers object for the ligand strand
	//int sysStrNum = -1;	// the strand number of the system
	//int ligStrNum = -1;	// the strand number of the ligand
	double overlapThresh = -10000.0f;	// hard overlap threshold used for checking sterics (this should be used when the atom positions will not be allowed to change after the steric check)
	double softOverlapThresh = -10000.0f; //soft overlap threshold used for checking sterics (this should be used when the atom positions may be allowed to change after the steric check)
	//	int curAANum[] = null;	// for each residue in the system strand, the
	// index of the current amino acid type; if the residue is not
	// rotamerizable (ie. it's not flexible) then the curAANum entry
	// should be -1
	//int curLigAANum = -1;	// the index of the current ligand type
	boolean computeEVEnergy = false;
	// do we compute EV energies during a conformation search
	boolean doMinimization = false;
	// do we some EV minimization steps during a conformation search
	boolean hElect = true;
	// should hydrogens be used in electrostatic energy calculations
	boolean hVDW = true;
	// should hydrogens be used in vdw energy calculations
	boolean hSteric = false; //should hydrogens be used in steric checks
	double vdwMultiplier = 1.0f;
	// vdw multiplier used in energy evaluation
	boolean addHydrogens = true;
	// during a mutation, should hydrogens be included when
	//  changing residue type
	boolean connectResidues = true;
	// during a mutation, should a new residue be bonded to
	//  the prior and subsequent resiudes if the numbering
	//  is sequential
	int curConfNum = 0;
	int numMinSteps = 35; //140
	// number of minimization steps to perform by simpmin

	boolean doPerturbations;//Indicates DEEPer
	boolean minimizePerturbations;//When doing DEEPer, indicates perturbation-minimizing mode (i.e. DEE/A* step)
	String pertFile;//File with the perturbation information

	double switchedTemplateE[] = null;//The changes in template energies caused by full structure switch perturbations

	AStar MSAStarSearch;

	private Emat arpMatrix = null;
	// all rotamer pairs lower min energy bound matrix, created with
	//  simplePairwiseMutationAllRotamerSearch and loaded with
	//  loadPairwiseEnergyMatrices()
	//private PairwiseEnergyMatrix arpMatrixMax = null;
	// all rotamer pairs lower max energy bound matrix, created with
	//  simplePairwiseMutationAllRotamerSearchMax and loaded with
	//  loadPairwiseEnergyMatricesMax()

	int ASAANums[] = null;
	// integer array containing the index for each AS residues
	//  that can be used with rotamerIndexOffset and for the arpMatrix
	int curStrRotNum[] = null;
	// integer array containing the currently assumed rotamer for
	//  each amino acid in the active site
	// is allocated during a rotamer search.
	// note that it is _not_ the same size as curAANum
	//int curLigRotNum = 0;
	// the current rotamer number of the ligand
	double bestEMin = 9999999.0f; //this should ONLY be accessed/modified using the synchronized methods below
	double bestEUnMin = 9999999.0f;
	// the best minimized and unminimized energy found thus far
	double lowestOverallBound = Double.POSITIVE_INFINITY;
	BigInteger numConfsTotal = new BigInteger("0");
	// the number of total conformations for the current configuration
	// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsLeft = new BigInteger("0");
	// the number of remaining conformations for the current configuration
	//  updated as the search progresses
	BigInteger numConfsBelowLevel[] = null;
	// the number of conformations below the specified level, level 0
	//  refers to the ligand level, if there's no ligand then level 0
	//  and level 1 have the same value
	// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsAboveLevel[] = null; //at level i, num confs from level i+1 to the last level
	BigInteger numConfsPrunedByE = new BigInteger("0");
	// the number of conformations not minimized because their 'best'
	//  energy (as computed from the arpMatrix) was too unfavorable
	//  based on the accuracy threshold below
	BigInteger numConfsPrunedByS = new BigInteger("0");
	// number of conformations pruned due to a steric clash
	BigInteger numConfsPrunedByMinDEE = new BigInteger("0");	//the number of confs pruned by MinDEE
	BigInteger numConfsEvaluated = new BigInteger("0");
	// number of conformations that got all the way down to the energy
	//  evaluation
	// Note that numConfsPrunedByE + numConfsPrunedByS + numConfsEvaluated
	//  should equal the total number of conformations
	double KSepsilon = 0.03;
	// the accuracy for computing energies for K*
	// a value of 0.03 means the energies computed
	//  will allow for a calculation of K*_approx
	//  that's within 3% of the true K*
	BigDecimal partial_q = new BigDecimal(0.0);
	// the partially computed partition function (updated as we go)
	BigDecimal partial_p = new BigDecimal(0.0);
	// the bound on the partition function of the pruned conformations
	BigDecimal initial_q = new BigDecimal(0.0);
	// used in mutation search as an initial partial_q if we're
	//  bootstrapping the search
	StrandRotamers[] strandRot = null;
	int numberOfStrands = -1;
	MutableResParams strandMut;
	//int mutRes2Strand[] = null;
	//int mutRes2StrandMutIndex[] = null;
	int numberMutable = 0;

//	boolean isTemplateOn = false;	
//	boolean minimizePairwise = true;//Minimize pairwise energies instead of minimizing energy of pair and subtracting out intra terms

	//Replacing isTemplateOn and minimizePairwise with an enum that controls both parameters
	public enum MINIMIZATIONSCHEME {
		PAIRWISE, WITHRES, WITHTEMPL
	}

	//these thresholds will be used to detect steric clashes for intra+template
	//and pairwise energies respectively, when we need to reject clashing voxels for the CETM and 
	//pruning information isn't available (i.e. for K* before splitting into single-mutant calculations)
	double pairSt = 100;
	double templateSt = 30; 

	HBondSettings hbonds;
	

	EPICSettings es = new EPICSettings();//turn off EPIC by default 
	boolean compCETM = false;//Indicates we are currently computing EPIC fits
	CETMatrix cetm = null;

	// Note that the mutation search functions in this class are relatively
	//  messy as a result of changing the algorithms multiple times. They
	//  could be rewritten to be much tighter and more elegant.

	//Variables for Improved Bounds algorithms (tuples and partitioned rotamers)
	HashMap<ArrayList<Index3>,EnergyTuple> energyTuples; //HashMap of all the tuples (3 or more rotamers) that have been calculated so far
	ArrayList<GlobalRCConf> generatedConfs; //Keep track of the conformations generated and their energies per residue
	ArrayList<ArrayList<Index3wVal>> worstRots; //Keep track of the worst rotamers from the generated conformations;
	int contractPos = -1;
	int contractPos2 = -1;
	boolean partitionedRotamers;
	boolean tuples;
	boolean superRotamers;
	boolean improvedBounds;
	
	PrintStream outPS = System.out;
	
	// the constructor if you also have a ligand
	RotamerSearch(Molecule theMolec, int numMut,int strandsPresent, boolean hE, 
			boolean hV, boolean hS, boolean addH, boolean conRes, double eps, 
			double stericThresh, double softStericThresh,	boolean ddDielect, 
			double dielectC, boolean doDihedral, boolean doSolv,double solvScFactor, 
			double vdwMult,	boolean doPerts, String pFile, boolean minPerts, boolean useTrips, boolean flagsAStar,
			EPICSettings epics,HBondSettings hbonds, MutableResParams strandMut) {

		hElect = hE;
		hVDW = hV;
		hSteric = hS;
		addHydrogens = addH;
		connectResidues = conRes;
		KSepsilon = eps;
		overlapThresh = stericThresh;
		softOverlapThresh = softStericThresh;
		vdwMultiplier = vdwMult;
		distDepDielect = ddDielect;
		dielectConst = dielectC;
		doDihedE = doDihedral;
		doSolvationE = doSolv;
		solvScale = solvScFactor;
		numberMutable = numMut;
		setBestE(stericE);

		doPerturbations=doPerts;
		pertFile = pFile;
		minimizePerturbations = minPerts;
		useTriples = useTrips;
		useFlagsAStar = flagsAStar;

		m=theMolec;
		a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier,hbonds);
		a96ff.calculateTypesWithTemplates(); //KER: Needed so that the N and C terminus flags are set correctly

		numberOfStrands = strandsPresent;
		strandRot = new StrandRotamers[numberOfStrands];

		//		curAANum = new int[m.numberOfResidues]; //

		//		if(doPerturbations){
		simpMin = new PMinimizer(minimizePerturbations);
		for(int i=0; i<numberOfStrands;i++){
			strandRot[i] = new StrandRCs(m.rotLibForStrand(i),m.strand[i]);
		}
		//		}
		//		else{                    
		//			simpMin = new SimpleMinimizer();
		//			bbMin = new BBMinimizer();
		//			brMin = new BackrubMinimizer();
		//			for(int i=0; i<numberOfStrands;i++){
		//				strandRot[i] = new StrandRotamers(m.rotLibForStrand(i),m.strand[i]);
		//			}
		//		}

		es = epics;
		this.hbonds = hbonds; 
		this.strandMut = strandMut;
		//sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
		/*if (ligStrNum>=0) { //there is a ligand
			m.strand[ligStrNum].residue[0].flexible = true;
			ligROT = new StrandRotamers(grl,m.strand[ligStrNum]);		
			ligROT.setAllowable(0,m.strand[ligStrNum].residue[0].name);
		}*/
	}


	//Set up perturbations and residue conformations (RCs) for DEEPer
	//This should be called after the allowable AAs are set
	public void setupRCs(boolean doPerturbations){

		if(doPerturbations)
			PertFileHandler.readPertFile(pertFile, m, strandRot,true);

		for (int str=0;str<numberOfStrands;str++){

			StrandRCs sr = (StrandRCs)strandRot[str];

			//Could be a problem here!  Ligand needs storage regardless because of rotation/translation!!
			//Should handle this by reversing translation/rotation
			//			if(addWTRot)
			//				sr.storeWTRotamers(m);

			sr.addUnperturbedRCs(m);
			sr.countRCs();
		}

	}


	// This function adds the AA type named name to the list
	//  of allowable types for residue number resNum in
	//  the system strand (resNum is strand based numbering)
	public void setAllowable(int resNum, String name, int strNum) {
		strandRot[strNum].setAllowable(resNum,name);
		m.strand[strNum].residue[resNum].flexible = true;
	}


	//	public boolean[][][][][][] getSplitFlags(){
	//		return splitFlags;
	//	}

	public boolean[][][][][][][][][] getTripleFlags(){
		return tripleFlags;
	}

	public synchronized double getBestE(){
		return bestEMin;
	}

	public synchronized void setBestE(double newBestE){
		bestEMin = newBestE;
	}

	//updates bestEMin only if (ue<bestEMin)
	public synchronized void updateBestE(double ue){
		bestEMin = Math.min(bestEMin, ue);
	}

	// This function clears the list of allowable AA types
	//  for residue number resNum in the system strand
	public void clearAllowable(int resNum, int strandNum) {
		strandRot[strandNum].clearAllowable(resNum);
	}

	// Refreshes the system strand
	public void refreshStrand(int str){
		strandRot[str] = strandRot[str].reInit(m.strand[str]);
	}

	public double getEc_const(){
		return Ec_const;
	}

	public double getBoundForPartition(){
		return boundForPartition;
	}

	private double[] calcTotalSnapshotEnergyParts(){
		return (a96ff.calculateTotalEnergy(m.actualCoordinates,-1)); //compute the energy
	}

	// This function computes one energy 
	private double calcTotalSnapshotEnergy(){

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy

		return (double)energyTerms[0]; //the total energy is in energyTerms[0]
	}
	
	/**
	 * This function computer the energy of the molecule. An additional array is added
	 * that contains all of the mutable residues. The energy function stores the single
	 * and pairwise energies associated with each mutable residue 
	 * @param mutRes
	 * @return
	 */
	private double calcTotalSnapshotEnergy(ArrayList<Integer> mutRes){

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1,mutRes); //compute the energy
		
		return energyTerms[0]; //the total energy is in energyTerms[0]
	}

	//// BEGIN CHECK_STERICS CODE SECTION
	//This version checks all residues against residue resNum (strand-relative numbering) of strand strandNum;
	//This function is only called *before* minimization for the pairwise matrix energy computation;
	//The residue numbers (strand-relative numbering) in excludeRes[] (from the system strand) are not included in the steric check
	private boolean RS_CheckAllSterics(int strandNum, int resNum, int excludeList[]) {

		ProbeStericCheck psc = new ProbeStericCheck();

		Residue res = m.strand[strandNum].residue[resNum];

		for(int i=0;i<res.numberOfAtoms;i++) {
			Atom a1 = res.atom[i];
			if ( hSteric || (!a1.elementType.equalsIgnoreCase("H")) ) {
				for(int q=0;q<m.numberOfResidues;q++) {
					if(!(q==res.moleculeResidueNumber)) {
						boolean resInList = false;
						//if (q==sysStrNum) //check for excluded residues in the system strand
						resInList = isInList(q,excludeList);
						for(int t=0;t<m.residue[q].numberOfAtoms;t++) {
							Atom a2 = m.residue[q].atom[t];
							if ( !resInList || a2.getIsBBatom() ) { //only check the backbone atoms for the sysStrand residues that are in excludeRes[]
								if ( hSteric || (!a2.elementType.equalsIgnoreCase("H"))) {
									if (!psc.isAllowedSteric(m, a1, a2, softOverlapThresh))
										return false;
								}
							}
						}
					}
				}
			}		
		}

		// If you got here then everything passed
		return true;
	}
	//Checks if a is in list[]
	private boolean isInList(int a, int list[]){
		for (int i=0; i<list.length; i++){
			if (list[i]==a)
				return true;
		}
		return false;
	}

	//// END CHECK_STERICS CODE SECTION

	//// BEGIN HELPER FUNCTION SECTION

	// Loads the min (minMatrix==true) or max (minMatrix==false) pairwise energy matrix
	public void loadPairwiseEnergyMatrices(String allRotamerPairsEnergyName) {

		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(allRotamerPairsEnergyName));
			arpMatrix = (Emat)in.readObject();
			in.close();
		}
		catch (Exception e){
			e.printStackTrace();
			return;
		}		

		//arpMatrix.reconnect(m);

	}


	//Load a scaled continuous energy term matrix
	public void loadCETMatrix(String fileName) {
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
			cetm = (CETMatrix)in.readObject();
			in.close();
		}
		catch (Exception e){
			e.printStackTrace();
			System.out.println("Could not load CET Matrix: "+fileName);
		}
	}

	//	public PairwiseEnergyMatrix getMinMatrix(){
	//		return arpMatrix;
	//	}
	//
	//	public PairwiseEnergyMatrix getMaxMatrix(){
	//		return arpMatrixMax;
	//	}

	public Emat getMinMatrix(){
		return arpMatrix;
	}

	public void setMinMatrix(Emat emat){
		arpMatrix = emat;
	}




	//Adds the reference energies to the intra-energies in arpMatrix;
	//If doMinimize is true, then arpMatrixMax is also updated appropriately
	//Adds the reference energies to the intra-energies in arpMatrix;
	//If doMinimize is true, then arpMatrixMax is also updated appropriately
	//	public void addEref(double eRef[][], boolean doMinimize, int strandMut[][]){
	//
	//		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
	//
	//		int ctr=0;
	//		for (int str=0;str<numberOfStrands;str++){
	//			for(int i=0;i<strandMut[str].length;i++){
	//				for(int j=0; j<strandRot[str].getNumAllowable(strandMut[str][i]);j++){
	//					int aaInd = strandRot[str].getIndexOfNthAllowable(strandMut[str][i],j);
	//					int numRot = getNumRot(str, strandMut[str][i], aaInd);
	//
	//					for (int k=0; k<numRot; k++){
	//						arpMatrix.addToIntraE(ctr,aaInd,k, -eRef[ctr][aaInd]);
	//						if (doMinimize)
	//							arpMatrixMax.addToIntraE(ctr,aaInd,k, -eRef[ctr][aaInd]);
	//						ind++;
	//					}
	//				}			
	//				ctr++;
	//			}		
	//		}
	//	}

	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	//	public double computeBestRotEnergyBound(/*int numTotalRotamers, int rotamerIndexOffset[]*/) {
	//
	//		double bestE = arpMatrix.getShellShellE(); // Add shell-shell energy
	//
	//		for(int i=0;i<ASAANums.length;i++) {
	//			bestE += arpMatrix.getIntraAndShellE(i, ASAANums[i], curStrRotNum[i]); //Add the intra-rotamer and shell energies for each rotamer
	//			for(int j=i+1;j<ASAANums.length;j++){ // Add the pairwise energies
	//				bestE += arpMatrix.getPairwiseE( i, ASAANums[i], curStrRotNum[i], j, ASAANums[j], curStrRotNum[j] );
	//			}
	//		}
	//
	//		return bestE;
	//	}

	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	public double computeBestRotEnergyBound(Index3[] rots) {

		double bestE = arpMatrix.getTemplMinE(); // Add shell-shell energy

		for(int i=0;i<rots.length;i++) {

			double tmpSingleE = arpMatrix.getSingleMinE(rots[i]); //Add the intra-rotamer energy
			bestE += tmpSingleE;
			for(int j=i+1;j<rots.length;j++){ // Add the pairwise energies
				if(arpMatrix.areNeighbors(i, j)){
					double tmpPairE = arpMatrix.getPairMinE(rots[i], rots[j]);
					bestE += tmpPairE;
				}
			}
		}

		return bestE;
	}


	//// END HELPER FUNCTION SECTION


	//// BEGIN MASTER MUTATION SEARCH SECTION

	// This function is the model for a mutation search
	// This is the function used by the master node to generate a list of
	//  mutations that it wishes to consider.
	// Utilizes a number of helper functions
	public int simpleMasterMutationSearch(MutableResParams strandMut, int numMutable,
			int theCurConfNum,
			Set<OneMutation> mutArray, double minVol, double maxVol) {

		curConfNum = theCurConfNum;

		AARotamerType curAAtypes[] = new AARotamerType[numMutable];
		
		masterMutationSearchHelper(0, numMutable, strandMut, mutArray, minVol, 
				maxVol, curAAtypes);

		return curConfNum;
	}
	// This function is similar to mutationSearchHelper
	//  the only difference is that we only compute volumes and an amino acid
	//  level energy approximation. I could have modified that function, but
	//  decided not to so as to keep that function fast (ie. this way the
	//  execution of a bunch of conditionals is saved in the normal search)
	public void masterMutationSearchHelper(int depth, int maxDepth,
			MutableResParams strandMut, Set<OneMutation> mutSet, double minVol, double maxVol, AARotamerType curAAtypes[]) {

		if (depth >= maxDepth) {
			// If we've arrived here then we're ready to
			//  compute a volume and approximate energy
			if(debug){
				System.out.print(".");
			}
			double curVolume = 0.0f;
			for (int i=0; i<maxDepth; i++){
				curVolume += curAAtypes[i].volume;
			}
			if ((curVolume > minVol) && (curVolume < maxVol)) {
				// Add mutation to mutation array
				OneMutation tMut = new OneMutation();
				assignAANums(strandMut);
				tMut.score = new BigDecimal("0.0");  // Added when aap removed 6/23/03
				tMut.resTypes = new int[maxDepth];
				for(int q=0;q<maxDepth;q++) {
					tMut.resTypes[q] = curAAtypes[q].index;
				}
				tMut.vol = curVolume;
				/*if (curConfNum >= mutArray.length) {
					// If there's no space left, make space in mutArray
					OneMutation newArray[] = new OneMutation[mutArray.length + 5000];
					System.arraycopy(mutArray, 0, newArray, 0, mutArray.length);
					mutArray = newArray;
				}*/
				getNumConfForMut(tMut);
				mutSet.add(tMut);
				curConfNum++;
			}
			return;
		}

		// Check with allowed AAs
		Residue r = m.residue[strandMut.allMut[depth]];

		for(AARotamerType aa: r.AATypesAllowed()) {
			curAAtypes[depth] = aa;
			masterMutationSearchHelper(depth+1,maxDepth,strandMut,mutSet,minVol,maxVol,curAAtypes);
		}
	}
	// Assigns elements of the ASAANums[] array
	private void assignAANums(MutableResParams strandMut) {

		ASAANums = new int[numberMutable];
		int ctr=0;
		for(int i=0;i<strandMut.allMut.length;i++){
			Residue r = m.residue[strandMut.allMut[i]];
			ASAANums[ctr] = strandRot[r.strandNumber].rl.getAARotamerIndex(r.name);
			ctr++;
		}
	}
	//Compute the number of conformations for the given mutation sequence
	private void getNumConfForMut(OneMutation tMut){
		tMut.numConfUB = BigInteger.ONE;
		tMut.numConfB = BigInteger.ONE;
		for (int i=0; i<tMut.resTypes.length; i++){
			int str = m.residue[strandMut.allMut[i]].strandNumber;
			int numRot = strandRot[str].rl.getNumRotForAAtype(tMut.resTypes[i]);
			if (numRot==0)
				numRot = 1;
			tMut.numConfUB = tMut.numConfUB.multiply(BigInteger.valueOf(numRot));
			tMut.numConfB = tMut.numConfUB;
		}
		/*int numLigRot = grl.getNumRotamers(ligROT.getCurRotType(0));
		if (numLigRot==0)
			numLigRot = 1;
		tMut.numConfB = tMut.numConfB.multiply(BigInteger.valueOf(numLigRot));*/
	}
	//// END MASTER MUTATION SEARCH SECTION
	///////////////////////////////////////////////////////////////////////////////////	

	///////////////////////////////////////////////////////////////////////////////////
	//// BEGIN PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION

	// This function helps to compute the min/max pairwise energy matrix, two
	//  residues are allowed to mutate (as specified in residueMutatable)
	//  a steric check is done for each rotamer pair if the steric check
	//  passes then the min/max energy is computed for the pair and is
	//  saved. ALL rotamer pairs that pass the steric threshold are
	//  saved. If a pair doesn't pass the threshold then an energy of
	//  10^38 is assigned.
	// To compute the pairwise interactions of each rotameric position with
	//  the ligand only one residue should be "allowed" in residueMutatable[]
	//  and the ligPresent should be set to true.
	// A shellRun computes the energy of the "allowed" residue with all other
	//  residues that are NOT in strandMut
	// Utilizes a number of helper functions
	public void simplePairwiseMutationAllRotamerSearch(MutableResParams strandMut, int mutableSpots, 
			boolean searchDoMinimize, boolean shellRun, boolean intraRun, 
			int residueMutatable[], boolean minimizeBB, boolean doBackrubs, boolean templateOnly, String backrubFile, MINIMIZATIONSCHEME minScheme,
			EmatCalcParams runParams, boolean doCompCETM) {
		//If doCompCETM is true we are computing a continuous energy term matrix; otherwise we're computing the min/max pairwise energy matrices

		compCETM = doCompCETM;

		doMinimization = searchDoMinimize;
		computeEVEnergy = true;

		if(doPerturbations)
			((PMinimizer)simpMin).minimizePerturbations=minimizePerturbations;//This field needs to be set properly when this function is called
		if(useCCD)//Need a list of degrees of freedom to use the CCD minimizer
			m.DOFs = DegreeOfFreedom.makeDOFArray(strandRot, strandMut, m);

		// Prepare Amber
		if(computeEVEnergy){
			// Amber should already be loaded
			// First turn off energy evaluation for all residues
			//   since we're only computing pairwise energies
			//   we only want specific residues on (will be turned on later)
			// If we're doing a shell run then turn the shell on
			switch(minScheme){
			case PAIRWISE: //If we don't want the template on during minimization turn it off
			case WITHRES: //if(!templateAlwaysOn){
				if (shellRun) {
					for(int i=0;i<m.numberOfResidues;i++){
						m.residue[i].setEnergyEval(true, true);
						m.residue[i].flexible = false;
					}
					//Need to turn off mutable res for the template calculation
					for(int i=0;i<strandMut.allMut.length;i++){ 
						m.residue[strandMut.allMut[i]].setEnergyEval(false, false);
						m.residue[strandMut.allMut[i]].flexible = false;
					}
				}
				else {
					for(int i=0;i<m.numberOfResidues;i++){
						m.residue[i].setEnergyEval(false, false);
						m.residue[i].flexible = false;
					}
				}
			break;
			case WITHTEMPL: //if we want the template on during minimization make sure it's on
				System.out.println("Template on during minimization of bounds (this might take longer)");
				// 2010: If templateAlwaysOn then always have the shell on for pairwise calculations.
				for(int i=0;i<m.numberOfResidues;i++){
					m.residue[i].setEnergyEval(true,true);
					m.residue[i].flexible = false;
				}
				for(int i=0;i<strandMut.allMut.length;i++){
					if(doPerturbations || minimizeBB || doBackrubs) //If the backbone moves it isn't part of the template
						m.residue[strandMut.allMut[i]].setEnergyEval(false,false);
					else //if the backbone can't move it can be part of the template
						//Theoretically the bb could be turned on, but in the current implementation
						//this would cause the bb to be counted as part of the actual E
						m.residue[strandMut.allMut[i]].setEnergyEval(false,false);
					m.residue[strandMut.allMut[i]].flexible = false;

				}
				break;
			default:
				System.out.println("Don't recognize minimization scheme: "+minScheme);
				System.exit(0);
				break;
			}


			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Error: simpMin not allocated and you are attempting minimization");
						System.exit(1);
					}
					bbMin = null;
					brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){ // phi/psi minimization
						if (bbMin == null) {
							System.out.println("Error: bbMin not allocated and you are attempting backbone minimization");
							System.exit(1);
						}
						simpMin = null;
						brMin = null;
					}
					else { //minimization with backrubs
						if (brMin == null) {
							System.out.println("Error: brMin not allocated and you are attempting backrub minimization");
							System.exit(1);
						}
						simpMin = null;
						bbMin = null;
					}
				}
			}
		}


		if (shellRun) { //compute the template energies
			if ((!minimizeBB)&&m.perts.length>0) {//side-chain minimization, so the template is fixed

				//The full structure switch perturbation can change the template
				//If we have one, it'll be in m.perts[0], and we need to calculate a template energy for each of its states
				if( m.perts[0].type.equalsIgnoreCase("FULL STRUCTURE SWITCH") ){

					int numTemplates = ((FullStructureSwitch)m.perts[0]).numStructs;
					switchedTemplateE = new double[numTemplates];

					for(int tempNum=0; tempNum<numTemplates; tempNum++){

						if(tempNum>0)
							m.perts[0].applyPerturbation(tempNum);

						a96ff.calculateTypesWithTemplates();
						a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);
						double minE = calcTotalSnapshotEnergy();

						if(tempNum==0){
							arpMatrix.templ_E = minE;
							//							retEMatrixMax.setShellShellE(minE);
						}
						else
							switchedTemplateE[tempNum] = minE - arpMatrix.templ_E;
					}
				}
			}
		}

		if (templateOnly) { //the template energies are computed only once	

			if (!minimizeBB && !doMinimization) {//No minimization, so the template is fixed

				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);

				double minE = calcTotalSnapshotEnergy();
				arpMatrix.setShellShellE( minE );
				//				retEMatrixMax.setShellShellE( minE );
			}
			if (!minimizeBB) {//side-chain minimization, so the template is fixed	
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);

				simpMin.initialize(m,numberOfStrands,a96ff,strandRot,doDihedE);

				double beginE = calcTotalSnapshotEnergy();
				simpMin.minimize(35);
				double curE = calcTotalSnapshotEnergy();

				double minE = Math.min(beginE, curE);

				arpMatrix.setShellShellE(minE);
				//				retEMatrixMax.setShellShellE(minE);
				m.updateCoordinates();
				m.revertPertParamsToCurState();
			}
			else if (minimizeBB){ //the template energies for backbone minimization are computed only once

				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);
				if (!doBackrubs){
					bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
					pairwiseRotamerEnergyBackboneHelper(-1,-1,-1,-1,shellRun, arpMatrix, true,-1,-1,-1,-1);					
					bbMin = new BBMinimizer(); //reset the backbone minimizer
				}
				else {
					brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true); //no ligand for template energy computation					
					pairwiseRotamerEnergyBackrubsHelper(-1,-1,-1,-1, shellRun, arpMatrix, true, -1,-1,-1,-1);					
					brMin = new BackrubMinimizer(); //reset the backbone minimizer
				}
			}
		}

		else if (intraRun) { //INTRA run
			computeIntraRotEnergies(strandMut,arpMatrix,residueMutatable,minimizeBB,doBackrubs, backrubFile);			
			return;
		}

		else {	//pairwise rotamer or rot-shell run	
			// Initialize curAANum array, we have to do this because we only recurse through the 9 core residues of the active site
			//  and not all 40 residues, so we use a residueMap and we have to do some preinitialization.
			//			for(int i=0;i<m.numberOfResidues;i++)
			//				curAANum[i] = -1;

			// Note: In this search we only search over the key active site  residues rather than all residues in the molecule, 
			//		thus maxDepth should be numInAS

			if(compCETM)
				cetm = new CETMatrix(arpMatrix,m);

			pairwiseEnergyComp (residueMutatable, shellRun, minimizeBB, doBackrubs, backrubFile,runParams, minScheme);

			if( ( switchedTemplateE != null ) && ( residueMutatable[0] != 0 ) )//This transfers energies from template changes in full structure switch perturbations to the first residue shell energies
				correctSwitchedTemplates(strandMut, arpMatrix);

			//this is now obsolete because we have CETMatrix.tryMolecCompression
			//if(compCETM)
			//    cetm.cleanupSVE();

			return;
		}
	}

	private void computeIntraRotEnergies(MutableResParams strandMut, Emat retEMatrixMin,
			int residueMutatable[],
			boolean minimizeBB, boolean doBackrubs, String backrubFile) {

		a96ff.refEnergy = true;
		for(int i=0;i<m.residue.length;i++)
			m.residue[i].setEnergyEval(false,false);

		// Go through each active site residue, each AA type they could be and all their rotamers, 
		// 	saving the computed energies to the appropriate place.
		int prevPos = -1;
		boolean firstTime = true;
		Iterator<EMatrixEntryWIndex> iter = arpMatrix.singlesIterator();
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			EMatrixEntry re = emeWI.eme;

			if(prevPos != emeWI.pos1())
				firstTime = true;

			System.out.print(".");

			boolean neededMutation = re.applyMutation(m, arpMatrix.resByPos,addHydrogens, connectResidues);
			re.flexible(m,arpMatrix.resByPos,true);
			re.applyRC(arpMatrix.resByPos,m);

			re.setEnergyEval(m,arpMatrix.resByPos,true,true);


			// Setup Amber, Setup Minimizer
			if (computeEVEnergy && (neededMutation || firstTime)){
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);
				if (doMinimization){
					if (!minimizeBB) //side-chain minimization
						simpMin.initialize(m,numberOfStrands,a96ff,strandRot,doDihedE);
					else { //backbone minimization
						if (!doBackrubs)
							bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
						else{
							brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
						}
					}
				}
			}




			computeIntraRotEnergiesHelper(emeWI, minimizeBB, doBackrubs);


			re.flexible(m,arpMatrix.resByPos,false);
			re.setEnergyEval(m,arpMatrix.resByPos,false,false);

			firstTime = false;	
			prevPos = emeWI.pos1();
		}


	}
	//Computes the intra-residue min and max energies for a given rotamer
	private void computeIntraRotEnergiesHelper(EMatrixEntryWIndex emeWI, /*int totRotForCur,*/ 
			boolean minimizeBB, boolean doBackrubs){

		double curEnergy = 0.0f;	
		double beginE = 0.0f;
		double minEnergy = (double)Math.pow(10,30);
		double maxEnergy = -(double)Math.pow(10,30);

		//Compute min and max energies
		if ( doMinimization ){								

			//Iniitialize Amber for the current residue
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);

			if ( (!minimizeBB) && (!doBackrubs) ) { //side-chain minimization

				//first, compute the initial energy at the new position
				beginE = calcTotalSnapshotEnergy();

				//minimize, starting at the initial position
				simpMin.minimize(numMinSteps);
				curEnergy = calcTotalSnapshotEnergy();

				//Compare to the min and max energies found so far and update, if necessary
				if (beginE<curEnergy)
					curEnergy = beginE;
				double lE = Math.min(beginE, curEnergy);
				double hE = Math.max(beginE, curEnergy);
				minEnergy = Math.min(minEnergy,lE);									
				maxEnergy = Math.max(maxEnergy,hE);

				m.updateCoordinates();//restore the actualCoordinates array to the initial values
				m.revertPertParamsToCurState();
			}
			else if (!doBackrubs) { //phi/psi backbone minimization, so only rotate the backbone O and HN		

				int at[] = new int[5];
				int numAtoms = 0;
				Residue r1 = null;
				int resID = arpMatrix.resByPos.get(emeWI.pos1()).get(0);//0 index assumes no super-rotamers
				r1 =  m.residue[resID]; 

				numAtoms = r1.numberOfAtoms;

				//get the atoms
				for (int i=0; i<numAtoms; i++){
					if (r1.atom[i].name.equalsIgnoreCase("CA"))
						at[0] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("C"))
						at[1] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("O"))
						at[2] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("N"))
						at[3] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("H"))
						at[4] = r1.atom[i].moleculeAtomNumber;
				}

				//get the CA-C bond
				double dx = m.actualCoordinates[at[1]*3] - m.actualCoordinates[at[0]*3];
				double dy = m.actualCoordinates[at[1]*3+1] - m.actualCoordinates[at[0]*3+1];
				double dz = m.actualCoordinates[at[1]*3+2] - m.actualCoordinates[at[0]*3+2];

				//get the N-CA bond
				double dxH = m.actualCoordinates[at[0]*3] - m.actualCoordinates[at[3]*3];
				double dyH = m.actualCoordinates[at[0]*3+1] - m.actualCoordinates[at[3]*3+1];
				double dzH = m.actualCoordinates[at[0]*3+2] - m.actualCoordinates[at[3]*3+2];

				//get the center of rotation for O (the actualCoordinates[] of C)
				double center[] = new double[3];
				center[0] = m.actualCoordinates[at[1]*3];
				center[1] = m.actualCoordinates[at[1]*3+1];
				center[2] = m.actualCoordinates[at[1]*3+2];

				//get the center of rotation for H (the actualCoordinates[] of N)
				double centerH[] = new double[3];
				centerH[0] = m.actualCoordinates[at[3]*3];
				centerH[1] = m.actualCoordinates[at[3]*3+1];
				centerH[2] = m.actualCoordinates[at[3]*3+2];

				double rotForInitPos = bbMin.getMaxDihedRot(); //get the max phi/psi rotation

				//Do the sampling and minimization
				for (int curSample=0; curSample<samplesBB; curSample++){

					//randomly generate the rotation angle
					double rotChange[] = new double[2];
					Random randNum = new Random();

					if (curSample!=0) {
						for (int i=0; i<2; i++)
							rotChange[i] = (randNum.nextDouble()-0.5f)*rotForInitPos*2.0f;
					}
					else {
						for (int i=0; i<2; i++)
							rotChange[i] = 0.0f;
					}

					//Compute the energy corresponding to the new positions
					if (curSample!=0){ //do not apply a change for the initial position
						m.rotateAtom(at[2], dx, dy, dz, center[0], center[1], center[2], rotChange[0], false);
						m.rotateAtom(at[4], dxH, dyH, dzH, centerH[0], centerH[1], centerH[2], rotChange[1], false);
					}

					//compute the initial energy at the new position
					curEnergy = calcTotalSnapshotEnergy();

					//Compare to the min and max energies found so far and update, if necessary;
					//For intra-energies with BB flexibility, the initial point is taken as the max energy,
					//		since the O and H positions are only sampled without minimization
					minEnergy = Math.min(minEnergy,curEnergy);
					if (curSample==0)
						maxEnergy = Math.max(maxEnergy,curEnergy);

					m.updateCoordinates();//restore the actualCoordinates array to the initial values
				}
			}
			else { //backrub minimization
				if (emeWI.pos1()>=0){ //not the ligand
					beginE = calcTotalSnapshotEnergy();
					double e[] = brMin.getMinMaxIntraEnergyBR(emeWI.pos1());
					double lE = Math.min(beginE, e[0]);
					double hE = Math.min(beginE, e[1]);
					minEnergy = Math.min(minEnergy,lE);									
					maxEnergy = Math.max(maxEnergy,hE);
				}
				else { //the ligand (backrubs are not applied and have no effect on the intra-energy)
					minEnergy = calcTotalSnapshotEnergy();
					maxEnergy = minEnergy;
				}
				m.updateCoordinates();//restore the actualCoordinates array to the initial values
			}
		}
		else if (computeEVEnergy){
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			curEnergy = calcTotalSnapshotEnergy();
			m.updateCoordinates();
			minEnergy = curEnergy;
			maxEnergy = curEnergy;
		}

		// Store result
		int pos = -1;
		int aa = -1;

		arpMatrix.setSingleE(emeWI.index, minEnergy);
		//		retEMatrixMax.setIntraE(pos, aa, curRot, maxEnergy);
	}
	// Sets up the computation for residue-to-template (SHL-AS, LIG-SHL) and rotamer-rotamer (AS-AS, LIG-AS) energies;
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable);
	//		Sets all non-mutatable residues in residueMap to either Gly or Ala (and leaves all Pro)
	public void pairwiseEnergyComp (int residueMutatable[],
			boolean shellRun, boolean minimizeBB, boolean doBackrubs, String backrubFile, EmatCalcParams runParams,
			MINIMIZATIONSCHEME minScheme) {

		//KER: this needs to be called to initialize the cterm and nterm flags before
		//KER: anything is mutated.
		a96ff.calculateTypesWithTemplates();

		for (int depth=0; depth<residueMutatable.length; depth++){				
			if (residueMutatable[depth]==0){ //not a mutatable residue, but in residueMap[]
				// Make this residue a "gly", if it is not already GLY or PRO;
				//		If minimizeBB, then change to ALA, if not already GLY or PRO
				for(Integer resID:arpMatrix.resByPos.get(depth)){
					m.residue[resID].flexible = false;
					switch(minScheme){
						case PAIRWISE:
						case WITHRES:
							m.residue[resID].setEnergyEval(false, false);
							break;
						case WITHTEMPL:
							if(doPerturbations || minimizeBB || doBackrubs) //If the backbone moves it isn't part of the template
								m.residue[resID].setEnergyEval(false,false);
							else //if the backbone can't move it can be part of the template
								//Theoretically the bb could be turned on, but in the current implementation
								//this would cause the bb to be counted as part of the actual E
								m.residue[resID].setEnergyEval(false,false);
							break;
					}
				}
			}
		}

		pairwiseEnergyCompAllMutatedResHelper(residueMutatable, shellRun, 
				minimizeBB, doBackrubs, runParams, 0, backrubFile, minScheme);
	}
	// Helper to pairwiseEnergyComp();
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable[]);
	//		the non-mutatable residues in residueMutatable[] have already been set to either Gly or Ala (all Pro remain)
	private void pairwiseEnergyCompAllMutatedResHelper(int residueMutatable[],
			boolean shellRun, boolean minimizeBB, boolean doBackrubs, 
			EmatCalcParams runParams, int curMut, String backrubFile, MINIMIZATIONSCHEME minScheme){

		Iterator<EMatrixEntryWIndex> rotamerEntries; 
		if(shellRun){
			rotamerEntries = arpMatrix.singlesIterator(runParams.pos1);
			a96ff.onlySingle = true;
			a96ff.onlyPair = false;
		}
		else if(runParams.AAs1 == null || runParams.AAs2 == null){
			rotamerEntries = arpMatrix.pairsIterator(runParams.pos1,runParams.pos2);
			a96ff.onlySingle = false;
			a96ff.onlyPair = true;
			a96ff.pair1 = arpMatrix.resByPos.get(runParams.pos1);
			a96ff.pair2 = arpMatrix.resByPos.get(runParams.pos2);
		}
		else{
			rotamerEntries = arpMatrix.pairsIterator(runParams.pos1,runParams.pos2,runParams.AAs1,runParams.AAs2);
			a96ff.onlySingle = false;
			a96ff.onlyPair = true;
			a96ff.pair1 = arpMatrix.resByPos.get(runParams.pos1);
			a96ff.pair2 = arpMatrix.resByPos.get(runParams.pos2);
		}

		if(minScheme.equals(MINIMIZATIONSCHEME.PAIRWISE)){ //The minimizer will only have the pairwise terms on (same as the energy eval)
			a96ffmin = a96ff;
		}else{
			a96ffmin = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier,hbonds);
			a96ffmin.onlySingle = true;
			a96ffmin.onlyPair = false;
		}

		boolean firstTime = true;
		int ctr = 0;
		while(rotamerEntries.hasNext()){

			ctr++;
			EMatrixEntryWIndex reWi= rotamerEntries.next();
			
			EMatrixEntry re = reWi.eme;

			if(ctr % 50 == 0){
				System.out.print(".");
				//outPS.flush();
			}

			//skip if pruned and CETM calculation
			if(compCETM && re.isPruned())
				continue;

			//KER: only calculate the energies if it is in provided collection
			if(runParams.rotamers == null || runParams.rotamers.contains(new Index3(reWi.index[0],reWi.index[1],reWi.index[2]))
					|| (!shellRun && runParams.rotamers.contains(new Index3(reWi.index[3],reWi.index[4],reWi.index[5])))  ){

				//ApplyMutation
				boolean neededMut = re.applyMutation(m,arpMatrix.resByPos,addHydrogens,connectResidues);
				re.flexible(m,arpMatrix.resByPos,true);
				re.setEnergyEval(m,arpMatrix.resByPos,true,true);
				re.applyRC(arpMatrix.resByPos,m);


				if (computeEVEnergy){
					if(neededMut || firstTime){
						a96ff.calculateTypesWithTemplates();
						a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);

						if(a96ffmin != a96ff){
							a96ffmin.calculateTypesWithTemplates();
							a96ffmin.initializeCalculation();
							a96ffmin.setNBEval(hElect,hVDW);
						}
						firstTime = false;
					}


					if (doMinimization){
						if (!minimizeBB){ //side-chain minimization

							//Setup energy function
							efunc = new ForceFieldEnergy(m,a96ff); //Just normal pair-wise only terms
							
							ContSCObjFunction ofmin = null;
							//Setup minimizers
							if(useCCD){
								//Setup energy function for minimization
								ForceFieldEnergy efuncmin = new ForceFieldEnergy(m,a96ffmin); //Has terms turned on for minimization
								boolean transRotStrands[] = re.transRotStrands(m,arpMatrix.resByPos, strandMut);
								ofmin = new ContSCObjFunction(m,numberOfStrands,efuncmin,strandRot,(doDihedE&&shellRun),transRotStrands);
								ccdMin = new CCDMinimizer(ofmin,true);
							}
							else
								simpMin.initialize(m,numberOfStrands,a96ffmin,strandRot,doDihedE);
							
							//Turn on doDihedE in the energy function
							if(doDihedE && shellRun){
								if(useCCD)
									efunc = ofmin.efunc;//ef will now include dihedral energies if we need them (i.e. if doDihedE && shellRun)
								else{
									DihedralEnergy de = new DihedralEnergy(m, numberOfStrands, simpMin.strDihedralAtNums, simpMin.numStrDihedrals, efunc.getAmber96ext());
									efunc = efunc.addTerm(de);
								}
							}
							
						}
						else { //backbone minimization
							if (!doBackrubs){
								bbMin.initialize(m, a96ffmin, strandMut, numberOfStrands);
							}
							else {
								brMin.initialize(m, a96ffmin, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
							}
						}
					}
				}


				pairwiseMutationAllRotamerSearch(reWi,minimizeBB, doBackrubs, runParams, 0,shellRun);
			}
		}
	}

	// Called by pairwiseEnergyCompAllMutatedResHelper();
	// Computes the energies among all rotamer combinations for residues res1 and res2 (if res2==1, then computes the res1-to-template energies)
	private void pairwiseMutationAllRotamerSearch(EMatrixEntryWIndex reWi,
			boolean minimizeBB, boolean doBackrubs, EmatCalcParams runParams, 
			int curMut, boolean shellRun) {


		//First, if using DEEPer and calculating a pairwise energy,
		//check if the pair is parametrically incompatible
		//and set the pairwise energy to infinity if it is
		//Do the same for residues with non-pre-Pro compatible BB dihedrals
		//right before a proline
		//and for RCs whose application failed
		boolean validConf = true;
		if(doPerturbations){

			Residue firstRes = m.residue[arpMatrix.resByPos.get(reWi.pos1()).get(0)], secondRes; //0 index assumes no super rotamers

			//For a shell run, if the next residue is a shell proline, make sure the RC's backbone dihedrals are OK for pre-pro
			//Also, if this residue is a proline and the previous residue is in the shell, make sure that shell residue's backbone dihedrals are OK for pre-pro
			if( shellRun ){

				if(!firstRes.validConf)//Pro ring not closed
					validConf = false;

				if( m.checkCBonded(firstRes.moleculeResidueNumber) ){

					Residue nextRes = m.strand[firstRes.strandNumber].residue[firstRes.strandResidueNumber+1];
					if( nextRes.getEnergyEvalSC() && nextRes.name.equalsIgnoreCase("PRO") ){//We are only calculating nextRes' energy if it's in the shell
						if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) ){
							validConf = false;
						}
					}
				}


				if( m.checkNBonded(firstRes.moleculeResidueNumber) && firstRes.name.equalsIgnoreCase("PRO") ){

					Residue prevRes = m.strand[firstRes.strandNumber].residue[firstRes.strandResidueNumber-1];
					if( prevRes.getEnergyEvalSC() ){//We are only calculating prevRes' energy if it's in the shell
						if( ! RamachandranChecker.getInstance().checkPrePro(m, prevRes.moleculeResidueNumber) ){
							validConf = false;
						}
					}
				}

				if( !validConf ){
					arpMatrix.setShellRotE( reWi.rot1index(), Double.POSITIVE_INFINITY );
					//						retEMatrixMax.setShellRotE( res1, res1AANum, res1RotNum, Double.POSITIVE_INFINITY );
					return;
				}
			}
			else if (!shellRun){

				boolean incompatible = !validConf;

				secondRes = m.residue[arpMatrix.resByPos.get(reWi.pos2()).get(0)];

				if( firstRes.strandNumber == secondRes.strandNumber ){
					if( (secondRes.strandResidueNumber == firstRes.strandResidueNumber + 1) && ( secondRes.name.equalsIgnoreCase("PRO") ) ){//firstRes is right before secondRes, which is proline
						if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) )
							incompatible = true;//the pair of RCs is incompatible, because secondRes has AA type proline and firstRes' backbone dihedrals are not OK for pre-pro
					}
					else if( (secondRes.strandResidueNumber == firstRes.strandResidueNumber - 1) && ( firstRes.name.equalsIgnoreCase("PRO") ) ){//secondRes if right before firstRes, which is proline
						if( ! RamachandranChecker.getInstance().checkPrePro(m, secondRes.moleculeResidueNumber) )
							incompatible = true;
					}
				}

				if( ! (firstRes.validConf&&secondRes.validConf) )
					incompatible = true;

				int globalRC1 = arpMatrix.singles.getRot(reWi.rot1index())[0]; //0 index implies no super-rotamers
				int globalRC2 = arpMatrix.singles.getRot(reWi.rot2index())[0]; //0 index implies no super-rotamers

				ResidueConformation rc1 = m.strand[firstRes.strandNumber].rcl.getRC(globalRC1);
				ResidueConformation rc2 = m.strand[secondRes.strandNumber].rcl.getRC(globalRC2);


				if( isParametricallyIncompatible(firstRes,rc1,secondRes,rc2) )
					incompatible = true;


				if(incompatible){
					arpMatrix.setPairwiseE(reWi.index, Double.POSITIVE_INFINITY);
					//						retEMatrixMax.setPairwiseE(res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, Double.POSITIVE_INFINITY);
					return;
				}
			}



		}


		//If we've gotten here then check the sterics to make sure we're sterically allowable;
		//Exclude the other residues from residueMap[], since they have been mutated to Gly/Ala
		int excludeRes[] = new int[strandMut.allMut.length];


		for(int j=0;j<excludeRes.length;j++){
			if(!posContainsMolRes(arpMatrix.resByPos.get(runParams.pos1), strandMut.allMut[j])){
				if(shellRun)
					excludeRes[j] = strandMut.allMut[j];
				else if(!posContainsMolRes(arpMatrix.resByPos.get(runParams.pos2), strandMut.allMut[j]))
					excludeRes[j] = strandMut.allMut[j];
				else{
					excludeRes[j] = -1;
				}
			}
			else
				excludeRes[j] = -1;
		}

		//		m.saveMolecule("beforeStericOrig1.pdb", 0.0f);
		//System.out.println("DELETE ME!");
		//KER: TODO: Fix this to sterically check.
		boolean stericallyGood = true;
		stericallyGood = a96ff.checkSterics(m.actualCoordinates);
		//			if (RS_CheckAllSterics(runParams.pos1,excludeRes)) {
		//				if (!shellRun) {
		//					if (RS_CheckAllSterics(runParams.pos2,excludeRes))
		//						stericallyGood = true;
		//				}
		//				else 
		//					stericallyGood = true;
		//			}

		if ((stericallyGood)) { //good steric found or doing backbone minimization

			// After minimization do a m.updateCoordinates() to resync the actualCoordinates which were changed
			//  in the minimization procedure									

			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					pairwiseRotamerEnergySidechainHelper (reWi,
							shellRun, runParams);
				}
				else { //backbone minimization
					if (!doBackrubs) // phi/psi minimization
						pairwiseRotamerEnergyBackboneHelper (reWi,shellRun, false,runParams);
					else { //backrubs
						pairwiseRotamerEnergyBackrubsHelper (reWi,shellRun, false, runParams);
					}
				}
			}
			else if (computeEVEnergy){

				double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
				double curEnergy = energyTerms[0];//calcTotalSnapshotEnergy();

				arpMatrix.setE(reWi.index,curEnergy);

			}
			else {
				System.out.println("This should not happen. No energy evaluation specified");
				System.exit(1);
			}					

			//restore to the coordinates before the energy computation
			m.updateCoordinates();
		}
		else {
			arpMatrix.setE(reWi.index, stericE);
			if(arpMatrix.doDih()){
				double[][] dihedrals = new double[2][0];
				arpMatrix.setDihedrals(reWi.index,dihedrals);
				arpMatrix.setMaxE(reWi.index,stericE);
			}
				
		}

		return;
	}



	//This method helps compute energy bounds for a given pair of rotamers for backbone phi/psi minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackboneHelper (EMatrixEntryWIndex emeWI,
			boolean shellRun, boolean templateOnly, EmatCalcParams runParams){		

		//Initialize Amber for the current pair
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);

		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);

		float beginE = 0.0f;

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		beginE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		//		beginE -= removeExtraE(templateOnly, shellRun,runParams);

		bbMin.minimizeFull(true); //minimize

		energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		//		curEnergy -= removeExtraE(templateOnly, shellRun,runParams);


		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	

		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		/////////////////////////////////////////////

		if (templateOnly){
			arpMatrix.templ_E = minEnergy;
		}
		else{
			arpMatrix.setE(emeWI.index, minEnergy);
		}

		return;
	}

	//This method helps compute energy bounds for a given pair of rotamers for backbone phi/psi minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackrubsHelper (EMatrixEntryWIndex emeWI,
			boolean shellRun, boolean templateOnly, EmatCalcParams runParams){		

		//Initialize Amber for the current pair
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);

		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);

		float beginE = 0.0f;

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		beginE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		//		beginE -= removeExtraE(templateOnly, shellRun,runParams);

		brMin.minimizeFull(shellRun,templateOnly); //minimize

		energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		//		curEnergy -= removeExtraE(templateOnly, shellRun,runParams);


		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	

		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		/////////////////////////////////////////////

		if (templateOnly){
			arpMatrix.templ_E = minEnergy;
		}
		else{
			arpMatrix.setE(emeWI.index, minEnergy);
		}

		return;
	}


	//This method computes the minimized energy for a given pair of rotamers (or rot-to-template)
	//Called by pairwiseMutationAllRotamerSearch(.)
	private void pairwiseRotamerEnergySidechainHelper (EMatrixEntryWIndex emeWI,
			boolean shellRun, EmatCalcParams runParams){

		double minEnergy = Math.pow(10,30);
		double maxEnergy = -Math.pow(10,30);

		double curEnergy = 0.0f;
		double beginE = 0.0f;

		if(useCCD)//The molecule has ideal values for the dihedrals now; record these
			((ContSCObjFunction)ccdMin.objFcn).updateIdealDihedrals();

		double[][] dihedrals = null;


		beginE = (double)efunc.getEnergy();

		//Minimize
		if(useCCD){
			ccdMin.minimize();
			if(arpMatrix.doDih()){
				Residue[] res1 = new Residue[arpMatrix.resByPos.get(runParams.pos1).size()];
				int ctr=0;
				for(int resID: arpMatrix.resByPos.get(runParams.pos1))
					res1[ctr++] = m.residue[resID];
				Residue[] res2 = new Residue[0];
				if(!shellRun){
					res2 = new Residue[arpMatrix.resByPos.get(runParams.pos2).size()];
					ctr=0;
					for(int resID: arpMatrix.resByPos.get(runParams.pos2))
						res2[ctr++] = m.residue[resID];
				}
					
				dihedrals = ccdMin.getCurSCDihedrals(res1, res2);
			}
		}
		else{
			simpMin.minimize(numMinSteps);
			if(arpMatrix.doDih()){
				dihedrals = parseDihedrals(simpMin, runParams); 
			}
		}

		curEnergy = (double)efunc.getEnergy();


//		if (shellRun && doDihedE){ //add dihedral energies
//			if(useCCD)
//				curEnergy += ((ContSCObjFunction)ccdMin.objFcn).de.getEnergy();
//			else
//				curEnergy += simpMin.computeDihedEnergy();
//		}

		if(compCETM){

			boolean stericsOK = true;
			//				TODO: Figure out if I need to replace the eliminatedRotAtRes check. 
			//				if(eliminatedRotAtRes==null){
			//no pruning information
			//so we'll do an extra steric check here instead
			if(shellRun){
				if( curEnergy >= templateSt )
					stericsOK = false;
			}
			else if( curEnergy >= pairSt  )
				stericsOK = false;
			//				}


			if(stericsOK){
				Residue firstRes = m.residue[arpMatrix.resByPos.get(emeWI.pos1()).get(0)];
				Residue secondRes = null;
				if(!shellRun)
					secondRes = m.residue[arpMatrix.resByPos.get(emeWI.pos2()).get(0)];

				ArrayList<ArrayList<AARotamerType>> aaTypes = arpMatrix.getAATypes(m, emeWI.index); 
				AARotamerType aaType1 = aaTypes.get(0).get(0); //0 index assumes no super-rotamers
				AARotamerType aaType2 = null;
				if(!shellRun)
					aaType2 = aaTypes.get(1).get(0); //0 index assumes no super-rotamers
				compEPICFit( firstRes, secondRes, aaType1, aaType2, emeWI, curEnergy, beginE, shellRun );
			}
		}


		minEnergy = Math.min(minEnergy,curEnergy);
		maxEnergy = Math.max(maxEnergy,beginE);		

		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		m.revertPertParamsToCurState();
		//			System.out.println(simpMin.a96ff.numEVevals+" E: "+beginE+" "+curEnergy+" "+(beginE-curEnergy)+" "+ emeWI.index[0]+" "+ emeWI.index[1]+" "+ emeWI.index[2]+" "+ emeWI.index[3]+" "+ emeWI.index[4]+" "+ emeWI.index[5]);
		//if (shellRun) {
		arpMatrix.setE(emeWI.index, minEnergy);
		if(arpMatrix.doDih()){
			arpMatrix.setDihedrals(emeWI.index,dihedrals);
			arpMatrix.setMaxE(emeWI.index,beginE);
		}

	}


	private double[][] parseDihedrals(SimpleMinimizer simpMin, EmatCalcParams runParams) {

		double[][] dihedrals = new double[2][];

		int[] mutPos = {runParams.pos1,runParams.pos2};

		for(int i=0; i<2;i++){
			int pos = mutPos[i];
			if(pos >=0){
				ArrayList<Double> dihForPos = new ArrayList<Double>();
				for(int res:arpMatrix.resByPos.get(pos)){
					for(int str=0; str<simpMin.strDihedToResNum.length;str++){
						for(int j=0; j<simpMin.strDihedToResNum[str].length;j++){
							if(simpMin.strDihedToResNum[str][j] == res)
								dihForPos.add(simpMin.strCumulativeDihedStep[str][j]);	
						}

					}
				}
				dihedrals[i] = new double[dihForPos.size()];
				int ctr=0;
				for(double d: dihForPos)
					dihedrals[i][ctr++] = d;

			}
		}

		return dihedrals;

	}


	private boolean posContainsMolRes(ArrayList<Integer> residueList,
			int moleculeResidueNumber) {
		for(Integer resID:residueList){
			if(m.residue[resID].moleculeResidueNumber == moleculeResidueNumber)
				return true;
		}
		return false;
	}

	// Called by pairwiseEnergyCompAllMutatedResHelper();
	// Computes the energies among all rotamer combinations for residues res1 and res2 (if res2==1, then computes the res1-to-template energies)
	//	private void pairwiseMutationAllRotamerSearch(int maxDepth, MutableResParams strandMut,
	//			int res1, int res2, int res1AANum, int res2AANum, int res1RotNum, int res2RotNum,
	//			Emat retEMatrixMin, 
	//			boolean minimizeBB, boolean doBackrubs, final int numMut, int mutDepth[], int curMut, boolean validConf) {
	//		//validConf indicates that all RC applications were successful (always true if not DEEPer)
	//		if (curMut>=numMut) {
	//
	//			int res1Strand = -1;
	//			int res1Resnum = -1;
	//			int res2Strand = -1;
	//			int res2Resnum = -1;
	//			int res1MutIndex = -1;
	//			int res2MutIndex = -1;
	//
	//			boolean shellRun = false;
	//			if (res2 == -1)
	//				shellRun = true;
	//
	//			//			res1Strand = strandMut.resStrand[res1];
	//			//			res1Resnum = strandMut[res1Strand][mutRes2StrandMutIndex[res1]];
	//			//			res1MutIndex = mutRes2StrandMutIndex[res1];
	//			//			if(res2!=-1){
	//			//				res2Strand = mutRes2Strand[res2];
	//			//				res2Resnum = strandMut[res2Strand][mutRes2StrandMutIndex[res2]];
	//			//				res2MutIndex = mutRes2StrandMutIndex[res2];
	//			//			}
	//
	//
	//			//First, if using DEEPer and calculating a pairwise energy,
	//			//check if the pair is parametrically incompatible
	//			//and set the pairwise energy to infinity if it is
	//			//Do the same for residues with non-pre-Pro compatible BB dihedrals
	//			//right before a proline
	//			//and for RCs whose application failed
	//			if(doPerturbations){
	//
	//				Residue firstRes = m.residue[strandMut.allMut[res1]], secondRes;
	//
	//				//For a shell run, if the next residue is a shell proline, make sure the RC's backbone dihedrals are OK for pre-pro
	//				//Also, if this residue is a proline and the previous residue is in the shell, make sure that shell residue's backbone dihedrals are OK for pre-pro
	//				if( shellRun ){
	//
	//					if(!firstRes.validConf)//Pro ring not closed
	//						validConf = false;
	//
	//					if( m.checkCBonded(firstRes.moleculeResidueNumber) ){
	//
	//						Residue nextRes = m.strand[res1Strand].residue[res1Resnum+1];
	//						if( nextRes.getEnergyEvalSC() && nextRes.name.equalsIgnoreCase("PRO") ){//We are only calculating nextRes' energy if it's in the shell
	//							if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) ){
	//								validConf = false;
	//							}
	//						}
	//					}
	//
	//
	//					if( m.checkNBonded(firstRes.moleculeResidueNumber) && firstRes.name.equalsIgnoreCase("PRO") ){
	//
	//						Residue prevRes = m.strand[res1Strand].residue[res1Resnum-1];
	//						if( prevRes.getEnergyEvalSC() ){//We are only calculating prevRes' energy if it's in the shell
	//							if( ! RamachandranChecker.getInstance().checkPrePro(m, prevRes.moleculeResidueNumber) ){
	//								validConf = false;
	//							}
	//						}
	//					}
	//
	//					if( !validConf ){
	//						retEMatrixMin.setShellRotE( res1, res1AANum, res1RotNum, Double.POSITIVE_INFINITY );
	//						//						retEMatrixMax.setShellRotE( res1, res1AANum, res1RotNum, Double.POSITIVE_INFINITY );
	//						return;
	//					}
	//				}
	//				else if (!shellRun){
	//
	//					boolean incompatible = !validConf;
	//
	//					secondRes = m.strand[res2Strand].residue[res2Resnum];
	//
	//					if( res1Strand == res2Strand ){
	//						if( (res2Resnum == res1Resnum + 1) && ( secondRes.name.equalsIgnoreCase("PRO") ) ){//firstRes is right before secondRes, which is proline
	//							if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) )
	//								incompatible = true;//the pair of RCs is incompatible, because secondRes has AA type proline and firstRes' backbone dihedrals are not OK for pre-pro
	//						}
	//						else if( (res2Resnum == res1Resnum - 1) && ( firstRes.name.equalsIgnoreCase("PRO") ) ){//secondRes if right before firstRes, which is proline
	//							if( ! RamachandranChecker.getInstance().checkPrePro(m, secondRes.moleculeResidueNumber) )
	//								incompatible = true;
	//						}
	//					}
	//
	//					if( ! (firstRes.validConf&&secondRes.validConf) )
	//						incompatible = true;
	//
	//					int globalRC1 = arpMatrix.singles.getRot(res1, res1AANum, res1RotNum)[0]; //0 index implies no super-rotamers
	//					int globalRC2 = arpMatrix.singles.getRot(res2, res2AANum, res2RotNum)[0]; //0 index implies no super-rotamers
	//
	//					ResidueConformation rc1 = strandRot[res1Strand].rcl.getRC(globalRC1);
	//					ResidueConformation rc2 = strandRot[res1Strand].rcl.getRC(globalRC2);
	//
	//
	//					if( isParametricallyIncompatible(firstRes,rc1,secondRes,rc2) )
	//						incompatible = true;
	//
	//
	//					if(incompatible){
	//						retEMatrixMin.setPairwiseE(res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, Double.POSITIVE_INFINITY);
	//						//						retEMatrixMax.setPairwiseE(res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, Double.POSITIVE_INFINITY);
	//						return;
	//					}
	//				}
	//
	//
	//
	//			}
	//
	//
	//			//If we've gotten here then check the sterics to make sure we're sterically allowable;
	//			//Exclude the other residues from residueMap[], since they have been mutated to Gly/Ala
	//			int excludeRes[] = new int[strandMut.allMut.length];
	//			System.arraycopy(strandMut.allMut, 0, excludeRes, 0, strandMut.allMut.length);
	//
	//			excludeRes[res1] = -1; //do not exclude res1, since it is fixed here
	//			if (res2!=-1)
	//				excludeRes[res2] = -1; //do not exclude res2, since it is fixed here
	//			boolean stericallyGood = false;
	//			if (RS_CheckAllSterics(res1Strand,res1Resnum,excludeRes)) {
	//				if (res2 != -1) {
	//					if (RS_CheckAllSterics(res2Strand,res2Resnum,excludeRes))
	//						stericallyGood = true;
	//				}
	//				else 
	//					stericallyGood = true;
	//			}
	//
	//			if ((stericallyGood)) { //good steric found or doing backbone minimization
	//
	//				// After minimization do a m.updateCoordinates() to resync the actualCoordinates which were changed
	//				//  in the minimization procedure									
	//
	//				if (doMinimization){
	//					if (!minimizeBB){ //side-chain minimization
	//						//We need this computation for the normal pairwise energy matrices and for unpruned RCs or pairs in the continuous energy matrix
	//						boolean needE = checkNeedE(res1,res1AANum,res1RotNum,res2,res2AANum,res2RotNum);
	//						if(needE)
	//							pairwiseRotamerEnergySidechainHelper (res1,res2,res1Strand,res1Resnum,res2Strand,res2Resnum,
	//									shellRun, retEMatrixMin, res1AANum, res1RotNum, res2AANum, res2RotNum);
	//					}
	//					else { //backbone minimization
	//						if (!doBackrubs) // phi/psi minimization
	//							pairwiseRotamerEnergyBackboneHelper (res1,res2,res1Strand,res1Resnum,
	//									shellRun, retEMatrixMin, false, res1AANum, res1RotNum, res2AANum, res2RotNum);
	//						else { //backrubs
	//							pairwiseRotamerEnergyBackrubsHelper (res1,res2,res1Strand,res1Resnum,
	//									shellRun, retEMatrixMin, false, res1AANum, res1RotNum, res2AANum, res2RotNum);
	//						}
	//					}
	//				}
	//				else if (computeEVEnergy){
	//
	//					double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
	//					double curEnergy = (double)energyTerms[0];//calcTotalSnapshotEnergy();
	//
	//					//Store the computed energy
	//					if (shellRun) {
	//						retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, curEnergy);
	//						//						retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, curEnergy);
	//					}
	//					else {
	//						retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, curEnergy );
	//						//						retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, curEnergy );
	//					}
	//				}
	//				else {
	//					System.out.println("This should not happen. No energy evaluation specified");
	//					System.exit(1);
	//				}					
	//
	//				m.updateCoordinates();
	//				m.revertPertParamsToCurState();
	//			}
	//			else {
	//				if (shellRun) {
	//					retEMatrixMin.setShellRotE( res1, res1AANum, res1RotNum, stericE );
	//					//					retEMatrixMax.setShellRotE( res1, res1AANum, res1RotNum, stericE );
	//				}
	//				else {
	//					retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, stericE );
	//					//					retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, stericE );
	//				}
	//			}
	//
	//			return;
	//		}
	//
	//		else {
	//			// Next check with different rotamers
	//			int depth = mutDepth[curMut];
	//			int str = strandMut.resStrand[depth];
	//			int strResNum = strandMut.resStrandNum[depth];
	//			int molResNum = strandMut.allMut[depth];
	//
	//			int numRot;
	//			if(doPerturbations)
	//				numRot = ((StrandRCs)strandRot[str]).getNumRCs(strResNum, curAANum[molResNum]);
	//			else
	//				numRot = strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum]);
	//
	//			if (numRot == 0) {//If there are no rotamers for this AA then allow the default conformation; this will only happen with Ala and Gly
	//				if (depth==res1)
	//					pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
	//							res1AANum, res2AANum, 0, res2RotNum, retEMatrixMin, minimizeBB, doBackrubs,
	//							numMut, mutDepth, curMut+1, validConf);
	//				else
	//					pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
	//							res1AANum, res2AANum, res1RotNum, 0, retEMatrixMin, minimizeBB, doBackrubs,
	//							numMut, mutDepth, curMut+1, validConf);
	//			}
	//			else {		
	//				for(int w=0;w<numRot;w++) {
	//
	//
	//					if(doPerturbations)
	//						validConf = validConf && ((StrandRCs)strandRot[str]).applyRC(m, strandMut.resStrandNum[depth], w);
	//					else
	//						strandRot[str].applyRotamer(m, strandMut.resStrandNum[depth], w);
	//
	//					if (depth==res1)
	//						pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
	//								res1AANum, res2AANum, w, res2RotNum, retEMatrixMin, minimizeBB, doBackrubs,
	//								numMut, mutDepth, curMut+1, validConf);
	//					else
	//						pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
	//								res1AANum, res2AANum, res1RotNum, w, retEMatrixMin, minimizeBB, doBackrubs,
	//								numMut, mutDepth, curMut+1, validConf);
	//				}
	//			}
	//		}
	//	}	
	//This method computes the minimized energy for a given pair of rotamers (or rot-to-template)
	//Called by pairwiseMutationAllRotamerSearch(.)
	//	private void pairwiseRotamerEnergySidechainHelper (int res1,int res2,int res1Strand,int res1Resnum,int res2Strand,int res2Resnum,
	//			boolean shellRun, Emat retEMatrixMin, int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){
	//
	//		double minEnergy = (double)Math.pow(10,30);
	//		double maxEnergy = -(double)Math.pow(10,30);
	//
	//		double curEnergy = 0.0f;
	//		double beginE = 0.0f;
	//
	//
	//		if(useCCD)//The molecule has ideal values for the dihedrals now; record these
	//			((ContSCObjFunction)ccdMin.objFcn).updateIdealDihedrals();
	//
	//
	//		if(minimizePairwise){//ef will be the total energy we are minimizing (e.g. pairwise = pair energy - intra energies)
	//
	//			beginE = (double)efunc.getEnergy();
	//			ccdMin.minimize();
	//			curEnergy = (double)efunc.getEnergy();
	//
	//			if(compCETM){
	//
	//				boolean stericsOK = true;
	//				if(eliminatedRotAtRes==null){
	//					//no pruning information
	//					//so we'll do an extra steric check here instead
	//					if(res2==-1){
	//						if( curEnergy >= templateSt )
	//							stericsOK = false;
	//					}
	//					else if( curEnergy >= pairSt  )
	//						stericsOK = false;
	//				}
	//
	//
	//				if(stericsOK)
	//					compEPICFit( res1Strand, res1Resnum, res1, res1AANum, res1RotNum,
	//							res2Strand, res2Resnum, res2, res2AANum, res2RotNum, curEnergy, beginE, shellRun );
	//			}
	//		}
	//		else{
	//			/////////////////////////////
	//			//formally making sure that when minimization is performed,
	//			//		the minimized value is never greater than the initial value
	//			double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
	//			beginE = (double)energyTerms[0];//calcTotalSnapshotEnergy();
	//
	//			//Minimize
	//			if(useCCD)
	//				ccdMin.minimize();
	//			else
	//				simpMin.minimize(numMinSteps);
	//
	//			energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
	//			curEnergy = (double)energyTerms[0];//calcTotalSnapshotEnergy();
	//
	//			if (shellRun && doDihedE){ //add dihedral energies
	//				if(useCCD)
	//					curEnergy += ((ContSCObjFunction)ccdMin.objFcn).de.getEnergy();
	//				else
	//					curEnergy += simpMin.computeDihedEnergy();
	//			}
	//		}
	//
	//		/////////////////////////////
	//
	//		//Compare to the min and max energies found so far and update, if necessary
	//		//With the dihedral energies included during minimization this is unnecessary
	//		if(doDihedE){
	//			if (beginE<curEnergy)
	//				curEnergy = beginE;
	//			double lE = Math.min(beginE, curEnergy);
	//			double hE = Math.max(beginE, curEnergy);
	//		}
	//		minEnergy = Math.min(minEnergy,curEnergy);
	//		maxEnergy = Math.max(maxEnergy,beginE);		
	//
	//
	//		m.updateCoordinates();//restore the actualCoordinates array to the initial values
	//		m.revertPertParamsToCurState();
	//
	//		if (shellRun) {
	//			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
	//			//			retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
	//		}
	//		else {
	//			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
	//			//			retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
	//		}
	//
	//		return;
	//	}

	private boolean checkNeedE( int res1, int res1AANum, int res1RotNum, int res2, int res2AANum, int res2RotNum ){
		//Returns if we need to compute the given shell or pairwise energy
		if(!compCETM)//Need all energies for normal energy matrices
			return true;


		//		if (eliminatedRotAtRes == null){//no pruning, so try to compute everything
		//			//(we'll do an extra steric check in pairwiseRotamerEnergySidechainHelper in this case)s
		//			return true;
		//		}


		//Now check that the RC(s) and, if applicable, pair are not pruned
		//If compCETM then we should have eliminatedRotAtRes and splitFlags available
		if( arpMatrix.getSinglePruned(res1,res1AANum,res1RotNum) )
			return false;

		if(res2!=-1){//Pairwise energy
			if( arpMatrix.getSinglePruned(res2,res2AANum,res2RotNum) )
				return false;
			if(compCETM) {
				//if we've done pairs pruning and are considering CET computation
				if( arpMatrix.getPairPruned(res1,res1AANum,res1RotNum,res2,res2AANum,res2RotNum) )
					return false;
			}
		}

		return true;
	}



	private void compEPICFit(Residue firstRes, Residue secondRes, AARotamerType aaType1, AARotamerType aaType2, EMatrixEntryWIndex emeWI, 
			double curEnergy, double beginE, boolean shellRun){
		//Do the EPIC series fits for this rotamer pair


		ContSCObjFunction of = (ContSCObjFunction)ccdMin.objFcn;
		DoubleMatrix1D startVec = of.curDOFVals.copy();

		double baseE = of.getValue(startVec);

		EnergyFunction oldEF = of.efunc;
		if(es.PBTest||es.quantumTest){//we'll fit to a Poisson-Boltzmann or quantum energy function for test purposes
			if(es.PBTest)
				of.efunc = new PBChangeEFcn(m,firstRes,secondRes);
			else{
				int res2 = -1;
				if(!shellRun)
					res2 = emeWI.pos2();
				of.efunc = new QuantumChangeEFcn(of.curDOFVals,emeWI.pos1(),emeWI.pos2(),oldEF);
			}
			baseE = 0;
		}

		EPICFitter fitter = new EPICFitter(ccdMin,es,startVec,baseE);

		EPoly bestSer = null;

		if(ccdMin.numDOFs>0){//need to actually compute polynomial fits

			//baseE should match curEnergy
			//(difference in calculations is efunc.getEnergy(0) vs. getEnergy()--should agree 
			//for all pairwise or shell E's, unless no flexibility and thus partial energy doesn't include
			//the residue or pair of interest)
			if( (Math.abs(baseE-curEnergy)>0.001) && !(es.PBTest||es.quantumTest) ){
				System.out.println("Warning: baseE="+baseE+" but curEnergy="+curEnergy+" in RotamerSearch.compEPICFit!");
			}

			FitParams curFitParams = FitParams.quadratic(of.getNumDOFs(),false);
			double bestResid = Double.POSITIVE_INFINITY;
			int bestBound = -1;
			ArrayList<EPoly> series = new ArrayList<>();//list of fit series

			for( int boundCount=0; bestResid>es.EPICGoalResid && curFitParams!=null; boundCount++ ){

				System.out.println("FIT NUMBER "+boundCount+": "+curFitParams.getDescription());

				EPoly curSeries;
				try{
					curSeries = fitter.doFit(curFitParams);
				}
				catch(Exception e){//sometimes singular matrices will arise during fitting
					//if so we skip that order.  More SVE etc. might help
					System.err.println("Fit failed: "+e.getMessage());
					e.printStackTrace();
					series.add(null);
					continue;
				}

				double meanResid = fitter.crossValidateSeries(curSeries,curFitParams);
				series.add(curSeries);

				if( curFitParams.order==2 && curFitParams.PCOrder<=2 && 
						curFitParams.explicitVDWCutoff==0 )//pure quadratic--to use as template for PCs
					fitter.PCTemplate = curSeries;

				if(meanResid<bestResid){
					bestResid = meanResid;
					bestBound = boundCount;
				}

				curFitParams = fitter.raiseFitOrder(curFitParams);
			}

			if(bestResid > es.EPICGoalResid){
				System.out.println("Failed to reach goal residual.");
			}
			System.out.println("Best residual: "+bestResid+" for bound number "+bestBound);



			int numDOFs = fitter.numDOFs;

			//TESTING FITS
			//KER: Printing similar information as before, could be removed to reduce array indexing
			int aaInd1=-1;int rotInd1=-1;int aaInd2=-1;int rotInd2=-1;int pos2=-1;
			if(shellRun){
				aaInd1 = aaType1.index;
				rotInd1 = m.strand[firstRes.strandNumber].rcl.getRC(((RotamerEntry)emeWI.eme).r.rotamers[0]).rot.aaIndex; //0 index assumes no supRot
			}else{
				aaInd1 = aaType1.index;
				rotInd1 = m.strand[firstRes.strandNumber].rcl.getRC(((RotamerPairEntry)emeWI.eme).r1.rotamers[0]).rot.aaIndex; //0 index assumes no supRot
				pos2 = emeWI.pos2();
				aaInd2 = aaType2.index;
				rotInd2 = m.strand[secondRes.strandNumber].rcl.getRC(((RotamerPairEntry)emeWI.eme).r2.rotamers[0]).rot.aaIndex; //0 index assumes no supRot
			}

			System.out.println("ROT:"+emeWI.pos1()+" "+aaInd1+" "+rotInd1+" "+pos2+" "+aaInd2+" "+rotInd2 + "(" + emeWI.toString()+ ")");
			System.out.println("curEnergy: "+curEnergy +" beginE: "+beginE);

			double testScales[] = new double[] { 0.01, 0.5, 5, 100 };//100
			int samplesPerScale = 3;



			double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
			double relMin[] = new double[numDOFs];
			for(int dof=0; dof<numDOFs; dof++){
				relMax[dof] = ccdMin.DOFmax.get(dof) - startVec.get(dof);
				relMin[dof] = ccdMin.DOFmin.get(dof) - startVec.get(dof);
			}


			for(double scale : testScales){
				for(int s=0; s<samplesPerScale; s++){

					//Generate vector relative to minimum
					double dx[] = new double[numDOFs];
					//and absolute
					DoubleMatrix1D sampAbs = DoubleFactory1D.dense.make(numDOFs);
					for(int dof=0; dof<numDOFs; dof++){
						double top = Math.min(relMax[dof], scale);
						double bottom = Math.max(relMin[dof], -scale);
						dx[dof] = bottom + fitter.rand.nextDouble()*(top-bottom);
						sampAbs.set(dof, startVec.get(dof)+dx[dof]);
					}

					double trueVal = of.getValue(sampAbs) - baseE;

					System.out.print("TEST: S="+scale+" dx=");
					for(int dof=0; dof<numDOFs; dof++)
						System.out.print(dx[dof]+" ");

					System.out.print("TRUE="+trueVal+" SER=");

					for(EPoly b : series){
						if(b!=null)
							System.out.print(b.evaluate(sampAbs,false,false)+" ");
					}

					System.out.println();
				}
			}

			bestSer = series.get(bestBound);
		}
		else
			bestSer = fitter.blank();

		bestSer.minE = curEnergy;//store voxel-minimum energy.  curEnergy not baseE (which =0 w/o flexibility)

		int pos2 = -1;
		if(!shellRun)
			pos2 = emeWI.pos2();

		storeDiscreteDOFs(emeWI.pos1(),aaType1,pos2,aaType2,bestSer);
		
		
		if (shellRun)
			cetm.setShellRotE(emeWI.index, bestSer);
		else
			cetm.setPairwiseE( emeWI.index, bestSer );    
		
		
		if(es.PBTest||es.quantumTest)//revert to regular energy function
			((ContSCObjFunction)ccdMin.objFcn).efunc = oldEF;
	}



	private void storeDiscreteDOFs(int res1, AARotamerType res1AA, int res2, AARotamerType res2AA, EPoly bestSer){
		//Store discrete DOFs in a continuous energy matrix
		//We know that res1 and res2 (if not -1) are in RC-assigned states
		//so we store the values of discrete DOFs affecting them
		//(for sequence we do mutation rather than AA type DOFs though)
		//(so the discrete DOFs can be either perturbations or mutations)
		bestSer.discreteDOFNums = new ArrayList<>();
		bestSer.discreteDOFVals = new ArrayList<>();

		for(DegreeOfFreedom dof : m.DOFs){
			if( arrayContains( dof.flexResAffected, res1 )
					|| arrayContains( dof.flexResAffected, res2 ) ){

				if(dof.type==DegreeOfFreedom.PERTURBATION){
					if(!dof.isContinuous()){//this is the opposite of the condition for inclusion in a ContSCObjFunction...
						bestSer.discreteDOFNums.add(dof.DOFNum);
						bestSer.discreteDOFVals.add(m.perts[dof.pertNum].curParam);//dof.pert.curParam);
					}
				}
				else if(dof.type==DegreeOfFreedom.MUTATION){
					//always discrete
					bestSer.discreteDOFNums.add(dof.DOFNum);

					//see if the mutation is active, record value accordingly
					if(  ( dof.flexResAffected[0]==res1 && res1AA.index==dof.mutAANum.index )
							|| ( dof.flexResAffected[0]==res2 && res2AA.index==dof.mutAANum.index ) ){

						bestSer.discreteDOFVals.add(1.);
					}
					else
						bestSer.discreteDOFVals.add(0.);
				}
			}
		}
	}


	public static boolean arrayContains(int[] a, int b){//Does a contain b?
		for(int ae : a){
			if( b == ae )
				return true;
		}
		return false;
	}


	//Adjusting one or a few elements of a vector
	//(functions made for numerical differentiation)
	public static DoubleMatrix1D adjustVec(DoubleMatrix1D a, int b, double c){
		DoubleMatrix1D ans = a.copy();
		ans.set( b, ans.get(b) + c );
		return ans;
	}

	public static DoubleMatrix1D adjustVec(DoubleMatrix1D a, int b1, double c1, int b2, double c2){
		DoubleMatrix1D ans = a.copy();
		ans.set( b1, ans.get(b1) + c1 );
		ans.set( b2, ans.get(b2) + c2 );
		return ans;
	}

	public static DoubleMatrix1D adjustVec(DoubleMatrix1D a, int b1, double c1, int b2, double c2, int b3, double c3){
		DoubleMatrix1D ans = a.copy();
		ans.set( b1, ans.get(b1) + c1 );
		ans.set( b2, ans.get(b2) + c2 );
		ans.set( b3, ans.get(b3) + c3 );
		return ans;
	}



	//This method helps compute energy bounds for a given pair of rotamers for backbone phi/psi minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackboneHelper (int res1,int res2,int res1Strand,int res1Resnum,
			boolean shellRun, Emat retEMatrixMin, boolean templateOnly,
			int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){		


		double minEnergy = (double)Math.pow(10,30);
		double maxEnergy = -(double)Math.pow(10,30);

		double beginE = 0.0f;

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		beginE = (double)energyTerms[0];//calcTotalSnapshotEnergy();


		bbMin.minimizeFull(true); //minimize

		energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		double curEnergy = (double)energyTerms[0];//calcTotalSnapshotEnergy();



		if (beginE<curEnergy)
			curEnergy = beginE;
		double lE = Math.min(beginE, curEnergy);
		double hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	

		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		/////////////////////////////////////////////

		if (templateOnly){
			retEMatrixMin.setShellShellE(minEnergy);
			//			retEMatrixMax.setShellShellE(maxEnergy);
		}
		else if (shellRun) {
			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
			//			retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
		}
		else {
			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
			//			retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
		}

		return;
	}
	//This method helps perform sampling for a given pair of rotamers for backrubs minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackrubsHelper (int res1,int res2,int res1Strand,int res1Resnum, boolean shellRun, Emat retEMatrixMin, boolean templateOnly,
			int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){


		double minEnergy = (double)Math.pow(10,30);
		double maxEnergy = -(double)Math.pow(10,30);

		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		double beginE = (double)energyTerms[0];//calcTotalSnapshotEnergy();

		double initActualCoords[] = m.getActualCoords(); //this is necessary for performing brMin.applyMaxBackrub() after brMin.minimizeFull()
		brMin.minimizeFull(shellRun,templateOnly); //minimize

		energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		double curEnergy = (double)energyTerms[0];//calcTotalSnapshotEnergy();
		m.setActualCoords(initActualCoords);

		brMin.applyMaxBackrub();
		energyTerms = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		double maxE = (double)energyTerms[0];//calcTotalSnapshotEnergy();
		beginE = (double)Math.min(beginE,maxE); //max energy will be the min of (initial energy) and (max sterically-allowed energy from brMin)

		if (beginE<curEnergy)
			curEnergy = beginE;
		double lE = Math.min(beginE, curEnergy);
		double hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	

		m.updateCoordinates();//restore the actualCoordinates array to the initial values from the atoms coord[]
		/////////////////////////////////////////////

		if (templateOnly){
			retEMatrixMin.setShellShellE(minEnergy);
			//			retEMatrixMax.setShellShellE(maxEnergy);
		}
		else if (shellRun) {
			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
			//			retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
		}
		else {
			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
			//			retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
		}

		return;
	}


	private void correctSwitchedTemplates(MutableResParams strandMut, Emat retEMatrixMin){
		//This function takes the changes in template energy due to a full structure switch's template changes,
		//which are stored in switchedTemplateE, and transfers them to the intra energies of RCs at the first flexible residue
		//that are in the corresponding perturbation states.  There, they are used for pruning and can be efficiently
		//incorporated into A*.  The two formulations of the energy matrix are biophysically equivalent.


		if( ! m.perts[0].type.equalsIgnoreCase("FULL STRUCTURE SWITCH") )
			return;

		SinglesIterator iter = retEMatrixMin.singlesIterator(0); //Add the energy to the first position for every RC
		Residue r = m.residue[retEMatrixMin.resByPos.get(0).get(0)]; // Second 0 index assumes no super-rotamers
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			int globalRC = retEMatrixMin.singles.getRot(emeWI.index)[0]; //0 index assumes no super-rotamers
			ResidueConformation rc = ((StrandRCs)strandRot[r.strandNumber]).rcl.getRC(globalRC);
			int resPertState = rc.pertState;
			int pertState = r.pertStates[resPertState][0];//Corresponding state of the full structure switch
			double switchE = switchedTemplateE[pertState];
			retEMatrixMin.addToShellRotE(emeWI.index, switchE);
		}

	}


	//// END PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION
	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//// BEGIN ROTAMER SEARCH SECTION

	//For the given mutation sequence,
	//Count the total number of conformations, the number of confs pruned by MinDEE, the
	//	number of remaining conformations; the number of conformations above each level
	//	(from level i+1 to the top level, which is the ligand); the number of rotamers for each flexible residue;
	//	the number of rotamers (total and non-pruned by MinDEE) for the flexible residues only;
	//The MinDEE matrix eliminatedRotAtRes[] should already be computed;
	//For each pruned rotamer at curLevel, count the number of *new* conformations that 
	//	are pruned as a result, then decrease the number of rotamers for curLevel to
	//	make sure we do not overcount the pruned conformations (the next pruned rotamer
	//	should not count the conformations that include the first rotamer, and therefore
	//	have already been pruned)
	private int countConfs(int numRes, int[] prunedStericByRes, int numRotForRes[], int numRotForResNonPruned[]){

		int numTotalRotRed = 0;
		numConfsTotal = BigInteger.valueOf(1);
		//numConfsAboveLevel = new BigInteger[numRes];

		int curNumRot;

		//Store the number of rotamers for each AA in the current mutation
		int[] tmpRotNonPruned = arpMatrix.remainingRot();
		int[] tmpRotPerPos = arpMatrix.numRotPerPos();
		for(int i=0; i<tmpRotNonPruned.length;i++){
			numRotForResNonPruned[i] = tmpRotNonPruned[i];
			numRotForRes[i] = tmpRotPerPos[i];
		}

		//Set the number of rotamers for each position in the Emat

		//int[] numRotForResPruned = new int[numRes];
		int[] numRotAfterSteric = new int[numRes];
		BigInteger numConfsAfterSteric = BigInteger.ONE;
		numConfsLeft = BigInteger.ONE;
		numConfsTotal = BigInteger.ONE;
		for(int i=0; i<numRotForRes.length;i++){
			numTotalRotRed += numRotForRes[i];

			numRotAfterSteric[i] = numRotForRes[i] - prunedStericByRes[i];



			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[i]));
			numConfsLeft = numConfsLeft.multiply(BigInteger.valueOf(numRotForResNonPruned[i]));
			numConfsAfterSteric = numConfsAfterSteric.multiply(BigInteger.valueOf(numRotAfterSteric[i]));

		}

		numConfsPrunedByS = numConfsTotal.subtract(numConfsAfterSteric);

		numConfsPrunedByMinDEE = numConfsTotal.subtract(numConfsLeft).subtract(numConfsPrunedByS);


		return numTotalRotRed;

	}

	//Computes the number of conformations pruned by a rotamer (eliminated by MinDEE) at curLevel:
	//	multiply the number of rotamers for each level different from curLevel
	private BigInteger compPrunedConfsByRot(int numRotForRes[], int numRes, int curLevel){

		BigInteger numPruned = BigInteger.ONE;

		for (int i=0; i<numRes; i++){
			if (i!=curLevel)
				numPruned = numPruned.multiply(BigInteger.valueOf(numRotForRes[i]));
		}

		return numPruned;
	}

	//Computes the number of conformations that are pruned by MinDEE due to steric clashes:
	//	for each level, count the number of rotamers that are not pruned by MinDEE due to steric clash, then 
	//	multiply for all levels; the result is the number of conformations that do not include
	//	pruned rotamers; return (totalNumConfs - this number) as the number of conformations pruned by MinDEE
	//	private BigInteger countPrunedByMinDEESteric(int numMutable, int numRes, /*boolean ligPresent,*/ 
	//			MutableResParams strandMut, int numRotForRes[], Emat emat){
	//
	//		BigInteger numConfNotPruned = BigInteger.ONE;
	//		int numNonPruned[] = new int[numRes];
	//		for (int i=0; i<numRotForRes.length; i++)
	//			numNonPruned[i] = numRotForRes[i];
	//
	//		for (int curLevel=0; curLevel<numRes; curLevel++){
	//			int str = strandMut.resStrand[curLevel];
	//			int strResNum = strandMut.resStrandNum[curLevel];
	//			int molResNum = m.residue[strandMut.allMut[curLevel]].moleculeResidueNumber;
	//			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
	//				if ((eliminatedRotAtRes.get(curLevel,curAANum[molResNum],curRot))&&(prunedIsSteric.get(curLevel,curAANum[molResNum],curRot))){
	//					numNonPruned[curLevel]--;
	//				}
	//			}
	//		}
	//
	//		for (int curLevel=0; curLevel<numRes; curLevel++)
	//			numConfNotPruned = numConfNotPruned.multiply(BigInteger.valueOf(numNonPruned[curLevel]));
	//
	//		return (numConfsTotal.subtract(numConfNotPruned));
	//	}
	//Sets-up the repeat mutation search if the desired accuracy has not been achived, although
	//	all non-pruned conformations have been generated by A*;
	//Updates eliminatedRotAtRes[] to not prune a subset of the pruned rotamers, such that the
	//	number of conformations not pruned by this is at least numPrunedConfsToAllow
	//	private void setupRepeatRun(BigInteger numPrunedConfsToAllow, int numRotForResNonPruned[], 
	//			int numLevels, int numMutable){
	//
	//		BigInteger totalNumConfsUnpruned = new BigInteger("0");
	//		Iterator<RotInfo<Boolean>> elimRotIter = eliminatedRotAtRes.iterator();
	//		Iterator<RotInfo<Boolean>> prunedStericIter = prunedIsSteric.iterator();
	//		while(elimRotIter.hasNext()){
	//			//for (int i=0; i<eliminatedRotAtRes.length; i++){
	//			RotInfo<Boolean> elimRotInfo = elimRotIter.next();
	//			RotInfo<Boolean> prunedStericInfo = prunedStericIter.next();
	//			//check to make sure we are always checking same index
	//			assert (elimRotInfo.curPos == prunedStericInfo.curPos && elimRotInfo.curAA == prunedStericInfo.curAA &&
	//					elimRotInfo.curRot == prunedStericInfo.curRot) : "ElimRot and PrunedSteric indexes don't match";
	//
	//			if ((elimRotInfo.state)&&(!prunedStericInfo.state)){ //pruned non-steric rotamer
	//				//if ((eliminatedRotAtRes[i])&&(!prunedIsSteric[i])){ //pruned non-steric rotamer
	//
	//				int curLevel = elimRotInfo.curPos;//(int)Math.floor(i/numTotalRotamers); //the residue number for the cur rot index
	//				//TODO: fix that all positions are assumed to have "numTotalRotamers" number of rotamers
	//				/*if (curLevel>=numInAS) //ligand level (the ligand may have more than numTotalRotamers rotamers)
	//					curLevel = numInAS; //the residue number for the cur rot index*/
	//				numRotForResNonPruned[curLevel]++;
	//				eliminatedRotAtRes.set(elimRotInfo,false);
	//
	//				BigInteger numConfsAtLevelUnpruned = new BigInteger("1"); //count the unpruned confs by unpruning the cur rot index
	//				for (int j=0; j<numLevels; j++){
	//					if (j!=curLevel){
	//						numConfsAtLevelUnpruned = numConfsAtLevelUnpruned.multiply(BigInteger.valueOf(numRotForResNonPruned[j]));
	//					}
	//				}
	//
	//				totalNumConfsUnpruned = totalNumConfsUnpruned.add(numConfsAtLevelUnpruned);
	//				if (totalNumConfsUnpruned.compareTo(numPrunedConfsToAllow)>=0) // num pruned confs reduced by the necessary number
	//					break;
	//			}
	//		}
	//	}
	// This function performs a rotamer search to compute
	//  a partition function for the slave node.
	// Utilizes a number of helper functions
	public AStarResults slaveDoRotamerSearch(int runNum, boolean searchComputeEVEnergy, boolean searchDoMinimization, int numInAS, 
			MutableResParams strandMut, boolean usingInitialBest, BigDecimal initialBest,
			CommucObj cObj, boolean minimizeBB, boolean doBackrubs, String backrubFile,
			SaveConfsParams saveConfsParams, int curMut,boolean useMaxKSconfs, BigInteger maxKSconfs,
			int[] prunedStericPerPos, double Ival, Settings.Enum enumSettings) {

		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer

		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
//		ASAANums = new int[numberMutable];
//		curStrRotNum = new int[numberMutable];
		//curASRotNum = new int[numInAS];
		//int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
		// This map maps the system strand residues back to the AS numbering
		// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		//curLigRotNum = 0;
		numConfsPrunedByE = BigInteger.ZERO;
		numConfsPrunedByS = BigInteger.ZERO;
		numConfsEvaluated = BigInteger.ZERO;
		numConfsPrunedByMinDEE = BigInteger.ZERO;
		allPruned = false;
		//confEnergies = new Vector(2048,1024);
		setBestE(9999999.0f);
		bestEUnMin = 9999999.0f;
		if (usingInitialBest)
			initial_q = initialBest.multiply(new BigDecimal((double)(1-KSepsilon)));
		else
			initial_q = new BigDecimal(0.0);
		partial_q = new BigDecimal(0.0);
		partial_p = new BigDecimal(0.0);

		for(int resID:strandMut.allMut){
			m.residue[resID].flexible = true;
		}


		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return null;
		}

		// Prepare Amber
		if(searchComputeEVEnergy){
			// Amber should already be loaded
			/*if(ligPresent) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}*/
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization or DEEPer
					if (simpMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return null;
					}
					bbMin = null;
					//brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){
						if (bbMin == null) {
							System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
							return null;
						}
						simpMin = null;
						//brMin = null;
					}
					/*else {
						if (brMin == null) {
							System.out.println("Warning: Attempting minimization run but brMin not allocated, RotamerSearch aborting");
							return;
						}
						simpMin = null;
						bbMin = null;
					}*/
				}
			}
		}

		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null){
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return null;
		}

		//		if (eliminatedRotAtRes == null){
		//			System.out.println("Warning: MinDEE matrix not computed");
		//			return null;
		//		}


		//if(useCCD || es.useEPIC)
		//    m.DOFs = DegreeOfFreedom.makeDOFArray(strandRot, strandMut, mutRes2Strand, mutRes2StrandMutIndex, m, null);
		if(es.useEPIC)//For EPIC we will refer to DOFNums in the CETMatrix, so DOFNums must match cetm.DOFList
			//if we recomputed m.DOFs here there might be fewer mutations available so the DOFNums might come out different
			m.DOFs = cetm.DOFList;
		else if(useCCD)//we have no cetm and we need a DOF array, but we don't need to match DOFNums, so just recompute the DOF array
			m.DOFs = DegreeOfFreedom.makeDOFArray(strandRot, strandMut, m);


		// Setup the residue number to AS number map
		/*for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			curResToASMap[i] = -1;
		}
		for(int i=0;i<residueMap.length;i++){
			curResToASMap[residueMap[i]] = i;
		}*/

		/*if (ligPresent)
			slaveMutationRotamerSearch(-1, numInAS, strandMut, AAList, numAAAllowed,
				numTotalRotamers, rotamerIndexOffset,	curResToASMap, ligPresent, minimizeBB, saveConfs, fName, doBackrubs, backrubFile);
		else*/
		AStarResults asr = slaveMutationRotamerSearch(runNum, 0, numberMutable, strandMut, minimizeBB, 
				doBackrubs, backrubFile,saveConfsParams, 
				curMut,useMaxKSconfs, maxKSconfs,prunedStericPerPos, Ival,enumSettings);

		// Store results to the communication object
		if (cObj!=null){
			
			cObj.searchNumConfsTotal[runNum] = numConfsTotal.intValue();
			if (allPruned) { //all of the conformations were pruned by MinDEE, as there is a residue with no remaining rotamers
				cObj.searchNumConfsPrunedByS[runNum] = 0;
				cObj.searchNumPrunedMinDEE[runNum] = numConfsTotal.intValue();
				cObj.searchNumConfsEvaluated[runNum] = 0;
				cObj.searchNumConfsLeft[runNum] = 0;
			}
			else {
				cObj.searchNumConfsPrunedByS[runNum] = numConfsPrunedByS.intValue();
				cObj.searchNumPrunedMinDEE[runNum] = numConfsPrunedByMinDEE.intValue();
				cObj.searchNumConfsEvaluated[runNum] = numConfsEvaluated.intValue();
				cObj.searchNumConfsLeft[runNum] = numConfsLeft.intValue();
			}
			
			cObj.allPruned[runNum] = allPruned;

			// Compute q_X
			cObj.q[runNum] = partial_q;
			cObj.bestE[runNum] = (double)bestEUnMin;
			cObj.bestEMin[runNum] = (double)getBestE();
		}
		else {
			System.out.println("Statistics (unbound):");
			System.out.println("Best Energy:  " + (double)getBestE());
			System.out.println("partial_q: " + partial_q);
			System.out.println("partial_p: " + partial_p);
			System.out.println("NumConfsTotal:      	" + numConfsTotal);
			System.out.println("NumConfsPrunedByMinDEE: " + numConfsPrunedByMinDEE);
			System.out.println("NumConfsPrunedByS:  	" + numConfsPrunedByS);
			System.out.println("NumConfsEvaluated:  	" + numConfsEvaluated);
			System.out.println("NumConfsLeft:       	" + numConfsLeft);
		}

		return asr;

	}

	// This function does a mutation search then for each allowable mutation
	//  a simple rotamer search is performed
	private AStarResults slaveMutationRotamerSearch(int runNum, int depth, int maxDepth,  
			MutableResParams strandMut, boolean minimizeBB,
			boolean doBackrubs, String backrubFile, 
			SaveConfsParams saveConfsParams, int curMut,
			boolean useMaxKSconfs, BigInteger maxKSconfs, int[] prunedStericByPos, double Ival,
			Settings.Enum enumSettings) {


		// If we've arrived here then we're ready to
		//  do a rotamerSearch
		if(debug)
			System.out.println("One Mutation Conf Found and is being tested");
		if (computeEVEnergy){
			for(Residue r:m.residue){
				r.setEnergyEval(true, true);
			}
			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);

			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					//if(ligPresent)
					if(useCCD){
						efunc = new ForceFieldEnergy(m, a96ff);
						ContSCObjFunction of = new ContSCObjFunction(m,numberOfStrands,efunc,strandRot,doDihedE,null);
						efunc = of.efunc;//ef will now include dihedral energies if appropriate
						ccdMin = new CCDMinimizer(of,false);
					}
					else
						simpMin.initialize(m,numberOfStrands,a96ff,strandRot,doDihedE);
					//else
					//	simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
				}
				else { //backbone minimization
					if (!doBackrubs){
						//if (ligPresent)
						//	bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
						//else
						bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
					}
					else {
						brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
					}
				}
			}
		}

		// Determine the following
		// -AANums of all AS residues (so that we can index properly into rotamerIndexOffset, others)
		//   residue i (i=0..8) has 3 letter code AAList(AANums[i])
		// -total number of conformations for this mutation -> assign to numConfsTotal, numConfsLeft
		// -number of conformations below each level (ie. if we prune one of the 3rd AAs rotamers
		//   because of steric overlap, how many confs are we pruning)
		assignAANums(strandMut);
		//int numConfs = computeTotalNumConfs(residueMap,curResToASMap,ligPresent);

		// also computes num of confs at each level
		//numConfsLeft = (numConfsTotal-numConfsPrunedByMinDEE); //MinDEE should have already been applied

		//if doing EPIC and checking thresholds need the lowest bound without polynomial fits
		es.lowestBound = 0;
		if(es.useEPIC && es.enforceEPICThresh){
			System.out.println("Getting lowest-bound energy...");
			es.gettingLowestBound = true;
			es.lowestBound = Double.NaN;//this value will be used if no conformations are available
			slaveRotamerSearchAStar(maxDepth, strandMut, minimizeBB, doBackrubs,
					saveConfsParams, useMaxKSconfs, maxKSconfs,prunedStericByPos, Ival, enumSettings);
			es.gettingLowestBound = false;
			System.out.println("Got lowest bound: "+es.lowestBound);
		}
		AStarResults asr = null;
		if(!Double.isNaN(es.lowestBound)){//if no conformations available for lower bound, no point in doing full A*

			long AStarStartTime = System.currentTimeMillis();

			//Perform A* search: compute the partial partition function q*
			asr = slaveRotamerSearchAStar(maxDepth, strandMut, minimizeBB, doBackrubs,
					saveConfsParams, useMaxKSconfs, maxKSconfs,prunedStericByPos, Ival,enumSettings);

			System.out.println("Partition function time (ms): "+(System.currentTimeMillis()-AStarStartTime));

		}
		MSAStarSearch = null;
		return asr;

	}
	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by slaveMutationRotamerSearch(.)
	private AStarResults slaveRotamerSearchAStar(int numMutable, MutableResParams strandMut, boolean minimizeBB,
			boolean doBackrubs, SaveConfsParams saveConfsParams, boolean useMaxKSconfs, BigInteger maxKSconfs,
			int[] prunedStericByPos, double Ival,Settings.Enum enumSettings){


		if(doPerturbations)//Make sure minimizer is set properly
			((PMinimizer)simpMin).minimizePerturbations = minimizePerturbations;

		if(saveConfsParams.saveTopConfs || saveConfsParams.printTopConfs)
			topConfs = new PriorityQueue<ConfPair>(saveConfsParams.numTopConfs);

		ExpFunction ef = new ExpFunction();


		int treeLevels; // total num levels in the conformation tree
		/*if (ligPresent) //determine num tree levels: if ligPresent, then numInAS+1
			treeLevels = numInAS+1;
		else*/
		treeLevels = numMutable;

		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues


		//Count the total number of conformations, the number of conformations pruned by MinDEE,
		//	the remaining conformations; the total num rotamers for the flexible residues, the num
		//	rotamers (total and non-pruned by MinDEE) for each fllexible residue;
		numTotalRotRed = countConfs(numMutable, prunedStericByPos, numRotForRes, numRotForResNonPruned);

		BigInteger k_const = numConfsPrunedByMinDEE;

		//KER: numConfsPrunedMinDEE is now correct so we don't have to subtract out the Steric part
		//		BigInteger numConfsPrunedMinDEESteric = countPrunedByMinDEESteric(numMutable, treeLevels, /*ligPresent,*/ 
		//				strandMut, numRotForRes, prunedIsSteric);
		//
		//		k_const = k_const.subtract(numConfsPrunedMinDEESteric); //only the non-steric prunings are used in the computation of k_const

		//Bound the contribution of the conformations pruned by MinDEE
		BigDecimal pStar = ef.exp(-Ec_const/constRT).multiply(new BigDecimal(k_const),ExpFunction.mc);

		final double ro = (double)KSepsilon /(double)(1-KSepsilon);

		//		System.out.println("k_const: "+k_const+" pStar: "+printBigNum(pStar,5)+" numConfsPrunedMinDEESteric: "+numConfsPrunedMinDEESteric);

		//Count the number of non-pruned rotamers
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				allPruned = true;
				BigDecimal e = ef.exp(-Ec_const/constRT);
				if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ) { //some non-sterics pruned but accuracy not achieved, so search must be repeated
					BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro),ExpFunction.mc));

					BigDecimal f = psi.divide(e,ExpFunction.mc); //rounding is ROUND_HALF_UP
					BigInteger l_const = k_const.subtract( BigInteger.valueOf( (long)Math.ceil(f.doubleValue()) ) );
					//					setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
					//					repeatSearch = true;
				}

				return null;
			}
			else
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];
		}

		//logPS.println(Ec_const+" "+numConfsPrunedByMinDEE+" "+numConfsPrunedMinDEESteric+" "+k_const+" "+pStar+" "+numTotalRotRedNonPruned);


		boolean[] eliminatedRotAtPosRed = new boolean[numTotalRotRed];//reduced MinDEE matrix


		boolean splitFlagsRed[][] = null, tripleFlagsRed[][][] = null;
		if( useFlagsAStar ){
			splitFlagsRed = new boolean[numTotalRotRedNonPruned][numTotalRotRedNonPruned];
			if( useTriples )
				tripleFlagsRed = new boolean[numTotalRotRedNonPruned][][];//As with the normal tripleFlags, tripleFlagsRed[a][b][c] is only defined for a>b>c
		}



		//Get the reduced min energy matrix
		//		ReducedEnergyMatrix arpMatrixRed = arpMatrix.  reduceMatrix( eliminatedRotAtPosRed,
		//				numRotForRes, numRotForResNonPruned, treeLevels,
		//				numTotalRotRedNonPruned, numMutable, strandMut,
		//				numTotalRotRed, this, true, splitFlagsRed, tripleFlagsRed );

		Emat arpMatrixRed = arpMatrix.unprunedMatrix();

		//Set-up the A* search
//		MSAStarSearch = new PGAStar(treeLevels, numRotForResNonPruned, arpMatrix,false, es, doPerturbations,m, strandRot, strandMut,cetm); //false = no reordering
		//		MSAStar AStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,null,
		//				splitFlagsRed,tripleFlagsRed,doPerturbations,es,
		//				m,strandRot,strandMut);

		if(MSAStarSearch == null){
			switch(enumSettings.asMethod){
			case ORIG:
				MSAStarSearch = new PGAStar(treeLevels, numRotForResNonPruned, arpMatrix,false,es,doPerturbations,m, strandRot, strandMut, cetm); //false = no reordering
				break;
			case ASORIG:
			case ASLP:
			case ASWCSP:
			case ASMPLP:
				MSAStarSearch = new PGgurobiAStar(treeLevels,numRotForResNonPruned,arpMatrix,enumSettings.asMethod,enumSettings.varOrder,es,doPerturbations,m, strandRot, strandMut, cetm);
				break;
//			case MIN:
//				MSAStarSearch = new PGgurobiMinAStar(treeLevels,numRotForResNonPruned,arpMatrix,m,arpMatrix.resByPos,a96ff,simpMin,null,doDihedE,false,false);
//				break;
//			case ASMPLP:
//				MSAStarSearch = new PGMPLPAStar(treeLevels,numRotForResNonPruned,arpMatrix);
//				break;
			case BYSEQ:
				System.out.println("AS Method BYSEQ isn't applicable for K* search");
				break;
			case BYSUBROT:
				MSAStarSearch = new PGgurobiAStarBySubRot(treeLevels,numRotForResNonPruned,arpMatrix,bestEMin+Ival,m);
				break;
			case WCSP:
				MSAStarSearch = new WCSPSearch(treeLevels,numRotForResNonPruned,arpMatrix,Ival);
				break;
			default:
				System.out.println("Don't Recognize AStar method... using ASWCSP with reordering");
				MSAStarSearch = new PGgurobiAStar(treeLevels,numRotForResNonPruned,arpMatrix,Settings.ASTARMETHOD.ASWCSP,Settings.VARIABLEORDER.MINFSCORE,es,doPerturbations,m, strandRot, strandMut, cetm);
				break;
			}
		}

		if(es.useEPIC){//set up fit series for A* search
			if(es.gettingLowestBound)//run A* without EPIC
				MSAStarSearch.es = new EPICSettings();
			else{
				MSAStarSearch.DOFList = m.DOFs;
				MSAStarSearch.CETM = cetm;
			}
		}


		EMatrixEntryWIndex conf[] = new EMatrixEntryWIndex[treeLevels]; //the conformation with the actual rotamer numbers
		EMatrixEntryWIndex parentConf[] = null; //For subrotamer search				
		boolean run1 = true;
		double lowestBound = stericE;
		double minELowerBound = stericE;

		PGQueueNode curNode = new PGQueueNode(1, new int[1], 0.0, 0, 0);
		while (numConfsLeft.compareTo(BigInteger.ZERO)==1){

			///Run A*
			curNode = MSAStarSearch.doAStar(run1); //the current rotamer sequence
			conf = curNode.actualConf;
			parentConf = createParentConf(conf);

			if (conf == null){ // no valid conformations remaining
				if (partial_q.multiply(new BigDecimal(ro)).compareTo(pStar)<0){ //approximation accuracy not achieved
					BigDecimal e = ef.exp(-Ec_const/constRT);
					if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ){ //some non-sterics pruned but accuracy not achieved, so search must be repeated
						BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));

						BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP
						BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
						//							setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
						//							repeatSearch = true;
					}
				}

				return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
			}


			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			//			conf = getActualConf(curConf,arpMatrix,treeLevels,numRotForResNonPruned);

			m.backupAtomCoord();

			int[] curAANums = new int[treeLevels];
			int[] curRotNums = new int[treeLevels];
			for (int i=0; i<treeLevels; i++){

				curAANums[i] = conf[i].index[1];
				curRotNums[i] = conf[i].index[2];

				//TODO: allow all rotamer libraries the ability to mutate
				if(parentConf != null){ //Will only not be null if we want to compute it 
					parentConf[i].eme.applyRC(arpMatrix.resByPos, m);
					outPS.print(parentConf[i].eme.printRes(m,arpMatrix.resByPos));
				}
				else{
					conf[i].eme.applyRC(arpMatrix.resByPos, m);
					outPS.print(conf[i].eme.printRes(m,arpMatrix.resByPos));
				}
			}
			
			for(int i=0;i<treeLevels;i++)System.out.print(conf[i]+" ");System.out.println();

			//Check the energy of the conformation and compute the score if necessary
			minELowerBound = computeBestRotEnergyBound(curAANums,curRotNums,null);/*numTotalRotamers,rotamerIndexOffset*///));

			if(run1){
				lowestBound = minELowerBound;
				if(lowestOverallBound == Double.POSITIVE_INFINITY)// || subRotamers)
					lowestOverallBound = lowestBound;
				run1 = false;
			}

			if(es.useEPIC){

				if(es.gettingLowestBound){//just want to get the lowest bound,
					//without polynomial fits, from this A* tree and return
					es.lowestBound = minELowerBound;
					return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				}


				ContSCObjFunction cof = null;
				if(es.useSVE)//will need cof to set DOFs for sve
					cof = (ContSCObjFunction)ccdMin.objFcn;

				//Setup indices
				Index3[] indices = new Index3[conf.length];
				for(int i=0; i<conf.length;i++){
					indices[i] = conf[i].rot1index3();
				}

				CETObjFunction lof = cetm.getObjFunc(indices, false, false, cof);
				CCDMinimizer lsbMin = new CCDMinimizer(lof,false);
				DoubleMatrix1D bestVals = lsbMin.minimize();


				double LSBE = lof.getValue(bestVals);//compute scaled polynomial fits' contribution to energy bound
				minELowerBound += LSBE;//and this will be used as the energy too!

				//check validity of thresholds
				es.checkThresh1Validity(lof,bestVals);
				if(es.enforceEPICThresh){
					if(minELowerBound > es.lowestBound+es.EPICThresh2){
						System.err.println("ERROR: EPICThresh2="+es.EPICThresh2+" insufficient for "
								+" energy="+minELowerBound+" lowestBound="+es.lowestBound);
						System.exit(1);
					}
				}

				if(es.useSVE){
					//cof minimized m...revert to unminimized
					m.updateCoordinates();
					m.revertPertParamsToCurState();
				}
			}
			
			//double psi = Math.max(initial_q,partial_q);
			BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));

			//double curThreshold = -constRT * (Math.log(psi)+Math.log(ro/numConfsLeft));
			double curThreshold = stericE;
			BigDecimal diff_qp = psi.subtract(pStar,ExpFunction.mc);
			if (diff_qp.compareTo(new BigDecimal(0.0))<0) { //the contribution of the pruned confs is bigger than ro*partial_q, so the search cannot be halted

				BigDecimal qBound = partial_q.add(ef.exp(-minELowerBound/constRT).multiply(new BigDecimal(numConfsLeft),ExpFunction.mc),ExpFunction.mc); //an upper bound on what partial_q can be

				if (pStar.compareTo(initial_q.max(qBound.multiply(new BigDecimal(ro),ExpFunction.mc)))>0 && 
						(!useMaxKSconfs || (useMaxKSconfs && numConfsLeft.add(numConfsEvaluated).compareTo(maxKSconfs) < 0))){ //approximation accuracy cannot be achieved

					BigDecimal e = ef.exp(-Ec_const/constRT);
					if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ){ //some non-sterics pruned but accuracy not achieved, so search must be repeated						

						BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP
						BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
						//						setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
						//						repeatSearch = true;
						AStarResults asr = new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
						asr.status = AStarResults.INCREASEIVAL;
						return asr;
					}
					else { //approximation accuracy not achieved and nothing to un-prune
						System.out.println("ERROR: Approximation accuracy not achieved, but no rotamers to unprune..");
						System.exit(1);
					}
				}
				else //it may be possible to achieve approximation accuracy
					curThreshold = stericE;
			}
			else
				curThreshold = -constRT * (ef.log(diff_qp)-Math.log(numConfsLeft.doubleValue()));

			System.out.println("conf: "+numConfsEvaluated.add(BigInteger.ONE)+" minELowerBound: "+minELowerBound+" curThreshold: "+curThreshold);
			System.out.println("pStar: "+printBigNum(pStar,3)+" qStar: "+printBigNum(partial_q,3)+" rho*qStar: "+printBigNum(partial_q.multiply(new BigDecimal(ro)),3));

			//Check if we are done
			if(lowestOverallBound+Ival < minELowerBound && numConfsEvaluated.compareTo(BigInteger.ONE) >= 0){ 
				//&& !subRotamers){//&& !useTopKHeuristic){  
				// We are not done and we pruned too much, repeat search
				AStarResults asr = new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				asr.status = AStarResults.INCREASEIVAL;
				return asr;
			}
			else if (minELowerBound > curThreshold && numConfsEvaluated.compareTo(BigInteger.ZERO)>0 ){
				//MH: always evaluate at least one conformation 
				AStarResults asr = new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				asr.status = AStarResults.DONE;
				return asr;
			}
			else if(useMaxKSconfs && numConfsEvaluated.compareTo(maxKSconfs) >= 0){
				AStarResults asr = new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				asr.status = AStarResults.DONE;
				return asr;
			}

			else { //the energy of the new conformation is still below the threshold

				// After minimization do a m.updateCoordinates() to
				//  resync the actualCoordinates which were changed
				//  in the minimization procedure
				//				numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
				//				numConfsLeft = numConfsLeft.subtract(BigInteger.ONE);
				//				if (stericF!=null)
				//					stericF.setNumConfsLeft(numConfsLeft); //update the number of conformations still to examine
				double[] myEnergy = new double[6]; //dim0: minE, dim1: unMinE, dim2-5 energy parts
				if (doMinimization){
					if(es.useEPIC){
						myEnergy[0] = minELowerBound;
						myEnergy[1] = myEnergy[0];//we don't generate an unminimized energy in this case
						//but this isn't really needed for the K* so we'll go with this
					}
					else if (!minimizeBB){ //side-chain minimization
						myEnergy = KSenergy(minimizeBB,doBackrubs,false, useCCD);	
					}
				}
				else if (computeEVEnergy){
					myEnergy[0] = (double)minELowerBound;
					myEnergy[1] = myEnergy[0];
				}

				m.updateCoordinates();
				m.revertPertParamsToCurState();

				if(parentConf != null) //If we have the parent conf we want to save it
					updateAll(saveConfsParams, parentConf, myEnergy,minELowerBound);
				else
					updateAll(saveConfsParams, conf, myEnergy,minELowerBound);

				System.out.println("energy: "+myEnergy[0]);
			}
		}
		if ((numConfsLeft.equals(BigInteger.ZERO))&&(!k_const.equals(BigInteger.ZERO))){ //no conformations remaining, non-sterics pruned			
			if (partial_q.multiply(new BigDecimal(ro),ExpFunction.mc).compareTo(pStar)<0){ //approximation accuracy not achieved, so repeat search
				BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
				BigDecimal e = ef.exp(-Ec_const/constRT);
				if (e.compareTo(BigDecimal.ZERO)!=0){
					BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP				
					BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
					//					setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE

					repeatSearch = true;
				}
			}
		}
		return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
	}


	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	public double computeBestRotEnergyBound(int[] AAnums, int[] ROTnums, Index3wVal[] boundPerPos) {

		Index3wVal pos1Index=null;
		
		double bestE = arpMatrix.getTemplMinE(); // Add shell-shell energy

		for (int i=0; i< AAnums.length; i++){ //compute using the formula

			double tmpSingleE = arpMatrix.getSingleMinE(i,AAnums[i],ROTnums[i]); //Add the intra-rotamer energy 
//			System.out.println("s_"+i+"_"+AAnums[i]+"_"+ROTnums[i]+" "+tmpSingleE);
			bestE += tmpSingleE;

			if(improvedBounds)
				pos1Index = updateBoundPerPos(i,tmpSingleE,boundPerPos);
			
			for(int j=i+1;j<AAnums.length;j++){
				if(!arpMatrix.areNeighbors(i, j))
					continue;
				double tmpPairE = arpMatrix.getPairwiseE(i,AAnums[i],ROTnums[i],j,AAnums[j],ROTnums[j]);
				bestE += tmpPairE;
//				System.out.println("p_"+i+"_"+AAnums[i]+"_"+ROTnums[i]+"_"+j+"_"+AAnums[j]+"_"+ROTnums[j]+" "+tmpPairE);
				
				if(improvedBounds)
					updateBoundPerPosPair(pos1Index,j,tmpPairE,boundPerPos);
			}


		}

		return bestE;

	}

	private Index3wVal updateBoundPerPos(int pos, double E, Index3wVal[] boundPerPos){
		Index3wVal iVal = null;
		if(superRotamers){
			iVal = boundPerPos[pos];
			iVal.val += E;
		}
		else if(tuples || partitionedRotamers){
			for(Index3wVal i3wV:boundPerPos)
				if(i3wV.containsPos(pos)){
					iVal = i3wV;
					iVal.val += E;
				}
		}
		return iVal;
	}

	private void updateBoundPerPosPair(Index3wVal pos1Index, int pos2, double E, Index3wVal[] boundPerPos){
		
			if(superRotamers){
				boundPerPos[pos2].val += E/2;
				pos1Index.val += E/2;
			}
			else if(tuples || partitionedRotamers){
				pos1Index.val += E/2;
				for(Index3wVal i3wV:boundPerPos)
					if(i3wV.containsPos(pos2))
						i3wV.val += E/2;
			}
		
	}

	public double[] KSenergy(boolean minimizeBB, boolean doBackrubs, boolean storeEnergyPerRes, boolean useCCD){

		double[] energy = new double[6];
		//energy[0] is minimized energy, energy [1] unminimized energy
		energy[0] = 0.0f;
		energy[1] = 0.0f;
		//energy[2-5] are the energy parts of the minimized energy
		energy[2] = 0.0f; energy[3] = 0.0f; energy[4] = 0.0f; energy[5] = 0.0f;

		if (doMinimization){

			if (!minimizeBB){ //side-chain minimization
				//unMin energy
				energy[1] = calcTotalSnapshotEnergy();			
				if(useCCD){//The ideal dihedrals are already recorded since we just reinitialized ccdMin
					((ContSCObjFunction)ccdMin.objFcn).updateIdealDihedrals(); //The molecule has ideal values for the dihedrals now; record these
					ccdMin.minimize();
					energy[0] = (double)this.efunc.getEnergy();
				}else{
					simpMin.minimize(numMinSteps);
					//minimized energy
					double tmpE[];
					//					if(!storeEnergyPerRes)
					tmpE = calcTotalSnapshotEnergyParts();
					//					else
					//						tmpE = calcTotalSnapshotEnergyParts(arpMatrix.allMutRes()); 
					energy[0] =  tmpE[0];
					energy[2] =  tmpE[1];energy[3] =  tmpE[2];energy[4] =  tmpE[3];energy[5] =  tmpE[4];

					if (doDihedE) //add dihedral energies
						energy[0] += simpMin.computeDihedEnergy();
				}


			}
			else { //backbone minimization
				energy[1] = calcTotalSnapshotEnergy();

				if (!doBackrubs)
					bbMin.minimizeFull(false);
				else
					brMin.minimizeFull();
				double tmpE[];
				tmpE = calcTotalSnapshotEnergyParts();
				energy[0] =  tmpE[0];
				energy[2] =  tmpE[1];energy[3] =  tmpE[2];energy[4] =  tmpE[3];energy[5] =  tmpE[4];
			}

			m.updateCoordinates();
		}
		else if (computeEVEnergy){
			//unMin energy
			energy[0] = calcTotalSnapshotEnergy();
			energy[1] = energy[0];

			m.restoreAtomCoord();
			m.updateCoordinates();

		}

		if(energy[1] < energy[0]){
			energy[0] = energy[1];
		}


		System.out.println("energy: "+energy[0]+" ");
		return energy;

	}



	/**
	 * @param lastRun
	 * @param saveTopConfs
	 * @param printTopConfs
	 * @param numTopConfs
	 * @param conf
	 * @param myEnergy
	 * 
	 */
	public void updateAll(SaveConfsParams scp, EMatrixEntryWIndex[] conf, double[] myEnergy,
			double minELowerBound ) {
		updateNumConfs();
		updatePartialQ(myEnergy[0]);
		updateEnergy(myEnergy);

		if ((scp.saveTopConfs || scp.printTopConfs)){
			ConfPair cp = new ConfPair(conf, myEnergy);
			ConfPair head = topConfs.peek();
			if(topConfs.size() >= scp.numTopConfs){
				if(cp.energy[0] < head.energy[0]){
					topConfs.add(cp);
					topConfs.remove(); //Will remove head of queue, so want to traverse
					//Backwards when generating conformations
				}
			}
			else
				topConfs.add(cp);
		}
	}

	//updates
	public synchronized void updateEnergy(double[] e){
		bestEMin = Math.min(bestEMin, e[0]);
		bestEUnMin = Math.min(bestEUnMin, e[1]);
	}

	public synchronized void updateNumConfs(){
		numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
		numConfsLeft = numConfsLeft.subtract(BigInteger.ONE);

	}

	public synchronized void updatePartialQ(double energy){
		ExpFunction ef = new ExpFunction();
		partial_q = partial_q.add(ef.exp(-((energy)) / constRT),ExpFunction.mc);
	}


	//	public void applyRotamers(MutableResParams strandMut, int[] conf) {
	//		//Apply the rotamers of the current conformation
	//		int curAS = 0;	
	//
	//		for(int i=0; i<strandMut.allMut.length;i++){
	//			int str=strandMut.resStrand[i];
	//			int strResNum = strandMut.resStrandNum[i];
	//			if(doPerturbations){
	//				boolean validRC = ((StrandRCs)strandRot[str]).applyRC(m, strResNum, conf[curAS]);
	//				if(!validRC)
	//					System.err.println("Error: invalid RC " + conf[curAS] + " at residue " + 
	//							strResNum + " of strand " + str );
	//
	//				curStrRotNum[curAS] = conf[curAS];
	//			}
	//			else{
	//
	//				int molResNum = strandMut.allMut[i];
	//
	//				if (strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum])!=0){//not GLY or ALA
	//					strandRot[str].applyRotamer(m, strResNum, conf[curAS]);
	//					curStrRotNum[curAS] = conf[curAS];
	//				}
	//				else { //GLY or ALA
	//					curStrRotNum[i] = 0;
	//				}
	//			}
	//
	//		}
	//
	//	}



	//// END ROTAMER SEARCH SECTION
	/////////////////////////////////////////////////////////////////////////////


	/*
	 * 
	 * DEE section
	 * 
	 */		
	//Compute the two interval terms in the summation of the MinDEE criteria
	//	public void doCompMinDEEIntervals(int numMutable, int strandMut[][], PrunedRotamers<Boolean> prunedRotAtRes,
	//			boolean scaleInt, double maxIntScale){
	//
	//		System.out.print("Computing MinDEE interval terms..");	
	//
	//		double dist[][] = null;
	//		if (scaleInt)
	//			dist = doCompDistSC(numMutable,strandMut);
	//
	//		MinDEEIntervals compInt = new MinDEEIntervals(arpMatrix, arpMatrixMax, numMutable, 	strandMut, 
	//				strandRot, prunedRotAtRes, scaleInt, dist, maxIntScale, mutRes2Strand,mutRes2StrandMutIndex, doPerturbations);
	//
	//		compInt.compMinDEEIntervals();
	//		indIntMinDEE = compInt.getIndInt();
	//		pairIntMinDEE = compInt.getPairInt();
	//
	//		System.out.println("done.");
	//
	//		if (debug){
	//			System.out.println();
	//			System.out.print(" ind: ");
	//			for (int curPos=0; curPos<indIntMinDEE.length; curPos++){
	//				System.out.print(indIntMinDEE[curPos]+" ");
	//			}
	//			System.out.println();
	//			System.out.print(" pair: ");
	//			for (int curPos=0; curPos<pairIntMinDEE.length; curPos++){
	//				System.out.print(pairIntMinDEE[curPos]+" ");
	//			}
	//			System.out.println();
	//		}
	//	}

	//Compute the distance between the side-chains for each active site residue pair (the ligand is not considered here);
	//	Returns the minimum distance between a pair of atoms in the two side-chains, for each side-chain pair
	private double[][] doCompDistSC(int numMutable, int strandMut[][]){

		double dist[][] = new double[numMutable][numMutable];

		for (int str1=0;str1<numberOfStrands;str1++){
			for(int i=0;i<strandMut[str1].length;i++){
				Residue r1 = m.strand[str1].residue[strandMut[str1][i]];
				for(int str2=str1;str2<numberOfStrands;str2++){
					for (int j=i+1; j<strandMut[str2].length; j++){
						Residue r2 = m.strand[str2].residue[strandMut[str2][j]];
						dist[i][j] = r1.getDist(r2,false);
						dist[j][i] = dist[i][j];
					}
				}
			}
		}

		/*for (int i=0; i<numResInActiveSite; i++){
			Residue r1 = m.strand[sysStrNum].residue[residueMap[i]];
			for (int j=i+1; j<numResInActiveSite; j++){
				Residue r2 = m.strand[sysStrNum].residue[residueMap[j]];
				dist[i][j] = r1.getDist(r2,false);
				dist[j][i] = dist[i][j];
			}
		}*/
		return dist;
	}

	//Prune all rotamers that are incompatible with the template (intra E + res-to-template E >= stericE)
	//Return number of rotamers pruned at each position
	public static int[] DoPruneStericTemplate(Emat arpMatrix, double stericE, boolean doPairs, PrintStream outPS){

		int[] numPrunedPerPos = new int[arpMatrix.numMutPos()];
		int numPruned = 0;
		int numPairPruned = 0;
		//Compute for the AS residues first
		Iterator<EMatrixEntryWIndex> iter = arpMatrix.singlesIterator(); 
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI= iter.next();
			EMatrixEntry re = emeWI.eme;
			if (re.minE()>=stericE && !re.isPruned()){
				arpMatrix.setSinglePruned(emeWI.index, true);
				numPruned++;
				numPrunedPerPos[emeWI.pos1()]++;
				outPS.println(emeWI.pos1()+" "+emeWI.aa1()+" "+emeWI.rot1());
			}
		}	

		if(doPairs){
			iter = arpMatrix.pairsIterator();
			while(iter.hasNext()){
				EMatrixEntryWIndex emeWI= iter.next();
				EMatrixEntry re = emeWI.eme;
				if (re.minE()>=stericE && !re.isPruned()){

					//int index_r = curPos*totalNumRotamers + rotamerIndexOffset[curAA] + curRot;
					arpMatrix.setPruned(emeWI.index, true);
					//re.setPrunedIsSteric(true);
					numPairPruned++;
					//System.out.println(emeWI.pos1()+" "+emeWI.aa1()+" "+emeWI.rot1()+" "+emeWI.pos2()+" "+emeWI.aa2()+" "+emeWI.rot2());
				}
			}
		}
		outPS.println("Number of rotamers pruned due to incompatibility with the template: "+numPruned);		
		outPS.println("Number of pairs pruned due to incompatibility with the template: "+numPairPruned);

		return numPrunedPerPos;
	}



	//Prune rotamer or residue-conformation pairs that are definitely impossible due to a steric clash or to parametric incompatibility
	//as measured by their interaction energy being greater than the cutoff
	//Particularly useful for DEEPer
	public void pruneRidiculousPairs(int numMutable, MutableResParams strandMut,
			Emat emat, double cutoff){

		//		int numAAtypes[] = new int[numMutable];


		//		for(int ctr=0; ctr<strandMut.allMut.length;ctr++){
		//			//the number of AAs allowed for each AS residue
		//			numAAtypes[ctr] = strandRot[strandMut.resStrand[ctr]].getNumAllowable(strandMut.resStrandNum[ctr]);
		//		}

		int numPruned = 0;

		PairsIterator iter = emat.pairsIterator();
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(!arpMatrix.getSinglePruned(emeWI.rot1index())){
				if(!arpMatrix.getSinglePruned(emeWI.rot2index())){
					if( emeWI.eme.minE() >= cutoff ){

						arpMatrix.setPairPruned(emeWI.index, true);
						arpMatrix.setSymmetricPairPruned(emeWI.index, true);

						numPruned++;
					}
				}
			}
		}

		//		for (int curPos=0; curPos<numMutable; curPos++){
		//			int str=strandMut.resStrand[curPos];
		//			int strResNum=strandMut.resStrandNum[curPos];
		//			for (int AA=0; AA<numAAtypes[curPos]; AA++){
		//				int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
		//
		//				//find how many rotamers are allowed for the current AA type at the given residue;
		//				//note that ala and gly have 0 possible rotamers
		//				int numRotForCurAAatPos = getNumRot( str, strResNum, curAA );
		//
		//
		//				for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){
		//
		//					if( !arpMatrix.getSinglePruned(curPos, curAA, curRot)){
		//
		//						for (int altPos=0; altPos<curPos; altPos++){
		//							int str2=strandMut.resStrand[altPos];
		//							int strResNum2=strandMut.resStrandNum[altPos];
		//							for (int AA2=0; AA2<numAAtypes[altPos]; AA2++){
		//								int altAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
		//
		//								//find how many rotamers are allowed for the current AA type at the given residue;
		//								//note that ala and gly have 0 possible rotamers
		//								int numRotForCurAAatPos2 = getNumRot( str2, strResNum2, altAA );
		//
		//
		//								for(int altRot=0; altRot<numRotForCurAAatPos2; altRot++){
		//
		//
		//									if( !arpMatrix.getSinglePruned(altPos, altAA, altRot)){
		//
		//										if( arpMatrix.getPairwiseE(curPos, curAA, curRot, altPos, altAA, altRot) >= cutoff ){
		//
		//											arpMatrix.setPairPruned(curPos,curAA,curRot,altPos,altAA,altRot, true);
		//											arpMatrix.setPairPruned(altPos,altAA,altRot,curPos,curAA,curRot, true);
		//
		//											numPruned++;
		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}


		if(doPerturbations)
			System.out.println("Number of pairs pruned due to steric clashes or parametric incompatibility: "+numPruned);
		else
			System.out.println("Number of pairs pruned due to steric clashes: "+numPruned);
	}


	//Remove RCs that are parametrically incompatible with all RCs at some other position (i.e. are impossible)
	//This can be done before energy precomputation to avoid precomputing useless energies
	//The RCs are actually removed from sysLR and ligROT
	public void removeImpossibleRCs(int numMutable, MutableResParams strandMut){

		boolean done = false;

		while( !done ){//We iterate until no more RCs can be removed

			done = true;

			int numAAtypes[] = new int[numMutable];

			boolean prunedRCs[][] = new boolean[numberOfStrands][];//Which RCs should be removed
			for(int str=0; str<numberOfStrands; str++)
				prunedRCs[str] = new boolean[((StrandRCs)strandRot[str]).rcl.allRCs.size()];

			for (int curPos=0; curPos<numMutable; curPos++){
				Residue res = m.residue[strandMut.allMut[curPos]];
				int str = res.strandNumber;
				int strResNum = res.strandResidueNumber;
				for(ResidueConformation rc1: res.allowedRCs()){

					for (int altPos=0; altPos<numMutable; altPos++){

						if( altPos != curPos ){

							int str2=strandMut.resStrand[altPos];
							int strResNum2=strandMut.resStrandNum[altPos];

							boolean prune = true;
							Residue res2 = m.strand[str2].residue[strResNum2];

							for (ResidueConformation altRC: res2.allowedRCs()){
								if( ! isParametricallyIncompatible(res, rc1, res2, altRC))
									prune = false;
							}

							prunedRCs[str][rc1.id] = prunedRCs[str][rc1.id] || prune;
						}
					}

					if(prunedRCs[str][rc1.id])
						done = false;
				}
			}

			for(int str=0; str<numberOfStrands; str++)
				((StrandRCs)strandRot[str]).removeRCs(prunedRCs[str],m);
		}
	}


	//Do simple Goldstein DEE
	public static void DoDEEGoldstein(Emat emat,double initEw, boolean useSF, 
			boolean typeDep, boolean useMinDEEPruningEw, double Ival,
			boolean distrDEE, boolean[] resInMut,int[] singleStartEnd, boolean removeRot){

		//			System.out.println("Starting pruning with DEE (simple Goldstein)");
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldsteinPierceNoIter DEERun = new DEEGoldsteinPierceNoIter(emat, initEw, 
				useSF, typeDep, useMinDEEPruningEw, Ival,distrDEE,resInMut,singleStartEnd,removeRot);

		DEERun.ComputeEliminatedRotConf();

		DEERun = null;

	}

	//Do Goldstein DEE pairs
	public static void DoDEEPairs(Emat emat, double initEw, boolean resInPair[], boolean doMinimize, 
			boolean useSF, boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean scaleInt, 
			double maxScale, boolean typeDep, boolean doIMinDEE, double Ival, int[] singleStartEnd){

		//			if(magicBullet)
		//				System.out.println("Starting pruning with DEE (mb pairs)");
		//			else
		//				System.out.println("Starting pruning with DEE (full pairs)");
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldsteinPairsNoIter DEERunPairs = new DEEGoldsteinPairsNoIter(emat, initEw, 
				resInPair, doMinimize, useSF, magicBullet, distrDEE, minimizeBB, scaleInt, 
				maxScale, typeDep, doIMinDEE, Ival, singleStartEnd);

		DEERunPairs.ComputeEliminatedRotConf();		

		DEERunPairs = null;
		//			System.out.println();
	}

	//SplitDEE (conformational splitting) with 1 or 2 split positions
	public static void DoDEEConfSplitting(Emat emat,double initEw, boolean resInMut[],
			boolean doMinimize, boolean useSF, int numSplits, 
			boolean distrDEE, boolean minimizeBB, boolean typeDep, boolean doIMinDEE, double Ival,
			int[] singleStartEnd){

		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		if (numSplits==1){ //1 split position
			//				System.out.println("Starting pruning with DEE (1-sp split-DEE)");
			DEESplit1fPierceNoIter DEERunConfSplitting = new DEESplit1fPierceNoIter(emat, initEw, 
					resInMut, doMinimize,  
					useSF, distrDEE, minimizeBB, typeDep,
					doIMinDEE, Ival,singleStartEnd);

			DEERunConfSplitting.ComputeEliminatedRotConf();
			//splitFlags = DEERunConfSplitting.getSplitFlags();

			DEERunConfSplitting = null;
			//				System.out.println();
		}
		else{ // 2 split positions
			//				System.out.println("Starting pruning with DEE (2-sp split-DEE)");
			DEESplit2fPierceNoIter DEERunConfSplitting = new DEESplit2fPierceNoIter(emat, initEw, 
					resInMut, doMinimize,  
					useSF, distrDEE, minimizeBB, typeDep, doIMinDEE, Ival,singleStartEnd);

			DEERunConfSplitting.ComputeEliminatedRotConf();
			//splitFlags = DEERunConfSplitting.getSplitFlags();

			DEERunConfSplitting = null;
			//				System.out.println();
		}

		//return eliminatedRotAtRes;
	}

	public void DoMinBounds(double pruningE, double initEw, boolean useSF, 
			boolean boundKS, boolean onlyBounds){

		if (boundKS && onlyBounds){ //cannot be both set at the same time
			System.out.println("ERROR (RotamerSearch, DoMinBounds): boundKS and onlyBounds cannot be simultaneously 'true'.");
			System.exit(1);
		}

		MSMinBounds minBoundsRun = new MSMinBounds(arpMatrix,pruningE,useSF,initEw,boundKS,onlyBounds);

		//		if ( (!boundKS) && (!onlyBounds) ) { //use Bounds to prune new rotamers
		//			eliminatedRotAtRes = minBoundsRun.ComputeEliminatedRotConf();
		//		}
		if (boundKS) { //compute Ec
			minBoundsRun.ComputeEliminatedRotConf();
			//prunedIsSteric = minBoundsRun.getPrunedSteric();
			Ec_const = minBoundsRun.getEc();
		}
		//		else { //(onlyBounds==true)
		//			minBoundsRun.ComputeEliminatedRotConf();
		//			boundForPartition = minBoundsRun.getBoundForPartition();
		//		}

		minBoundsRun = null;

		//return eliminatedRotAtRes;
	}

	//	TODO: Implement DoDEEIndirect Pruning
	public static void DoDEEIndirect(Emat emat, MutableResParams strandMut,
			double initEw, boolean resInPair[], boolean doMinimize,
			boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean typeDep,
			boolean doIMinDEE, double Ival, boolean doPerturbations, Molecule m, StrandRotamers[] strandRot){


		DEEIndirect DEERun;

		int numPos = strandMut.numMutPos();

		ArrayList<boolean[]> mpz = new ArrayList<boolean[]>();
		if(doPerturbations)
			mpz = getMinimalPruningZones(m,strandMut);

		//		eliminatedRotAtRes = prunedRotAtRes;

		for(int zoneNum=-1; zoneNum<mpz.size(); zoneNum++){

			boolean inZ[];//Indicates if a residue is in Z

			if(zoneNum == -1){//All active-site residues and the ligand are in the pruning zone
				inZ = new boolean[numPos];
				for(int pos=0; pos<numPos; pos++)
					inZ[pos] = true;
			}
			else
				inZ = mpz.get(zoneNum);

			//TODO: Implement triples pruning
			//				if(useTriples)
			//					DEERun = new DEEIndirect(arpMatrix, numMutable, strandMut, initEw,
			//							strandRot, resInPair, doMinimize,
			//							splitFlags, magicBullet, distrDEE, minimizeBB, 
			//							typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations, inZ);
			//	
			//				else
			DEERun = new DEEIndirect(emat, strandMut, initEw,
					strandRot, resInPair, doMinimize,
					magicBullet, distrDEE, minimizeBB, 
					typeDep, doIMinDEE, Ival, null, doPerturbations, inZ);


			DEERun.ComputeEliminatedRotConf();
			//				splitFlags = DEERun.getSplitFlags();

			DEERun = null;

		}

		//		return eliminatedRotAtRes;
	}


	public void DoDEETriples(int numMutable, MutableResParams strandMut,
			double initEw, boolean resInTriple[], boolean doMinimize,
			boolean magicBullet, int magicBulletNum, boolean distrDEE, boolean minimizeBB, boolean typeDep,
			boolean doIMinDEE, double Ival){

		DEEGoldsteinTriples DEERun = new DEEGoldsteinTriples(arpMatrix, numMutable, strandMut, initEw, strandRot,  
				resInTriple, doMinimize, arpMatrix.pairs.pruned, magicBullet, magicBulletNum, distrDEE, minimizeBB,
				typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations );

		if(tripleFlags == null)//At this point tripleFlags == DEERun.tripleFlags
			DEERun.initializeTripleFlags();

		DEERun.ComputeEliminatedRotConf();
		tripleFlags = DEERun.getTripleFlags();

		DEERun = null;
	}


	//SplitDEE (conformational splitting) with 1 or 2 split positions
	//	public PrunedRotamers<Boolean> DoDEEConfSplitting(int numMutable, 
	//			int strandMut[][], double initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInMut[],
	//			boolean doMinimize, boolean useSF, int numSplits, 
	//			boolean distrDEE, boolean minimizeBB, boolean typeDep, boolean doIMinDEE, double Ival){
	//
	//		//arpmatrix has a row/column for the backbone energies, so we just need
	//		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
	//		if (numSplits==1){ //1 split position
	//			DEESplit1f DEERunConfSplitting = new DEESplit1f(arpMatrix, arpMatrixMax, 
	//					numMutable, strandMut, initEw, 
	//					strandRot, prunedRotAtRes, resInMut, doMinimize, indIntMinDEE, pairIntMinDEE, 
	//					splitFlags, useSF, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex, typeDep,
	//					doIMinDEE, Ival, doPerturbations);
	//
	//			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
	//			splitFlags = DEERunConfSplitting.getSplitFlags();
	//
	//			DEERunConfSplitting = null;
	//		}
	//		else{ // 2 split positions
	//			DEESplit2f DEERunConfSplitting = new DEESplit2f(arpMatrix, arpMatrixMax, 
	//					numMutable, strandMut, initEw, 
	//					strandRot, prunedRotAtRes, resInMut, doMinimize, indIntMinDEE, pairIntMinDEE, 
	//					splitFlags, useSF, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex, typeDep,
	//					doIMinDEE, Ival, doPerturbations);
	//
	//			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
	//			splitFlags = DEERunConfSplitting.getSplitFlags();
	//
	//			DEERunConfSplitting = null;
	//		}
	//
	//		return eliminatedRotAtRes;
	//	}
	//	//Do Goldstein DEE pairs
	//	public void DoDEEPairs(int numMutable, int strandMut[][], 
	//			double initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInPair[], boolean doMinimize, 
	//			boolean useSF, boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean scaleInt, 
	//			double maxScale, boolean typeDep, boolean doIMinDEE, double Ival){
	//
	//		//arpmatrix has a row/column for the backbone energies, so we just need
	//		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
	//		DEEGoldsteinPairs DEERunPairs = new DEEGoldsteinPairs(arpMatrix, arpMatrixMax, numMutable, strandMut, initEw, strandRot, prunedRotAtRes, resInPair, doMinimize, 
	//				splitFlags, useSF, magicBullet, distrDEE, minimizeBB, scaleInt, maxScale, 
	//				mutRes2Strand, mutRes2StrandMutIndex, typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations );
	//
	//		DEERunPairs.ComputeEliminatedRotConf();		
	//		splitFlags = DEERunPairs.getSplitFlags();
	//
	//		DEERunPairs = null;
	//	}


	public static ArrayList<boolean[]> getMinimalPruningZones(Molecule m, MutableResParams strandMut){//Get the minimal pruning zones for indirect pruning
		//This is for DEEPer
		//answer.get(a)[b] indicates if flexible residue b is in pruning zone a


		ArrayList<boolean[]> mpz = new ArrayList<boolean[]>();

		if(m.perts.length==0)
			return mpz;//Return an empty ArrayList because there are no minimal pruning zones

		int numMutable = strandMut.numMutPos();
		int start = 0;
		if(m.perts[0].resDirectlyAffected.length == numMutable)//If there is an initial perturbation affecting everything (e.g. a full structure switch) don't count it
			start = 1;

		for(int pertNum=start;pertNum<m.perts.length;pertNum++){

			Perturbation pert = m.perts[pertNum];

			int flexResAffected[] = new int[pert.resDirectlyAffected.length];//Residues affected by the perturbation (numbering with respect to all flexible residues)
			for(int a=0;a<pert.resDirectlyAffected.length;a++){
				int molResNum = pert.resDirectlyAffected[a];//The molecule residue number
				int strResNum = m.residue[molResNum].strandResidueNumber;
				int str = m.residue[molResNum].strandNumber;
				for(int b=0;b<numMutable;b++){
					int strb = strandMut.resStrand[b];

					//If the b'th flexible residue has the indicated molResNum
					if( (strb==str) && (strandMut.resStrandNum[b]==strResNum) )
						flexResAffected[a] = b;
				}
			}

			ArrayList<Integer> zones = new ArrayList<Integer>();//Existing pruning zones that residues affected by this perturbation are found in
			for( int zoneNum=0; zoneNum<mpz.size(); zoneNum++){
				for(int a=0;a<flexResAffected.length;a++){
					if( mpz.get(zoneNum)[flexResAffected[a]] && (!zones.contains(zoneNum)) )
						zones.add(zoneNum);
				}
			}

			if( zones.isEmpty() ){//Need to create a new pruning zone
				boolean[] new_mpz = new boolean[numMutable];
				for(int a=0;a<flexResAffected.length;a++)
					new_mpz[flexResAffected[a]] = true;
				mpz.add(new_mpz);
			}
			else if ( zones.size() == 1 ){//Put this perturbation's affected residues in an existing pruning zone
				for(int a=0;a<flexResAffected.length;a++)
					mpz.get(zones.get(0))[flexResAffected[a]] = true;
			}
			else{//zones.size() > 1: Merge the pruning zones and put this perturbation's affected residues in it
				boolean[] new_mpz = new boolean[numMutable];
				for(int a=0;a<zones.size();a++){
					for(int b=0;b<numMutable;b++){
						if( mpz.get(zones.get(a))[b] )
							new_mpz[b] = true;
					}
				}
				for(int a=0;a<flexResAffected.length;a++)
					new_mpz[flexResAffected[a]] = true;

				mpz.add(new_mpz);
				for(int a=zones.size()-1;a>=0;a--)//zones is in ascending order...delete in descending order to avoid messing up the numbering
					mpz.remove(zones.get(a).intValue());
			}
		}

		return mpz;
	}


	/*
	 * 
	 * End of DEE section
	 * 
	 */	

	//////////////////////////////////////////////////////////////////////////
	//	Compute Min GMEC section
	//////////////////////////////////////////////////////////////////////////
	public AStarResults doAStarGMEC(String fileName, boolean searchComputeEVEnergy, 
			boolean searchDoMinimization,int numMutable, 
			MutableResParams strandMut, double Ew, double bestScore, 
			CommucObj cObj, boolean approxMinGMEC, double lambda, boolean minimizeBB, boolean useEref, 
			boolean doBackrubs, String backrubFile, boolean useMinDEEPruningEw, double Ival,Settings.Enum enumSettings,Settings.Output outSettings) {

		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer
		// If we're doing a mutation search then residues
		//  are allowed to mutate


		if(doPerturbations)//Make sure minimizer is set properly
			((PMinimizer)simpMin).minimizePerturbations = minimizePerturbations;//which should generally be true


		numConfsEvaluated = BigInteger.ZERO;
		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
		//ASAANums = new int[numberMutable];
		//curStrRotNum = new int[numberMutable];
		//curASRotNum = new int[numInAS];
		//int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
		// This map maps the system strand residues back to the AS numbering
		// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		//curLigRotNum = 0;

		//boolean ligPresent = (ligStrNum>=0); //determine if a ligand is present

		for(int i=0; i<strandMut.allMut.length;i++){
			m.residue[strandMut.allMut[i]].flexible = true;
		}

		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return null;
		}

		// Prepare Amber
		if(searchComputeEVEnergy){
			
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null&&ccdMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return null;
					}
					bbMin = null;
					brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){
						if (bbMin == null) {
							System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
							return null;
						}
						simpMin = null;
						brMin = null;
					}
					else {
						if (brMin == null) {
							System.out.println("Warning: Attempting minimization run but brMin not allocated, RotamerSearch aborting");
							return null;
						}
						simpMin = null;
						bbMin = null;
					}
				}
			}
		}

		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null) {
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return null;
		}

		long startEnum = System.currentTimeMillis();
		AStarResults asr = doAStarGMECHelper(numMutable, strandMut, fileName, Ew, bestScore, cObj, 
				approxMinGMEC, lambda, minimizeBB, useEref, doBackrubs, backrubFile, useMinDEEPruningEw, Ival, enumSettings,outSettings);
		long endEnum = System.currentTimeMillis();
		KSParser.metrics.EnumTime += (endEnum - startEnum);
		KSParser.metrics.totalNumConfs += asr.numConfsEvaluated;
		
		if(MSAStarSearch != null){
			MSAStarSearch.stopSlaves();
			KSParser.metrics.numExpanded += MSAStarSearch.numExpanded;
			KSParser.metrics.totNumNodes += MSAStarSearch.curExpansion.totalNodes;
			MSAStarSearch = null;
		}
		return asr;
	}
	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by mutationRotamerSearch(.)
	private AStarResults doAStarGMECHelper(int numMutable, MutableResParams strandMut, String fileName, 
			double Ew, double bestScore, CommucObj cObj, 
			boolean approxMinGMEC, double lambda, boolean minimizeBB, boolean useEref,  
			boolean doBackrubs, String backrubFile, boolean useMinDEEPruningEw, double Ival, Settings.Enum enumSettings, Settings.Output outSettings){

		boolean outputFile = (fileName!=null); //output to file

		int treeLevels; // total num levels in the conformation tree
		treeLevels = numMutable;

		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues

		int numTotalConf = 1;
		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){ 
			curNumRot = 0;
			for(int i=0; i<arpMatrix.singles.E[curLevel].length;i++){
				curNumRot += arpMatrix.singles.E[curLevel][i].length; 
			}
			numRotForRes[curLevel] = curNumRot;
			numRotForResNonPruned[curLevel] = numRotForRes[curLevel];
			numTotalRotRed += numRotForRes[curLevel];

			numTotalConf *= numRotForRes[curLevel];
		}

		int numPrunedThisLevel;

		for (int curLevel=0; curLevel<treeLevels; curLevel++){ //for each residue
			//int str = mutRes2Strand[curLevel];
			//int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			//int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;

			int curIndex;
			System.out.println("curLevel "+curLevel);
			numPrunedThisLevel = 0;
			int numWTrots = 0;
			Iterator<EMatrixEntryWIndex> iter = arpMatrix.singlesIterator(curLevel);
			while(iter.hasNext()){
				EMatrixEntry re = iter.next().eme;
				if (re.isPruned())
					numPrunedThisLevel++;
			}
			numRotForResNonPruned[curLevel] -= numPrunedThisLevel;
			//System.out.println("NumWTRots: "+numWTrots);
		}

		int numConfNonPrunedMinDEE = 1;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				System.out.println("MinDEE has pruned all possible conformations: position "+curLevel); //at least the minGMEC should be remaining
				if (cObj!=null) {//output num conf info to cObj
					//TODO: confirm program never gets here
					System.out.println("PROGRAM SHOULD NEVER REACH HERE");
					System.exit(0);
				}

				return null;
			}
			else {
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];

				numConfNonPrunedMinDEE *= numRotForResNonPruned[curLevel];
			}
		}

		if (cObj!=null) {//output num conf info to cObj
			System.out.println("PROGRAM SHOULD NEVER REACH HERE");
			System.exit(0);
			/*cObj.EL_searchNumConfsTotal = numTotalConf;
			cObj.EL_searchNumPrunedMinDEE = (numTotalConf - numConfNonPrunedMinDEE);*/
		}

		System.out.println("ASTAR PRUNING INFO: Total Rots before pruning for each residue: ");
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForRes[i]+" ");System.out.println();
		System.out.println("Total number of rotamers before pruning: "+numTotalRotRed+" ");
		System.out.println("ASTAR PRUNING INFO: Total Rots non-pruned for each residue ");
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForResNonPruned[i]+" ");System.out.println();
		System.out.println("Total number of rotamers after pruning: "+numTotalRotRedNonPruned+" ");

		//Allocate this so reduceMatrix can fill it in
		boolean[] eliminatedRotAtPosRed = new boolean[numTotalRotRed];//reduced MinDEE matrix

		boolean splitFlagsRed[][] = null, tripleFlagsRed[][][] = null;
		if( useFlagsAStar ){
			splitFlagsRed = new boolean[numTotalRotRedNonPruned][numTotalRotRedNonPruned];
			if( useTriples )
				tripleFlagsRed = new boolean[numTotalRotRedNonPruned][][];//As with the normal tripleFlags, tripleFlagsRed[a][b][c] is only defined for a>b>c
		}



		//Get the reduced min energy matrix
		Emat arpMatrixRed = arpMatrix.unprunedMatrix();


		//Declaring the logPS output here prevents opening an empty file and returning
		//	(for example, if all conformations are pruned by MinDEE above)
		PrintStream logPS = null; //the output file for conf info
		if (outputFile){
			try {			
				FileOutputStream fileOutputStream = new FileOutputStream(fileName,true); //append file if more than 1 partition
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				logPS = new PrintStream( bufferedOutputStream );
			}
			catch (Exception ex) {
				System.out.println("ERROR: An exception occured while opening log file");
			}
		}


		if(useCCD || es.useEPIC )
			m.DOFs = DegreeOfFreedom.makeDOFArray(strandRot, strandMut, m);

		//Set-up the A* search
		//		MSAStarSearch = new PGAStar(treeLevels, numRotForResNonPruned, arpMatrix,false,es,doPerturbations,m, strandRot, strandMut, cetm); //false = no reordering
		//		MSAStar MSAStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,null,
		//				splitFlagsRed,tripleFlagsRed,doPerturbations,es,
		//				m,strandRot,strandMut);
		boolean run1 = false;
		if(MSAStarSearch == null){
			switch(enumSettings.asMethod){
			case ORIG:
				MSAStarSearch = new PGAStar(treeLevels, numRotForResNonPruned, arpMatrix,false,es,doPerturbations,m, strandRot, strandMut, cetm); //false = no reordering
				//				MSAStarSearch = new PGgurobiAStar(treeLevels, numRotForResNonPruned, arpMatrix,asMethod);
				break;
			case ASORIG:
			case ASLP:
			case ASWCSP:
			case ASMPLP:
				MSAStarSearch = new PGgurobiAStar(treeLevels,numRotForResNonPruned,arpMatrix,enumSettings.asMethod,enumSettings.varOrder,es,doPerturbations,m, strandRot, strandMut, cetm);
				break;
			//			case MIN:
			//				MSAStarSearch = new PGgurobiMinAStar(treeLevels,numRotForResNonPruned,arpMatrix,m,arpMatrix.resByPos,a96ff,simpMin,arpMatrix.eRef,doDihedE,useEref,EnvironmentVars.useEntropy);
			//				break;
		
			//			case ASMPLP:
			//				MSAStarSearch = new PGMPLPAStar(treeLevels,numRotForResNonPruned,arpMatrix);
			//				break;
			case BYSEQ:
				MSAStarSearch = new PGgurobiAStarBySeq(treeLevels,numRotForResNonPruned,arpMatrix,bestScore+Ew,enumSettings.varOrder);
				break;
//			case BYSEQREORDER:
//				MSAStarSearch = new PGgurobiAStarBySeq(treeLevels,numRotForResNonPruned,arpMatrix,bestScore+Ew,true);
//				break;
			case BYSUBROT:
				MSAStarSearch = new PGgurobiAStarBySubRot(treeLevels,numRotForResNonPruned,arpMatrix,bestScore+Ew,m);
				break;
			case WCSP:
				MSAStarSearch = new WCSPSearch(treeLevels,numRotForResNonPruned,arpMatrix,Ew);
				break;
			case ILP:
				MSAStarSearch = new ILPSearch(treeLevels,numRotForResNonPruned,arpMatrix,Ew);
				break;
			default:
				System.out.println("Don't Recognize AStar method... using WCSP");
				MSAStarSearch = new WCSPSearch(treeLevels,numRotForResNonPruned,arpMatrix,Ew);
				break;
			}
			//			MSAStarSearch.rotIndexes = indicesEMatrix;
			HashMap<ArrayList<Index3>,EnergyTuple> energyTuples = new HashMap<ArrayList<Index3>,EnergyTuple>();
			MSAStarSearch.setEnergyTuples(energyTuples);
			run1 = true;
		}

		EMatrixEntryWIndex conf[] = new EMatrixEntryWIndex[treeLevels];//the rotamer sequence
		EMatrixEntryWIndex parentConf[] = null;//the parent rotamer sequence
		int numConfsOutput = 0;//num confs output to file
		double lowestBound = stericE;
		double minELowerBound = stericE;

		if(partitionedRotamers || superRotamers || tuples){
			generatedConfs = new ArrayList<GlobalRCConf>();
			worstRots = new ArrayList<ArrayList<Index3wVal>>();
		}

		if(es.useEPIC){
			if(es.gettingLowestBound)//need to run non-EPIC A*
				MSAStarSearch.es = new EPICSettings();
			else{//set up the MSAStarSearch to handle fit series
				MSAStarSearch.DOFList = m.DOFs;
				MSAStarSearch.CETM = cetm;
				lowestBound = es.lowestBound;
			}
		}


		/*if (!doMinimization)
			approxMinGMEC = false; //reset approxMinGMEC, since it is valid only for MinDEE, and not for traditional DEE*/

		updateBestE((double)bestScore);


		ArrayList<double[]> CETRecord = new ArrayList<double[]>();
		//this records data on the polynomial fits for use in analysis
		//if we are running checkEPIC

		improvedBounds = tuples || partitionedRotamers || superRotamers;
		double minIval = Double.POSITIVE_INFINITY;
		PGQueueNode curNode = new PGQueueNode(1, new int[1], 0.0, 0, 0);
		while (true){

			//debugPS.println("curMinE: "+curMinE);debugPS.flush();

			if (cObj!=null){				
				/*if (numConfsEvaluated>=cObj.confSeq.length){
					CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[2*cObj.confSeq.length];
					System.arraycopy(cObj.confSeq,0,tempConf,0,cObj.confSeq.length);
					cObj.confSeq = tempConf;
				}
				cObj.confSeq[numConfsEvaluated] = cObj.new ConfInfo(treeLevels);*/
				System.out.println("PROGRAM SHOULD NEVER REACH HERE");
				System.exit(0);
				/*cObj.EL_searchNumConfsEvaluated = numConfsEvaluated.intValue();*/
			}

			long startAS = System.currentTimeMillis();
			curNode = MSAStarSearch.doAStar(run1); //the current rotamer sequence); //the current rotamer sequence
			long endAS = System.currentTimeMillis();
			KSParser.metrics.AStime += (endAS - startAS);
			
			//check if the conformation is valid
			if (curNode == null || curNode.actualConf == null){ // no valid conformations remaining

				if (cObj!=null)
					cObj.bestScore = new BigDecimal(getBestE()); //update the best score so far to supply to the next partition

				//					if(!keepAStree){
				MSAStarSearch.stopSlaves();
				KSParser.metrics.numExpanded += MSAStarSearch.numExpanded;
				KSParser.metrics.totNumNodes += MSAStarSearch.curExpansion.totalNodes;
				MSAStarSearch = null;
				//					}

				if (outputFile){
					logPS.flush(); //there may still be results to output
				}
				return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
			}
			
			conf = curNode.actualConf;
			
			parentConf = createParentConf(conf);
			
			
			System.out.println("confNum: "+(numConfsEvaluated.add(BigInteger.ONE)));

			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			//EMatrixEntryWIndex conf[] = new EMatrixEntryWIndex[curConf.length]; //the conformation with the actual rotamer numbers

			//			conf = getActualConf(curConf,arpMatrix,treeLevels,numRotForResNonPruned);


			//KER: Need to calculate Types with template before we mutate anything to set the
			//KER: nterm and cterm flags for the residues
			if(doMinimization) //Will only mutate residues when minimization is on
				a96ff.calculateTypesWithTemplates();
			//Extract the AA numbers for the current conformation and appply the corresponding AA types
			int[] curAANums = new int[treeLevels];
			int[] curRotNums = new int[treeLevels];
			for (int i=0; i<treeLevels; i++){

				curAANums[i] = conf[i].index[1];
				curRotNums[i] = conf[i].index[2];

				//TODO: allow all rotamer libraries the ability to mutate
				//KER: only apply mutation if we actually need to calcualate an energy
				//  or if we are saving PDBs right now
				if(doMinimization || (outSettings!= null && outSettings.savePDBs)){
					conf[i].eme.applyMutation(m, arpMatrix.resByPos, addHydrogens,connectResidues );
					if(parentConf != null){ //Will only not be null if we want to compute it
						parentConf[i].eme.applyRC(arpMatrix.resByPos, m);
						System.out.print(parentConf[i].eme.printRes(m,arpMatrix.resByPos));
					}else{
						System.out.print(conf[i].eme.printRes(m,arpMatrix.resByPos));
						conf[i].eme.applyRC(arpMatrix.resByPos, m);
					}
				}

			}

			// After minimization do a m.updateCoordinates() to
			//  resync the actualCoordinates which were changed
			//  in the minimization procedure

			//////////////////////////////////////////////////////////////////////////////////////////
			// Energy computation
			double unMinE = 0.0f;
			double minE = 0.0f;
			
			LinkedList<EnergyTuple> curTuples;
			Index3wVal[] boundPerPos = null;
			Index3wVal[] energyPerPos = null;
			if(partitionedRotamers || tuples){
				curTuples = curNode.curTuples;
				EnergiesPerPos epp = initializeBoundAndEnergyPerPos(curAANums, curRotNums, curTuples);
				boundPerPos = epp.boundPerPos;
				energyPerPos = epp.energyPerPos;
			}

			
			minELowerBound = computeBestRotEnergyBound(curAANums, curRotNums,boundPerPos);///*numTotalRotamers,rotamerIndexOffset*/);

			double minTime = 0;
			double minTimeEPIC = 0;

			if(es.gettingLowestBound){
				es.lowestBound = minELowerBound;
				return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
//				return 0;//not going to use the interval in this case, just wanted the lowestBound
			}


			//debugPS.println("minELowerBound: "+minELowerBound);debugPS.flush();


			if (run1 && !es.useEPIC) //this is the first extracted conformation, so it has the lowest energy lower bound, so store it
				//for EPIC the lowestBound has been precomputed
				lowestBound = minELowerBound;


			boolean done = false;

			double polyE = 0;//energy, calculated with polynomial fits (0 if not using EPIC)

			DoubleMatrix1D EPICMinVals = null;
			//DOF values obtained by minimizing EPIC, listed as an input vector for ContSCObjFunction

			if(es.useEPIC){
				//compute energy from polynomial fits

				ContSCObjFunction cof = null;

				if(es.useSVE)//will need cof to set DOFs for sve
					cof = new ContSCObjFunction(m,numberOfStrands,null,strandRot,false,null);

				//Setup indices
				Index3[] indices = new Index3[conf.length];
				for(int i=0; i<conf.length;i++){
					indices[i] = conf[i].rot1index3();
				}

				CETObjFunction lof = cetm.getObjFunc(indices,false,false,cof);
				CCDMinimizer lmin = new CCDMinimizer(lof,false);

				DoubleMatrix1D lminVals = lmin.minimize();
				minTimeEPIC = lmin.minTime;
				double LSBE = lof.getValue( lminVals );//comput continuous energy term contribution to energy bound

				if(es.useSVE){//we have cof, lof.sveOFMap
					EPICMinVals = DoubleFactory1D.dense.make(lof.sveOFMap.length);
					for(int dof=0; dof<lof.sveOFMap.length; dof++)
						EPICMinVals.set(dof, lminVals.get(lof.sveOFMap[dof]));
				}

				polyE = minELowerBound + LSBE;
				double noShellShell = polyE - arpMatrix.getShellShellE();
				System.out.println("CET term: "+LSBE+" Full bound: "+polyE+" Without shell-shell E:"+noShellShell);

				if(es.checkEPIC){
					System.out.print("CET components: ");
					lof.printTerms(lminVals);
					System.out.println();
				}

				minE = polyE;

				es.checkThresh1Validity(lof,lminVals);

				//we know we're enumerating conformations in order of energy
				//so if minE is more than Ew above the GMEC so far, then we've got all the conformations we need
				if( minE > (getBestE()+Ew) )
					done = true;


				if(es.useSVE){
					//cof minimized m...revert to unminimized state
					m.updateCoordinates();
					m.revertPertParamsToCurState();
				}
			}
			else if ((minELowerBound>(getBestE()+Ew)) && (!run1)){ //we already have all confs within Ew of the minGMEC
				done = true;
			}
			else if (approxMinGMEC){ //running the heuristic halting condition

				if ((minELowerBound>(lowestBound+lambda)) && (!run1)) //compare the current bound to the lowest bound
					done = true;
			}
			// 2010: This is a debatable issue: Our "new" useMinDEEPruningEw method doesn't prune anything
			//  withing Ew of the lowestBound, while the old MinDEE doesn't prune anything within Ew of
			//   the minGMEC.  In any case, if we pruned too much, returned the new energy window so 
			//   the process can be repeated.  

			/* 2013: iMinDEE pruning ensures all rotamers are present for conformations whose lower bounds are within
			 * Ival+Ew of the lowestBound.  So if we enumerate in order of lower bound and finish with all lower bounds
			 * below lowestBound+Ival+Ew, then we know that less iMinDEE pruning would not add any conformations to our enumeration.  
			 * But with EPIC, we are enumerating in order of minimized E.  We have precomputed the lowestBound.  
			 * If we finish with all lower bounds below lowestBound+Ival+Ew, then we still may have missed some
			 * but if we finish with all minimized E's below lowestBound+Ival+Ew, then we know any missed conf
			 * would have to have minimized E and thus lower bound below lowestBound+Ival+Ew.  This is a contradiction.  
			 * Either way, an Ival of (GMEC so far - lowestBound) is returned and will def be sufficient for the repeated run
			 */
			if(useMinDEEPruningEw){  
				if(es.useEPIC){
					if(lowestBound+Ival+Ew < minE){
						if(es.checkEPIC){
							EPICFitter.analyzeLSBRecord(CETRecord);
							cetm.analyzeFitTypes();
						}
						double IvalTol = RotamerSearch.constRT;//we want to raise the Ival to make sure that
						//we have all conformations up to minE, where minE may vary a little bit from run to run
						//depending on the EpicFitter samples (hence we need a tolerance).  
						//But the variation should be (generally much) less than thermal energy, so we use that as tolerance
						if(run1)//if first conformation, haven't updated best E yet, so just use our single conformational energy so far 
							return new AStarResults(minE,lowestBound,numConfsEvaluated.longValue(),minELowerBound);
//							return (minE-lowestBound+IvalTol);
						else//use best so far
							return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
//							return (getBestE()-lowestBound+IvalTol);
					}
				}
				else if(lowestBound+Ival+Ew < minELowerBound || numConfsEvaluated.longValue() >= enumSettings.numToEnumerate /*|| getBestE()+Ew < minELowerBound*/){ // We are not done and we pruned too much, repeat search
					return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				}			  
			}


			if (done){ //we already have all required conformations

				if (cObj!=null)
					cObj.bestScore = new BigDecimal(getBestE()); //update the best score so far to supply to the next partition

				if(es.checkEPIC){
					EPICFitter.analyzeLSBRecord(CETRecord);
					cetm.analyzeFitTypes();
				}

				/*CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[numConfsEvaluated];				
				System.arraycopy(cObj.confSeq,0,tempConf,0,numConfsEvaluated);				
				cObj.confSeq = tempConf;*/

				if (outputFile){
					logPS.flush();
				}
				if(useMinDEEPruningEw && !approxMinGMEC)
					return new AStarResults(getBestE(),lowestBound,numConfsEvaluated.longValue(),minELowerBound);
				else
					return new AStarResults(lowestBound,lowestBound,numConfsEvaluated.longValue(),minELowerBound); //stop the search
			}

			else{
				if(es.checkEPIC || !es.useEPIC){

					double checkminE = 0;//Check EPIC energy using real energy function

					if (computeEVEnergy){
						if (doMinimization){
							a96ff.calculateTypesWithTemplates();
							a96ff.initializeCalculation();
							a96ff.setNBEval(hElect,hVDW);
							if (!minimizeBB){
								
								if(useCCD){
									efunc = new ForceFieldEnergy(m, a96ff);
									ContSCObjFunction of = new ContSCObjFunction(m,numberOfStrands,efunc,strandRot,doDihedE,null);
									efunc = of.efunc;//ef will now include dihedral energies if appropriate
									ccdMin = new CCDMinimizer(of,false);
								}
								else
									simpMin.initialize(m,numberOfStrands,a96ff,strandRot,doDihedE);
							}
							else { //backbone minimization
								if (!doBackrubs){
									/*if (ligPresent)
                                                                            bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
                                                                            else*/
									bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
								}
								else {
									brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
								}
							}
						}
					}
					if (doMinimization){
						if (!minimizeBB) {//side-chain minimization
							unMinE = calcTotalSnapshotEnergy();
							if(useCCD){//The ideal dihedrals are already recorded since we just reinitialized ccdMin
								ccdMin.minimize();
								minTime = ccdMin.minTime;
								minE = (double)efunc.getEnergy();


								if(es.useEPIC && es.useSVE){
									ccdMin.objFcn.setDOFs(EPICMinVals);
									checkminE = (double)efunc.getEnergy();
								}
							}
							else{
								simpMin.minimize(numMinSteps);
								if(tuples || partitionedRotamers)
									minE = calcTotalSnapshotEnergy(arpMatrix.allMutRes());
								else
									minE = calcTotalSnapshotEnergy();
								if (doDihedE) //add dihedral energies
									minE += simpMin.computeDihedEnergy();
							}
						}
						else {//backbone minimization
							unMinE = calcTotalSnapshotEnergy();
							if (!doBackrubs)
								bbMin.minimizeFull(false);
							else{
								brMin.minimizeFull();
							}
							minE = calcTotalSnapshotEnergy();
						}
						m.updateCoordinates();
						m.revertPertParamsToCurState();

						double totEref = 0.0f;
						double totEntropy = 0.0f;
						if (useEref)
							totEref = getTotSeqEref(arpMatrix.eRef,conf,energyPerPos);
						if (EnvironmentVars.useEntropy)
							totEntropy = getTotSeqEntropy(strandMut,energyPerPos);
						unMinE += totEntropy - totEref;
						minE += totEntropy - totEref;

						if(es.checkEPIC)//record different conformational energy estimates
							CETRecord.add(new double[]{minELowerBound,polyE,minE,minTime,minTimeEPIC});

						checkminE += totEntropy - totEref;

					}
					else if (computeEVEnergy){ //no minimization, so traditional DEE
						//minE = calcTotalSnapshotEnergy(); //the sum of the precomputed energy terms
						minE = minELowerBound;
						unMinE = minE;
						//						m.updateCoordinates();
						//						m.revertPertParamsToCurState();

					}

					if(es.useEPIC && es.useSVE)
						System.out.println("Full energy function evaluated at EPIC minimum: "+checkminE);
				}

				//If using partitioned rotamers store the conformation information
				if(doMinimization){
					computeEnergyPerPos(energyPerPos);
					minIval = storeGeneratedConformation(useEref, curAANums, curRotNums,
							boundPerPos, energyPerPos, arpMatrix.resByPos,minIval);
				}
				

				updateBestE(minE); //for the halting condition
				////////////////////////////////////////////////////////////////////////////////////////////

				
				
				System.out.println(minELowerBound+" "+minE+" "+getBestE());				

				//Since we only need to save the information for conformations whose energy is within Ew of
				//	the minGMEC, we can compare the minimized energy for the current conformation to the
				//	lowest minimized energy in the search so far and output only the conformations that
				//	are within Ew. This optimization is important, since writing to the output file is
				//	relatively very exomputationally expensive; moreover, the output file can become very
				//	big, so the optimization reduces the storage requirement. This approach is most beneficial
				//	when low minimized energies are returned early in thse serach, so there will be only a small
				//	number of extra output conformations.
				if (outputFile){//Output to file

					if ((approxMinGMEC)||(minE<=(getBestE()+Ew))){ //heuristic stopping condition or minE within Ew of the current lowest energy

						numConfsOutput++;
						logPS.print(numConfsOutput+" ");
						for (int i=0; i<treeLevels; i++){
							if(parentConf != null) //will only not be null when using subRot or partitioned rotamers
								logPS.print(parentConf[i].eme.printRes(m,arpMatrix.resByPos));
							else
								logPS.print(conf[i].eme.printRes(m,arpMatrix.resByPos));	
						}

						logPS.print("unMinE: "+unMinE+" ");
						logPS.print("minE: "+minE+" ");
						if (doMinimization)
							logPS.print("minBound: "+minELowerBound+" ");
						logPS.print("bestE: "+getBestE());
						logPS.print(" timeToConf: "+ (System.currentTimeMillis() - KSParser.metrics.loopStart)+" ");
						logPS.print("numConfs: "+(KSParser.metrics.totalNumConfs+numConfsEvaluated.longValue()+1)+" ");
						if(MSAStarSearch != null){
							logPS.print("numExpanded: "+(KSParser.metrics.numExpanded + MSAStarSearch.numExpanded)+" ");
							logPS.print("totNumNodes: "+(KSParser.metrics.totNumNodes + MSAStarSearch.curExpansion.numNodes())+" ");
						}
						logPS.println();
						logPS.flush();
						
						if(outSettings != null && outSettings.savePDBs){
							
							String filename = String.format(outSettings.pdbOutDir+"/%1$s_%2$03d_min.pdb",outSettings.runName,numConfsOutput);
							m.saveMolecule(filename, minE);
						}
						
					}
				}
			}

			numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);

			run1 = false;
		}
	}


	private double storeGeneratedConformation(boolean useEref, int[] curAANums,
			int[] curRotNums, Index3wVal[] boundPerPos,
			Index3wVal[] energyPerPos,ArrayList<ArrayList<Integer>> resByPos, double minIval) {
		if(tuples || partitionedRotamers || superRotamers){

			//KER: Also, generate Global Confs to be returned
			int numMutRes = arpMatrix.numMutRes();
			int[] globalRots = new int[numMutRes];
			int[] resNum = new int[numMutRes];
			double[] EforRes = new double[numMutRes];

			int resCtr = 0;
			int indexCtr = 0;
//					System.out.println("EforRes: ");
			for(ArrayList<Integer> resAtPos: resByPos){
				int AA = curAANums[resCtr];
				int ROT = curRotNums[resCtr];
				int[] rotIndex = {resCtr,AA,ROT};	
				int[] rotamers = arpMatrix.singles.getRot(rotIndex);
				for(int j=0; j<resAtPos.size();j++){
					resNum[indexCtr] = resAtPos.get(j);
					globalRots[indexCtr] = rotamers[j];

					EforRes[indexCtr] += a96ff.singleE.get(resNum[indexCtr]);
					
					if(useEref || EnvironmentVars.useEntropy){
						Residue r = m.residue[resNum[indexCtr]];
						String pdbNum = r.getResNumberString();
						Strand s = m.strand[r.strandNumber];
						ResidueConformation rc = s.rcl.getRC(globalRots[indexCtr]);
						if(useEref)
							EforRes[indexCtr] -= arpMatrix.eRef.get(pdbNum)[rc.rot.aaType.index];
						if(EnvironmentVars.useEntropy)
							EforRes[indexCtr] += rc.rot.aaType.entropy;
					}
					
					for(Entry<Amber96ext.Pair, Double> entry: a96ff.pairsE.entrySet()){
						if(entry.getKey().pos1 == resNum[indexCtr])
							EforRes[indexCtr] += entry.getValue()/2;
					}

					indexCtr++;
				}
				resCtr++;
			}
			generatedConfs.add(new GlobalRCConf(globalRots,resNum,EforRes));

			ArrayList<Index3wVal> badRots = new ArrayList<Index3wVal>();

			System.out.print("Diffs: " );
			double diffTotal = 0;
			for(int i=0; i<energyPerPos.length;i++){
				double diff = energyPerPos[i].val - boundPerPos[i].val;
				System.out.print(diff+" ");
				diffTotal += diff;
				badRots.add(new Index3wVal(boundPerPos[i].i3s, diff,i));
			}
			System.out.println("Total: "+diffTotal);


			Collections.sort(badRots);
			if(diffTotal < minIval){
				minIval = diffTotal;
				contractPos = badRots.get(0).pos;
				if(badRots.size() > 1)
					contractPos2 = badRots.get(1).pos;
				else
					contractPos2 = -1;
			}
			worstRots.add(badRots);

		}
		
		return minIval;
	}


	private EMatrixEntryWIndex[] createParentConf(EMatrixEntryWIndex[] conf) {
		EMatrixEntryWIndex[] parentConf = null;
		if(MSAStarSearch instanceof PGgurobiAStarBySubRot) {
			parentConf = ((PGgurobiAStarBySubRot)MSAStarSearch).parentConf;
		}
		else if(MSAStarSearch instanceof WCSPSearch){
			if(conf == null)
				return null;
			if(partitionedRotamers){ //Build Parent Conf
				parentConf = new EMatrixEntryWIndex[conf.length];
				for(int i=0; i<parentConf.length;i++){
					EMatrixEntryWIndex curRot = conf[i];
					Strand s = m.strand[m.residue[arpMatrix.resByPos.get(curRot.pos1()).get(0)].strandNumber];
					int parentGlobalID = s.rcl.getRC(arpMatrix.singles.getRot(curRot.index)[0]).parent; //0 index assumes no superrotamers
					SuperRotamer r = new SuperRotamer(parentGlobalID);
					RotamerEntry rE = new RotamerEntry(curRot.pos1(),r);
					EMatrixEntryWIndex emeWI = new EMatrixEntryWIndex(rE, curRot.index);
					parentConf[i] = emeWI;
				}
			}
		}
		return parentConf;
	}


	private void computeEnergyPerPos(Index3wVal[] energyPerPos) {
		//KER: if superRotamer we have to walk through resByPos
		//KER: but if tuples we have to walk though energyPerPos
//					if(superRotamers){
//						for(int p1=0; p1<resByPos.size();p1++){
//
//							for(int molResNum: resByPos.get(p1)){
//
//								energyPerPos[p1].val += a96ff.singleE.get(molResNum);
//								//System.out.println(p1+" "+p1+" "+(a96ff.singleE.get(j)));
//								for(int p2 = p1+1; p2<resByPos.size();p2++){
//									for(int molResNum2: resByPos.get(p2)){
//
//										Amber96ext.Pair p = a96ff.new Pair(molResNum,molResNum2); 
//										energyPerPos[p1].val += a96ff.pairsE.get(p)/2;
//										energyPerPos[p2].val += a96ff.pairsE.get(p)/2;
//										//System.out.println(p1+" "+p2+" "+a96ff.pairsE.get(p));
//									}
//
//
//								}
//							}
//						}
//					}
		if(tuples || partitionedRotamers){
			ArrayList<Integer> mutRes = arpMatrix.allMutRes();
			double eFunctMinE = 0.0;
			for(double e: a96ff.pairsE.values()){
				eFunctMinE += e;
			}
			//KER: Divide by two because pairs are listed twice
			eFunctMinE /= 2;

			for(double e: a96ff.singleE.values()){
				eFunctMinE += e;
			}

			//System.out.println("Printing Energies: ");
			for(int p1=0; p1<mutRes.size();p1++){
				Index3wVal iVal = null;
				for(Index3wVal i3wV: energyPerPos)
					if(i3wV.containsPos(p1))
						iVal = i3wV;

				iVal.val += a96ff.singleE.get(mutRes.get(p1));
				//System.out.println(p1+" "+p1+" "+(a96ff.singleE.get(mutRes.get(p1))));
				for(int p2 = p1+1; p2<mutRes.size();p2++){
					Index3wVal p2Val = null;
					for(Index3wVal i3wV_2: energyPerPos)
						if(i3wV_2.containsPos(p2))
							p2Val = i3wV_2;

					Amber96ext.Pair p = a96ff.new Pair(mutRes.get(p1),mutRes.get(p2)); 
					iVal.val += a96ff.pairsE.get(p)/2;
					p2Val.val += a96ff.pairsE.get(p)/2;
					//System.out.println(p1+" "+p2+" "+a96ff.pairsE.get(p));
				}


			}
		}
	}


	private EnergiesPerPos initializeBoundAndEnergyPerPos(int[] curAANums,
			int[] curRotNums, LinkedList<EnergyTuple> curTuples) {
		
		EnergiesPerPos epp = new EnergiesPerPos();
		
		int numPos = arpMatrix.numMutPos();
		
		//If there are tuples we need to subtract the number of mut pos
		if(curTuples != null){
			for(EnergyTuple tup: curTuples){
				numPos -= tup.rots.length-1;
			}
		}
		
		Index3wVal[] boundPerPos = new Index3wVal[numPos];
		Index3wVal[] energyPerPos = new Index3wVal[numPos];
		//Add the tuples first
		int ctr=0;
		boolean[] alreadyVisited = new boolean[arpMatrix.numMutPos()];
		if(curTuples != null){
			for(EnergyTuple tup: curTuples){
				boundPerPos[ctr] = new Index3wVal(tup.rots,0);
				energyPerPos[ctr] = new Index3wVal(tup.rots,0);
				ctr++;
				for(int q=0; q<tup.rots.length;q++)
					alreadyVisited[tup.rots[q].pos] = true;
			}
		}

		//Add the remaining residues
		for(int i=0; i<alreadyVisited.length;i++){
			if(!alreadyVisited[i]){
				Index3[] rots = {new Index3(i,curAANums[i],curRotNums[i])};
				boundPerPos[ctr] = new Index3wVal(rots,0);
				energyPerPos[ctr] = new Index3wVal(rots,0);
				ctr++;
			}

		}
		
		epp.boundPerPos = boundPerPos;
		epp.energyPerPos = energyPerPos;
		
		return epp;
	}


	private double getTotSeqEref(HashMap<String, double[]> eRef, EMatrixEntryWIndex[] conf,Index3wVal[] energyPerPos) {
		double totEref = 0;
		for(EMatrixEntryWIndex emeWI: conf){
			double tmpE = ((RotamerEntry)emeWI.eme).getEref(eRef, m, arpMatrix.resByPos);
			totEref += tmpE;
//			if(superRotamers)
//				energyPerPos[emeWI.pos1()].val -= tmpE;
			if(tuples || partitionedRotamers)
				for(Index3wVal i3wV: energyPerPos)
					if(i3wV.containsPos(emeWI.pos1()))
						i3wV.val -= tmpE;
		}
		return totEref;
	}


	//Returns the reference energy for the current amino acid sequence assignment (for the mutatable positions only)
	public double getTotSeqEntropy(MutableResParams strandMut,Index3wVal[] energyPerPos){

		double totEref = 0.0f;

		for(int i=0; i<strandMut.numMutPos();i++){
			int str = strandMut.resStrand[i];
			if(m.strand[str].isProtein){
				int strResNum = strandMut.resStrandNum[i];
				String resName = m.strand[str].residue[strResNum].name;
				double tmpE = EnvironmentVars.getEntropyTerm(resName); 
				totEref += tmpE;
//				if(superRotamers)
//					energyPerPos[emeWI.pos1()].val += tmpE;
				if(tuples || partitionedRotamers){
					//Find position in matrix
					int molResNum = strandMut.allMut[i];
					int pos1 = 0;
					for(int curPos=0; curPos<arpMatrix.resByPos.size();curPos++){
						for(int res: arpMatrix.resByPos.get(curPos))
							if(res == molResNum)
								pos1 = i;
						
					}
					for(Index3wVal i3wV: energyPerPos)
						if(i3wV.containsPos(pos1))
							i3wV.val += tmpE;
				}
			}
		}


		return totEref;
	}


	//	private int getAAIndex(int rotIndex, int curRes, String strandDefault[][], int strandMut[][]){
	//
	//		int rotSum = 0;
	//
	//		int str = mutRes2Strand[curRes];
	//		int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
	//		int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
	//
	//		for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){
	//			int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,i);
	//			int curRot = getNumRot( str, strResNum, curAA );
	//
	//			rotSum += curRot;
	//			if (rotSum>rotIndex)
	//				return curAA;
	//		}
	//		return -1;
	//	}

	//	private int getRotSum(int rotIndex, int curRes, String strandDefault[][], MutableResParams strandMut){
	//
	//		int rotSum = 0;
	//		int str = strandMut.resStrand[curRes];
	//		int strResNum = strandMut.resStrandNum[curRes];
	//		for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){
	//			int curRot = getNumRot( str, strResNum, strandRot[str].getIndexOfNthAllowable(strResNum,i) );
	//
	//			if ((rotSum+curRot)>rotIndex)
	//				return rotSum;
	//			else
	//				rotSum += curRot;
	//		}
	//		return -1;
	//	}


	//find how many rotamers or RCs are allowed for the given AA type index at the given
	//residue of the given strand
	// giving 1 for gly or ala
	//	public int getNumRot(int str, int strResNum, int AAind){
	//
	//		if(doPerturbations){
	//			return ((StrandRCs)strandRot[str]).getNumRCs(strResNum, AAind);
	//		}
	//		else{
	//			int ans = strandRot[str].rl.getNumRotForAAtype(AAind);
	//			if (ans==0)	//ala or gly
	//				return 1;
	//			else
	//				return ans;
	//		}
	//	}

	//////////////////////////////////////////////////////////////////////////
	//	End of Compute Min GMEC section
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//	Generate Backbones section
	//////////////////////////////////////////////////////////////////////////
	public void doGenBackbones(String runName, int numMutable, MutableResParams strandMut, double theta, double alpha,
			int numSamples, boolean systematicSampling){

		//setup the log file to store information about the generated backbones
		PrintStream logPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream((runName+".bb.all"));
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
			System.exit(1);
		}

		Backbone bb = new Backbone();

		//the initial phi and psi angles for the residues with flexible backbones
		double initFi[] = new double[numMutable];
		double initPsi[] = new double[numMutable];

		for (int i=0; i<numMutable; i++){
			Residue r = m.residue[strandMut.allMut[i]];
			int str = r.strandNumber;
			int strResNum = r.strandResidueNumber;
			int molResNum = r.moleculeResidueNumber;

			initFi[i] = bb.getFiPsi(m, str, strResNum, 0);
			initPsi[i] = bb.getFiPsi(m, str, strResNum, 1);

			System.out.println("AS residue: "+i+" ("+m.strand[str].residue[strResNum].getResNumber()+") phi: "+initFi[i]+" psi: "+initPsi[i]);

			if ((initFi[i]==0.0)||(initPsi[i]==0.0)){
				System.out.println("ERROR: Strand "+str+" Residue "+strResNum+" does not have a valid (phi,psi) pair.");
				System.exit(1);
			}
		}

		if (systematicSampling){ //systematic sampling (generates a large number of backbones)
			//the number of steps (from -theta to +theta)
			int numSteps = 1;
			if (alpha!=0.0)
				numSteps = 2*(int)(theta/alpha) + 1; //count the initial phi/psi values as well

			doGenBackbonesHelper(runName, numMutable, strandMut, theta, alpha, numSteps, bb, initFi, initPsi, logPS, 0);
		}
		else { //random sampling (generates a small number of backbones)
			outputBB(logPS, numMutable, runName, strandMut, bb); //output the initial backbone first
			doGenBackbonesRandHelper(runName, numMutable, strandMut, theta, numSamples, bb, logPS);
		}
	}

	//Generates up to numSamples backbones by applying random phi/psi changes within theta degrees 
	//		This version is used to generate a small number of backbones (the systematic version below
	//		will generate a very large number of backbones, even for a very small number of steps)
	private void doGenBackbonesRandHelper(String runName, int numMutable, MutableResParams strandMut, double theta,
			int numSamples, Backbone bb, PrintStream logPS){

		Random randNum = new Random();

		for (int curSample=0; curSample<numSamples; curSample++){ //generate up to numSamples backbones

			double curFiChange[] = new double[numMutable]; //the phi changes for each residue
			double curPsiChange[] = new double[numMutable]; //the psi changes for each residue

			for(int curRes=0; curRes<numMutable; curRes++){ //apply the random (phi,psi) changes for each residue
				int str = strandMut.resStrand[curRes];
				int strResNum = strandMut.resStrandNum[curRes];
				int molResNum = strandMut.allMut[curRes];
				//Get the random (phi,psi) change
				curFiChange[curRes] = 2*(randNum.nextDouble()-0.5f)*(double)theta; //phi change within (-theta,theta)
				curPsiChange[curRes] = 2*(randNum.nextDouble()-0.5f)*(double)theta; //psi change within (-theta,theta)

				//Apply the (phi,psi) change
				bb.applyFiPsi(m,str,strResNum,curFiChange[curRes],0);
				bb.applyFiPsi(m,str,strResNum,curPsiChange[curRes],1);
			}

			//we have a full backbone conformation, so output
			outputBB(logPS, numMutable, runName, strandMut, bb);

			for(int curRes=0; curRes<numMutable; curRes++){ //restore the original (phi,psi) for each residue
				int str = strandMut.resStrand[curRes];
				int strResNum = strandMut.resStrandNum[curRes];
				bb.applyFiPsi(m,str,strResNum,-curFiChange[curRes],0);
				bb.applyFiPsi(m,str,strResNum,-curPsiChange[curRes],1);
			}
		}
	}

	//Recursively generates backbones by applying phi/psi changes within theta degrees 
	//		at a step of alpha degrees to each residue with a flexible backbone (systematic sampling)
	private void doGenBackbonesHelper(String runName, int numMutable, MutableResParams strandMut, double theta, double alpha, 
			int numSteps, Backbone bb, double initFi[], double initPsi[], PrintStream logPS, int curRes){

		if (curRes==numMutable){//we have a full backbone conformation, so output
			outputBB(logPS, numMutable, runName, strandMut, bb);
		}
		else { //we only have a partial backbone, so apply changes to the current residue

			if (curRes==0) 
				System.out.println("Starting..");

			//First, move the phi/psi angle to -theta, then apply the changes to +theta, at alpha steps

			Residue r = m.residue[strandMut.allMut[curRes]];
			int str = r.strandNumber;
			int strResNum = r.strandResidueNumber;

			//move phi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
			bb.applyFiPsi(m,str,strResNum,-(theta+alpha),0);

			//apply phi changes up to +theta
			for (int curStepFi=0; curStepFi<numSteps; curStepFi++){ //apply each alpha step up to a displacement of +theta from the initial phi

				if (curRes==0)
					System.out.print("*");

				bb.applyFiPsi(m,str,strResNum,alpha,0);

				//move psi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
				bb.applyFiPsi(m,str,strResNum,-(theta+alpha),1);

				for (int curStepPsi=0; curStepPsi<numSteps; curStepPsi++){ //apply each alpha step up to a displacement of +theta from the initial psi

					if (curRes==0)
						System.out.print(".");

					bb.applyFiPsi(m,str,strResNum,alpha,1);
					if (checkStericsBBonly(str,strResNum)) {//passed steric test, so move to the next residue
						doGenBackbonesHelper(runName,numMutable,strandMut,theta,alpha,numSteps,bb,initFi,initPsi,logPS,curRes+1);
					}
				}

				//restore initial psi
				bb.applyFiPsi(m,str,strResNum,-theta,1);

				if (curRes==0)
					System.out.println();
			}

			//restore initial phi
			bb.applyFiPsi(m,str,strResNum,-theta,0);

			if (curRes==0)
				System.out.println("done");
		}
	}

	//Checks if the current backbone is sterically allowed and outputs the pdb and log information
	private void outputBB(PrintStream logPS, int numMutable, String runName, MutableResParams strandMut, Backbone bb){

		//Check all residues for sterics against the whole strand since all residues
		//		are already assigned (up to this point, we have checked for sterics only
		//		against the residues up to a given AS residue);
		//The backbone movement may have actually caused unallowed sterics between rigid
		//		parts of the molecule, so we check for this possibility
		boolean stericAllowed = true;
		for(int str=0; str<numberOfStrands;str++){
			for (int i=0; i<m.strand[str].numberOfResidues; i++){
				if (!checkAllStericsBBonly(str,i)){
					stericAllowed = false;
					break;
				}
			}
		}

		if (stericAllowed){ //all sterics are allowed

			String fileName = (runName+System.currentTimeMillis()+".pdb");

			m.saveMolecule(fileName, 0.0f);//save the molecule

			//output the molecule file name and all (phi,psi) pairs for the residues with flexible backbones
			logPS.print(fileName+" ");
			for (int i=0; i<strandMut.allMut.length; i++){
				Residue r = m.residue[strandMut.allMut[i]];
				logPS.print("( "+bb.getFiPsi(m, r.strandNumber, r.strandResidueNumber, 0)+" , "+bb.getFiPsi(m, r.strandNumber, r.strandResidueNumber, 1)+" ) ");
			}

			logPS.println();
			logPS.flush();
		}
	}
	// This function checks the sterics of the given conformation;
	//  it checks backbone atoms of the given residue resNum against all
	//  backbone atoms that are in residues 0..resNum-1 for the given strand only
	//  If any two atoms overlap by more than overlapThresh then
	//  false is returned
	// strandNum is the number of the strand containing resNum
	// Hydrogens are NOT used in checking sterics in this function
	private boolean checkStericsBBonly(int strandNum, int resNum) {

		Residue res = m.strand[strandNum].residue[resNum];

		Atom tmpAtm = null;
		int resToCheck = 0;

		for(int i=0;i<res.numberOfAtoms;i++) {

			if (isBBatom(res.atom[i])){ //backbone atom

				resToCheck = resNum;
				for(int w=0;w<resToCheck;w++) {
					for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
						tmpAtm = m.strand[strandNum].residue[w].atom[t];

						if (isBBatom(tmpAtm)){
							if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
								if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - softOverlapThresh)){
									if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
										return false;
									}
								}
							}
						}
					}
				}	
			}
		}

		return true; //if we are here, then everything passed the steric test
	}

	//This function is similar to checkStericsBBonly(), but it checks all residues in strandNum against resNum, 
	//		instead of just checking the residues up to resNum
	private boolean checkAllStericsBBonly(int strandNum, int resNum) {

		Residue res = m.strand[strandNum].residue[resNum];

		Atom tmpAtm = null;
		int resToCheck = 0;

		for(int i=0;i<res.numberOfAtoms;i++) {

			if (isBBatom(res.atom[i])){ //backbone atom

				resToCheck = m.strand[strandNum].numberOfResidues;
				for(int w=0;w<resToCheck;w++) {
					if (w!=resNum){
						for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
							tmpAtm = m.strand[strandNum].residue[w].atom[t];

							if (isBBatom(tmpAtm)){
								if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
									if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - softOverlapThresh)){
										if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
											return false;
										}
									}
								}
							}
						}
					}
				}	
			}
		}

		return true; //if we are here, then everything passed the steric test
	}

	//Determines if the given atom is a backbone atom
	private boolean isBBatom(Atom at){

		return ((at.name.equalsIgnoreCase("N"))||(at.name.equalsIgnoreCase("CA"))
				||(at.name.equalsIgnoreCase("C"))||(at.name.equalsIgnoreCase("O")));
	}
	//////////////////////////////////////////////////////////////////////////
	//	End of Generate Backbones section
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////

	/*//Determines how many residue positions in the system strand (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(int pos, double dist, boolean as[]){

		Residue r1 = m.strand[sysStrNum].residue[pos];
		for (int i=0; i<m.strand[sysStrNum].numberOfResidues; i++){
			if (i!=pos){
				Residue r2 = m.strand[sysStrNum].residue[i];
				if (r1.getDistSC(r2)<=dist)
					as[i] = true;
				else
					as[i] = false;
			}
		}

		return as;
	}*/


	//////////////////////////////////////////////////////////////////////////////////


	//Computes and stores the pairwise energies for the current molecule; this can be used to compare
	//		the actual energies for each pair of residues (computed here) against the respective lower bounds computed
	//		during the pairwise energy minimization
	/*private void pairwiseEnergySidechainMolecule (int residueMap[]){

		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		boolean savedFlexible[] = new boolean[m.numberOfResidues];

		// Save the energy eval flag, clear them at the same time
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			savedFlexible[i] = m.residue[i].flexible;
		}

		PrintStream logPS = setupOutputFile("energies.out");

		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(false, false);
			m.residue[i].flexible = false;
		}

		for (int i=0; i<residueMap.length; i++){

			m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(true, true);
			m.strand[sysStrNum].residue[residueMap[i]].flexible = true;

			double intraEi = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[i]].moleculeResidueNumber);

			for (int j=i+1; j<residueMap.length; j++){

				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(true, true);
				m.strand[sysStrNum].residue[residueMap[j]].flexible = true;

				double intraEj = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[j]].moleculeResidueNumber);

				double pairE = calcTotalSnapshotEnergy() - intraEi - intraEj;

				logPS.println(i+" "+j+" "+pairE+" "+intraEj);

				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(false, false);
				m.strand[sysStrNum].residue[residueMap[j]].flexible = false;
			}

			for(int j=0;j<m.numberOfResidues;j++){
				m.residue[j].setEnergyEval(true, true);
				m.residue[j].flexible = false;
			}
			if (ligStrNum != -1) {
				for(int j=0;j<m.strand[ligStrNum].numberOfResidues;j++)
					m.strand[ligStrNum].residue[j].setEnergyEval(false, false);
			}
			for(int j=0;j<residueMap.length;j++)
				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(false, false);
			m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(true, true);
			m.strand[sysStrNum].residue[residueMap[i]].flexible = true;

			double shellE = computeEnergyOfOnlyTemplate(residueMap);

			double resShellE = calcTotalSnapshotEnergy() - intraEi - shellE;

			for(int j=0;j<m.numberOfResidues;j++){
				m.residue[j].setEnergyEval(false, false);
				m.residue[j].flexible = false;
			}

			logPS.println(i+" "+intraEi+" "+resShellE+" "+shellE);
			logPS.flush();
		}
		logPS.close();

		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			m.residue[i].flexible = savedFlexible[i];
		}
	}*/

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


	public String printBigNum(BigDecimal bd, int sigDigits){
		int numDigits = (bd.toPlainString()).indexOf('.');
		if(numDigits == -1)
			numDigits = (bd.toPlainString().length());

		int newScale = -1*(numDigits - sigDigits);

		return bd.setScale(newScale,BigDecimal.ROUND_HALF_UP).toString();
	}

	public void resetMatrices() {
		arpMatrix = null;
		cetm = null;
	}

	public static double[][] loadErefMatrix(String fileName){
		double[][] eref = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
			eref = (double [][])in.readObject();
			in.close();
		}
		catch (Exception e){}
		return eref;
	}

	//	public void addEntropyTerm(boolean doMinimize, MutableResParams strandMut){
	//
	//		//NumMutRes (Don't count non-protein strands)
	//		int numMutRes = strandMut.allMut.length;
	//
	//		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
	//		for (int i=0; i<numMutRes; i++){
	//			int str = strandMut.resStrand[i];
	//			int strResNum = strandMut.resStrandNum[i];
	//			if(m.strand[str].isProtein){
	//				for (int j=0; j<strandRot[str].getNumAllowable(strResNum); j++){
	//					int aaInd = strandRot[str].getIndexOfNthAllowable(strResNum,j);
	//					String AAname = strandRot[str].rl.getAAName(aaInd);
	//					int numRot = getNumRot(str, strResNum, aaInd);
	//					for (int k=0; k<numRot; k++){
	//						arpMatrix.addToIntraE( i, aaInd, k, (double) EnvironmentVars.getEntropyTerm(AAname) );
	//						//						if (doMinimize)
	//						//							arpMatrixMax.addToIntraE( i, aaInd, k, (double) EnvironmentVars.getEntropyTerm(AAname) );
	//						ind++;
	//					}
	//				}
	//			}
	//			else{
	//				System.out.println("Entropy term not added for non-amino acid residue: "+m.strand[str].residue[strResNum].fullName);
	//			}
	//		}
	//
	//	}

	//	private void saveConf(int[] conf, double energy, String filename, MutableResParams strandMut,
	//			boolean minimizeBB, boolean doBackrubs) {
	//
	//		applyRotamers(strandMut, conf);
	//		m.backupAtomCoord();//The backup coordinates will now be in the unminimized state so we can return to it
	//
	//		if(doMinimization){
	//			if(!minimizeBB)
	//				simpMin.minimize(numMinSteps);
	//			else{
	//				if(!doBackrubs)
	//					bbMin.minimizeFull(false);
	//				else{
	//					System.out.println("BACKRUBS NOT IMPLEMENTED YET");
	//					//brMin.minimizeFull();
	//				}
	//			}
	//		}
	//		m.saveMolecule(filename, energy);
	//		m.revertPertParamsToCurState();
	//		m.restoreAtomCoord();
	//		m.updateCoordinates();
	//
	//	}




	//Check parametric incompatibility of a pair of RCs, with their residues and AA types also given
	public boolean isParametricallyIncompatible(Residue firstRes, ResidueConformation RC1, Residue secondRes, ResidueConformation RC2){

		//Get the residues' residue perturbations states
		int pertState1 = RC1.pertState;//((StrandRCs)strandRot[firstRes.strandNumber]).RCPertStates[firstRes.strandResidueNumber][curAA1][RC1];
		int pertState2 = RC2.pertState;//((StrandRCs)strandRot[secondRes.strandNumber]).RCPertStates[secondRes.strandResidueNumber][curAA2][RC2];
		//Now see if the residue perturbation states are compatible

		for(int p1=0;p1<firstRes.perts.length;p1++){//Check all the perturbations of firstRes in this RC for compatibility
			for(int p2=0;p2<secondRes.perts.length;p2++){
				if(firstRes.perts[p1] == secondRes.perts[p2]){//The residues share a perturbation
					if(firstRes.pertStates[pertState1][p1] != secondRes.pertStates[pertState2][p2])//The states of this perturbation don't match for the given RCs
						return true;//So the RCs are parametrically incompatible
				}
			}
		}

		return false;//No incompatibility found, so we're good
	}

	/**
	 * The eRef matrix is based on the initial mutRes indexes. Therefore,
	 * this function must be called before any positions are merged.
	 * @param eRef
	 * @param pos
	 */
	public static void addEref(Emat emat,Molecule m,HashMap<String,double[]> eRef, int pos, ArrayList<Index3> rotamersToUpdate) {
		//Loop through emat rotamers and add eRef
		SinglesIterator iter;
		if(pos < 0)
			iter = emat.singlesIterator();
		else
			iter = emat.singlesIterator(pos);
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(rotamersToUpdate == null || rotamersToUpdate.contains(new Index3(emeWI.index[0],emeWI.index[1],emeWI.index[2])))
				emat.addE(emeWI.index,-((RotamerEntry)emeWI.eme).getEref(eRef,m,emat.resByPos));
		}

		emat.hasEref = true;
	}

	public static void addEref(Emat emat, Molecule mol, HashMap<String,double[]> eRef,int pos){
		addEref(emat,mol,eRef,pos,null);
	}

	public void addEref(HashMap<String,double[]> eRef,int pos){
		addEref(arpMatrix,m,eRef,pos,null);
	}

	public void addEref(HashMap<String,double[]> eRef){
		addEref(arpMatrix,m,eRef,-1,null);
	}

	public void addEref(Emat emat, HashMap<String,double[]> eRef){
		addEref(emat,m,eRef,-1,null);
	}

	/**
	 * The eRef matrix is based on the initial mutRes indexes. Therefore,
	 * this function must be called before any positions are merged.
	 * @param eRef
	 * @param pos
	 */
	/**
	 * The eRef matrix is based on the initial mutRes indexes. Therefore,
	 * this function must be called before any positions are merged.
	 * @param eRef
	 * @param pos
	 */
	public static void addEntropyTerm(Emat emat,Molecule m, int pos, ArrayList<Index3> rotamersToUpdate) {

		//Loop through emat rotamers and add eRef
		SinglesIterator iter;
		if(pos < 0)
			iter = emat.singlesIterator();
		else
			iter = emat.singlesIterator(pos);
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(rotamersToUpdate == null || rotamersToUpdate.contains(new Index3(emeWI.index[0],emeWI.index[1],emeWI.index[2])))
				emat.addE(emeWI.index,((RotamerEntry)emeWI.eme).getEntropy(m,emat.resByPos));
		}

		emat.hasEntropy = true;
	}

	public static void addEntropyTerm(Emat emat, Molecule mol, int pos){
		addEntropyTerm(emat,mol,pos,null);
	}

	public void addEntropyTerm(int pos){
		addEntropyTerm(arpMatrix,m,pos,null);
	}

	public void addEntropyTerm(){
		addEntropyTerm(arpMatrix,m,-1,null);
	}

	public void addEntropyTerm(Emat emat){
		addEntropyTerm(emat,m,-1,null);
	}


	public void printTopConfs(int runNum, int curMut, boolean saveTopConfs, boolean printTopConfs,
			boolean minimizeBB, boolean doBackrubs, String outDir) {
		//Print the top structures
		if((saveTopConfs || printTopConfs) && runNum >= 0){
			ConfPair confPairs[] = new ConfPair[topConfs.size()];
			confPairs = topConfs.toArray(confPairs);
			Arrays.sort(confPairs);
			int ctr = 0;
			PrintStream printStream = null;
			String outName = ""+curMut+"_";

			for(Residue r: m.residue){
				if(r.isMutable)
					outName += RotamerLibrary.getOneLet(r.name);
			}

			outName += "_"+runNum;
			outName = outDir+"/"+outName;
			File tmpDir = new File(outDir);
			if(tmpDir != null && !tmpDir.exists()){
				tmpDir.mkdirs();
			}
			if(printTopConfs){
				try{
					FileOutputStream fileOutputStream = new FileOutputStream(outName);
					BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( 
							fileOutputStream);
					printStream = new PrintStream(bufferedOutputStream);
				}
				catch(Exception e){
					e.printStackTrace();
					System.out.println("Couldn't write confOut file.");
				}
			}
			for(int i = confPairs.length-1; i>=0; i--){
				String filename = String.format("%1$s_%2$03d.pdb",outName,ctr++);
				if(saveTopConfs){
					saveConf(confPairs[i].conf, confPairs[i].energy[0], filename,minimizeBB, doBackrubs);
					System.out.println("Saved Conf " + i);
				}
				if(printStream != null){ //printTopConfs
					printStream.print(""+ctr+" "+confPairs[i].energy[0]+" ");
					for(Residue r: m.residue){
						if(r.isMutable)
							printStream.print(r.name+" ");
					}
					for (int j=0; j<confPairs[i].conf.length; j++){
						printStream.print(confPairs[i].conf[j].eme.printRes(m,arpMatrix.resByPos));	
					}

					printStream.print("unMin: "+confPairs[i].energy[1]+" eParts: " + confPairs[i].energy[2]+" "+confPairs[i].energy[3]+" "+confPairs[i].energy[4]+" "+confPairs[i].energy[5]);

					printStream.println("");
				}
			}
			if(printTopConfs || printStream != null){
				printStream.close();
			}	
		}

	}


	private void saveConf(EMatrixEntryWIndex[] conf, double energy, String filename, 
			boolean minimizeBB, boolean doBackrubs) {
		m.backupAtomCoord();
		applyRCs(conf);
		if(doMinimization){
			if(!minimizeBB){
				if(useCCD)
					ccdMin.minimize();
				else
					simpMin.minimize(numMinSteps);
			}
			else{
				if(!doBackrubs)
					bbMin.minimizeFull(false);
				else{
					System.out.println("BACKRUBS NOT IMPLEMENTED YET");
					//brMin.minimizeFull();
				}
			}
		}
		m.saveMolecule(filename, energy);
		m.restoreAtomCoord();
		m.updateCoordinates();
		m.revertPertParamsToCurState();
	}


	/**
	 * @param conf
	 */
	public void applyRCs(EMatrixEntryWIndex[] conf) {
		//Apply the rotamers of the current conformation


		for(int i=0; i<arpMatrix.resByPos.size();i++){
			for(int j=0; j<arpMatrix.resByPos.get(i).size();j++){
				Residue r = m.residue[arpMatrix.resByPos.get(i).get(j)];
				MutUtils.applyRC(m, r, m.strand[r.strandNumber].rcl.getRC(arpMatrix.singles.getRot(conf[i].index)[j]));

			}
		}
	}
	
	private class EnergiesPerPos{
		Index3wVal[] boundPerPos;
		Index3wVal[] energyPerPos;
	}

}//end of RotamerSearch class
