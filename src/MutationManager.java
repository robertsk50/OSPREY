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

////////////////////////////////////////////////////////////////////////////////////////////
// MutationManager.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   --------------      ------------------------     ------------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2002-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
import java.math.*;

import mpi.*;

/**
 * The MutationManager class maintains a list of mutations to be tested, maintains
 *  their scores, prints a log file, and generally manages the mutations to test.
 */
public class MutationManager
{

	//the algorithm options that define what pruning criteria will be applied
	//	NOTE!!! These must be the same as in KSParser.java
	final int optSimple = 1;
	final int optBounds = 2;
	final int optSplit = 3;
	final int optPairs = 4;

	final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)

    Queue<CommucObj> waitingRuns;
    PartitionMessage[][] completedRuns = null;
	
	// Information needed by all mutations
	CommucObj cObjArray[] = null;
	//int residueMap[] = null;
	//String resDefault[] = null;
	MutableResParams strandMut = null;
	String[][] strandDefault = null;
	boolean[] strandPresent = null;
	String[][] strandLimits = null;
	int strandsPresent = 0;
	boolean typeDep = false;
	//String ligType = null;
	boolean ligPresent = false;
	int numMutations = 0;
	String arpFilenameMin = null;
	String arpFilenameMax = null;
	boolean resMutatable[][] = null;
	int algOption = 0;
	int numSplits = 0;
	String AAallowed[][] = null;
	String minDEEfile = null;
	double initEw = 0.0f;
	double pruningE = (double)Math.pow(10,38);
	double gamma = 0.01;  // The gamma used in inter-mutation pruning
	double epsilon = 0.03f;  // The epsilon used in intra-mutation pruning
	double stericThresh = -10000.0f;
	double softStericThresh = -10000.0f;
	//int numInAS = 0;
	int mutableSpots = 0;
	boolean computeEVEnergy = true;
	boolean doMinimization = true;
	boolean minimizeBB = false;
	boolean doBackrubs = false;
	String backrubFile = null;
	boolean repeatSearch = true;
	boolean calculateVolumes = true;
	BigDecimal bestScore = new BigDecimal("0.0");
	BigDecimal q_L = new BigDecimal("0.0");
	ParamSet sParams = null;
	boolean approxMinGMEC = false;
	double lambda = (double)Math.pow(10,38);
	double stericE = Math.pow(10,38);
	boolean distDepDielect = true;
	double dielectConst = 1.0;
	boolean doDihedE = false;
	boolean doSolvationE = false;
	double solvScale = 1.0;
	double vdwMult = 1.0;
	boolean scaleInt = false;
	double maxIntScale = 1.0f;
	boolean useEref = false;
	boolean entropyComp = false; //this *must* be false for the pairwise matrix energy computation
	double asasE[][][][] = null;
	boolean compASdist = false;
	boolean asDist[][] = null;
	double dist = 0.0f;
	int numberOfStrands = 0;
	PrintStream logPS = null;
	OneMutation mutArray[] = null;
	int mutRes2Strand[] = null;
	int mutRes2StrandMutIndex[] = null;	
	Molecule m;

	private boolean saveTopConfs;
	private boolean printTopConfs;
	private int numTopConfs;
	private KSParser.ASTARMETHOD asMethod;
	
	//Variables specific to PEM computation	
	Emat pairEMatrixMin = null;
	//PairwiseEnergyMatrix pairEMatrixMax = null;
	//double[][] eRefMatrix = null;
	double curMaxE = -(double)Math.pow(10,30);
	//int numLigRotamers = 0;

	double pairEMatrixMinEntropy[][] = null;
	//KER: This is only used for the SCMF entropy run.
	RotamerLibrary rotLib = null;

	//Variables specific to distributed DACS and distributed DEE computations
	//PrunedRotamers<Boolean> prunedRot = null;
	String rotFile = null;
	boolean useSF = false;
//	boolean splitFlags[][][][][][] = null;
//	String sfFile = null;
	boolean distrDACS = false;
	boolean distrDEE = false;
	int numSpPos = -1;
	int msp[] = null;
	int typeDEE = -1;
	int initDepth = -1;
	int subDepth = -1;
	int diffFact = -1;
	double minRatioDiff = 0.0;
	BigInteger numInitUnprunedConf = null;
	String outputPruneInfo = null;
	String outputConfInfo = null;

	boolean neighborList;
	double distCutoff;
	private String pdbOutDir = "";
	private DEEsettings deeSettings;
	
	String pdbName;

	String eRefMatrix;
	boolean PEMcomp = false; //true if PEM computation is performed; false if mut search is performed
	private boolean templateAlwaysOn;
	private boolean addOrigRots = false;
	private boolean useMaxKSconfs;
	private BigInteger numKSconfs;
	private int curStrForMatrix;

	boolean useTriples;
	boolean useFlagsAStar;
	boolean magicBulletTriples;
	int magicBulletNumTriples;
	//DEEPer
	boolean doPerturbations;
	boolean minimizePerts;
	String pertFile;
	boolean addWTRot;
	boolean idealizeSC;

	boolean useCCD;
	boolean minimizePairwise = true;

	double EConvTol;

	boolean compCETM;//we are computing the level-set-based bound matrix
	CETMatrix cetm;
	EPICSettings es;
	private double Ival;
	private KSParser.DEEMETHOD deeMethod;
	boolean doDih = false;
	private HashMap<String, double[]> eRef;

	// Generic constructor
	MutationManager(String logName, OneMutation mArray[], boolean PEMcomputation) {
		pdbName = logName;
		if (logName!=null){ //open log file for writing
			try {
				FileOutputStream fileOutputStream = new FileOutputStream(logName);
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				logPS = new PrintStream( bufferedOutputStream );
			}
			catch (Exception ex) {
				System.out.println("ERROR: An exception occured while opening log file");
			}
		}

		mutArray = mArray;
		cObjArray = new CommucObj[mutArray.length];

		PEMcomp = PEMcomputation;

		//This is only used for K* runs
		waitingRuns = new LinkedList<CommucObj>();
		completedRuns = new PartitionMessage[mutArray.length][];
		
		distrDACS = false; //if this is a distributed DACS computation, the flag will be set with setDistrDACS()
		distrDEE = false; //if this is a distributed DEE computation, the flag will be set with setDistrDEE()
	}

	// Returns the next mutation packaged in a communication object
	public synchronized CommucObj getNextComObj(int curMutIndex) {

		//to debug a given residue pair's energy precomputation
		//curMutIndex = 34;

		CommucObj cObj = new CommucObj();
		cObj.numberOfStrands = numberOfStrands;
		cObj.strandMut = strandMut;
		cObj.strandDefault = strandDefault;
		cObj.typeDep = typeDep;
		cObj.addOrigRots = addOrigRots;
		//cObj.ligPresent = ligPresent;
		//cObj.ligType = ligType;
		cObj.distrDEE = distrDEE;
		cObj.distrDACS = distrDACS;
		cObj.arpFilenameMin = arpFilenameMin;
		cObj.arpFilenameMax = arpFilenameMax;
		cObj.params = sParams;
		cObj.stericThresh = stericThresh;
		cObj.softStericThresh = softStericThresh;
		//cObj.numInAS = numInAS;
		cObj.saveTopConfs = saveTopConfs;
		cObj.printTopConfs = printTopConfs;
		cObj.numTopConfs = numTopConfs;
		cObj.rl = rotLib;
		cObj.computeEVEnergy = computeEVEnergy;
		cObj.doMinimization = doMinimization;
		cObj.minimizeBB = minimizeBB;
		cObj.doBackrubs = doBackrubs;
		cObj.backrubFile = backrubFile;
		cObj.calculateVolumes = calculateVolumes;
		cObj.distDepDielect = distDepDielect;
		cObj.dielectConst = dielectConst;
		cObj.doDihedE = doDihedE;
		cObj.doSolvationE = doSolvationE;
		cObj.solvScale = solvScale;
		cObj.vdwMult = vdwMult;
		cObj.PEMcomp = PEMcomp;
		cObj.mutableSpots = mutableSpots;
		cObj.mutRes2Strand = mutRes2Strand;
		cObj.mutRes2StrandMutIndex = mutRes2StrandMutIndex;
		cObj.strandPresent = strandPresent;
		cObj.curMut = curMutIndex;
		cObj.strandLimits = strandLimits;
		cObj.strandsPresent = strandsPresent;
		cObj.templateAlwaysOn = templateAlwaysOn;
		cObj.doDih = doDih;
		cObj.neighborList = neighborList;
		cObj.distCutoff = distCutoff;
		cObj.pdbOutDir = pdbOutDir;

		cObj.useFlagsAStar = useFlagsAStar;
		cObj.useTriples = useTriples;
		cObj.doPerturbations = doPerturbations;
		cObj.magicBulletTriples = magicBulletTriples;
		cObj.magicBulletNumTriples = magicBulletNumTriples;
		if(doPerturbations){
			cObj.minimizePerts = minimizePerts;
			cObj.pertFile = pertFile;
			cObj.addWTRot = addWTRot;
			cObj.idealizeSC = idealizeSC;
		}

		cObj.useCCD = useCCD;
		cObj.minimizePairwise = minimizePairwise;
		cObj.es = es;

		cObj.EConvTol = EConvTol;

		if (PEMcomp) {//PEM computation

			if (!entropyComp){
				//Only send the molecule if it is the first job to the processor
//				if(curMutIndex < KSParser.numProc-1)
					cObj.m = m;
//				else
//					cObj.m = null;

				cObj.flagMutType = mutArray[curMutIndex].flagMutType;
				cObj.compCETM = compCETM;
				if(compCETM){
					if(cObj.flagMutType.equals("AS-AS"))
						cObj.emat = Emat.dualPosMat(pairEMatrixMin, false, mutArray[curMutIndex].runParams.pos1,mutArray[curMutIndex].runParams.pos2);
					else{
						int pos1 = -1;
						if(mutArray[curMutIndex].runParams != null)
							pos1 = mutArray[curMutIndex].runParams.pos1;
						cObj.emat = new Emat(pairEMatrixMin, pos1);
					}
				}else{
					if(cObj.flagMutType.equals("AS-AS"))
						cObj.emat = new Emat(pairEMatrixMin, false, mutArray[curMutIndex].runParams);
					else
						cObj.emat = new Emat(pairEMatrixMin,true,mutArray[curMutIndex].runParams);
				
				}
				
				cObj.curStrForMatrix = curStrForMatrix;
				cObj.flagMutType = mutArray[curMutIndex].flagMutType;
				cObj.curMut = mutArray[curMutIndex].mutNum;
				cObj.resMut = new int[pairEMatrixMin.singles.E.length];
				cObj.runParams = mutArray[curMutIndex].runParams;
				for(int i=0;i<cObj.resMut.length;i++) {
					cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
				
				}

			}
			else { //entropy energy computation run
				cObj.entropyComp = entropyComp;
				if (compASdist){ //AS-AS distance computation
					cObj.compASdist = compASdist;
					cObj.asDist = new boolean[mutableSpots];
					cObj.dist = dist;
				}
				else { //AS-AS or SHL-AS energy matrix computation
					//TODO: fix this code so it makes sense for multiple strands
					assert false==true;
					cObj.flagMutType = mutArray[curMutIndex].flagMutType;
					cObj.resMut = new int[pairEMatrixMin.singles.E.length];
					cObj.runParams = mutArray[curMutIndex].runParams;
					for(int i=0;i<cObj.resMut.length;i++) {
						cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
					
					}
				}
			}
		}			
		else {//mutation search

			if ((!distrDEE)&&(!distrDACS)){ //mutation search run, not (distributed DACS or distributed DEE)

				//cObj.q_L = q_L;
				cObj.pdbName = pdbName;
				cObj.numMutations = numMutations;
				cObj.repeatSearch = repeatSearch;
				cObj.algOption = algOption;
				cObj.numSplits = numSplits;
				cObj.initEw = initEw;
				cObj.scaleInt = scaleInt;
				cObj.maxIntScale = maxIntScale;
				cObj.pruningE = pruningE;
				cObj.stericE = stericE;
				cObj.gamma = gamma;
				cObj.epsilon = epsilon;
				cObj.useMaxKSconfs = useMaxKSconfs;
				cObj.numKSconfs = numKSconfs;
				cObj.asMethod = asMethod;
				cObj.m = m;
				cObj.deeSettings = deeSettings;

				cObj.currentMutation = new int[mutableSpots];
				for(int i=0;i<mutableSpots;i++) {
					cObj.currentMutation[i] = mutArray[curMutIndex].resTypes[i];
				}
				
				if(mutArray[curMutIndex].duplicateMut!=null){
					cObj.duplicateMut = new int[mutArray[curMutIndex].duplicateMut.length];
					//cObj.computedPartFuns = new PartitionMessage[mutArray[curMutIndex].duplicateMut.length];
					for(int i=0;i<cObj.duplicateMut.length;i++){
						cObj.duplicateMut[i] = mutArray[curMutIndex].duplicateMut[i];
					}
					completedRuns[curMutIndex] = new PartitionMessage[mutArray[curMutIndex].duplicateMut.length];
					
					//Check if the duplicate run has been completed
					for(int i=0;i<cObj.duplicateMut.length;i++){
						if(cObj.duplicateMut[i]>=0){
							if(completedRuns[cObj.duplicateMut[i]] != null && completedRuns[cObj.duplicateMut[i]][i]!=null){
								if(cObj.computedPartFuns == null)
									cObj.computedPartFuns = new PartitionMessage[cObj.duplicateMut.length]; 
								cObj.computedPartFuns[i] = new PartitionMessage(completedRuns[cObj.duplicateMut[i]][i]);
							}
						}
					}
				}
				
				cObj.bestScore = bestScore;
			}
			else {//distributed DACS or distributed DEE

				cObj.initEw = initEw;
				cObj.scaleInt = scaleInt;
				cObj.maxIntScale = maxIntScale;
				cObj.useSF = useSF;
				cObj.useEref = useEref;
				cObj.curStrForMatrix = curStrForMatrix;

				if (distrDACS){ //distributed DACS
					cObj.numMutations = numMutations;
					cObj.rotFileIn = rotFile;
					cObj.pruningE = pruningE;
					cObj.approxMinGMEC = approxMinGMEC;
					cObj.lambda = lambda;
					cObj.algOption = algOption;
					cObj.initDepth = initDepth;
					cObj.subDepth = subDepth;
					cObj.diffFact = diffFact;
					cObj.minRatioDiff = minRatioDiff;
					cObj.msp = msp;
					cObj.numInitUnprunedConf = numInitUnprunedConf;
					cObj.currentMutation = new int[mutableSpots];
					cObj.outputPruneInfo = outputPruneInfo;
					cObj.outputConfInfo = outputConfInfo;
					cObj.partIndex = new Index3[initDepth];
					cObj.asMethod = asMethod;
					for (int i=0; i<initDepth; i++)
						cObj.partIndex[i] = mutArray[curMutIndex].index[i];
					int ctr = 0;
					for(int i=0;i<strandDefault.length;i++)
						for(int j=0;j<strandDefault[i].length;j++){
							//KER: DACS Not Working right now
							//TODO: FIX DACS!!!
//							cObj.currentMutation[ctr] = strandDefault[i][j];
							ctr++;
						}
					cObj.bestScore = bestScore;
				}
				else { //distributed DEE
					//DistrDEE Specific
					int firstPos = -1;
					int secondPos = -1;
					
					cObj.resMut = new int[mutArray[curMutIndex].resMut.length];
					for (int i=0; i<cObj.resMut.length; i++){
						cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
						if(cObj.resMut[i] == 1 && firstPos == -1)
							firstPos = i;
						else if(cObj.resMut[i] == 1 && secondPos == -1)
							secondPos = i;
					}
					
					if(deeMethod.equals(KSParser.DEEMETHOD.GOLDSTEIN))
						cObj.emat = new Emat(pairEMatrixMin, firstPos);
					else{
//						File ematDir = new File("PEM");
//						if(ematDir != null && !ematDir.exists()){
//							ematDir.mkdirs();
//						}
						cObj.emat = Emat.dualPosMat(pairEMatrixMin, false, firstPos,secondPos);
//						emat.save("PEM/PEM_distr_COM+.dat",aaRotLib);
					}
					
					cObj.initEw = initEw;
					cObj.scaleInt = scaleInt;
					cObj.maxIntScale = maxIntScale;
					cObj.useSF = useSF;
					cObj.Ival = Ival;
					cObj.deeMethod = deeMethod;
					
					if(curMutIndex < KSParser.numProc-1)
						cObj.loadEmat = true;
					else
						cObj.loadEmat = false;
					
					
					cObj.pairStartEnd = mutArray[curMutIndex].pairStartEnd;

					cObj.numSpPos = numSpPos;
					cObj.typeDEE = typeDEE;
				}
			}
		}

		cObj.mutationNumber = curMutIndex;
		curMutIndex++;

		return(cObj);
	}


	// Output a finished mutation to the results file
	public synchronized void processFinishedMutation(CommucObj cObj) {

		if (PEMcomp){ //energy matrix computation
			if(!entropyComp){
			if(cObj.flagMutType.equals("TEMPL")){
				if(compCETM)
					cetm.mergeIn(cObj.cetm,cObj.resMut);
				else
					pairEMatrixMin.setTemplMinE(cObj.compEE.get(0).minE);
			}
			else if(cObj.flagMutType.equals("INTRA")){
				if(!compCETM){
					//KER: initialize eref to big E
					Set<String> keys = eRef.keySet();
					for(String pdbNum: keys){
						double[] eRefs = eRef.get(pdbNum);
						for(int j=0; j<eRefs.length;j++)
							eRefs[j]=Double.POSITIVE_INFINITY;
					}
					//KER: I now output a separate Eref matrix
					for(EMatrixEntrySlim re: cObj.compEE){
						RotamerEntry eme = (RotamerEntry)pairEMatrixMin.singles.getTerm(re.index);
						Residue mutRes = m.residue[pairEMatrixMin.resByPos.get(re.index[0]).get(0)];
						int aaInd = m.strand[mutRes.strandNumber].rcl.getRC(eme.r.rotamers[0]).rot.aaType.index;
						eRef.get(mutRes.getResNumberString())[aaInd] = Math.min(eRef.get(mutRes.getResNumberString())[aaInd],re.minE);	
					}
					//KER: Remove all of the bigE
					for(String pdbNum: keys){
						double[] eRefs = eRef.get(pdbNum);
						for(int j=0; j<eRefs.length;j++)
							if(eRef.get(pdbNum)[j]==Double.POSITIVE_INFINITY)
								eRef.get(pdbNum)[j] = 0.0f;
					}
					pairEMatrixMin.eRef = eRef;//outputObject(eRef,eRefMatrix+".dat");
				}
			}
			else{
				if(compCETM)
					cetm.mergeIn(cObj.cetm,cObj.resMut);
				else{
					for(EMatrixEntrySlim re : cObj.compEE){
						if(cObj.flagMutType.equals("SHL-AS")){
							pairEMatrixMin.singles.setE(re);
							if(re.rotDih1 != null){
								pairEMatrixMin.singles.setMaxE(re);
								pairEMatrixMin.singles.setDihed(re);
							}
						}
						else{
							pairEMatrixMin.pairs.setE(re);
							if(re.rotDih1 != null){
								pairEMatrixMin.pairs.setDihed(re);
								pairEMatrixMin.pairs.setMaxE(re);
							}
	
						}
					}
				}
			}
		}
		else { //entropy E matrix computation
			if (compASdist){ //AS-AS distance computation
				asDist[cObj.mutationNumber] = cObj.asDist;
			}
			else {
				//TODO: Fix entropy Ematrix
				System.out.println("Need to fix Entropy");
//				if (cObj.flagMutType.equalsIgnoreCase("INTRA")){
//					for (int i=0; i<countNewEntries; i++){
//						int index1 = 1 + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i2] + cObj.compEE[i].i3;
//						pairEMatrixMinEntropy[cObj.mutationNumber*rotLib.getTotalNumRotamers()+index1][0] = cObj.compEE[i].minE;
//					}
//				}
//				else { //AS-AS run
//					for (int i=0; i<countNewEntries; i++){
//						int ind1 = -1;
//						int ind2 = -1;
//						int index1 = cObj.compEE[i].i1*rotLib.getTotalNumRotamers() + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i2] + cObj.compEE[i].i3;
//						int index2 = cObj.compEE[i].i4*rotLib.getTotalNumRotamers() + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i5] + cObj.compEE[i].i6;
//						if (index1<rotLib.getTotalNumRotamers()){
//							ind1 = index1;
//							ind2 = index2-rotLib.getTotalNumRotamers();
//						}
//						else {
//							ind1 = index2;
//							ind2 = index1-rotLib.getTotalNumRotamers();
//						}
//						asasE[cObj.strandMut[0][0]][cObj.strandMut[0][1]][ind1][ind2] = cObj.compEE[i].minE;
//					}
//				}
			}

		}	
		//Output mutation information to results file (for resume)
		/*System.out.println("MutNUM: "+cObj.curMut+" produced "+countNewEntries+" new entries.");
			logPS.print("Completed mutation "+cObj.curMut);
				logPS.print(" SlaveNum "+cObj.slaveNum);
				logPS.print(" Time "+(cObj.elapsedTime/60.0));
				if (!entropyComp){
					for(int i=0;i<cObj.mutableSpots;i++)
						logPS.print(" "+cObj.resMut[i]);
					logPS.print(" "+cObj.flagMutType);
				}
				logPS.println();
				logPS.flush();*/
	}			
	else { //mutation search
		if ((!distrDEE)&&(!distrDACS)){ //Hybrid MinDEE-K*, not (distributed DACS or distributed DEE)
			System.out.println("MutNUM: "+cObj.mutationNumber);
			logPS.print("Completed mutation "+cObj.mutationNumber);
			BigDecimal score = new BigDecimal("0.0");
			BigDecimal denom = new BigDecimal("1.0");

			ExpFunction e = new ExpFunction();

			for(int i=0; i<cObj.numComplexes-1;i++){
				denom = denom.multiply(cObj.q[i],ExpFunction.mc);
			}
			if (denom.compareTo(new BigDecimal("0.0")) != 0)
				score = cObj.q[cObj.numComplexes-1].divide(denom,ExpFunction.mc);
			logPS.print(" Score "+score);

			logPS.print(" Volume "+mutArray[cObj.mutationNumber].vol);
			logPS.print(" SlaveNum "+cObj.slaveNum);
			logPS.print(" Time ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print((cObj.q_Time[i]/60.0)+" ");
			logPS.print(" InitBest "+cObj.bestScore);
			BigDecimal bs = cObj.bestScore;
			if (score.compareTo(cObj.bestScore) >0)
				bs = score;
			logPS.print(" FinalBest "+bs);
			for(int i=0;i<cObj.mutableSpots;i++)
				logPS.print(" "+cObj.currentMutation[i]);

			for(int i=0;i<cObj.numComplexes;i++){
				logPS.print(" "+i+"ConfInfo "+cObj.searchNumConfsEvaluated[i]+" "+cObj.searchNumPrunedMinDEE[i]+" "
						+cObj.searchNumConfsPrunedByS[i]+" "+cObj.searchNumConfsLeft[i]);
			}
			logPS.print(" MinEMinimized ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.bestEMin[i]+" ");
			logPS.print(" MinEUnMinimized ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.bestE[i]+" ");
			logPS.print(" EffectiveEpsilon: ");
			/*for(int i=0;i<cObj.numComplexes;i++)
					logPS.print(cObj.effEpsilon[i]+" ");*/
			logPS.print(" Partial_q_E ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.q[i]+" ");
			logPS.print(" E_total ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.searchNumConfsTotal[i]+" ");	
			logPS.print(" SecondEw ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.repeatEW[i]+" ");
			logPS.print(" E_allPruned ");
			for(int i=0;i<cObj.numComplexes;i++)
				logPS.print(cObj.allPruned[i]+" ");
			logPS.println();
			if (score.compareTo(bestScore) >0){
				logPS.println("BestScoreChange "+bestScore+" to "+score);
				bestScore = score;
			}
			logPS.flush();	
		}
		else if (distrDACS){ //distributed DACS
			bestScore = bestScore.min(cObj.bestScore);
			pruningE = bestScore.doubleValue();
			logPS.print("Completed mutation "+cObj.mutationNumber);
			logPS.print(" Score "+cObj.bestScore);
			logPS.print(" BestScore "+bestScore);
			logPS.print(" PartitionIndices");
			for (int i=0; i<initDepth; i++)
				logPS.print(" "+cObj.partIndex[i]);
			logPS.print(" Time "+(cObj.elapsedTime/60.0));
			logPS.println();
			logPS.flush();
			System.out.println("Partition "+cObj.mutationNumber+" done; best energy: "+cObj.bestScore);
		}
		else {//distributed DEE

			Emat emat = cObj.emat;
			
			int mutPos[] = {-1,-1};
			
			for(int i=0; i<cObj.resMut.length;i++){
				if(cObj.resMut[i] == 1){
					if(mutPos[0] == -1)
						mutPos[0] = i;
					else
						mutPos[1] = i;
				}
			}
			
			if(mutPos[1] == -1){//Singles
				Iterator<EMatrixEntryWIndex> iter = emat.singlesIterator(mutPos[0]);
				while(iter.hasNext()){
					EMatrixEntryWIndex single = iter.next();
					if(single.eme.isPruned()){
						pairEMatrixMin.setPruned(single.index, true);
					}
				}	
			}else{ //Pairs
				Iterator<EMatrixEntryWIndex> iter = emat.pairsIterator(mutPos[0], mutPos[1]);
				while(iter.hasNext()){
					EMatrixEntryWIndex pair = iter.next();
					if(pair.eme.isPruned()){
						pairEMatrixMin.setPruned(pair.index, true);
						pairEMatrixMin.setSymmetricPairPruned(pair.index, true);
					}
				}
			}
			}
//			else { //singles DEE
				//TODO:Fix singles DEE
				
//				Iterator<RotInfo<Boolean>> iter = prunedRot.iterator();
//				while(iter.hasNext()){
//					//for (int i=0; i<prunedRot.length; i++){
//					RotInfo<Boolean> ri = iter.next();
//					if (cObj.prunedRot.get(ri)) //get the pruned rotamers from this computation
//						prunedRot.set(ri, true);
//				}
//				outputObject(prunedRot,rotFile); //must be output after each result read
			}
		}

private synchronized Object readFromFile(Object inObj, String inFile){
	try{
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
		inObj = in.readObject();
		in.close();
	}
	catch (Exception e){
		System.out.println(e.toString());
		System.out.println("ERROR: An exception occurred while reading from object file");
		System.exit(0);
	}
	return inObj;
}

private void outputObject(Object outObj, String outFile){
	FileOutputStream fout = null;
	FileChannel ch = null;
	FileLock lock = null;
	try{
		fout = new FileOutputStream(outFile);
		ch = fout.getChannel();
		lock = ch.lock();
		ObjectOutputStream out = new ObjectOutputStream(fout);
		out.writeObject(outObj);
		//out.close();
	}
	catch (Exception e){
		System.out.println(e.toString());
		System.out.println("ERROR: An exception occurred while writing object file");
		System.exit(0);
	}
	finally {
		try {
			lock.release();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: unable to release lock on file "+outFile);
			System.exit(1);
		}
	}
}

private void deleteFile(String file) {

	Runtime rt = Runtime.getRuntime();
	try {
		File tmpSFMat = new File(file);
		if(tmpSFMat.exists())
			tmpSFMat.delete();
	}
	catch(Exception ex) {
		System.out.println("Exception: runtime");
		System.out.println(ex.getMessage());
	}	
}

public void closeLog() {
	if (logPS!=null){
		logPS.flush();
		logPS.close();
	}
}

public void setNumOfStrands(int s) {
	numberOfStrands = s;
}
public void setStrandMut(MutableResParams sm) {
	strandMut = sm;
}
public void setStrandDefault(String sd[][]) {
	strandDefault = sd;
}
public void setStrandPresent(boolean sp[]) {
	strandPresent = sp;
}
public void setStrandLimits(String sl[][]) {
	strandLimits = sl;
}
public void setStrandsPresent(int sp) {
	strandsPresent = sp;
}

/*public void setLigType(String lt) {
		ligType = lt;
	}*/
public void setLigPresent(boolean lp) {
	ligPresent = lp;
}
public void setNumMutations(int nm){
	numMutations = nm;
}
public void setResMutatable (int resMut[][]){
	resMutatable = new boolean[resMut.length][resMut[0].length];
	for (int i=0; i<resMutatable.length; i++){
		for (int j=0; j<resMutatable[0].length; j++){
			if (resMut[i][j]==1)
				resMutatable[i][j] = true;
			else
				resMutatable[i][j] = false;
		}
	}
}
public void setAAallowed(String aal[][]){
	AAallowed = aal;
}
public void setarpFilenameMin(String afnm) {
	arpFilenameMin = afnm;
}
public void setarpFilenameMax(String afnm) {
	arpFilenameMax = afnm;
}
public void setAlgOption(int ao){
	algOption = ao;
}
public void setNumSplits(int ns){
	numSplits = ns;
}
public void setMinDEEFileName(String mdf) {
	minDEEfile = mdf;
}
public void setInitEw(double iew){
	initEw = iew;
}
public void setPruningE(double pe){
	pruningE = pe;
}
public double getPruningE(){
	return pruningE;
}
public void setGamma(double g) {
	gamma = g;
}
public void setEpsilon(double g) {
	epsilon = g;
}
/*public void setEpsilon(double g) {
		epsilon = (new Double(g)).doubleValue();
	}*/
public void setParams(ParamSet theParams) {
	sParams = theParams;
}
public void setStericThresh(double st) {
	stericThresh = st;
}
public void setSoftStericThresh(double st){
	softStericThresh = st;
}
/*public void setNumInAS(int nas) {
		numInAS = nas;
	}*/

public void setComputeEVEnergy(boolean ceve) {
	computeEVEnergy = ceve;
}
public void setDoMinimization(boolean dm) {
	doMinimization = dm;
}
public void setMinimizeBB(boolean mbb){
	minimizeBB = mbb;
}
public void setDoBackrubs(boolean br){
	doBackrubs = br;
}
public void setBackrubFile(String brf){
	backrubFile = brf;
}
public void setRepeatSearch(boolean rs){
	repeatSearch = rs;
}
public void setCalculateVolumes(boolean cv) {
	calculateVolumes = cv;
}
/*public void setnumLigRotamers(int nlr) {
		numLigRotamers = nlr;
	}*/
public HashMap<String,double[]> getErefMatrix(){
	return eRef;
}

public void setPairEMatrixMin(Emat emat){
	pairEMatrixMin = emat;
}

public void setRotFile(String rf){
	rotFile = rf;
}
public void setUseSF(boolean usf){
	useSF = usf;
}
//public void setSpFlags(boolean spFlags[][][][][][]){
//	splitFlags = spFlags;
//}
//public void setSfFile(String sff){
//	sfFile = sff;
//}
public void setDistrDACS(boolean dDACS){
	distrDACS = dDACS;
}
public boolean getDistrDACS(){
	return distrDACS;
}
public void setDistrDEE(boolean dDEE){
	distrDEE = dDEE;
}
public void setBestScore(BigDecimal bs){
	bestScore = bs;
}
public void setNumSpPos(int spp){
	numSpPos = spp;
}
public void setMSP(int m[]){
	msp = m;
}
public void setTypeDEE(int t){
	typeDEE = t;
}
public void setInitDepth(int id){
	initDepth = id;
}
public void setSubDepth(int sd){
	subDepth = sd;
}
public void setDiffFact(int df){
	diffFact = df;
}
public void setMinRatioDiff(double mrd){
	minRatioDiff = mrd;
}
public void setNumInitUnprunedConf(BigInteger niuc){
	numInitUnprunedConf = niuc;
}
public void setOutputPruneInfo(String opi){
	outputPruneInfo = opi;
}
public void setOutputConfInfo(String oci){
	outputConfInfo = oci;
}
public void setApproxMinGMEC(boolean amg){
	approxMinGMEC = amg;
}
public void setLambda(double l){
	lambda = l;
}
public void setStericE(double se){
	stericE = se;
}
public void setDistDepDielect(boolean ddd){
	distDepDielect = ddd;
}
public void setDielectConst(double dc){
	dielectConst = dc;
}
public void setDoDihedE(boolean dde){
	doDihedE = dde;
}
public void setDoSolvationE(boolean dse){
	doSolvationE = dse;
}
public void setSolvScale(double ss){
	solvScale = ss;
}
public void setVdwMult(double vm){
	vdwMult = vm;
}
public void setScaleInt(boolean si){
	scaleInt = si;
}
public void setMaxIntScale(double is){
	maxIntScale = is;
}
public void setUseEref(boolean uer){
	useEref = uer;
}
public void setLigPartFn(BigDecimal ql){
	q_L = ql;
}
public void setEntropyComp(boolean ec){
	entropyComp = ec;
}
public Emat getMinEmatrix(){
	return pairEMatrixMin;
}
public void setPairEntropyMatrix(double aae[][][][]){
	asasE = aae;
}
public double [][][][] getPairEntropyEmatrix(){
	return asasE;
}
public void setIntraEntropyMatrixMin(double pemMin[][]){
	pairEMatrixMinEntropy = pemMin;
}
public void setASdistMatrix(boolean ad[][]){
	asDist = ad;
}
public void setASdist(double d){
	dist = d;
}
public void setCompASdist(boolean ad){
	compASdist = ad;
}
public boolean [][] getASdistMatrix(){
	return asDist;
}
public double [][] getMinEmatrixEntropy(){
	return pairEMatrixMinEntropy;
}
public void setMutableSpots(int ms) {
	mutableSpots = ms;
}

public void setTypeDep(boolean tD) {
	// TODO Auto-generated method stub
	typeDep = tD;
}
public void setRotamerLibrary(RotamerLibrary rl){
	rotLib = rl;
}
public void setErefMatrixName(String erm) {
	eRefMatrix  = erm;
}

public void setTemplateAlwaysOn(boolean templAlwaysOn) {
	templateAlwaysOn = templAlwaysOn;
}

public void setAddOrigRots(boolean aor) {
	addOrigRots   = aor;
}

public void setSaveTopConfs(boolean stc) {
	// TODO Auto-generated method stub
	saveTopConfs = stc;
}

public void setPrintTopConfs(boolean ptc) {
	// TODO Auto-generated method stub
	printTopConfs = ptc;
}

public void setNumTopConfs(int ntc) {
	// TODO Auto-generated method stub
	numTopConfs = ntc;
}

public void setUseMaxKSconfs(boolean useKSconfs) {
	useMaxKSconfs = useKSconfs;
}

public void setNumKSconfs(BigInteger nkc) {
	numKSconfs = nkc;
}

public void setcurStrForMatrix(int cSFM) {
	curStrForMatrix = cSFM;
}

public void setUseFlagsAStar(boolean useFlagsAStar) {
	this.useFlagsAStar = useFlagsAStar;
}
public void setUseTriples(boolean useTriples) {
	this.useTriples = useTriples;
}
public void setDoPerturbations(boolean dp){
	doPerturbations=dp;
}
public void setMinimizePerts(boolean mp){
	minimizePerts=mp;
}
public void setPertFile(String pf){
	pertFile = pf;
}
public void setAddWTRot(boolean awr){
	addWTRot = awr;
}
public void setIdealizeSC(boolean idealizeSC) {
	this.idealizeSC = idealizeSC;
}
public void setMagicBulletNumTriples(int magicBulletNumTriples) {
	this.magicBulletNumTriples = magicBulletNumTriples;
}
public void setMagicBulletTriples(boolean magicBulletTriples) {
	this.magicBulletTriples = magicBulletTriples;
}

public void setMinimizePairwise(boolean minimizePairwise) {
	this.minimizePairwise = minimizePairwise;
}

public void setUseCCD(boolean useCCD) {
	this.useCCD = useCCD;
}

public void setEConvTol(double EConvTol) {
	this.EConvTol = EConvTol;
}

public void setCompCETM(boolean compCETM) {
	this.compCETM = compCETM;
}

public void setCETM(CETMatrix cetm) {
	this.cetm = cetm;
}

public void setES(EPICSettings es) {
	this.es = es;
}

public void setMolecule(Molecule mol){
	this.m = mol;
}

public void setIval(double ival) {
	this.Ival = ival;
}

public void setDEEMethod(KSParser.DEEMETHOD deeMethod) {
	this.deeMethod = deeMethod;
}

public void setDoDih(boolean doDih2) {
	doDih = doDih2;
}

public void setNeighborList(boolean nl) {
	this.neighborList= nl;
	
}

public void setDistCutoff(double d) {
	this.distCutoff = d;
	
}

public void setErefMatrix(HashMap<String,double[]> eRef){
		this.eRef = eRef;
}

public void setASMethod(KSParser.ASTARMETHOD asMethod) {
	this.asMethod = asMethod;
}

public void setDEEsettings(DEEsettings deeSettings) {
	this.deeSettings = deeSettings;	
}

}
