import java.io.ObjectInputStream.GetField;
import java.util.Iterator;

/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.0
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
//	DEEGoldstein.java
//
//	Version:           2.0
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ivelin Georgiev (2004-2009)
 * 
 */

/**
 * Performs simple Goldstein DEE rotamer pruning
 * 
 */
public class DEEGoldsteinPierceNoIter {

	//two pairwise energy matrices: one for the min energies and one for the max
	private Emat pairwiseMinEnergyMatrix = null;
	//private float pairwiseMaxEnergyMatrix[][][][][][] = null;

	//eliminated rotamers at position i, for all positions
	//private PrunedRotamers<Boolean> eliminatedRotAtPos = null;

	//number of residues under consideration
	//private int numSiteResidues;

	//for each residue, number of possible amino acids
	//private int numTotalRot;

	//number of possible rotamers for the ligand
	//int numLigRot;

	//offset of the given rotamer in the total rotamer set (?152?)
	//int rotIndOffset[];

	//the number of AA types allowed for each AS residue
	//int numAAtypes[] = null;

	//number of rotamers for the current AA type at the given residue
	//int numRotForAAtypeAtRes[];

	//this value depends on the particular value specified in the pairwise energy matrices;
	//		in KSParser, this value is 10^38;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private float bigE = (float)Math.pow(10,38);

	//steric energy that determines incompatibility of a rotamer with the template
	float stericE = bigE;

	private double curEw = 0.0;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)
	

	//the minimum difference in the checkSum when a rotamer cannot be pruned
	private double minDiff = -(float)Math.pow(10,30);

	//the rotamer library
	//RotamerLibrary rl = null;

	//The system rotamer handler
	//StrandRotamers sysLR = null;

	//The mapping from AS position to actual residue numbers
	//int residueMap[] = null;

	//boolean for turning on type dependent DEE i.e. only rotamers of the same type can prune each other
	boolean typeDependent = false;

	//the number of runs
	int numRuns = 1;

	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
//	boolean doMinimize = false;

	//the single and pair interval terms in the MinDEE criterion
//	double indIntMinDEE[] = null;
//	double pairIntMinDEE[] = null;

	//split flags for all rotamer pairs
	//boolean splitFlags[][][][][][] = null;

	//determines if split flags are used
	boolean useFlags = false;

	//determines if backbone minimization is performed
//	boolean minimizeBB = false;

	//the template interval energy (0.0 if fixed backbone)
//	double templateInt = 0.0f;

	// 2010: iMinDEE
	boolean doIMinDEE = false;
	double Ival = 0.0;
	double initEw = 0.0;

//	boolean getFalseMax = true;

	private boolean distrDEE;
	private boolean[] resInMut;
	private int[] singleStartEnd;

	boolean removeRot = false;

	//the current ligand amino acid index
	//int ligAANum = -1;

	//private int numMutable;
	//MutableSpot strandMut = null;

	//constructor
	DEEGoldsteinPierceNoIter(Emat arpMatrix,
			double initEw, boolean useSF, 
			boolean typeDep, boolean iMinDEE, double Ival, 
			boolean distrDEE, boolean[] resInMut, int[] singleStartEnd, boolean removeRot) {

//		doMinimize = doMin;
		typeDependent = typeDep;
		doIMinDEE = iMinDEE;
		this.removeRot = removeRot; 
//		getFalseMax = !doMinimize || doIMinDEE;

		pairwiseMinEnergyMatrix = arpMatrix;
		// 2010: No max matrix if doIMinDEE set
		/*if (doMinimize && !doIMinDEE) //max matrix is different
			pairwiseMaxEnergyMatrix = arpMatrixMax;
		else //no minimization, so the same matrix // 2010: if doIMinDEE is set to true then it is the same as DEE
			pairwiseMaxEnergyMatrix = pairwiseMinEnergyMatrix;*/

		this.distrDEE = distrDEE;
		this.resInMut = resInMut;
		this.singleStartEnd = singleStartEnd;

		//splitFlags = spFlags;
		//eliminatedRotAtPos = prunedRotAtRes;
		//rotIndOffset = rotamerIndexOffset;		
		//residueMap = resMap;
//		indIntMinDEE = indInt;
//		pairIntMinDEE = pairInt;
		//sysLR = systemLRot;
		//rl = rlP;
		useFlags = useSF;
//		minimizeBB = minBB;


		//numMutable = numResMutable;

		//numSiteResidues = numResInActiveSite;		// tested with 9
		//numTotalRot = numTotalRotamers;				// ?152?
		//numLigRot = numLigRotamers;					// 0 if no ligand

		/*numAAtypes = new int[strandMut.numMutPos()];

		for(int i=0; i<strandMut.numMutPos();i++){
			numAAtypes[i] = strandMut.getNumAllowable(i);
		}*/

		this.Ival = Ival;
		this.initEw = initEw;
		curEw = initEw + Ival;

		numRuns = 1;

//		templateInt = 0.0f;
//		if (minimizeBB){ //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)
//			templateInt = pairwiseMinEnergyMatrix.getTemplMaxE(doIMinDEE)-pairwiseMinEnergyMatrix.getTemplMinE();//pairwiseMaxEnergyMatrix[pairwiseMaxEnergyMatrix.length-1][0][0][0][0][0] - pairwiseMinEnergyMatrix[pairwiseMinEnergyMatrix.length-1][0][0][0][0][0];
//		}
	}

	//return the split flags for all rotamer pairs
	/*public boolean[][][][][][] getSplitFlags(){
		return splitFlags;
	}*/

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding r at i can be eliminated, and false otherwise
	public void ComputeEliminatedRotConf(){

		int numRotForCurAAatPos;
		int prunedTotal = 0;
		int prunedCurRun = 0;
		boolean done = false;
		numRuns = 1;
		while (!done){

			prunedCurRun = 0;
			//Compute for the AS residues first
			for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){
				if(distrDEE && !resInMut[curPos])
					continue;

				int singleCtr = -1;
				for(int curAA=0; curAA<pairwiseMinEnergyMatrix.singles.E[curPos].length;curAA++){
					for(int curRot=0; curRot<pairwiseMinEnergyMatrix.singles.E[curPos][curAA].length;curRot++){
				
					singleCtr++;
					if (pairwiseMinEnergyMatrix.getSinglePruned(curPos, curAA, curRot))//skip if pruned
						continue;

					
					if(distrDEE && (singleCtr < singleStartEnd[0] || singleCtr >= singleStartEnd[1]) )
						continue;

					if (CanEliminate(curPos,curAA,curRot)){
						pairwiseMinEnergyMatrix.singles.pruned[curPos][curAA][curRot] = true;

						prunedCurRun++;
						prunedTotal++;
						//if(prunedCurRun % 20 == 0)
						//	System.out.print(".");
						
//						if(!distrDEE && prunedCurRun>1000 && removeRot){
//							//if we've pruned a bunch, then shrink the array (should save time)
//							pairwiseMinEnergyMatrix.removePrunedRotReducedMem(false);
//	
//							//Reset search
//							curPos = 0;
//							curAA = 0;
//							curRot = 0;
//							
//							prunedCurRun = 0;
//	
//						}
					}

				}}
			}

			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
		}

//		System.out.println("DEEGoldstein pruned rotamers: "+prunedTotal);

	}

	//Called only by ComputeEliminatedRotConf(.)
	/*
	 * The logic is as follows:
	 * 	for every AS residue
	 * 		for every AA type at this residue
	 * 			for every possible rotamer at the given AA type
	 * 				check against every other possible rotamer for all AA types at the same residue;
	 * 				eliminate if provable
	 * 
	 * That is, for each residue: out of the 152 rotamer possibilities, choose 1
	 * and compare it to the other 151 until elimination can be proven or there
	 * are no more comparisons left. Repeat this for all 152 rotamers and all AS residues
	 * 
	 */
	private boolean CanEliminate (int pos,int i_r_aa,int i_r_rot){


		double minDiffPairVoxelE;

		double checkSum;

		
		//For all i_t
		for(int i_t_aa=0; i_t_aa<pairwiseMinEnergyMatrix.singles.E[pos].length;i_t_aa++){
			for(int i_t_rot=0; i_t_rot<pairwiseMinEnergyMatrix.singles.E[pos][i_t_aa].length;i_t_rot++){
				
		
			if(pairwiseMinEnergyMatrix.getSinglePruned(pos, i_t_aa, i_t_rot)) //skip if pruned
				continue;

			//If typeDependent we can set same AAs to have curEw of 0
			//Otherwise the curEw = initEw
			if(typeDependent && i_t_aa == i_r_aa)
				curEw = Ival;
			else{
				curEw = Ival+initEw;
			}


			//Remove if statement and do one more calculation that should never return true
//			if (!((i_r.aa1()==i_t.aa1())&&(i_r.rot1()==i_t.rot1()))){

				checkSum = pairwiseMinEnergyMatrix.singles.E[pos][i_r_aa][i_r_rot] - pairwiseMinEnergyMatrix.singles.E[pos][i_t_aa][i_t_rot];//pairwiseMinEnergyMatrix.getIntraE(altRe.index).maxE(doIMinDEE);//  [posNum][altAA][altRot][posNum][0][0];		//formula term 2
//				checkSum += -templateInt;

				minDiffPairVoxelE = SumMinDiffPVE(pos, i_r_aa,i_r_rot, i_t_aa,i_t_rot);	//formula term 5

				checkSum += minDiffPairVoxelE;

				if (checkSum > curEw){
					return true;}//this rotamer can be pruned/eliminated
				else {
					minDiff = Math.max(minDiff,checkSum);
				}
//			}
		}}


		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}

	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////

	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (int atPos, int i_r_aa,int i_r_rot,int i_t_aa,int i_t_rot){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			

			if (curPos != atPos && pairwiseMinEnergyMatrix.areNeighbors(atPos,curPos))
				sum += IndMinDiffPVE(atPos, i_r_aa,i_r_rot, i_t_aa,i_t_rot, curPos);
		}

		return sum;
	}

	//Called by SumMaxMaxPVE(.)
	//min_s (E_\ominus(i_r,j_s)-E_\oplus(i_t,j_s))
	//secondPos = j
	private double IndMinDiffPVE (int i, int i_r_aa,int i_r_rot, int i_t_aa,int i_t_rot, int j){

		double minE = bigE;
		double curEmin, curEmax;

		boolean found = false;
		boolean validK = false;

		for (int j_s_aa=0; j_s_aa<pairwiseMinEnergyMatrix.singles.pruned[j].length;j_s_aa++){
			for(int j_s_rot=0; j_s_rot<pairwiseMinEnergyMatrix.singles.pruned[j][j_s_aa].length;j_s_rot++){
				//s at j
				//int[] curAAindex = {secondPos,curAA,curRot};
				if ((pairwiseMinEnergyMatrix.getSinglePruned(j,j_s_aa,j_s_rot))) //skip rotamer if it is pruned
					continue;

				validK = true;

				if ((pairwiseMinEnergyMatrix.getPairPruned(i,i_r_aa,i_r_rot, j,j_s_aa,j_s_rot)))// skip pair if pruned
					continue; 

				//if firstAA2 is pruned then firstAA1 still might be better

				found = true;

				curEmin = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot];
				curEmax = pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][j][j_s_aa][j_s_rot];

				if ((curEmin-curEmax) < minE)
					minE = curEmin-curEmax;

			}

		}


		if(validK && !found)	
			minE = Double.POSITIVE_INFINITY;
		else if(!found)
			minE = Double.NEGATIVE_INFINITY; //KER: if there's not a valid pair then this cannot be pruned;


		return minE;
	}

	//Called by SumMaxMaxPVE(.)
	/*private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2){

		double minE = bigE;
		double curEmin, curEmax;

		int index1, index2, index3;

		//r at i
		index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;

		//t at i
		index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;

		boolean found = false;

		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index3]))){ //not pruned 

			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){

				//s at j (the ligand residue)
				index2 = numSiteResidues*numTotalRot + curLigPos;

				if ((!eliminatedRotAtPos[index2])){ //not pruned 

					if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged

						curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos];
						curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAA2][firstRot2][numSiteResidues][ligAANum][curLigPos];
						//if (/*(curEmin<=stericEThreshPair)&&*//*(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
						//}
						found = true;
					}
				}
			}

			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}

		if (!found)
			minE = 0.0; //contributes nothing to the sum

		return minE;
	}*/
	//////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////
	//Same as CanEliminate(), just checks the ligand rotamers for pruning
	//Called by ComputeEliminatedRotConf()
	/*private boolean CanEliminateLig (int curLigRot){

		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;

		int index_r, index_t;

		double checkSum;

		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r = numSiteResidues*numTotalRot + curLigRot;

		if ((!eliminatedRotAtPos[index_r])){ //not pruned 

			minIndVoxelE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][0]; //formula term 1
			minShellResE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][1];

			if ((minIndVoxelE + minShellResE)>=stericE) //rotamer incompatible with template, so prune
				return true;

			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE[numSiteResidues];							//formula term 3
				pairVoxelInterval = pairIntMinDEE[numSiteResidues];							//formula term 4
			}
			else { //traditional-DEE, so no interval terms
				indVoxelInterval = 0.0;
				pairVoxelInterval = 0.0;
			}

			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning
			for (int altRot=0; altRot<numLigRot; altRot++){

				//if t and r are not actually the same lig rotamer
				if (curLigRot!=altRot){

					//at this point, we know what r at i and t at i are

					index_t = numSiteResidues*numTotalRot + altRot;

					maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0];		//formula term 2
					maxShellResE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][1];

					//if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
					if ((!eliminatedRotAtPos[index_t])){ //not pruned 

						minDiffPairVoxelE = SumMinDiffPVELig(curLigRot, altRot);	//formula term 5

						checkSum = -templateInt + (minIndVoxelE + minShellResE) - (maxIndVoxelE + maxShellResE)
									- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

						if (checkSum > curEw){
							//System.out.println(checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
							return true;}//this rotamer can be pruned/eliminated
						else {
							minDiff = Math.max(minDiff,checkSum);
						}
					}
				}
			}
		}
		else //aready pruned
			return true;

		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}*/

	//Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
	//Called by CanEliminateLig()
	/*private double SumMinDiffPVELig (int withRot1, int withRot2){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			

				sum += IndMinDiffPVELig(withRot1, withRot2, curPos);
		}

		return sum;
	}*/

	//Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
	//Called by SumMinDiffPVELig()
	/*private double IndMinDiffPVELig (int firstRot1, int firstRot2, int secondPos){

		double minE = bigE;
		double curEmin, curEmax;

		int index1, index2, index3;
		int numRotForAAatPos;

		//r at i
		index1 = numSiteResidues*numTotalRot + firstRot1;

		//t at i
		index3 = numSiteResidues*numTotalRot + firstRot2;

		boolean found = false;

		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index3]))){ //not pruned 

			for (int AA=0; AA<numAAtypes[secondPos]; AA++){

				int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);

				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;

				for (int curRot=0; curRot<numRotForAAatPos; curRot++){

					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1

					//s at j
					index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;

					if ((!eliminatedRotAtPos[index2])){ //not pruned 

						if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged

							curEmin = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][firstRot1][secondPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][firstRot2][secondPos][curAA][curRot];
							//if (/*(curEmin<=stericEThreshPair)&&*//*(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
								if ((curEmin-curEmax) < minE)
									minE = curEmin-curEmax;
							//}
							found = true;
						}
					}					
				}
			}

			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}

		if (!found) //no possible pairs found
			minE = 0.0; //contributes nothing to the sum

		return minE;
	}*/
}