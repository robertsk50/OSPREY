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
//	DEESplit1f.java
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

import java.util.Iterator;
import java.util.Random;

/**
 * Performs full split-DEE (conformational splitting) with 1 plit position
 * 
 */
public class DEESplit2fPierceNoIter {

	//two pairwise energy matrices: one for the min energies and one for the max
	private Emat pairwiseMinEnergyMatrix = null;

	boolean typeDependent = false;

	//this value depends on the particular value specified in the pairwise energy matrices;
	//		in KSParser, this value is 10^38;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private double bigE = Double.POSITIVE_INFINITY;

	//steric energy that determines incompatibility of a rotamer with the template
	double stericE = bigE;

	//size of the pairwise energy matrix
	private int PEMsize;

	private double curEw = 0.0f;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)

	//the minimum difference in the checkSum when a rotamer cannot be pruned
	private double minDiff = Double.NEGATIVE_INFINITY;

	int numRotForRes[] = null;

	//the number of split positions
	int numSplits = 0;

	//split flags for all rotamer pairs
	boolean splitFlags[][][][][][] = null;

	//determines if split flags are used
	boolean useFlags = false;

	//the number of runs
	int numRuns = 1;

	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
	boolean doMinimize = false;

	// 2010: iMinDEEboolean 
	boolean doIMinDEE = false;
	double Ival = 0.0;
	double initEw = 0.0;

	//determines which residue is checked (only for the distributed DEE)
	boolean resInMut[] = null;

	//determines if distributed DEE is performed
	boolean distrDEE = false;

	//the single and pair interval terms in the MinDEE criterion
//	double indIntMinDEE[] = null;
//	double pairIntMinDEE[] = null;

	//determines if backbone minimization is performed
	boolean minimizeBB = false;

	//the template interval energy (0.0 if fixed backbone)
	double templateInt = 0.0f;


	boolean getFalseMax = false;
	//the current ligand amino acid index
	//int ligAANum = -1;

	private double[][] minDiffMat;

	private int[] singleStartEnd;

	//private int numMutable;
	//StrandRotamers strandRot[] = null;
	//int strandMut[][] = null;
	//int mutRes2Strand[] = null;
	//int mutRes2MutIndex[] = null;

	//constructor
	DEESplit2fPierceNoIter(Emat arpMatrix, double initEw, boolean[] residueMut,
			boolean doMin, boolean useSF, boolean dDEE, boolean minBB,
			boolean typeDep, boolean aIMinDEE, double aIval, int[] singleStartEnd ) {

		doMinimize = doMin;
		typeDependent = typeDep;
		doIMinDEE = aIMinDEE;
		this.Ival = aIval;
		this.initEw = initEw;

		getFalseMax = !doMinimize || doIMinDEE;

		pairwiseMinEnergyMatrix = arpMatrix;

//		indIntMinDEE = indInt;
//		pairIntMinDEE = pairInt;

		useFlags = useSF;
		resInMut = residueMut;
		distrDEE = dDEE;
		minimizeBB = minBB;
		
		this.singleStartEnd = singleStartEnd;


		compNumRotForRes(); //compute the number of rotamers for each residue position

		curEw = initEw+Ival;

		numRuns = 1;

		templateInt = 0.0f;
		if (minimizeBB) //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)			
			templateInt = pairwiseMinEnergyMatrix.getTemplMaxE(doIMinDEE)-pairwiseMinEnergyMatrix.getTemplMinE();//pairwiseMaxEnergyMatrix[pairwiseMaxEnergyMatrix.length-1][0][0][0][0][0] - pairwiseMinEnergyMatrix[pairwiseMinEnergyMatrix.length-1][0][0][0][0][0];
	}

	//return the split flags for all rotamer pairs
	public boolean[][][][][][] getSplitFlags(){
		return splitFlags;
	}

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding r at i can be eliminated, and false otherwise
	public void ComputeEliminatedRotConf(){

		int numRotForCurAAatPos;

		int prunedCurRun = 0;
		boolean done = false;
		numRuns = 1;

		while (!done){

			prunedCurRun = 0;

//			System.out.println("Current run: "+numRuns);

			//Compute for the AS residues first
			for (int i=0; i<pairwiseMinEnergyMatrix.numMutPos(); i++){ //For each position i

				if(distrDEE && !resInMut[i]) //Skip if distrDEE and not current position
					continue;
				
				int singleCtr = -1;

				for(int i_r_aa=0;i_r_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_r_aa++){
					for(int i_r_rot=0;i_r_rot < pairwiseMinEnergyMatrix.singles.E[i][i_r_aa].length;i_r_rot++){

					singleCtr++;
					if (pairwiseMinEnergyMatrix.singles.pruned[i][i_r_aa][i_r_rot])//Skip if pruned
						continue;

					
					if(distrDEE && (singleCtr < singleStartEnd[0] || singleCtr >= singleStartEnd[1]) )
						continue;
					
					setupDiffMatrix(i,i_r_aa,i_r_rot);
					int[] mbSplitPos = getMBSplitPos(i,i_r_aa,i_r_rot);

					
					
					if (CanEliminateUsing(i,i_r_aa,i_r_rot,mbSplitPos[0],mbSplitPos[1])){
						pairwiseMinEnergyMatrix.singles.pruned[i][i_r_aa][i_r_rot] = true;
						prunedCurRun++;
					}
					
				}}


			}


//			System.out.println("Split2f RotPruned:"+prunedCurRun);


			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
			
			//Don't repeat a SplitFlags2 run
			done=true; 
		}

		//return eliminatedRotAtPos;
	}

	private int[] getMBSplitPos(int i, int i_r_aa, int i_r_rot) {

		double minDiff[] = new double[pairwiseMinEnergyMatrix.numMutPos()];

		//Find min_k min_t min_v [E(i_r,k_v) - E(i_t,k_v)
		for(int k=0;k<pairwiseMinEnergyMatrix.numMutPos();k++){
			if(k == i || !pairwiseMinEnergyMatrix.areNeighbors(i, k)) // k != i
				continue;

			double minInteraction = Double.POSITIVE_INFINITY; 
			
			int i_t_ctr = -1; 
			//For each candidate rotamer t at i
			for(int i_t_aa=0;i_t_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_t_aa++){
				for(int i_t_rot=0;i_t_rot < pairwiseMinEnergyMatrix.singles.E[i][i_t_aa].length;i_t_rot++){
				i_t_ctr++;
				
				if(pairwiseMinEnergyMatrix.singles.pruned[i][i_t_aa][i_t_rot] || 
						(i_r_aa == i_t_aa && i_r_rot == i_t_rot)) //skip if pruned or the same as i_r
					continue;

				for(int k_v_aa=0;k_v_aa < pairwiseMinEnergyMatrix.singles.E[k].length;k_v_aa++){
					for(int k_v_rot=0;k_v_rot < pairwiseMinEnergyMatrix.singles.E[k][k_v_aa].length;k_v_rot++){
				
					if(pairwiseMinEnergyMatrix.singles.pruned[k][k_v_aa][k_v_rot]) //skip if pruned
						continue;
					if(pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot])
						continue;
					
					//If i_t k_v are pruned i_r might still be better than it

					double diff = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot] - 
									pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][k][k_v_aa][k_v_rot];

					if(diff < minInteraction)
						minInteraction = diff;

				}}
			}}

			minDiff[k] = minInteraction;

		}

		double[] mins = {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY};
		int[] minInds = {-1,-1};

		for(int curPos=0; curPos<minDiff.length;curPos++){
			if(curPos == i || !pairwiseMinEnergyMatrix.areNeighbors(i, curPos))
				continue;
			if(minDiff[curPos] < mins[0]){
				minInds[0] = curPos;
				mins[0] = minDiff[curPos];
			}
			else if(minDiff[curPos] < mins[1]){
				minInds[1] = curPos;
				mins[1] = minDiff[curPos];

			}
		}



		return minInds;

	}

	private void setupDiffMatrix(int i,int i_r_aa,int i_r_rot) {

		minDiffMat = new double[pairwiseMinEnergyMatrix.numMutPos()][pairwiseMinEnergyMatrix.numRot(i)];
		
		
		for(int j=0;j<minDiffMat.length;j++){
			for(int k=0;k<minDiffMat[j].length;k++)
				minDiffMat[j][k] = Double.POSITIVE_INFINITY;
		}

		//For each candidate rotamer t at i
		int i_t_ctr = -1; 
		for(int i_t_aa=0;i_t_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_t_aa++){
			for(int i_t_rot=0;i_t_rot < pairwiseMinEnergyMatrix.singles.E[i][i_t_aa].length;i_t_rot++){
			i_t_ctr++;
			
			if(pairwiseMinEnergyMatrix.singles.pruned[i][i_t_aa][i_t_rot]) //skip if pruned
				continue;


			for (int j=0; j<pairwiseMinEnergyMatrix.numMutPos(); j++){//For all other positions j, j!=i 
				if(j == i || !pairwiseMinEnergyMatrix.areNeighbors(i, j)) //j != i
					continue;

				boolean foundValidJ = false;
				//For all rotamers at j
				for(int j_u_aa=0;j_u_aa < pairwiseMinEnergyMatrix.singles.E[j].length;j_u_aa++){
					for(int j_u_rot=0;j_u_rot < pairwiseMinEnergyMatrix.singles.E[j][j_u_aa].length;j_u_rot++){
						if(pairwiseMinEnergyMatrix.singles.pruned[j][j_u_aa][j_u_rot]) //skip if pruned
							continue;
						if(pairwiseMinEnergyMatrix.getPairPruned(i,i_r_aa,i_r_rot,j,j_u_aa,j_u_rot)) //skip if i_r,j_u pruned
							continue;
						//if i_t,j_u pruned then i_r might be better than it so we can't prune
	
						foundValidJ = true;
						//store Yjt = min_u [E(i_r,j_u) - E(i_t,j_u)]
						minDiffMat[j][i_t_ctr] = Math.min(minDiffMat[j][i_t_ctr],pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][j][j_u_aa][j_u_rot]-pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][j][j_u_aa][j_u_rot]);
				}}
				if(!foundValidJ) //KER: If there wasn't a valid j_u then we can't use i_t to prune
					minDiffMat[j][i_t_ctr] = Double.NEGATIVE_INFINITY;
			}

		}}

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
	private boolean CanEliminateUsing(int i, int i_r_aa, int i_r_rot,int k1, int k2){

		if(k1 == -1 || k2 == -1)
			return false;
		
		double checkSum = 0;
		
		boolean partitionPruned[][] = new boolean[numRotForRes[k1]][numRotForRes[k2]]; //prune r at i for each partition

		for (int partP=0; partP<partitionPruned.length; partP++)
			for(int j=0; j<partitionPruned[partP].length;j++)
				partitionPruned[partP][j] = false;

		//For each competitor t at i
		int i_t_ctr = -1; 
		for(int i_t_aa=0;i_t_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_t_aa++){
			for(int i_t_rot=0;i_t_rot < pairwiseMinEnergyMatrix.singles.E[i][i_t_aa].length;i_t_rot++){
			i_t_ctr++;
			
			if(pairwiseMinEnergyMatrix.singles.pruned[i][i_t_aa][i_t_rot]) //skip if pruned
				continue;

			//If typeDependent we can set same AAs to have curEw of 0
			//Otherwise the curEw = initEw
			if(typeDependent && i_t_aa == i_r_aa)
				curEw = Ival;
			else{
				curEw = Ival+initEw;
			}
			
			if(i_t_aa == i_r_aa && i_t_rot == i_r_rot) //Don't prune against self
				continue;

			checkSum = pairwiseMinEnergyMatrix.singles.E[i][i_r_aa][i_r_rot] - pairwiseMinEnergyMatrix.singles.E[i][i_t_aa][i_t_rot];


			for (int j=0; j<pairwiseMinEnergyMatrix.numMutPos(); j++){//for all other positions j, j!=i, j!=k 
				if(j == k1 || j == k2 || j == i || !pairwiseMinEnergyMatrix.areNeighbors(i, j))
					continue;

				checkSum += minDiffMat[j][i_t_ctr];
			}

			int k1_v_ctr = -1;
			for(int k1_v_aa=0;k1_v_aa < pairwiseMinEnergyMatrix.singles.E[k1].length;k1_v_aa++){//For each partition v at k1
				for(int k1_v_rot=0;k1_v_rot < pairwiseMinEnergyMatrix.singles.E[k1][k1_v_aa].length;k1_v_rot++){
				k1_v_ctr++;
				
				double splitDiff1 = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][k1][k1_v_aa][k1_v_rot] - 
										pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][k1][k1_v_aa][k1_v_rot]; 
				
				int k2_w_ctr = -1;
				for(int k2_w_aa=0;k2_w_aa < pairwiseMinEnergyMatrix.singles.E[k2].length;k2_w_aa++){
					for(int k2_w_rot=0;k2_w_rot < pairwiseMinEnergyMatrix.singles.E[k2][k2_w_aa].length;k2_w_rot++){
					k2_w_ctr++;
				
					double splitDiff2 = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][k2][k2_w_aa][k2_w_rot] - 
											pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][k2][k2_w_aa][k2_w_rot]; 
					
					if( (checkSum+splitDiff1 + splitDiff2) > curEw)
						partitionPruned[k1_v_ctr][k2_w_ctr] = true;
					/*else if(pairwiseMinEnergyMatrix.getPairPruned(i_r.index, k1_v.index) ||
							pairwiseMinEnergyMatrix.getPairPruned(i_r.index, k2_v.index))
						partitionPruned[k1_v_ctr][k2_v_ctr] = true;*/
				}}

			}}
		}}

		boolean allPruned = true;
		for(int part1=0; part1<partitionPruned.length;part1++)
			for(int part2=0; part2<partitionPruned[part1].length;part2++)
				if(!partitionPruned[part1][part2]){
					allPruned = false;
					break;
				}

		if(allPruned)
			return true;





		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
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
	//	private boolean CanEliminate (EMatrixEntryWIndex i_r){
	//
	//		double checkSum = 0;
	//
	//		int posNum = i_r.pos1();
	//
	//
	//
	//
	//		//For each splitting position k
	//		for (int splitPos=0; splitPos<pairwiseMinEnergyMatrix.numMutPos(); splitPos++){ //for each splitting position
	//
	//			if(splitPos == posNum) //For each split pos not equal to curPos
	//				continue;
	//			boolean partitionPruned[] = new boolean[numRotForRes[splitPos]]; //prune r at i for each partition
	//			for (int partP=0; partP<partitionPruned.length; partP++)
	//				partitionPruned[partP] = false;
	//
	//			//For each competitor t at i
	//			SinglesIterator i_t_iter = pairwiseMinEnergyMatrix.singlesIterator(i_r.pos1());
	//			int i_t_ctr = -1; 
	//			while(i_t_iter.hasNext()){
	//				i_t_ctr++;
	//				EMatrixEntryWIndex i_t = i_t_iter.next();
	//				if(i_t.eme.isPruned()) //skip if pruned
	//					continue;
	//				
	//				if(i_t.aa1() == i_r.aa1() && i_t.rot1() == i_r.rot1()) //Don't prune against self
	//					continue;
	//
	//				checkSum = i_r.eme.minE() - i_t.eme.minE();
	//
	//
	//				for (int altPos=0; altPos<pairwiseMinEnergyMatrix.numMutPos(); altPos++){//for all other positions j, j!=i, j!=k 
	//					if(altPos == splitPos || altPos == posNum)
	//						continue;
	//
	//					checkSum += minDiffMat[altPos][i_t_ctr];
	//				}
	//
	//				SinglesIterator k_v_iter = pairwiseMinEnergyMatrix.singlesIterator(splitPos);
	//				int k_v_ctr = -1;
	//				while(k_v_iter.hasNext()){//For each partition v at k
	//					k_v_ctr++;
	//					EMatrixEntryWIndex k_v = k_v_iter.next();
	//					double splitDiff = pairwiseMinEnergyMatrix.getPairMinE(i_r.index, k_v.index) - 
	//					pairwiseMinEnergyMatrix.getPairMinE(i_t.index, k_v.index); 
	//					if( (checkSum+splitDiff) > curEw)
	//						partitionPruned[k_v_ctr] = true;
	//
	//				}
	//			}
	//
	//			boolean allPruned = true;
	//			for(int i=0; i<partitionPruned.length;i++)
	//				if(!partitionPruned[i]){
	//					allPruned = false;
	//					break;
	//				}
	//
	//			if(allPruned)
	//				return true;
	//
	//		}
	//
	//
	//
	//
	//		//We have tried all of the other rotamers at the current position and none
	//		//of them is able to prune the given rotamer, so we return false
	//		return false;
	//	}

	//Compute the number of rotamers for each residue position (assign to numRotForRes[])
	private void compNumRotForRes(){

		//boolean ligPresent = (numLigRot==0); //ligand present
		int treeLevels = pairwiseMinEnergyMatrix.numMutPos();
		/*if (ligPresent)
			treeLevels++;*/

		numRotForRes = new int[treeLevels];

		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			numRotForRes[curLevel] = pairwiseMinEnergyMatrix.numRot(curLevel);
		}
	}

	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////

	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (EMatrixEntryWIndex curRot, EMatrixEntryWIndex altRot, int splitPos){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			

			if ((curPos != curRot.pos1())&&(curPos!=splitPos)) // j!=i and j!=k

				sum += IndMinDiffPVE(curRot, altRot, curPos);
		}

		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			sum += LigandIndMinDiffPVE (atPos, withAA1, withRot1, withAA2, withRot2);
		}*/

		return sum;
	}

	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (EMatrixEntryWIndex firstRot1, EMatrixEntryWIndex firstRot2, int secondPos){

		double minE = bigE;
		double curEmin, curEmax;

		int index1, index2, index3;
		int numRotForAAatPos;

		//r at i
		//index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;

		//int str2=mutRes2Strand[secondPos];
		//int strResNum2=strandMut[str2][mutRes2MutIndex[secondPos]];


		//t at i
		//index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;

		boolean found = false;

		if (((!firstRot1.eme.isPruned()))&&
				((!firstRot2.eme.isPruned()))){ //not pruned 


			for(int secondAA=0; secondAA<pairwiseMinEnergyMatrix.singles.pruned[secondPos].length;secondAA++){
				//SinglesIterator iter = pairwiseMinEnergyMatrix.singlesIterator(secondPos);
				//while(iter.hasNext()){
				//	EMatrixEntryWIndex curRot = iter.next();
				//for (int AA=0; AA<numAAtypes[secondPos]; AA++){

				//int curAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA);

				//numRotForAAatPos = strandRot[str2].rl.getNumRotForAAtype(curAA);
				//if (numRotForAAatPos==0)	//ala or gly
				//	numRotForAAatPos = 1;
				for(int secondRot=0; secondRot<pairwiseMinEnergyMatrix.singles.pruned[secondPos][secondAA].length;secondRot++){
					//for (int curRot=0; curRot<numRotForAAatPos; curRot++){

					//s at j
					//index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;
					int[] secondAAindex = {secondPos,secondAA,secondRot};
					if ((!pairwiseMinEnergyMatrix.getSinglePruned(secondAAindex))){ //not pruned 
						//EMatrixEntryWIndex oneCur = new EMatrixEntryWIndex(pairwiseMinEnergyMatrix,firstRot1,curRot);
						//EMatrixEntryWIndex twoCur = new EMatrixEntryWIndex(pairwiseMinEnergyMatrix,firstRot2,curRot);
						if ((!useFlags)||(!pairwiseMinEnergyMatrix.getPairPruned(firstRot1.index, secondAAindex))){ //not using split flags or not flagged

							curEmin = pairwiseMinEnergyMatrix.getPairMinE(firstRot1.index, secondAAindex);//pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][secondPos][curAA][curRot];
							curEmax = pairwiseMinEnergyMatrix.getPairMinE(firstRot2.index, secondAAindex);//pairwiseMaxEnergyMatrix[firstPos][firstAA2][firstRot2][secondPos][curAA][curRot];
							//if (/*(curEmin<=stericEThreshPair)&&*/(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
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

		if (!found){ //no possible pairs found
			//minE = 0.0; //contributes nothing to the sum
			minE = Double.NEGATIVE_INFINITY; //KER: if there's not a valid pair then this cannot be pruned;
		}


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

		if ((!eliminatedRotAtPos[index_r])){ //not already pruned

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

					maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0]; //formula term 2
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
