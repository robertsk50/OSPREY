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
//	DEEGoldsteinPairs.java
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
 * Performs DEE Goldstein pairs rotamer pruning
 * 
 */
public class DEEGoldsteinPairsNoIter {

	//two pairwise energy matrices: one for the min energies and one for the max
	private Emat pairwiseMinEnergyMatrix = null;

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
	//		in KSParser, this value is 99999.0;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private double bigE = Double.POSITIVE_INFINITY;

	private double curEw = 0.0f;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)

	//the minimum difference in the checkSum when a rotamer cannot be pruned
	private double minDiff = Double.NEGATIVE_INFINITY;

	//int count = 0;

	//boolean for turning on type dependent DEE i.e. only rotamers of the same type can prune each other
	boolean typeDependent = false;

	//the rotamer library
	//RotamerLibrary rl = null;

	//The system rotamer handler
	//StrandRotamers sysLR = null;

	//The mapping from AS position to actual residue numbers
	//int residueMap[] = null;

	//the number of runs
	int numRuns = 1;

	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
	boolean doMinimize = false;

	//the single and pair interval terms in the DE pairs MinDEE criterion
	double indIntMinDEE2Pos[][] = null;
	double pairIntMinDEE2Pos[][] = null;

	//determines if magic bullet or full pairs is used
	boolean magicBullet = false;

	//split flags for all rotamer pairs
	//boolean splitFlags[][][][][][] = null;

	//determines if split flags are used
	boolean useFlags = false;

	//determines which two residues are in the current pair (only for the distributed DEE)
	boolean resInPair[] = null;

	//determines if distributed DEE is performed
	boolean distrDEE = false;

	//determines if backbone minimization is performed
	boolean minimizeBB = false;

	//the template interval energy (0.0 if fixed backbone)
	double templateInt = 0.0;

	// 2010: iMinDEEboolean 	
	boolean doIMinDEE = false;
	double Ival = 0.0;
	double initEw = 0.0;

	boolean getFalseMax = true;
	//the max scaling factor for the interval terms
	double maxScale = 1.0f;
	
	int pairStartEnd[];
	

	//the current ligand amino acid index
	//int ligAANum = -1;

	//private int numMutable;
	//MutableSpot strandMut = null;

	//constructor
	DEEGoldsteinPairsNoIter(Emat arpMatrix, double initEw,	boolean residueMut[],
			boolean doMin, boolean useSF, boolean mb, boolean dDEE, boolean minBB, boolean scaleInt, double maxSc, 
			boolean typeDep, boolean aIMinDEE, double aIval, int[] pairStartEnd) {

		doMinimize = doMin;
		typeDependent = typeDep;
		doIMinDEE = aIMinDEE;
		this.Ival = aIval;
		this.initEw = initEw;

		getFalseMax = !doMinimize || doIMinDEE;

		pairwiseMinEnergyMatrix = arpMatrix;

		// 2010: No max matrix if useMinDEEPruningEw set
		//if (doMinimize && !doIMinDEE) //max matrix is different
		//	pairwiseMaxEnergyMatrix = arpMatrixMax;
		//else //no minimization, so the same matrix // 2010: if useMinDEEPruningEw is set to true then it is the same as DEE
		//	pairwiseMaxEnergyMatrix = pairwiseMinEnergyMatrix;

		//residueMap = resMap;
		//sysLR = systemLRot;
		//rl = rlP;
		useFlags = useSF;
		resInPair = residueMut;
		distrDEE = dDEE;
		magicBullet = mb;
		minimizeBB = minBB;

		//strandMut = strMut;
		this.pairStartEnd = pairStartEnd;
		
		

		//numSiteResidues = numResInActiveSite;		// tested with 9

		//numLigRot = numLigRotamers;					// 0 if no ligand
		//if (numLigRot>0)
		//	ligAANum = ligROT.getIndexOfNthAllowable(0,0);

		/*numAAtypes = new int[strandMut.numMutPos()];

		for(int i=0; i<strandMut.numMutPos();i++){
			numAAtypes[i] = strandMut.getNumAllowable(i);
		}*/

		curEw = initEw + Ival;

		numRuns = 1;

		templateInt = 0.0f;
		if (minimizeBB){ //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)			
			templateInt = pairwiseMinEnergyMatrix.getTemplMaxE(doIMinDEE)-pairwiseMinEnergyMatrix.getTemplMinE();//pairwiseMaxEnergyMatrix[pairwiseMaxEnergyMatrix.length-1][0][0][0][0][0] - pairwiseMinEnergyMatrix[pairwiseMinEnergyMatrix.length-1][0][0][0][0][0];
		}

		if (!scaleInt) //no scaling of the interval terms performed
			maxScale = 1.0f;
		else
			maxScale = maxSc;
	}

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding rotamer pair can be eliminated, and false otherwise
	public void ComputeEliminatedRotConf(){

		// 2010: No interval terms for useMinDEEPruningEw
		if (doMinimize && !doIMinDEE){ //compute the DE pairs MinDEE interval terms

			int numRes = pairwiseMinEnergyMatrix.numMutPos();		
			

			indIntMinDEE2Pos = new double[numRes][numRes];
			pairIntMinDEE2Pos = new double[numRes][numRes];

			for (int posNum1=0; posNum1<numRes; posNum1++){
				for (int posNum2=posNum1+1; posNum2<numRes; posNum2++){

					indIntMinDEE2Pos[posNum1][posNum2] = SumMaxIndInt(posNum1,posNum2);				//formula term 3
					pairIntMinDEE2Pos[posNum1][posNum2] = SumSumMaxPairInt(posNum1,posNum2);		//formula term 4
					indIntMinDEE2Pos[posNum2][posNum1] = indIntMinDEE2Pos[posNum1][posNum2];		//formula term 3
					pairIntMinDEE2Pos[posNum2][posNum1] = pairIntMinDEE2Pos[posNum1][posNum2];		//formula term 4
				}
			}
		}

		//Check for pairs pruning
		int numRotForCurAAatPos1;

		int prunedCurRun = 0;
		boolean done = false;
		int[] mbEntry = null;
		numRuns = 1;
		while (!done){

			prunedCurRun = 0;

//			System.out.println("Current run: "+numRuns);

	
			for (int i=0; i<pairwiseMinEnergyMatrix.numMutPos(); i++){

				if ((!distrDEE)||(resInPair[i])){ //not distrDEE or cur res is in distr pair

					for (int j=i+1; j<pairwiseMinEnergyMatrix.numMutPos(); j++){
						if(!pairwiseMinEnergyMatrix.areNeighbors(i, j))
							continue;
						

						if ((!distrDEE)||(resInPair[j])){ //not distrDEE or cur res is in distr pair

							if(magicBullet){
								mbEntry = findMagicBullets(i, j);
								if(mbEntry == null)
									continue;
							}
							
							int pairCtr = -1;
							
							//For all rot at i, i_r
							for(int i_r_aa=0;i_r_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_r_aa++){
								for(int i_r_rot=0;i_r_rot < pairwiseMinEnergyMatrix.singles.E[i][i_r_aa].length;i_r_rot++){
								
								//For all rot at j, j_s
								for(int j_s_aa=0;j_s_aa < pairwiseMinEnergyMatrix.singles.E[j].length;j_s_aa++){
									for(int j_s_rot=0;j_s_rot < pairwiseMinEnergyMatrix.singles.E[j][j_s_aa].length;j_s_rot++){	
							

									pairCtr++;
									//I moved the i_r check down here to get the pairsCtr correct
									if (pairwiseMinEnergyMatrix.singles.pruned[i][i_r_aa][i_r_rot] || 
											pairwiseMinEnergyMatrix.singles.pruned[j][j_s_aa][j_s_rot])//not already pruned
										continue;

									
									if(distrDEE && (pairCtr < pairStartEnd[0] || pairCtr >= pairStartEnd[1]) )
										continue;
									
									if (pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot]) //rotamer pair not already pruned
										continue;
										
										
									boolean pruned = false;

									if(magicBullet){
										//MBEntry holds i_t and j_u
										pruned = CanEliminateUsing(i,i_r_aa,i_r_rot,j,j_s_aa,j_s_rot,
												mbEntry[1],mbEntry[2],mbEntry[4],mbEntry[5]);
									} else {
										pruned = CanEliminate(i,i_r_aa,i_r_rot,j,j_s_aa,j_s_rot);
									}

									if (pruned){

										pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot]=true;
										pairwiseMinEnergyMatrix.pairs.pruned[j][j_s_aa][j_s_rot][i][i_r_aa][i_r_rot]=true;

										prunedCurRun++;
									}
									

								}}

							}}
						}
					}
					//System.out.println("done");
				}
			}


//			System.out.println();
//			System.out.println("Number of pairs pruned this run: "+prunedCurRun);
//			System.out.println("DEE: The minimum difference is "+minDiff);


			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;

			if ((!magicBullet)&&(!distrDEE)) //full non-distributed pairs, so do not repeat (computationally-expensive)
				done = true;
			
			//Never repeat a pairs run
			done=true;
		}

//		System.out.println("DEEGoldsteinPairs pruned rotamers: "+prunedCurRun);

	}


	//Check if the pair denoted by pairEntry can be pruned by mbEntry
	//KER: Debug method that ignores the mbEntry given to the function
	private boolean CanEliminateUsing(int i,int i_r_aa,int i_r_rot,int j,int j_s_aa,int j_s_rot,
			int i_t_aa,int i_t_rot,int j_u_aa,int j_u_rot) {

		double checkSum=0;


		double minIndVoxelE =pairwiseMinEnergyMatrix.singles.E[i][i_r_aa][i_r_rot] + pairwiseMinEnergyMatrix.singles.E[j][j_s_aa][j_s_rot]; //Individual energies
		double minPairE = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot]; //Pair energy

//		if(typeDependent && (i_r_aa != i_t_aa || j_s_aa != j_u_aa))
//			return false;
		
		//If typeDependent we can set same AAs to have curEw of 0
		//Otherwise the curEw = initEw
		if(typeDependent && i_t_aa == i_r_aa && j_u_aa == j_s_aa)
			curEw = Ival;
		else{
			curEw = Ival+initEw;
		}
		
		//skip if this is the exact same pair
		if (i_r_aa == i_t_aa && i_r_rot == i_t_rot && j_s_aa == j_u_aa && j_s_rot == j_u_rot) 
			return false;

		double maxIndVoxelE = pairwiseMinEnergyMatrix.singles.E[i][i_t_aa][i_t_rot] + 
								pairwiseMinEnergyMatrix.singles.E[j][j_u_aa][j_u_rot];
		double maxPairE = pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][j][j_u_aa][j_u_rot];//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2];

		
		//get the contribution from the active site residue rotamers
		double minDiffPairVoxelE = SumMinDiffPVE(i,i_r_aa,i_r_rot,j,j_s_aa,j_s_rot,i_t_aa,i_t_rot,j_u_aa,j_u_rot);	//formula term 5

		checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE) + minDiffPairVoxelE; 

		if (checkSum > curEw)
			return true;
		else{
			if(checkSum != 0.0)
				minDiff = Math.max(minDiff,checkSum);
		}

		return false;
	}



	//Check if the pair denoted by pairEntry can be pruned by mbEntry
	//KER: Debug method that ignores the mbEntry given to the function
	private boolean CanEliminateIgnoreMBMod(EMatrixEntryWIndex pos1entry,
			EMatrixEntryWIndex pos2entry, EMatrixEntryWIndex pairEntry,
			EMatrixEntryWIndex[] mbEntry) {

		double checkSum=0;

		EMatrixEntryWIndex mbPos1 = mbEntry[0];
		EMatrixEntryWIndex mbPos2 = mbEntry[1];

		//Individual energies
		double minIndVoxelE = pos1entry.eme.minE() + pos2entry.eme.minE();
		//Pair energy
		double minPairE = pairEntry.eme.minE();


		/*if(pos1entry.aa1() == mbEntry.aa1() && pos1entry.rot1() == mbEntry.rot1())
			return false;
		if(pos2entry.aa1() == mbEntry.aa2() && pos2entry.rot1() == mbEntry.rot2())
			return false;*/

		Iterator<EMatrixEntryWIndex> iter1 = pairwiseMinEnergyMatrix.singlesIterator(((RotamerEntry)pos1entry.eme).pos);
		//for (int altAA1=0; altAA1<numAAtypes[posNum1]; altAA1++){
		while(iter1.hasNext()){
			EMatrixEntryWIndex altAA1 = iter1.next();

			if(pos1entry.aa1() == altAA1.aa1() && pos1entry.rot1() == altAA1.rot1())
				continue;

			if(altAA1.eme.isPruned())
				continue;

			Iterator<EMatrixEntryWIndex> iter2 = pairwiseMinEnergyMatrix.singlesIterator(((RotamerEntry)pos2entry.eme).pos);
			while(iter2.hasNext()){
				EMatrixEntryWIndex altAA2 = iter2.next();

				if(pos2entry.aa1() == altAA2.aa1() && pos2entry.rot1() == altAA2.rot1())
					continue;

				if(altAA2.eme.isPruned())
					continue;

				/*if(mbEntry[2].aa1() == altAA1.aa1() && mbEntry[2].aa2() == altAA2.aa1() &&
						mbEntry[2].rot1() == altAA1.rot1() && mbEntry[2].rot2() == altAA2.rot1()){
					//System.out.println("Found MB");
				}else
					continue;*/

				//Energies of MB
				//double maxIndVoxelE = pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot1index()) + pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot2index());
				//double maxPairE     = mbEntry.eme.minE();

				double MBmaxIndVoxelE = pairwiseMinEnergyMatrix.getSingleMinE(mbPos1.rot1index()) + pairwiseMinEnergyMatrix.getSingleMinE(mbPos2.rot1index());
				double MBmaxPairE     = mbEntry[2].eme.minE();


				double maxIndVoxelE = altAA1.eme.minE() + altAA2.eme.minE();//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0];		//formula term 2
				//maxShellResE = RotamerSearch.getShellRotE(pairwiseMaxEnergyMatrix, posNum1, altAA1, altRot1)+RotamerSearch.getShellRotE(pairwiseMaxEnergyMatrix, posNum2, altAA2, altRot2);//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][1] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][1];
				EMatrixEntryWIndex altPair = new EMatrixEntryWIndex(pairwiseMinEnergyMatrix,altAA1,altAA2);
				double maxPairE = altPair.eme.minE();//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2];

				double MBminDiffPairVoxelE = 0;
				double minDiffPairVoxelE = 0;
				//get the contribution from the active site residue rotamers
				for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			

					if ((curPos != pos1entry.pos1())&&(curPos!=pos2entry.pos1())){
						minDiffPairVoxelE += IndMinDiffPVE(pairEntry,altPair,curPos);
						MBminDiffPairVoxelE += IndMinDiffPVE(pairEntry,mbEntry[2],curPos);
					}
				}

				checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE) + minDiffPairVoxelE; 

				if (checkSum > curEw)
					return true;
				else{
					if(checkSum != 0.0)
						minDiff = Math.max(minDiff,checkSum);

					if (magicBullet){ //magic bullet pairs, so no further checks
						return false;
					}
				}
			}
		}

		return false;
	}



	//Check if the pair denoted by pairEntry can be pruned by mbEntry
	//KER: Debug method that ignores the mbEntry given to the function
	private boolean CanEliminateIgnoreMB(EMatrixEntryWIndex pos1entry,
			EMatrixEntryWIndex pos2entry, EMatrixEntryWIndex pairEntry,
			EMatrixEntryWIndex mbEntry) {

		double checkSum=0;

		//Individual energies
		double minIndVoxelE = pos1entry.eme.minE() + pos2entry.eme.minE();
		//Pair energy
		double minPairE = pairEntry.eme.minE();


		/*if(pos1entry.aa1() == mbEntry.aa1() && pos1entry.rot1() == mbEntry.rot1())
			return false;
		if(pos2entry.aa1() == mbEntry.aa2() && pos2entry.rot1() == mbEntry.rot2())
			return false;*/

		Iterator<EMatrixEntryWIndex> iter1 = pairwiseMinEnergyMatrix.singlesIterator(((RotamerEntry)pos1entry.eme).pos);
		//for (int altAA1=0; altAA1<numAAtypes[posNum1]; altAA1++){
		while(iter1.hasNext()){
			EMatrixEntryWIndex altAA1 = iter1.next();

			if(pos1entry.aa1() == altAA1.aa1() && pos1entry.rot1() == altAA1.rot1())
				continue;

			if(altAA1.eme.isPruned())
				continue;

			Iterator<EMatrixEntryWIndex> iter2 = pairwiseMinEnergyMatrix.singlesIterator(((RotamerEntry)pos2entry.eme).pos);
			while(iter2.hasNext()){
				EMatrixEntryWIndex altAA2 = iter2.next();

				if(pos2entry.aa1() == altAA2.aa1() && pos2entry.rot1() == altAA2.rot1())
					continue;

				if(altAA2.eme.isPruned())
					continue;

				//Energies of MB
				//double maxIndVoxelE = pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot1index()) + pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot2index());
				//double maxPairE     = mbEntry.eme.minE();

				double maxIndVoxelE = altAA1.eme.minE() + altAA2.eme.minE();//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0];		//formula term 2
				//maxShellResE = RotamerSearch.getShellRotE(pairwiseMaxEnergyMatrix, posNum1, altAA1, altRot1)+RotamerSearch.getShellRotE(pairwiseMaxEnergyMatrix, posNum2, altAA2, altRot2);//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][1] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][1];
				EMatrixEntryWIndex altPair = new EMatrixEntryWIndex(pairwiseMinEnergyMatrix,altAA1,altAA2);
				double maxPairE = altPair.eme.minE();//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2];

				double minDiffPairVoxelE = 0;
				//get the contribution from the active site residue rotamers
				for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			

					if ((curPos != pos1entry.pos1())&&(curPos!=pos2entry.pos1())){
						minDiffPairVoxelE += IndMinDiffPVE(pairEntry,altPair,curPos);
					}
				}

				checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE) + minDiffPairVoxelE; 

				if (checkSum > curEw)
					return true;
				else{
					if(checkSum != 0.0)
						minDiff = Math.max(minDiff,checkSum);

					/*if (magicBullet) //magic bullet pairs, so no further checks
					return false;
				}*/
				}
			}
		}

		return false;
	}


	//Check if the pair denoted by pairEntry can be pruned by mbEntry
	//	private boolean CanEliminateUsing(EMatrixEntryWIndex pos1entry,
	//			EMatrixEntryWIndex pos2entry, EMatrixEntryWIndex pairEntry,
	//			EMatrixEntryWIndex mbEntry) {
	//		
	//		double checkSum=0;
	//		
	//		//Individual energies
	//		double minIndVoxelE = pos1entry.eme.minE() + pos2entry.eme.minE();
	//		//Pair energy
	//		double minPairE = pairEntry.eme.minE();
	//		//Energies of MB
	//		double maxIndVoxelE = pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot1index()) + pairwiseMinEnergyMatrix.getSingleMinE(mbEntry.rot2index());
	//		double maxPairE     = mbEntry.eme.minE();
	//		
	//		if(pos1entry.aa1() == mbEntry.aa1() && pos1entry.rot1() == mbEntry.rot1())
	//			return false;
	//		if(pos2entry.aa1() == mbEntry.aa2() && pos2entry.rot1() == mbEntry.rot2())
	//			return false;
	//		
	//		double minDiffPairVoxelE = 0;
	//		//get the contribution from the active site residue rotamers
	//		for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			
	//				
	//			if ((curPos != pos1entry.pos1())&&(curPos!=pos2entry.pos1())){
	//				minDiffPairVoxelE += IndMinDiffPVE(pairEntry,mbEntry,curPos);
	//			}
	//		}
	//		
	//		checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE) + minDiffPairVoxelE; 
	//		
	//		if (checkSum > curEw)
	//			return true;
	//		else
	//			minDiff = Math.max(minDiff,checkSum);
	//		
	//				
	//		return false;
	//	}

	private int[] findMagicBullets(int i, int j) {

		double bestEMax = Double.POSITIVE_INFINITY;
		int[] retEntry = null;
		
		//For all rot at i, i_r
		for(int i_r_aa=0;i_r_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_r_aa++){
			for(int i_r_rot=0;i_r_rot < pairwiseMinEnergyMatrix.singles.E[i][i_r_aa].length;i_r_rot++){
				
			if (pairwiseMinEnergyMatrix.singles.pruned[i][i_r_aa][i_r_rot]) //if pruned skip
				continue;

			double i_r_E = pairwiseMinEnergyMatrix.singles.E[i][i_r_aa][i_r_rot];
			
			//For all rot at j, j_s
			for(int j_s_aa=0;j_s_aa < pairwiseMinEnergyMatrix.singles.E[j].length;j_s_aa++){
				for(int j_s_rot=0;j_s_rot < pairwiseMinEnergyMatrix.singles.E[j][j_s_aa].length;j_s_rot++){
				
					if(pairwiseMinEnergyMatrix.singles.pruned[j][j_s_aa][j_s_rot])//if pruned skip
						continue;

					if (pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot]) //rotamer pair not already pruned
						continue;
					
					double checkSum = 0;
	
					//Energies of individual rotamers					
					checkSum += i_r_E + pairwiseMinEnergyMatrix.singles.E[j][j_s_aa][j_s_rot];
				
					//Energy of pair
					checkSum += pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot];
	
					for (int k=0; k<pairwiseMinEnergyMatrix.numMutPos();k++){ //for k
						if((k != i) && (k != j) 
								&& pairwiseMinEnergyMatrix.areNeighbors(i, k) 
								&& pairwiseMinEnergyMatrix.areNeighbors(j, k)){
							double maxTerm = Double.NEGATIVE_INFINITY;
	
							for(int k_v_aa=0; k_v_aa<pairwiseMinEnergyMatrix.singles.pruned[k].length;k_v_aa++){
								for(int k_v_rot=0; k_v_rot< pairwiseMinEnergyMatrix.singles.pruned[k][k_v_aa].length;k_v_rot++){ //for all rot at k
									
									if(pairwiseMinEnergyMatrix.singles.pruned[k][k_v_aa][k_v_rot])
										continue;
	
									if(pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot] ||  
											pairwiseMinEnergyMatrix.pairs.pruned[j][j_s_aa][j_s_rot][k][k_v_aa][k_v_rot])
										continue;
	
									double checkMax = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot]
											+ pairwiseMinEnergyMatrix.pairs.E[j][j_s_aa][j_s_rot][k][k_v_aa][k_v_rot];
	
									if( checkMax > maxTerm)
										maxTerm = checkMax;
								}
							}
							//Once we've gone through all t for given k add the maxTerm to total
							checkSum += maxTerm;
						}
					}
	
					if(checkSum < bestEMax){ //We want the least E_max(i_r,j_s) over all r and s
						if(retEntry == null)
							retEntry = new int[6];
						
						retEntry[0] = i;
						retEntry[1] = i_r_aa;
						retEntry[2] = i_r_rot;
						retEntry[3] = j;
						retEntry[4] = j_s_aa;
						retEntry[5] = j_s_rot;
						
						bestEMax = checkSum;
	
					}
			}}
		}}
		return retEntry;
	}


	//	private EMatrixEntryWIndex findMagicBullets(int curPos1, int curPos2) {
	//		
	//		
	//		double bestEMax = Double.POSITIVE_INFINITY;
	//		EMatrixEntryWIndex retEntry = null;
	//		Iterator<EMatrixEntryWIndex> iter1 = pairwiseMinEnergyMatrix.singlesIterator(curPos1);
	//		//for (int altAA1=0; altAA1<numAAtypes[posNum1]; altAA1++){
	//		while(iter1.hasNext()){
	//			EMatrixEntryWIndex rot1 = iter1.next();
	//			if(rot1.eme.isPruned()) //if pruned skip
	//				continue;
	//				
	//				
	//				Iterator<EMatrixEntryWIndex> iter2 = pairwiseMinEnergyMatrix.singlesIterator(curPos2);
	//				while(iter2.hasNext()){
	//				EMatrixEntryWIndex rot2 = iter2.next();
	//				if(rot2.eme.isPruned())//if pruned skip
	//					continue;
	//				
	//					double checkSum = 0;
	//					
	//					//Energies of individual rotamers					
	//					checkSum += rot1.eme.minE() + rot2.eme.minE();
	//					//Energy of pair
	//					EMatrixEntryWIndex pair = new EMatrixEntryWIndex(pairwiseMinEnergyMatrix,rot1,rot2);
	//					if(pair.eme.isPruned()) //skip if pruned
	//						continue;
	//					
	//					checkSum += pair.eme.minE();
	//					
	//					for (int posk=0; posk<pairwiseMinEnergyMatrix.numMutPos();posk++){
	//						if((posk != curPos1) && (posk != curPos2)){
	//							double maxTerm = Double.NEGATIVE_INFINITY;
	//							
	//							for(int poskAA=0; poskAA<pairwiseMinEnergyMatrix.singles.pruned[posk].length;poskAA++){
	//								for(int poskRot=0; poskRot< pairwiseMinEnergyMatrix.singles.pruned[posk][poskAA].length;poskRot++){
	//									int[] poskIndex = {posk,poskAA,poskRot};
	//									if(!pairwiseMinEnergyMatrix.getSinglePruned(poskIndex)){
	//										if((!pairwiseMinEnergyMatrix.getPairPruned(rot1.index, poskIndex)) &&  
	//												(!pairwiseMinEnergyMatrix.getPairPruned(rot2.index, poskIndex))){
	//											double checkMax = pairwiseMinEnergyMatrix.getPairMinE(rot1.index, poskIndex)
	//																+ pairwiseMinEnergyMatrix.getPairMinE(rot2.index, poskIndex);
	//											
	//											if( checkMax > maxTerm)
	//												maxTerm = checkMax;
	//										}
	//									}
	//								}
	//							}
	//							
	//							checkSum += maxTerm;
	//							
	//						}
	//					}
	//					
	//					if(checkSum < bestEMax){
	//						retEntry = pair;
	//						bestEMax = checkSum;
	//					
	//					}
	//					
	//				
	//			}
	//			
	//		}
	//		
	//		
	//		
	//		return retEntry;
	//	}



	//Called only by ComputeEliminatedRotConf(.)
	private boolean CanEliminate (int i, int i_r_aa, int i_r_rot,int j, int j_s_aa, int j_s_rot){

		double minIndVoxelE, maxIndVoxelE;
		double minPairE, maxPairE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;

		double checkSum;

		minIndVoxelE = pairwiseMinEnergyMatrix.singles.E[i][i_r_aa][i_r_rot] + pairwiseMinEnergyMatrix.singles.E[j][j_s_aa][j_s_rot];// pairwiseMinEnergyMatrix[posNum1][AANumAtPos1][rotNumAtPos1][posNum1][0][0] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][0]; //formula term 1
		minPairE = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][j][j_s_aa][j_s_rot];//pairwiseMinEnergyMatrix[posNum1][AANumAtPos1][rotNumAtPos1][posNum2][AANumAtPos2][rotNumAtPos2];

			// 2010: Don't compute the intervals if useMinDEEPruningEw is true since they are not neces
			//if (doMinimize && !doIMinDEE){ //MinDEE, so compute the interval terms
			//	indVoxelInterval = indIntMinDEE2Pos[posNum1][posNum2];							//formula term 3
			///	pairVoxelInterval = pairIntMinDEE2Pos[posNum1][posNum2];						//formula term 4
			//}
			//else { //traditional-DEE, so no interval terms // 2010: if useMinDEEPruningEw is set to true then it is the same as DEE
			indVoxelInterval = 0.0;
			pairVoxelInterval = 0.0;
			//}


			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning
			//int numRotForAAatPos1;
			for(int i_t_aa=0;i_t_aa < pairwiseMinEnergyMatrix.singles.E[i].length;i_t_aa++){
				for(int i_t_rot=0;i_t_rot < pairwiseMinEnergyMatrix.singles.E[i][i_t_aa].length;i_t_rot++){

				
					if(pairwiseMinEnergyMatrix.singles.pruned[i][i_t_aa][i_t_rot]) //skip if pruned
						continue;

							
//					if(typeDependent && i_t_aa != i_r_aa) //Type Dep Pruning
//						continue;
					
					//For all other rotamers at position j, j_u
					for(int j_u_aa=0;j_u_aa < pairwiseMinEnergyMatrix.singles.E[j].length;j_u_aa++){
						for(int j_u_rot=0;j_u_rot < pairwiseMinEnergyMatrix.singles.E[j][j_u_aa].length;j_u_rot++){

							if(pairwiseMinEnergyMatrix.singles.pruned[j][j_u_aa][j_u_rot]) //skip if j_u pruned
								continue;

							if(pairwiseMinEnergyMatrix.pairs.pruned[i][i_t_aa][i_t_rot][j][j_u_aa][j_u_rot]) //skip if pair pruned
								continue;
							
//							if(typeDependent && j_u_aa != j_s_aa)
//								continue;
							
							//If typeDependent we can set same AAs to have curEw of 0
							//Otherwise the curEw = initEw
							if(typeDependent && i_t_aa == i_r_aa && j_u_aa == j_s_aa)
								curEw = Ival;
							else{
								curEw = Ival+initEw;
							}

							//skip if this is the exact same pair
							if (i_r_aa == i_t_aa && i_r_rot == i_t_rot && j_s_aa == j_u_aa && j_s_rot == j_u_rot) 
								continue;



							maxIndVoxelE = pairwiseMinEnergyMatrix.singles.E[i][i_t_aa][i_t_rot] + 
											pairwiseMinEnergyMatrix.singles.E[j][j_u_aa][j_u_rot];//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0];		//formula term 2
							

							maxPairE = pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][j][j_u_aa][j_u_rot];//pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2];

							minDiffPairVoxelE = SumMinDiffPVE(i,i_r_aa,i_r_rot,j,j_s_aa,j_s_rot,i_t_aa,i_t_rot,j_u_aa,j_u_rot);	//formula term 5

							checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE)
									- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

							if (checkSum > curEw){
								//System.out.println(index_r+" "+index_t+" "+checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);

								return true;}//this rotamer can be pruned/eliminated
							else {
								if(checkSum != 0.0)
									minDiff = Math.max(minDiff,checkSum);

								if (magicBullet) //magic bullet pairs, so no further checks
									return false;
							}
							
					}}
				
			}}
		

		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}

	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////

	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (int i,int i_r_aa,int i_r_rot,int j,int j_s_aa,int j_s_rot,
			int i_t_aa,int i_t_rot,int j_u_aa,int j_u_rot){

		//int atPos1, int withAA1, int withRot1, int altAA1, int altRot1,
		//int atPos2, int withAA2, int withRot2, int altAA2, int altRot2){

		double sum = 0;

		//get the contribution from the active site residue rotamers at positions that are not i or j
		for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){			

			if (curPos == i || curPos == j)
				continue;
				
			sum += IndMinDiffPVE(i,i_r_aa,i_r_rot,j,j_s_aa,j_s_rot,i_t_aa,i_t_rot,j_u_aa,j_u_rot,curPos);
			
		}

		return sum;
	}

	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (int i,int i_r_aa,int i_r_rot,int j,int j_s_aa,int j_s_rot,
			int i_t_aa,int i_t_rot,int j_u_aa,int j_u_rot, int k){


		//int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
		//int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1, int k){

		double minE = bigE;
		double curEmin, curEmax;
		
		int numRotForAAatPos;

		boolean found = false;
		double min1,min2;
		double max1,max2;
		
		for (int k_v_aa=0; k_v_aa<pairwiseMinEnergyMatrix.singles.pruned[k].length;k_v_aa++){
			for (int k_v_rot=0; k_v_rot<pairwiseMinEnergyMatrix.singles.pruned[k][k_v_aa].length; k_v_rot++){
				
				if (pairwiseMinEnergyMatrix.singles.pruned[k][k_v_aa][k_v_rot]) //not pruned
					continue;

				if(pairwiseMinEnergyMatrix.areNeighbors(i, k)){
					if (pairwiseMinEnergyMatrix.pairs.pruned[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot]) //not using split flags or not flagged
						continue;
					min1 = pairwiseMinEnergyMatrix.pairs.E[i][i_r_aa][i_r_rot][k][k_v_aa][k_v_rot];
					max1 = pairwiseMinEnergyMatrix.pairs.E[i][i_t_aa][i_t_rot][k][k_v_aa][k_v_rot];
				}else{
					min1 = 0;
					max1 = 0;
				}
				if(pairwiseMinEnergyMatrix.areNeighbors(j, k)){
					if (pairwiseMinEnergyMatrix.pairs.pruned[j][j_s_aa][j_s_rot][k][k_v_aa][k_v_rot]) //not using split flags or not flagged
						continue;
					min2 = pairwiseMinEnergyMatrix.pairs.E[j][j_s_aa][j_s_rot][k][k_v_aa][k_v_rot];
					max2 = pairwiseMinEnergyMatrix.pairs.E[j][j_u_aa][j_u_rot][k][k_v_aa][k_v_rot];
				}else{
					min2 = 0;
					max2 = 0;
				}
				
				curEmin = min1+min2;//pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][k][curAA][curRot] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][k][curAA][curRot];
				curEmax = max1+max2;//pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][k][curAA][curRot] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][k][curAA][curRot];

				if ((curEmin-curEmax) < minE)
					minE = curEmin-curEmax;

				found = true;
			}
		}

			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		

		if (!found) //no possible pairs found
			minE = Double.NEGATIVE_INFINITY; //contributes nothing to the sum

		return minE;
	}


	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (EMatrixEntryWIndex pair, EMatrixEntryWIndex altPair,
			int thirdPos){


		double minE = bigE;
		double curEmin, curEmax;

		boolean found = false;
		boolean foundk = false;
		boolean foundValid = false;
		boolean found_i_r_j_s = false;
		boolean found_i_u_j_v = false;
		
		double min1,min2;
		double max1,max2;

		for (int thirdAA=0; thirdAA<pairwiseMinEnergyMatrix.singles.pruned[thirdPos].length;thirdAA++){

			for (int thirdRot=0; thirdRot<pairwiseMinEnergyMatrix.singles.pruned[thirdPos][thirdAA].length; thirdRot++){

				int[] thirdIndex = {thirdPos, thirdAA, thirdRot};
				if (pairwiseMinEnergyMatrix.getSinglePruned(thirdIndex)) //not pruned 
					continue;

				foundk = true;

				

				//if alt pair is pruned then first pair might be better at this position

				//KER: Need to check if the 
				
				if(pairwiseMinEnergyMatrix.areNeighbors(pair.pos1(), thirdPos)){
					if(pairwiseMinEnergyMatrix.getPairPruned(pair.rot1index(), thirdIndex))
						continue;
					min1 = pairwiseMinEnergyMatrix.getPairMinE(pair.rot1index(), thirdIndex);
					max1 = pairwiseMinEnergyMatrix.getPairMinE(altPair.rot1index(),thirdIndex);
				}else{
					min1 = 0;
					max1 = 0;
				}
				if(pairwiseMinEnergyMatrix.areNeighbors(pair.pos2(), thirdPos)){
					if(pairwiseMinEnergyMatrix.getPairPruned(pair.rot2index(), thirdIndex))
						continue;
					min2 = pairwiseMinEnergyMatrix.getPairMinE(pair.rot2index(), thirdIndex);
					max2 = pairwiseMinEnergyMatrix.getPairMinE(altPair.rot2index(), thirdIndex);
				}else{
					min2 = 0;
					max2 = 0;
				}
				
				
				curEmin = min1 + min2;//pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][thirdPos][curAA][curRot] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][thirdPos][curAA][curRot];
				curEmax = max1 + max2;//pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][thirdPos][curAA][curRot] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][thirdPos][curAA][curRot];

				if ((curEmin-curEmax) < minE)
					minE = curEmin-curEmax;

				found = true;
				//}
			}
		}

		//if(foundk && !found_i_r_j_s && found_i_u_j_v)
		if(foundk && !found)
			minE = Double.POSITIVE_INFINITY;
		if (!found){ //no possible pairs found
			//minE = 0.0; //contributes nothing to the sum
			minE = Double.NEGATIVE_INFINITY; //KER: if there's not a valid pair then this cannot be pruned;
		}

		return minE;
	}

	//Called by SumMaxMaxPVE(.)
	/*private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
			int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1){

		double minE = bigE;
		double curEmin, curEmax;

		int index_r1, index_r2, index_t1, index_t2, index2;

		//r at i
		index_r1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;

		//t at i
		index_t1 = firstPos*numTotalRot + rotIndOffset[firstAltAA1] + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;

		boolean found = false;

		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned

			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){

				//s at j (the ligand residue)
				index2 = numSiteResidues*numTotalRot + curLigPos;

				if ((!eliminatedRotAtPos[index2])){ //not pruned 

					if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged

						curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][numSiteResidues][ligAANum][curLigPos];
						curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][numSiteResidues][ligAANum][curLigPos];

						if ((curEmin-curEmax) < minE)
							minE = curEmin-curEmax;

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

	//////////////////////////////////////////////////////////////////////////////////		
	//Called only by CanEliminate() adn CanEliminateLig()
	private double SumMaxIndInt (int withoutPos1, int withoutPos2){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<pairwiseMinEnergyMatrix.numMutPos(); curPos++){

			if ((curPos != withoutPos1)&&(curPos != withoutPos2)){

				sum += maxScale*MaxIndInt(curPos);
			}
		}

		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)) //if we are not currently checking ligand rotamers for pruning
				sum += maxScale*LigandMaxIndInt();
		}*/

		return sum;
	}

	//Called by SumMaxIndInt(.)
	private double MaxIndInt (int atPos){

		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;

		int numRotForAAatPos;

		//int str=mutRes2Strand[atPos];
		//int strResNum=strandMut[str][mutRes2MutIndex[atPos]];
		Iterator<EMatrixEntryWIndex> iter1 = pairwiseMinEnergyMatrix.singlesIterator(atPos);
		while(iter1.hasNext()){
			EMatrixEntryWIndex atPosEntry = iter1.next(); 

			curEInt = IndInt(atPosEntry);
			if (curEInt > maxEInt){
				maxEInt = curEInt;
			}
		}

		return maxEInt;
	}

	//Called by MaxIndInt(.)
	private double IndInt (EMatrixEntryWIndex atPosEntry){

		//s at j
		//Index3 index1 = new Index3(atPos,atAA,atRot);//atPos*numTotalRot + rotIndOffset[atAA] + atRot;

		if (!atPosEntry.eme.isPruned()){ //not pruned 

			double maxE = atPosEntry.eme.minE();//pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][0];
			double minE = atPosEntry.eme.minE();//pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][0];

			//double maxShell = RotamerSearch.getShellRotE(pairwiseMaxEnergyMatrix, atPos, atAA, atRot);//pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][1];
			//double minShell = RotamerSearch.getShellRotE(pairwiseMinEnergyMatrix, atPos, atAA, atRot);//pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][1];

			//if ((maxE<=stericEThreshIntra)&&(minE<=stericEThreshIntra)
			//		&&(maxShell<=stericEThreshPair)&&(minShell<=stericEThreshPair))
			return ((maxE) - (minE));
		}
		else
			return 0.0;//contributes nothing
	}

	//Called by SumMaxIndInt(.)
	/*private double LigandMaxIndInt(){

		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;

		for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){

			curEInt = LigandIndInt(curLigPos);

			if (curEInt > maxEInt){
				maxEInt = curEInt;
			}
		}

		return maxEInt;
	}*/

	//Called by LigandMaxIndInt(.)
	/*private double LigandIndInt (int ligRot){

		//s at j (the ligand residue)
		int index1 = numSiteResidues*numTotalRot + ligRot;

		if ((!eliminatedRotAtPos[index1])){ //not pruned 

			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];

			double maxShell = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];
			double minShell = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];

			//if ((maxE<=stericEThreshIntra)&&(minE<=stericEThreshIntra)
			//		&&(maxShell<=stericEThreshPair)&&(minShell<=stericEThreshPair))
				return ((maxE+maxShell) - (minE+minShell));
		}
		else
			return 0.0;//contributes nothing
	}*/
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//Called only by CanEliminate() and CanEliminateLig()
	private double SumSumMaxPairInt(int withoutPos1, int withoutPos2){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos1=0; curPos1<pairwiseMinEnergyMatrix.numMutPos(); curPos1++){
			if ((curPos1 != withoutPos1)&&(curPos1 != withoutPos2)){
				for (int curPos2=0; curPos2<curPos1; curPos2++){
					if ((curPos2 != withoutPos1)&&(curPos2 != withoutPos2)){

						sum += maxScale*MaxPairInt(curPos1,curPos2);
					}
				}
			}
		}

		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position k here for which to add;
			//the range of j is the number of active site residues
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)){ //if we are not currently checking ligand rotamers for pruning
				for (int curPos=0; curPos<numSiteResidues; curPos++){
					if ((curPos != withoutPos1)&&(curPos != withoutPos2)){

						sum += maxScale*LigandMaxPairInt(curPos);
					}
				}
			}
		}*/

		return sum;
	}

	//Called by SumSumMaxPairInt(.)
	private double MaxPairInt (int atPos1, int atPos2){

		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;

		int numRotForAAatPos1;

		/*int str1=mutRes2Strand[atPos1];
		int strResNum1=strandMut[str1][mutRes2MutIndex[atPos1]];
		int str2=mutRes2Strand[atPos2];
		int strResNum2=strandMut[str2][mutRes2MutIndex[atPos2]];*/
		Iterator<EMatrixEntryWIndex> iter1 = pairwiseMinEnergyMatrix.singlesIterator(atPos1);
		while(iter1.hasNext()){
			EMatrixEntryWIndex atPos1Entry = iter1.next(); 
			//for (int curAA1=0; curAA1<numAAtypes[atPos1]; curAA1++){

			//int curAA1 = strandMut.getIndexOfNthAllowable(atPos1,AA1);

			//numRotForAAatPos1 = eliminatedRotAtPos.numRot(atPos1, curAA1);//strandMut.getNumRotForAAtype(atPos1,curAA1);
			//if (numRotForAAatPos1==0)	//ala or gly
			//	numRotForAAatPos1 = 1;

			//for (int curRot1=0; curRot1<numRotForAAatPos1; curRot1++){

			//int numRotForAAatPos2;

			Iterator<EMatrixEntryWIndex> iter2 = pairwiseMinEnergyMatrix.singlesIterator(atPos2);
			while(iter2.hasNext()){
				EMatrixEntryWIndex atPos2Entry = iter2.next(); 
				//for (int curAA2=0; curAA2<numAAtypes[atPos2]; curAA2++){

				//int curAA2 = strandMut.getIndexOfNthAllowable(atPos2,AA2);

				//numRotForAAatPos2 = eliminatedRotAtPos.numRot(atPos2, curAA2);//strandMut.getNumRotForAAtype(atPos2,curAA2);
				//if (numRotForAAatPos2==0)	//ala or gly
				//	numRotForAAatPos2 = 1;

				//for (int curRot2=0; curRot2<numRotForAAatPos2; curRot2++){

				curEInt = PairInt(atPos1Entry,atPos2Entry);
				if (curEInt > maxEInt){
					maxEInt = curEInt;
				}

			}
		}

		return maxEInt;
	}

	//Called by MaxPairInt(.)
	private double PairInt (EMatrixEntryWIndex atPos1Entry, EMatrixEntryWIndex atPos2Entry){

		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		//Index3 index1 = new Index3(atPos1,atAA1,atRot1);//atPos1*numTotalRot + rotIndOffset[atAA1] + atRot1;//u at k
		//Index3 index2 = new Index3(atPos2,atAA2,atRot2);//atPos2*numTotalRot + rotIndOffset[atAA2] + atRot2;//s at j

		if (((!atPos1Entry.eme.isPruned()))&&
				((!atPos2Entry.eme.isPruned()))){ //not pruned 

			double maxE = pairwiseMinEnergyMatrix.getPairMaxE(atPos1Entry.index, atPos2Entry.index,doIMinDEE);//[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];
			double minE = pairwiseMinEnergyMatrix.getPairMinE(atPos1Entry.index, atPos2Entry.index);//[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];

			//if ((maxE<=stericEThreshPair)&&(minE<=stericEThreshPair))
			return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
	}

	//Called by SumSumMaxPairInt(.)
	/*private double LigandMaxPairInt (int atPos){

		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;

		int numRotForAAatPos;

		for (int AA=0; AA<numAAtypes[atPos]; AA++){

			int curAA = sysLR.getIndexOfNthAllowable(residueMap[atPos],AA);

			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
			if (numRotForAAatPos==0)	//ala or gly
				numRotForAAatPos = 1;

			for (int curRot=0; curRot<numRotForAAatPos; curRot++){

				for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){

					curEInt = LigandPairInt(atPos, curAA, curRot, curLigPos);
					if (curEInt > maxEInt){
						maxEInt = curEInt;
					}
				}
			}
		}

		return maxEInt;
	}*/

	//Called by LigandMaxPairInt(.)
	/*private double LigandPairInt (int atPos, int atAA, int atRot, int ligRot){

		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = numSiteResidues*numTotalRot + ligRot;//u at k (the ligand residue)
		int index2 = atPos*numTotalRot + rotIndOffset[atAA] + atRot;//s at j

		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index2]))){ //not pruned 

			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];

			//if ((maxE<=stericEThreshPair)&&(minE<=stericEThreshPair))
				return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
	}*/
	///////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////
	//Same as CanEliminate(), just checks the ligand rotamers for pruning
	//Called by ComputeEliminatedRotConf()
	/*private boolean CanEliminateLig (int curLigRot,int posNum2, int AANumAtPos2, int rotNumAtPos2){

		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double minPairE, maxPairE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;

		int index_r1, index_r2, index_t1, index_t2;

		double checkSum;

		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r1 = numSiteResidues*numTotalRot + curLigRot;
		index_r2 = posNum2*numTotalRot + rotIndOffset[AANumAtPos2] + rotNumAtPos2;
		minIndVoxelE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][0] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][0]; //formula term 1
		minShellResE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][1] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][1];
		minPairE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][posNum2][AANumAtPos2][rotNumAtPos2];

		//if ((minIndVoxelE<=stericEThreshIntra)&&(minShellResE<=stericEThreshPair)){//check only if not an unallowed steric
		if ((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2])
				&&(!splitFlags[index_r1][index_r2])){ //not pruned

			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE2Pos[numSiteResidues][posNum2];							//formula term 3
				pairVoxelInterval = pairIntMinDEE2Pos[numSiteResidues][posNum2];						//formula term 4
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

					index_t1 = numSiteResidues*numTotalRot + altRot;

					if ((!eliminatedRotAtPos[index_t1])){ //not pruned 

						int numRotForAAatPos2;

						for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++){

							int altAA2 = sysLR.getIndexOfNthAllowable(residueMap[posNum2],AA2);

							numRotForAAatPos2 = rl.getNumRotForAAtype(altAA2);
							if (numRotForAAatPos2==0)	//ala or gly
								numRotForAAatPos2 = 1;

							for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++){

								//if t and r are not actually the same rotamer of the same AA
								if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){

									index_t2 = posNum2*numTotalRot + rotIndOffset[altAA2] + altRot2;

									if ((!eliminatedRotAtPos[index_t2])){ //not pruned 

										maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0]; //formula term 2
										maxShellResE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][1] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][1];	
										maxPairE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][posNum2][altAA2][altRot2];				

										minDiffPairVoxelE = SumMinDiffPVELig(curLigRot, altRot, posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2);	//formula term 5

										checkSum = -templateInt + (minIndVoxelE + minShellResE + minPairE) - (maxIndVoxelE + maxShellResE + maxPairE)
													- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

										if (checkSum > curEw){
											//System.out.println(checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);

											return true;}//this rotamer can be pruned/eliminated
										else {
											minDiff = Math.max(minDiff,checkSum);

											if (magicBullet) //magic bullet pairs, so no further checks
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

		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}*/

	//Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
	//Called by CanEliminateLig()
	/*private double SumMinDiffPVELig (int withRot1, int altRot1, int atPos2, int withAA2, int withRot2, int altAA2, int altRot2){

		double sum = 0;

		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){	

			if ((curPos!=atPos2))

				sum += IndMinDiffPVELig(withRot1, altRot1, atPos2, withAA2, withRot2, altAA2, altRot2, curPos);
		}

		return sum;
	}*/

	//Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
	//Called by SumMinDiffPVELig()
	/*private double IndMinDiffPVELig (int firstRot1, int firstAltRot1, int secondPos, int secondAA1, int secondRot1, 
			int secondAltAA1, int secondAltRot1, int thirdPos){

		double minE = bigE;
		double curEmin, curEmax;

		int index_r1, index_r2, index_t1, index_t2, index2;
		int numRotForAAatPos;

		//r at i
		index_r1 = numSiteResidues*numTotalRot + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;

		//t at i
		index_t1 = numSiteResidues*numTotalRot + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;

		boolean found = false;

		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned

			for (int AA=0; AA<numAAtypes[thirdPos]; AA++){

				int curAA = sysLR.getIndexOfNthAllowable(residueMap[thirdPos],AA);

				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;

				for (int curRot=0; curRot<numRotForAAatPos; curRot++){

					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1

					//s at j
					index2 = thirdPos*numTotalRot + rotIndOffset[curAA] + curRot;

					if ((!eliminatedRotAtPos[index2])){ //not pruned 

						if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged

							curEmin = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][firstRot1][thirdPos][curAA][curRot] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][thirdPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][firstAltRot1][thirdPos][curAA][curRot] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][thirdPos][curAA][curRot];

							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;

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
