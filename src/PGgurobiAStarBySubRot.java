import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

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
//	MSAStar.java
//
//	Version:           2.0
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.io.*;

import mpi.MPIException;

/**
 * Written by Pablo Gainza (2004-2012)
 */

/**
 * Uses A* search for single or multiple mutation sequences simultaneously to return the minimum-energy conformation; 
 * 		each consecutive run returns the next lowest-energy conformation, in order.
 * 
 */
public class PGgurobiAStarBySubRot extends AStar{

	//number of residues under consideration
//	private int numTreeLevels;

	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numNodesForLevel[] = null;

	//the total number of possible rotamers for the given mutation
	private int numTotalNodes;

	private int[] numParentRotPerLvl = null;
	private int[][] parentRotIndexPerLvl = null;
	private int[][] numSubRotPerParentRot = null;
	//First dim is the pos, second dim is the parent;
	private ArrayList<HashMap<Integer,ArrayList<Index3>>> subRotsPerLvlPerParent = null;

	//the current sequence: the number of the corresponding rotamer for each level assigned so far 
	private int curConf[] = null;

	//the reduced min pairwise energy matrix
//	private EMatrixEntryWIndex pairwiseMinEnergyMatrix [][] = null;
	
	Emat emat;
	EMatrixEntryWIndex[] bestConf; 
	EMatrixEntryWIndex[] parentConf;
	Molecule m;
	int topConf = 0;

	ArrayList<ArrayList<EnergyTuple>> tuplesPerPos;


	int topL = 0;
	int numTopL = 0;

	double upperE;



	//constructor
	/*
	 * We assume that the parameters supplied (energy and DEE information) have already been modified
	 * 		to consider only the residues and rotamers for the possible mutations; i.e., the matrices are of
	 * 		reduced size
	 */
	PGgurobiAStarBySubRot (int treeLevels, int numRotForRes[], Emat arpMatrix, double upperE, Molecule m){

		emat = arpMatrix;
		this.m = m;
		this.upperE = upperE; 
		numTreeLevels = treeLevels;

		numNodesForLevel = new int [numTreeLevels];
//		nodeIndexOffset = new int [numTreeLevels];
		numTotalNodes = 0;//arpMatrix.length-1;

		
		for (int i=0; i<numTreeLevels; i++){
//			nodeIndexOffset[i] = numTotalNodes;
			numNodesForLevel[i] = numRotForRes[i];
			numTotalNodes += numNodesForLevel[i];
		}
		
		numParentRotPerLvl = new int[numTreeLevels];
		parentRotIndexPerLvl = new int[numTreeLevels][];
		numSubRotPerParentRot = new int[numTreeLevels][];
		//First dim is the pos, second dim is the parent;
		subRotsPerLvlPerParent = new ArrayList<HashMap<Integer,ArrayList<Index3>>>();
		for(int i=0; i<numTreeLevels;i++)
			subRotsPerLvlPerParent.add(new HashMap<Integer,ArrayList<Index3>>());
		
		SinglesIterator iter = emat.singlesIterator();
		//This assumes we aren't doing super-rotamers
		//Since we are doing sub-rotamers, this should be a safe assumption
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(!emeWI.eme.isPruned()){
				int parent = m.getRC(emat.resByPos.get(emeWI.pos1()).get(0),emat.singles.getRot(emeWI.index)[0]).parent; 
				Index3 rot = new Index3(emeWI.index);
				if(!subRotsPerLvlPerParent.get(emeWI.pos1()).containsKey(parent))
					subRotsPerLvlPerParent.get(emeWI.pos1()).put(parent, new ArrayList<Index3>());
				
				subRotsPerLvlPerParent.get(emeWI.pos1()).get(parent).add(rot);
			}
		}
		
		for(int i=0; i<numTreeLevels;i++){
			numParentRotPerLvl[i] = subRotsPerLvlPerParent.get(i).size();
			parentRotIndexPerLvl[i] = new int[numParentRotPerLvl[i]];
			numSubRotPerParentRot[i] = new int[numParentRotPerLvl[i]];
			//Won't be sorted, but that's ok
			int ctr=0;
			for(Integer globalRot: subRotsPerLvlPerParent.get(i).keySet()){
				parentRotIndexPerLvl[i][ctr] = globalRot; 
				numSubRotPerParentRot[i][ctr] = subRotsPerLvlPerParent.get(i).get(globalRot).size();
				ctr++;
			}
		}
		
		
		//the min energy matrix: the last column contains the intra-energy for each rotamer; the last row
		//		contains the shell-residue energy for each rotamer
		
		


		/*for (int i=0; i<numTotalNodes+1; i++){
			for (int j=0; j<numTotalNodes+1; j++){
				pairwiseMinEnergyMatrix[i][j] = arpMatrixRed[i][j];
			}			
		}*/

		//the current expansion list
		curExpansion = new PGExpansionQueue();

		//the current conformation
		curConf = new int [numTreeLevels];
		for (int i=0; i<numTreeLevels; i++){
			curConf[i] = -1;
		}
		
		initialize2Dto3D(emat, numNodesForLevel);
		
		boolean useGurobi = false;
		boolean useWCSP = true;
		
		if(useWCSP || useGurobi){
		GurobiCalcInfo gci = new GurobiCalcInfo(arpMatrix,numNodesForLevel,numTotalNodes,false,upperE,
				numParentRotPerLvl,parentRotIndexPerLvl,numSubRotPerParentRot,subRotsPerLvlPerParent);
		
		CommucObj[] cObj = new CommucObj[1];
		cObj[0] = new CommucObj();
		if(useGurobi)
			cObj[0].gurobiCalc = true;
		else
			cObj[0].wcspCalc = true;
		cObj[0].gurobiCalcInfo = gci;
		
		
		//If running on the master node then use other nodes to compute energies
		//Set the slaves to the right mode now.
		try{
		if(MPItoThread.Rank() == 0 && KSParser.numProc >= 3){
			for (int curProc=1; curProc<KSParser.numProc; curProc++){ 
				try {
					MPItoThread.Send(cObj, 0, 1, ThreadMessage.OBJECT, curProc, KSParser.regTag);
				} catch (MPIException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}}
		catch(Exception E){
			System.out.println("Can't start slaves for AS");
		}
		}
		
		
		
		

	}

	//Find the lowest-energy conformation and return an array with the number of the corresponding
	//		chosen rotamer for each residue;
	//		the mapping to the rotamer information is done by the calling procedure;
	//		the chosen conformation should be marked so that the next call to AStar will return the
	//			conformation with the next lowest energy value, and so on
	/*
	 * Look at the minimum value in the expansion list; determine the node level corresponding to this value;
	 * 		expand the node into the next level; update the expansion list (adding the new nodes and deleting
	 * 		the expanded node) and determine the f(n)=g(n)+h(n) scores for the new nodes
	 * 
	 * To get the next lowest conformation, the state of the expansion queue is saved after each complete run
	 * 		of A*: after the lowest-energy conformation is found, the queue is saved, and A* returns this
	 * 		conformation. To find the second conformation, A* runs on the saved queue, and this is repeated
	 * 		for all subsequent conformations
	 */
	public PGQueueNode doAStar (boolean run1){

		int curLevelNum = 0;
		double hScore;
		double gScore;
		double fScore;
		PGQueueNode expNode = null;
		PGQueueNode newNode = null;

		int countNodes = 0;
		
		if (run1) {//if this is the first run of A*, then we need to set-up the empty head node
		
			//Expand nodes that only have 1 or 2 possibilities so we don't waste our time on them
			//First count how many there are
			int numInitialNodes = 1;
			for(int i=0; i<numParentRotPerLvl.length;i++)
				if(numParentRotPerLvl[i] <= 2)
					numInitialNodes *= numParentRotPerLvl[i];
			
			int[][] confs = new int[numInitialNodes][];
			makeConfs(0,numParentRotPerLvl,numTreeLevels,new ArrayList<Integer>(),confs); 
			
			
			for(int i=0; i<confs.length;i++){
				//create an empty PGQueueNode.  We do this by assigning a -1 at level 0.
				newNode = new PGQueueNode (numTreeLevels, confs[i], Double.NEGATIVE_INFINITY,0,confs[i][0]);
	
				//insert in the expansion list
				curExpansion.insert(newNode);
			}
		}							


		boolean done = false;		
		//start forming the conformation by expanding the lowest-valued node and updating the expansion queue
		/*
		 * While not at the last level
		 * 	For the current minimum node in the queue
		 * 		For all possible nodes at the next level
		 * 				Compute the f(n) scores; Add to the queue
		 * 		Delete the expanded node from the queue
		 */
		while (!done) {	

			for (int i=0; i<numTreeLevels; i++){ //reinitialize for each consecutive node to be expanded
				curConf[i] = -1;
			}

			expNode = (PGQueueNode) curExpansion.getMin();//get the current min node

			if (expNode==null){//the queue is empty
				return expNode; //so return a sequence of -1's to flag the stop of the search
			}

			else { //the queue is not empty

				printState(expNode);

				for (int i=0; i< numTreeLevels; i++){//get the corresponding conformation leading to this node
					curConf[i] = expNode.confSoFar[i];
				}

				//if the current node is fully assigned, we have found a full conformation
				if (expNode.emptyLevels.isEmpty()){
					//curExpansion.delete(expNode);//delete this node to set-up for the next min conformation (next run of A*)
					if(retConfE - expNode.fScore > 0.4 && retConfE < 0){
						System.out.println(retConfE +" "+expNode.fScore);
						if(expNode.curTuples != null)
							System.out.println(expNode.curTuples);
					}

					retConfE = expNode.fScore;
//					GurobiOptimization gOpt = gurobiFscoreOptimizer(expNode,true);
//					double finalScore = gOpt.getObjVal();

					WCSPOptimization optimizer = new WCSPOptimization(expNode,emat,numParentRotPerLvl,
							parentRotIndexPerLvl, numSubRotPerParentRot, subRotsPerLvlPerParent,
							numNodesForLevel,upperE);
					double finalScore = optimizer.optimize(null);
					optimizer.cleanUp();
					
					
					if(finalScore != retConfE){
						System.out.println("Fscore: "+retConfE +"  Gscore: "+finalScore);
						expNode.fScore = finalScore;
						curExpansion.insert(expNode);
					}
					else{
						//Convert seq to conf
//						expNode.confSoFar = gOpt.getConfFromSeq(expNode.confSoFar,seqIndexOffset);
						parentConf = new EMatrixEntryWIndex[expNode.confSoFar.length];
						for(int i=0; i<parentConf.length;i++){
							EMatrixEntryWIndex curRot = optimizer.bestConf.conf[i];
							int parentGlobalID = m.getRC(emat.resByPos.get(curRot.pos1()).get(0), emat.singles.getRot(curRot.index)[0]).parent;
							SuperRotamer r = new SuperRotamer(parentGlobalID);
							RotamerEntry rE = new RotamerEntry(curRot.pos1(),r);
							EMatrixEntryWIndex emeWI = new EMatrixEntryWIndex(rE, curRot.index);
							parentConf[i] = emeWI;
						}
						
						
						bestConf = optimizer.bestConf.conf;  
						done = true;
					}
				}
				else {// Queue not empty, we can continue.

					//					PGQueueNode[] nextLevelNodes = pickNextLevel(expNode);
					
					PGQueueNode[] nextLevelNodes = pickNextLevel(expNode);
					
					for(int rot = 0; rot < nextLevelNodes.length; rot++){

//						//Recompute better bound when inserting
//						if(energyTuples.size() > 0)
//							nextLevelNodes[rot].fScore = gurobiTupleFscore(nextLevelNodes[rot], expNode);
//						else
//							nextLevelNodes[rot].fScore = gurobiFscore(nextLevelNodes[rot]);

						if(expNode.fScore != 0 && nextLevelNodes[rot].fScore - expNode.fScore < -0.1)
							System.out.println("Something went wrong with fScores");

						curExpansion.insert(nextLevelNodes[rot]);	
						countNodes++;
					}

					//delete the expanded node from the queue, since it has already contributed
					//curExpansion.delete(expNode);
				}	
			}
		}
		System.out.println("Number of A* nodes inserted in the queue: "+countNodes);

		expNode.actualConf = bestConf;
		
		return expNode;
	}
	
	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.

	private void makeConfs(int level, int[] numSeqForLevel, int numTreeLevels,
			ArrayList<Integer> curConf, int[][] confs) {
		if(curConf.size() == numTreeLevels){
			int confNum = -1;
			for(int i=0; i<confs.length;i++)
				if(confs[i] == null){
					confNum = i;
					break;
				}
			
			confs[confNum] = new int[numTreeLevels];
			for(int i=0; i< curConf.size();i++)
				confs[confNum][i] = curConf.get(i);
		}else{ //Have to build the conformation
			
			if(numSeqForLevel[level] <= 2){
				for(int i=0; i<numSeqForLevel[level];i++){
					curConf.add(i);
					makeConfs(level+1,numSeqForLevel,numTreeLevels,curConf,confs);
					curConf.remove(curConf.size()-1);
				}
			}else{
				curConf.add(-1);
				makeConfs(level+1,numSeqForLevel,numTreeLevels,curConf,confs);
				curConf.remove(curConf.size()-1);
			}
			
			
		}
		
	}

	private PGQueueNode[] pickNextLevel(PGQueueNode dequeuedNode){		
		PGQueueNode[][] allLevelNodes = new PGQueueNode[dequeuedNode.confSoFar.length][];   
		
		int maxMinFScoreIndex = -1;
		double maxMinFScore = Double.NEGATIVE_INFINITY;
		int numRunning = 0;
		
		int numLevels = dequeuedNode.emptyLevels.size();
		boolean reorder = false;
		if(!reorder)
			numLevels = 1;
		// Iterate through all positions that have not been assigned in dequeuedNode.
		for(int iter = 0; iter < numLevels; iter++){
			int pos = dequeuedNode.emptyLevels.get(iter);
			double minAtThisLevel = Double.POSITIVE_INFINITY;
			int minIndexAtThisLevel = -1; 
			allLevelNodes[pos] = new PGQueueNode[numParentRotPerLvl[pos]];
	
			// And for each rotamer at that level
			for(int rot = 0; rot < numParentRotPerLvl[pos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
				//KER: Trying out gCompute as a heuristic
				//newNode.fScore = gCompute(newNode);//fCompute(newNode);
				try{
					if(MPItoThread.Rank() == 0 && KSParser.numProc >= 3){
						int curProc = numRunning+1;
						if(numRunning == KSParser.numProc-1){
							PGQueueNode[] node = new PGQueueNode[1];
							Object s = MPItoThread.Recv(node, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);
							allLevelNodes[node[0].level][node[0].nodeNum] = node[0];
							curProc = MPItoThread.getStatusSource(s);
							numRunning--;
						}
						PGQueueNode[] nodeToSend = new PGQueueNode[1];
						nodeToSend[0] = newNode;
						MPItoThread.Send(nodeToSend, 0, 1, ThreadMessage.OBJECT, curProc, KSParser.regTag);
						numRunning++;
							
					}else{
						newNode.fScore = toulbarFscore(newNode);
//						newNode.fScore = gurobiFscore(newNode,false);
						allLevelNodes[newNode.level][newNode.nodeNum] = newNode;
					}
				}catch(Exception E){
					E.printStackTrace();
				}
				/*if(newNode.fScore - dequeuedNode.fScore < -0.3 && dequeuedNode.fScore < 0){
					System.out.println("DELETE ME");
				}*/
					

			}			
		
		}
		
		try{
		while(numRunning > 0){
			PGQueueNode[] node = new PGQueueNode[1];
			Object s = MPItoThread.Recv(node, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);
			allLevelNodes[node[0].level][node[0].nodeNum] = node[0];
			numRunning--;
		}}catch(Exception E){
			E.printStackTrace();
		}
		
//		double[] diffAtLevel = new double[dequeuedNode.confSoFar.length];
		double[][] minsAtLevel = new double[dequeuedNode.confSoFar.length][2];
		
		
		
		//Try ordering by the largest gap between the best and second best
		//There should always be a second best because we expand all the
		//1 level nodes already
		for(int iter=0; iter<numLevels;iter++){
			int pos = dequeuedNode.emptyLevels.get(iter);
			minsAtLevel[pos][0] = Double.POSITIVE_INFINITY;
			minsAtLevel[pos][1] = Double.POSITIVE_INFINITY;
			
			for(int rot=0; rot<numParentRotPerLvl[pos]; rot++){	
				if(allLevelNodes[pos][rot].fScore < minsAtLevel[pos][0]){
					if(minsAtLevel[pos][0] < minsAtLevel[pos][1])
						minsAtLevel[pos][1] = minsAtLevel[pos][0];
					minsAtLevel[pos][0] = allLevelNodes[pos][rot].fScore;
					
				}
				else if(allLevelNodes[pos][rot].fScore < minsAtLevel[pos][1]){
					minsAtLevel[pos][1] = allLevelNodes[pos][rot].fScore;
				}
			}
			
			double diff = Math.abs(minsAtLevel[pos][0] - minsAtLevel[pos][1]);
			
			if( diff > maxMinFScore){
				maxMinFScore = diff;
				maxMinFScoreIndex = pos;
			}
		}
			
		// For now return the first level.  This should work the same as the old A*.... 
		if(!dequeuedNode.emptyLevels.isEmpty()){
			return allLevelNodes[maxMinFScoreIndex];			
		}
		else{
			return null;
		}
	}


	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
//	private PGQueueNode[] pickNextLevelByTuple(PGQueueNode dequeuedNode){		
//
//		PGQueueNode[] bestLevelNodes = null;
//
//		int maxMinFScoreIndex = -1;
//		double maxMinFScore = Double.NEGATIVE_INFINITY;
//
//		int maxTuples = 0;
//		int maxTuplePos = -1;
//
//		ArrayList<Index3> dequeuedNodeConf = new ArrayList<Index3>();
//		for(int pos: dequeuedNode.nonEmptyLevels){
//			int index2 = nodeIndexOffset[pos]+seqIndexOffset[pos][dequeuedNode.confSoFar[pos]] + dequeuedNode.confSoFar[pos];
//			dequeuedNodeConf.add(rotIndexes[index2]);
//		}
//
//
//		// Iterate through all positions that have not been assigned in dequeuedNode.
//		for(int iter = 0; iter < dequeuedNode.emptyLevels.size(); iter++){
//			int pos = dequeuedNode.emptyLevels.get(iter);
//			int numTuples = 0;
//
//			//Find the number of tuples that share at least one rotamer 
//			for(EnergyTuple tuple: tuplesPerPos.get(pos)){
//				if(tuple.shareRots(dequeuedNodeConf)){
//					numTuples++;
//				}
//			}
//
//
//			if(numTuples > maxTuples){
//				maxTuples = numTuples;
//				maxTuplePos = pos;
//			}
//		}
//
//		//If there isn't a spot with more tuples
//		//Use the position with the most tuples
//		if(maxTuplePos == -1){
//			for(int pos: dequeuedNode.emptyLevels){
//				if(tuplesPerPos.get(pos).size() > maxTuples){
//					maxTuples = tuplesPerPos.get(pos).size();
//					maxTuplePos = pos;
//				}
//			}
//		}
//
//
//		if(maxTuplePos >= 0){
//			//Generate the nodes for the best tuple position
//			bestLevelNodes = new PGQueueNode[numSeqForLevel[maxTuplePos]];
//			for(int aa = 0; aa < numSeqForLevel[maxTuplePos]; aa++){				
//				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, maxTuplePos, aa);  // No energy yet
//				bestLevelNodes[aa] = newNode;
//			}
//
//
//		}else{//Finally, if we still haven't found a position do the normal check
//			ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>(); 
//			for(int iter = 0; iter < dequeuedNode.emptyLevels.size(); iter++){
//				int pos = dequeuedNode.emptyLevels.get(iter);
//
//				double minAtThisLevel = Double.POSITIVE_INFINITY;
//				int minIndexAtThisLevel = -1; 
//				PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
//				// And for each rotamer at that level
//				for(int aa = 0; aa < numSeqForLevel[pos]; aa++){				
//					PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, aa);  // No energy yet
//					//KER: Trying out gCompute as a heuristic
//					//newNode.fScore = gCompute(newNode);//fCompute(newNode);
//					newNode.fScore = gurobiFscore(newNode,false);
//					/*if(newNode.fScore - dequeuedNode.fScore < -0.3 && dequeuedNode.fScore < 0){
//						System.out.println("DELETE ME");
//					}*/
//					curLevelPGQueueNodes[aa] = newNode;	
//
//					if(newNode.fScore < minAtThisLevel){
//						minAtThisLevel = newNode.fScore;
//						minIndexAtThisLevel = aa;
//					}
//				}			
//				allLevelNodes.add(curLevelPGQueueNodes);
//
//				if(minAtThisLevel > maxMinFScore){
//					maxMinFScore = minAtThisLevel;
//					maxMinFScoreIndex = iter;
//				}			
//			}
//
//			// For now return the first level.  This should work the same as the old A*....
//			if(!allLevelNodes.isEmpty()){
//				bestLevelNodes = allLevelNodes.get(maxMinFScoreIndex);			
//			}
//
//		}
//
//		return bestLevelNodes;
//
//	}


	//Updates and prints the state of the queue
	private void printState(PGQueueNode expNode){

		numExpanded++;

		if (expNode.nonEmptyLevels.size()>topL){
			topL = expNode.nonEmptyLevels.size();
			numTopL = 1;
		}
		else if (expNode.level+1==topL)
			numTopL++;


		if((numExpanded%2)==0){
			System.out.print(curExpansion.numNodes()+" "+expNode.fScore+" level:"+expNode.level+" numElem:"+expNode.nonEmptyLevels.size() + " elem:");
			for (int i=0;i<numTreeLevels;i++){
				System.out.print(expNode.confSoFar[i]+" ");
			}
			System.out.println();
			System.out.println("Top level:"+topL+" #"+numTopL);
		}
	}
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////


	
	//KER: Use gurobi LP to get a bound on the score so far
//	private double gurobiFscore(PGQueueNode node,boolean doInteger) {
//		GurobiOptimization optimizer = new GurobiOptimization(node,pairwiseMinEnergyMatrix,numNodesForLevel,nodeIndexOffset,numRotRemainingBySeq,seqIndexOffset,numTotalNodes,doInteger);
//
//		return optimizer.optimize();
//
//	}
	
	//KER: Use gurobi LP to get a bound on the score so far
	private double toulbarFscore(PGQueueNode node) {
		WCSPOptimization optimizer = new WCSPOptimization(node,emat,numParentRotPerLvl,
											parentRotIndexPerLvl, numSubRotPerParentRot, 
											subRotsPerLvlPerParent,	numNodesForLevel,upperE);
		System.out.print(".");
		//double E =optimizer.optimize(null);
		double E = optimizer.getBound(null);
		optimizer.cleanUp();
		return E;

	}
	

	//KER: Use gurobi LP to get a bound on the score so far
//	private GurobiOptimization gurobiFscoreOptimizer(PGQueueNode node,boolean doInteger) {
//		GurobiOptimization optimizer = new GurobiOptimization(node,pairwiseMinEnergyMatrix,numNodesForLevel,nodeIndexOffset,numRotRemainingBySeq,seqIndexOffset,numTotalNodes,doInteger);
//
//		optimizer.optimize();
//		
//		return optimizer;
//
//	}

	@Override
	public void stopSlaves(){
		try{
		if(MPItoThread.Rank() == 0 && KSParser.numProc >= 3){
			int[] nothing = new int[1];
			for (int curProc=1; curProc<KSParser.numProc; curProc++){ 
				try {
					MPItoThread.Send(nothing, 0, 1, ThreadMessage.INT, curProc, KSParser.doneTag);
				} catch (MPIException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}}catch(Exception E){
			System.out.println("Can't stop slaves for minAS");
		}
	}


}
