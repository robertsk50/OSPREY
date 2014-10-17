import java.util.ArrayList;
import java.util.Arrays;
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

import cern.colt.matrix.DoubleMatrix1D;
import mpi.MPIException;

/**
 * Written by Pablo Gainza (2004-2012)
 */

/**
 * Uses A* search for single or multiple mutation sequences simultaneously to return the minimum-energy conformation; 
 * 		each consecutive run returns the next lowest-energy conformation, in order.
 * 
 */
public class PGgurobiAStarBySeq extends AStar{

	//number of residues under consideration
//	private int numTreeLevels;

	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numNodesForLevel[] = null;

	//the total number of possible rotamers for the given mutation
	private int numTotalNodes;

	private int numSeqForLevel[] = null;
	private int seqIndexOffset[][] = null;
	private int[][] numRotRemainingBySeq = null;
	private int[][] seqIndicesPerLevel = null;
	//the offset in the array index for each level
//	private int nodeIndexOffset[] = null;

	//the current sequence: the number of the corresponding rotamer for each level assigned so far 
	private int curConf[] = null;

	//the reduced min pairwise energy matrix
//	private EMatrixEntryWIndex pairwiseMinEnergyMatrix [][] = null;
	
	Emat emat;
	EMatrixEntryWIndex[] bestConf; 

	int topConf = 0;

	ArrayList<ArrayList<EnergyTuple>> tuplesPerPos;


	int topL = 0;
	int numTopL = 0;
	int numFS = 0;
	
	int numRunning = 0; 
	double upperE;

	Settings.VARIABLEORDER variableOrder;

	//Array used to store order for med cost function ordering
	private Integer dom_cmed_order[] = null;

	//array to hold sorted order of number of nodes per level
	//Sorted in descending order
	private Integer sortedIndicesNumSeqsForLevel[] = null;
	
	//constructor
	/*
	 * We assume that the parameters supplied (energy and DEE information) have already been modified
	 * 		to consider only the residues and rotamers for the possible mutations; i.e., the matrices are of
	 * 		reduced size
	 */
	PGgurobiAStarBySeq (int treeLevels, int numRotForRes[], Emat arpMatrix, double upperE,Settings.VARIABLEORDER varOrder){
				
		emat = arpMatrix;
		this.upperE = upperE; 
		numTreeLevels = treeLevels;
		
		numNodesForLevel = new int [numTreeLevels];
//		nodeIndexOffset = new int [numTreeLevels];
		numTotalNodes = 0;//arpMatrix.length-1;

		
		for (int i=0; i<numTreeLevels; i++){
			numNodesForLevel[i] = numRotForRes[i];
			numTotalNodes += numNodesForLevel[i];
		}
		
		numSeqForLevel = emat.numAANotPruned();
		seqIndicesPerLevel = emat.seqIndicesPerLevel(numSeqForLevel);
		numRotRemainingBySeq = emat.numRotRemainingBySeq();
		seqIndexOffset = new int[numTreeLevels][];
		
		for (int i=0; i<numSeqForLevel.length;i++){
			seqIndexOffset[i] = new int[numSeqForLevel[i]];
			seqIndexOffset[i][0] = 0;
			for(int j=1; j<seqIndexOffset[i].length;j++){
				seqIndexOffset[i][j] = seqIndexOffset[i][j-1]+numRotRemainingBySeq[i][j-1];
			}
		}

		//the min energy matrix: the last column contains the intra-energy for each rotamer; the last row
		//		contains the shell-residue energy for each rotamer
		
		
		variableOrder = varOrder;
		if(variableOrder.compareTo(Settings.VARIABLEORDER.DOM_CMED) == 0)
			initDOMCMEDorder();
		//Sort indices for numNodesForLevel
		ArrayIndexComparator comparator = new ArrayIndexComparator(numSeqForLevel);
		sortedIndicesNumSeqsForLevel = comparator.createIndexArray();
		//Sorts in descending order
		Arrays.sort(sortedIndicesNumSeqsForLevel,comparator);

		//the current expansion list
		curExpansion = new PGExpansionQueue();

		//the current conformation
		curConf = new int [numTreeLevels];
		for (int i=0; i<numTreeLevels; i++){
			curConf[i] = -1;
		}
		
		
		boolean useGurobi = false;
		boolean useWCSP = true;
		
		if(useWCSP || useGurobi){
		GurobiCalcInfo gci = new GurobiCalcInfo(arpMatrix,numNodesForLevel,
				numRotRemainingBySeq,seqIndexOffset,numTotalNodes,false,seqIndicesPerLevel,upperE);
		
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
	
	/**
	 * Set the variable order based on minimizing the term |x_i|/med(cost functions of x_i)
	 * where x_i is the variable for residue position i.
	 */
	private void initDOMCMEDorder() {
		double minMedianCost = Double.POSITIVE_INFINITY;
		double[][] medianCosts = new double[numNodesForLevel.length][numNodesForLevel.length];
		for(int pos1=0; pos1<medianCosts.length;pos1++){
			for(int pos2=pos1; pos2<medianCosts.length;pos2++){
				Iterator<EMatrixEntryWIndex> iter;
				if(pos1 == pos2){
					iter = emat.singlesIterator(pos1);
				}else{
					iter = emat.pairsIterator(pos1, pos2);
				}
				
				ArrayList<Double> costs = new ArrayList<Double>();
				while(iter.hasNext()){
					EMatrixEntryWIndex emeWI = iter.next();
					if(!emeWI.eme.isPruned())
						costs.add(emeWI.eme.minE());
				}
				
				//Calculate Median
				Collections.sort(costs);//costs.sort(doubleComparator); Java7 does not have sort function
				int size = costs.size();
				double median;
				if(size % 2 == 0)
					median = (costs.get(size/2)+costs.get((size/2)-1))/2;
				else
					median = costs.get(size/2);
				
				medianCosts[pos1][pos2] = median;
				medianCosts[pos2][pos1] = median;
				
				if(median < minMedianCost)
					minMedianCost = median;
				
			}	
		}
		
		//We need all of the energies to be positive for the division to work
		//So we add all energies by the min value
		for(int i=0; i<medianCosts.length;i++)
			for(int j=0; j<medianCosts[i].length;j++)
				medianCosts[i][j] -= minMedianCost;
		//Calculate the dom/sum of median costs for each variable
		double[] dom_cmed = new double[numNodesForLevel.length];
		for(int i=0; i<dom_cmed.length;i++){
			double sumofcosts = 0;
			for(int j=0; j<medianCosts[i].length;j++)
				sumofcosts+= medianCosts[i][j];
			
			double val = (double)numNodesForLevel[i]/sumofcosts;
			dom_cmed[i] = val;
		}
		
		//Get the order
		ArrayIndexComparator comparator = new ArrayIndexComparator(dom_cmed);
		dom_cmed_order = comparator.createIndexArray();
		//Sorts in descending order
		Arrays.sort(dom_cmed_order,comparator);
			
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
			for(int i=0; i<numSeqForLevel.length;i++)
				if(numSeqForLevel[i] <= 2)
					numInitialNodes *= numSeqForLevel[i];
			
			int[][] confs = new int[numInitialNodes][];
			makeConfs(0,numSeqForLevel,numTreeLevels,new ArrayList<Integer>(),confs); 
			
			
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

			if (expNode==null || expNode.fScore > upperE){//the queue is empty
				return null; //so return a sequence of -1's to flag the stop of the search
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
					
					
					boolean foundBetter = false;
					try{
						while(numRunning > 0){
							PGQueueNode[] node = new PGQueueNode[1];
							Object s = MPItoThread.Recv(node, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);
							curExpansion.insert(node[0]);
							countNodes++;
							if(node[0].fScore < expNode.fScore)
								foundBetter = true;
							numRunning--;
						}
					}catch(Exception E){
						E.printStackTrace();
					}

					retConfE = expNode.fScore;
//					GurobiOptimization gOpt = gurobiFscoreOptimizer(expNode,true);
//					double finalScore = gOpt.getObjVal();

					WCSPOptimization optimizer = new WCSPOptimization(expNode,emat,seqIndicesPerLevel,numRotRemainingBySeq,
							numNodesForLevel,upperE);
					double finalScore = optimizer.optimize(null);
					optimizer.cleanUp();
					
					
					if(finalScore != retConfE || foundBetter){
						System.out.println("Fscore: "+retConfE +"  Gscore: "+finalScore);
						expNode.fScore = finalScore;
						curExpansion.insert(expNode);
						
					}
					else{
						//Convert seq to conf
//						expNode.confSoFar = gOpt.getConfFromSeq(expNode.confSoFar,seqIndexOffset);
						
						bestConf = optimizer.bestConf.conf;  
						done = true;
					}
				}
				else {// Queue not empty, we can continue.

					//					PGQueueNode[] nextLevelNodes = pickNextLevel(expNode);
					
					PGQueueNode[] nextLevelNodes=null;
					//Choose the next variable (residue position to expand)
					switch(variableOrder){
						case MINDOM: //choose the variable with the min domain size
							nextLevelNodes = pickNextLevelByDom(expNode,true);
							break;
						case MAXDOM: //choose the variable with the max domain size
							nextLevelNodes = pickNextLevelByDom(expNode,false);
							break;
						case DOM_CMED: //choose the variable that minimizes |dom|/median(cost functions)
							nextLevelNodes = pickNextLevelByDomCmed(expNode);
							break;
						case MINFSCORE: //choose the variable that has the max min fscore
							nextLevelNodes = pickNextLevelByMinFscore(expNode);
							break;
						case MEDFSCORE: //choose the variable that has the max med fscore
							nextLevelNodes = pickNextLevelByMedFscore(expNode);
							break;
						case BYTUPLE:
							nextLevelNodes = pickNextLevelByTuple(expNode);
							break;
						case SEQUENTIAL: //no ordering
						default:
							nextLevelNodes = pickNextLevelSequential(expNode);
							break;
					}
					for(int rot = 0; rot < nextLevelNodes.length; rot++){

						//Recompute better bound when inserting
						try{
							if(MPItoThread.Rank() == 0 && KSParser.numProc >= 3){
								int curProc = numRunning+1;
								if(numRunning == KSParser.numProc-1){
									PGQueueNode[] node = new PGQueueNode[1];
									Object s = MPItoThread.Recv(node, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);
									curExpansion.insert(node[0]);	
									countNodes++;
//									nextLevelNodes[node[0].nodeNum] = node[0];
									curProc = MPItoThread.getStatusSource(s);
									numRunning--;
								}
								PGQueueNode[] nodeToSend = new PGQueueNode[1];
								nodeToSend[0] = nextLevelNodes[rot];
								MPItoThread.Send(nodeToSend, 0, 1, ThreadMessage.OBJECT, curProc, KSParser.regTag);
								numRunning++;
									
							}else{
								nextLevelNodes[rot].fScore = wcspFscore(nextLevelNodes[rot]);
								curExpansion.insert(nextLevelNodes[rot]);	
								countNodes++;
							}
						}catch(Exception E){
							E.printStackTrace();
						}
						/*if(newNode.fScore - dequeuedNode.fScore < -0.3 && dequeuedNode.fScore < 0){
							System.out.println("DELETE ME");
						}*/
							

					}			
				
					
					
					
						
//					for(int rot = 0; rot < nextLevelNodes.length; rot++){
//						
//						if(expNode.fScore != 0 && nextLevelNodes[rot].fScore - expNode.fScore < -0.1)
//							System.out.println("Something went wrong with fScores");
//
//						curExpansion.insert(nextLevelNodes[rot]);	
//						countNodes++;
//					}

					//delete the expanded node from the queue, since it has already contributed
					//curExpansion.delete(expNode);
				}	
			}
		}
		

		//Get the emat indices for the current conformation to return
		expNode.actualConf = bestConf;

		System.out.println("Number of A* nodes inserted in the queue: "+countNodes);
				
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

//	private PGQueueNode[] pickNextLevel(PGQueueNode dequeuedNode){		
//		ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>();   
//		
//		int maxMinFScoreIndex = -1;
//		double maxMinFScore = Double.NEGATIVE_INFINITY;
//		int numRunning = 0;
//		
//		int numLevels = dequeuedNode.emptyLevels.size();
//		
//		if(!reorder)
//			numLevels = 1;
//		
//		// Iterate through all positions that have not been assigned in dequeuedNode.
//		for(int iter = 0; iter < numLevels; iter++){
//			int pos = dequeuedNode.emptyLevels.get(iter);
//			double minAtThisLevel = Double.POSITIVE_INFINITY;
//			int minIndexAtThisLevel = -1; 
//			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
//				
//			// And for each rotamer at that level
//			for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
//				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
//				//KER: Trying out gCompute as a heuristic
//				//newNode.fScore = gCompute(newNode);//fCompute(newNode);
//				newNode.fScore = fCompute(newNode);
//				curLevelPGQueueNodes[rot] = newNode;
//					
//				if(newNode.fScore < minAtThisLevel){
//					minAtThisLevel = newNode.fScore;
//					minIndexAtThisLevel = rot;
//				}
//			}			
//			allLevelNodes.add(curLevelPGQueueNodes);
//
//			if(minAtThisLevel > maxMinFScore){
//				maxMinFScore = minAtThisLevel;
//				maxMinFScoreIndex = iter;
//			}		
//		
//		}
//				
//		// For now return the first level.  This should work the same as the old A*.... 
//		if(!dequeuedNode.emptyLevels.isEmpty()){
//			return allLevelNodes.get(maxMinFScoreIndex);			
//		}
//		else{
//			return null;
//		}
//	}


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
	private double wcspFscore(PGQueueNode node) {
		WCSPOptimization optimizer = new WCSPOptimization(node,emat,seqIndicesPerLevel,numRotRemainingBySeq,numNodesForLevel,upperE);
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

	
	//////////////////////////////////////////////////////////////////////////
	// Now h(x) 
	//////////////////////////////////////////////////////////////////////////
	//Compute the h(n) score for the new node expanded by expNode
	//		called by doAStar(.)
	// dLevel is the current level that is being evaluated.
	private double fCompute (PGQueueNode node){

		double hn = 0.0f;

		for (int curALindex=0 ;curALindex<node.confSoFar.length;curALindex++){
			// For every level after the current one, we calculate the "heuristic" at that level
			hn += EnergyAtLevel(node, curALindex);
		}

		return hn;
	}
	//Called by hCompute(.)
	//  dLevel is called topLevel here.  
	//  At each level we compute the intra energy, the shell energy, the pairwise energy with respect to things that have already been assigned (e.g. interaction energies with nodes <= dLevel)
	//		and the energies with anything further ahead in the tree.
	private double EnergyAtLevel(PGQueueNode node, int curALindex){

		double minE = (double)Math.pow(10,30);
		double curE;		
		Index3 index1;

		double minShellResE;
		double minIndVoxE;			//formula term 1
		double sumMinPairE;			//formula term 2
		double sumMinMinPairE;		//formula term 3


		int curLevel = curALindex;// node.emptyLevels.get(curALindex);

		int startAA;
		int endAA;
		if(node.confSoFar[curLevel] >= 0){
			startAA = seqIndicesPerLevel[curLevel][node.confSoFar[curLevel]];
			endAA = seqIndicesPerLevel[curLevel][node.confSoFar[curLevel]];
		}else{
			startAA = 0;
			endAA = emat.singles.E[curLevel].length-1;
		}
		
		for (int a1=startAA; a1<=endAA;a1++){		//the rotamers at j
			for(int r1=0; r1< emat.singles.E[curLevel][a1].length;r1++){
				index1 = new Index3(curLevel,a1,r1);//twoDTo3D[curLevel][i1];	//the index of s at j
	
				//minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);//pairwiseMinEnergyMatrix[numTotalNodes][index1];//the shell-residue E is in the last row
				
				minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();//the intra-energy is in the last column
				
				
//				sumMinPairE = hSumMinPVE (node, index1);
				sumMinMinPairE = sumMinMinPVE(node, curALindex, index1);
	
				curE = minIndVoxE + sumMinMinPairE;
				if (curE<minE)		//compare to the min energy found so far
					minE = curE;
			}
		}

		return minE;

	}

	//Called by EnergyAtLevel(.)
	//  Here we compute the energy with the pairwise energy with respect to things that have already been assigned so far at this queue node (e.g. interaction energies with nodes inside node.nonEmptyNodes)
//	private double hSumMinPVE (PGQueueNode node, Index3 index1){
//
//		double sum = 0.0f;
//		Index3 index2;
//
//		
//		for (int levelALindex=0; levelALindex<node.nonEmptyLevels.size(); levelALindex++){
//			int level = node.nonEmptyLevels.get(levelALindex);
//			if(!excludeLevel[level]){
//				index2 = twoDTo3D[level][node.confSoFar[level]];//nodeIndexOffset[level] + node.confSoFar[level]; //the index of r at i
//
//				if(emat.areNeighbors(index1.pos, index2.pos))//[index2][index1] != null) //This should only happen when the two positions are not neighbors
//					sum += emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index2][index1].eme.minE(); //the pairwise energy between the two nodes
//			}
//		}
//
//		return sum;
//	}

	//Called by EnergyAtLevel(.)
	//  Here we compute the energy with the pairwise energy with respect to things that have not been assigned at this queue node (e.g. interaction energies with nodes <= dLevel)
	private double sumMinMinPVE(PGQueueNode node , int jALindex, Index3 index1){

		double sum = 0.0f;


		for (int curALindex=jALindex+1; curALindex<node.confSoFar.length; curALindex++){
			sum += indMinMinPVE(node, curALindex, index1);
		}

		return sum;
	}

	//Called by sumMinMinPVE(.)
	private double indMinMinPVE (PGQueueNode node, int kALindex, Index3 index1){

		double minEn = (double)Math.pow(10,30);
		double curEn;
		Index3 secondIndex;
		int kLevel = kALindex;//node.emptyLevels.get(kALindex);

		if(!emat.areNeighbors(index1.pos, kLevel))//pairwiseMinEnergyMatrix[index1][secondIndex] == null) //These positions aren't neighbors so return 0
			return 0.0;
		
		int startAA;
		int endAA;
		if(node.confSoFar[kLevel] >= 0){
			startAA = seqIndicesPerLevel[kLevel][node.confSoFar[kLevel]];
			endAA = seqIndicesPerLevel[kLevel][node.confSoFar[kLevel]];
		}else{
			startAA = 0;
			endAA = emat.singles.E[kLevel].length-1;
		}
		
		
		for (int a2=startAA; a2<=endAA; a2++){ //u at k
			
			for(int r2=0; r2<emat.singles.E[kLevel][a2].length;r2++){

//				secondIndex = twoDTo3D[kLevel][i2];//nodeIndexOffset[kLevel]+i2;
				
				
				curEn = emat.pairs.E[index1.pos][index1.aa][index1.rot][kLevel][a2][r2];//emat.getPairMinE(index1, secondIndex);//pairwiseMinEnergyMatrix[index1][secondIndex].eme.minE();	
				if (curEn<minEn){
					minEn = curEn;
				}
			}
		}

		return minEn;
	}

	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByMinFscore(PGQueueNode dequeuedNode){		

		ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>();   

		int maxMinFScoreIndex = -1;
		double maxMinFScore = Double.NEGATIVE_INFINITY;

		// Iterate through all positions that have not been assigned in dequeuedNode.
		int iterationLength = dequeuedNode.emptyLevels.size();
		for(int iter = 0; iter < iterationLength; iter++){
			int pos = dequeuedNode.emptyLevels.get(iter);
			double minAtThisLevel = Double.POSITIVE_INFINITY;
			int minIndexAtThisLevel = -1; 
			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
			// And for each rotamer at that level
			for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
				newNode.fScore = fCompute(newNode);
				curLevelPGQueueNodes[rot] = newNode;	

				if(newNode.fScore < minAtThisLevel){
					minAtThisLevel = newNode.fScore;
					minIndexAtThisLevel = rot;
				}
			}			
			allLevelNodes.add(curLevelPGQueueNodes);

			if(minAtThisLevel > maxMinFScore){
				maxMinFScore = minAtThisLevel;
				maxMinFScoreIndex = iter;
			}			
		}
		// For now return the first level.  This should work the same as the old A*.... 
		if(!allLevelNodes.isEmpty()){
			return allLevelNodes.get(maxMinFScoreIndex);			
		}
		else{
			return null;
		}
	}

	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByMedFscore(PGQueueNode dequeuedNode){		

		ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>();   

		int maxMedFScoreIndex = -1;
		double maxMedFScore = Double.NEGATIVE_INFINITY;

		// Iterate through all positions that have not been assigned in dequeuedNode.
		int iterationLength = dequeuedNode.emptyLevels.size();
		for(int iter = 0; iter < iterationLength; iter++){
			int pos = dequeuedNode.emptyLevels.get(iter);
			double[] fscoresAtThisLevel = new double[numSeqForLevel[pos]]; 
			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
			// And for each rotamer at that level
			for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
				newNode.fScore = fCompute(newNode);
				curLevelPGQueueNodes[rot] = newNode;	
				fscoresAtThisLevel[rot] = newNode.fScore;
			}			
			allLevelNodes.add(curLevelPGQueueNodes);

			Arrays.sort(fscoresAtThisLevel);
			int size = numSeqForLevel[pos];
			double medianAtCurLevel;
			if(size % 2 == 0)
				medianAtCurLevel = (fscoresAtThisLevel[size/2]+fscoresAtThisLevel[(size/2)-1])/2;
			else
				medianAtCurLevel = fscoresAtThisLevel[size/2];
			
			if(medianAtCurLevel > maxMedFScore){
				maxMedFScore = medianAtCurLevel;
				maxMedFScoreIndex = iter;
			}			
		}
		// For now return the first level.  This should work the same as the old A*.... 
		if(!allLevelNodes.isEmpty()){
			return allLevelNodes.get(maxMedFScoreIndex);			
		}
		else{
			return null;
		}
	}

	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByDom(PGQueueNode dequeuedNode, boolean minFirst){		

		//Find the level to expand next
		int levelToExpand;
		if(minFirst){
			//Sorted in descending order so we want the start from the far right if empty
			levelToExpand = sortedIndicesNumSeqsForLevel[dequeuedNode.emptyLevels.size()-1]; 
		}
		else{
			//Sorted in descending order so we want the start from the far left for max if empty
			levelToExpand = sortedIndicesNumSeqsForLevel[dequeuedNode.nonEmptyLevels.size()];
		}
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
			PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
			newNode.fScore = fCompute(newNode);
			curLevelPGQueueNodes[rot] = newNode;	
		}			
			
		return curLevelPGQueueNodes;
	}

	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelSequential(PGQueueNode dequeuedNode){		

		//Find the level to expand next
		int levelToExpand=dequeuedNode.nonEmptyLevels.size();
		
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
			PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
			newNode.fScore = fCompute(newNode);
			curLevelPGQueueNodes[rot] = newNode;	
		}			
			
		return curLevelPGQueueNodes;
	}
	
	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByDomCmed(PGQueueNode dequeuedNode){		

		//Find the level to expand next
		int levelToExpand=dom_cmed_order[dequeuedNode.emptyLevels.size()-1];
		
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
			PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
			newNode.fScore = fCompute(newNode);
			curLevelPGQueueNodes[rot] = newNode;	
		}			
			
		return curLevelPGQueueNodes;
	}
	
	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByTuple(PGQueueNode dequeuedNode){		

		PGQueueNode[] bestLevelNodes = null;

		int maxMinFScoreIndex = -1;
		double maxMinFScore = Double.NEGATIVE_INFINITY;

		int maxTuples = 0;
		int maxTuplePos = -1;

		ArrayList<Index3> dequeuedNodeConf = new ArrayList<Index3>();
		for(int pos: dequeuedNode.nonEmptyLevels){
			Index3 index2 = twoDTo3D[pos][dequeuedNode.confSoFar[pos]];
			dequeuedNodeConf.add(index2);
		}


		// Iterate through all positions that have not been assigned in dequeuedNode.
		for(int iter = 0; iter < dequeuedNode.emptyLevels.size(); iter++){
			int pos = dequeuedNode.emptyLevels.get(iter);
			int numTuples = 0;

			//Find the number of tuples that share at least one rotamer 
			for(EnergyTuple tuple: tuplesPerPos.get(pos)){
				if(tuple.shareRots(dequeuedNodeConf)){
					numTuples++;
				}
			}


			if(numTuples > maxTuples){
				maxTuples = numTuples;
				maxTuplePos = pos;
			}
		}

		//If there isn't a spot with more tuples
		//Use the position with the most tuples
		if(maxTuplePos == -1){
			for(int pos: dequeuedNode.emptyLevels){
				if(tuplesPerPos.get(pos).size() > maxTuples){
					maxTuples = tuplesPerPos.get(pos).size();
					maxTuplePos = pos;
				}
			}
		}


		if(maxTuplePos >= 0){
			//Generate the nodes for the best tuple position
			bestLevelNodes = new PGQueueNode[numSeqForLevel[maxTuplePos]];
			for(int rot = 0; rot < numSeqForLevel[maxTuplePos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, maxTuplePos, rot);  // No energy yet
				newNode.fScore = fCompute(newNode);
				bestLevelNodes[rot] = newNode;
			}


		}else{//Finally, if we still haven't found a position do the normal check
			ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>(); 
			
			int iterationLength = dequeuedNode.emptyLevels.size();
			for(int iter = 0; iter < iterationLength; iter++){
				int pos = dequeuedNode.emptyLevels.get(iter);

				double minAtThisLevel = Double.POSITIVE_INFINITY;
				int minIndexAtThisLevel = -1; 
				PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numSeqForLevel[pos]];
				// And for each rotamer at that level
				for(int rot = 0; rot < numSeqForLevel[pos]; rot++){				
					PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
					//KER: Trying out gCompute as a heuristic
					//newNode.fScore = gCompute(newNode);//fCompute(newNode);
					newNode.fScore = fCompute(newNode);
					/*if(newNode.fScore - dequeuedNode.fScore < -0.3 && dequeuedNode.fScore < 0){
						System.out.println("DELETE ME");
					}*/
					curLevelPGQueueNodes[rot] = newNode;	

					if(newNode.fScore < minAtThisLevel){
						minAtThisLevel = newNode.fScore;
						minIndexAtThisLevel = rot;
					}
				}			
				allLevelNodes.add(curLevelPGQueueNodes);

				if(minAtThisLevel > maxMinFScore){
					maxMinFScore = minAtThisLevel;
					maxMinFScoreIndex = iter;
				}			
			}

			// For now return the first level.  This should work the same as the old A*....
			if(!allLevelNodes.isEmpty()){
				bestLevelNodes = allLevelNodes.get(maxMinFScoreIndex);			
			}

		}

		return bestLevelNodes;

	}



}
