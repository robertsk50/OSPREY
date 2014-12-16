import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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

/**
 * Written by Pablo Gainza (2004-2012) and Kyle Roberts (2013-2014)
 *
 */

/**
 * Uses A* search for single or multiple mutation sequences simultaneously to return the minimum-energy conformation; 
 * 		each consecutive run returns the next lowest-energy conformation, in order.
 * 
 */
public class PGgurobiAStar extends AStar{

	//number of residues under consideration
//	private int numTreeLevels;

	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numNodesForLevel[] = null;

	//array to hold sorted order of number of nodes per level
	//Sorted in descending order
	private Integer sortedIndicesNumNodesForLevel[] = null;
	
	//Array used to store order for med cost function ordering
	private Integer dom_cmed_order[] = null;
	
	//Array used to store order for hmean cost function static ordering
	private Integer hmean_static_order[] = null;
	
	//the total number of possible rotamers for the given mutation
	private int numTotalNodes;

	
	
	//the offset in the array index for each level
//	private int nodeIndexOffset[] = null;

	//the current sequence: the number of the corresponding rotamer for each level assigned so far 
	private int curConf[] = null;

	//the reduced min pairwise energy matrix
//	private EMatrixEntryWIndex pairwiseMinEnergyMatrix [][] = null;
	private Emat emat;

	int topConf = 0;

	ArrayList<ArrayList<EnergyTuple>> tuplesPerPos;


	int topL = 0;
	int numTopL = 0;

	enum HEURISTIC{
		AS,GUROBI,WCSP,MPLP;
	}
 
	HEURISTIC heuristic = null;
	Settings.VARIABLEORDER variableOrder = null;

	Molecule m;
	StrandRotamers[] strandRot;
	MutableResParams strandMut;
	EPICSettings es;
	CETMatrix CETM;
	boolean doPerturbations;

	private int numFS = 0; 
	

	MSMPLP mpLP = null;

	private int numPreexpandedLevels = 0;

	private int preExpandDomainSize = 2;
	

	//constructor
	/*
	 * We assume that the parameters supplied (energy and DEE information) have already been modified
	 * 		to consider only the residues and rotamers for the possible mutations; i.e., the matrices are of
	 * 		reduced size
	 */
	PGgurobiAStar (int treeLevels, int numRotForRes[], Emat emat, Settings.ASTARMETHOD asMethod, Settings.VARIABLEORDER varOrder,
			EPICSettings es, boolean doPerturbations, Molecule m, StrandRotamers[] src, MutableResParams strandMut,
			CETMatrix CETM){

		this.m = m;
		this.strandRot  = src;
		this.strandMut = strandMut;
		this.CETM = CETM;
		this.es = es;
		this.doPerturbations = doPerturbations;
		this.emat = emat;
		
		switch(asMethod){
		case ORIG:
			heuristic = HEURISTIC.AS;
			break;
		case ASORIG:
			heuristic = HEURISTIC.AS;
			break;
		case ASLP:
			heuristic = HEURISTIC.GUROBI;
			break;
		case ASWCSP:
			heuristic = HEURISTIC.WCSP;
			break;
		case ASMPLP:
			heuristic = HEURISTIC.MPLP;
			break;
		default:
			System.out.println("ASMethod not recognized for PGgurobiAStar");
			System.exit(0);
			break;
		}

		numTreeLevels = treeLevels;

		numNodesForLevel = new int [numTreeLevels];
		numTotalNodes = 0;

		for (int i=0; i<numTreeLevels; i++){
			numNodesForLevel[i] = numRotForRes[i];
			numTotalNodes += numNodesForLevel[i];
		}
		
		variableOrder = varOrder;
		if(variableOrder.compareTo(Settings.VARIABLEORDER.DOM_CMED) == 0)
			initDOMCMEDorder();
		if(variableOrder.compareTo(Settings.VARIABLEORDER.HMEANSTATIC) == 0)
			initHMeanOrder();

		//Sort indices for numNodesForLevel
		ArrayIndexComparator comparator = new ArrayIndexComparator(numNodesForLevel);
		sortedIndicesNumNodesForLevel = comparator.createIndexArray();
		//Sorts in descending order
		Arrays.sort(sortedIndicesNumNodesForLevel,comparator);
		

		
		

		twoDTo3D = new Index3[numTreeLevels][];
		SinglesIterator iter = emat.singlesIterator();
		int ctr=0;
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(!emeWI.eme.isPruned()){
			if(twoDTo3D[emeWI.pos1()] == null){
				twoDTo3D[emeWI.pos1()] = new Index3[numNodesForLevel[emeWI.pos1()]];
				ctr=0;
			}
			
			twoDTo3D[emeWI.pos1()][ctr] = new Index3(emeWI.rot1index());
			ctr++;
			}
		}

		if(heuristic == HEURISTIC.MPLP)
			mpLP = new MSMPLP(numTreeLevels, numNodesForLevel, twoDTo3D,emat);

		//the current expansion list
		curExpansion = new PGExpansionQueue();

		//the current conformation
		curConf = new int [numTreeLevels];
		for (int i=0; i<numTreeLevels; i++){
			curConf[i] = -1;
		}

		KSParser.metrics.initASMetrics(numTreeLevels);
		
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
			if(numNodesForLevel[i] <= preExpandDomainSize) //Force low domain nodes to be low because they will be preexpanded anyway
				val = 0;
			dom_cmed[i] = val;
		}
		
		//Get the order
		ArrayIndexComparator comparator = new ArrayIndexComparator(dom_cmed);
		dom_cmed_order = comparator.createIndexArray();
		//Sorts in descending order
		Arrays.sort(dom_cmed_order,comparator);
			
	}
	
	/**
	 * Set the variable order based on minimizing the term for all x_i= sum_{j}(hmean(cost functions of x_i, x_j))
	 * where x_i is the variable for residue position i.
	 */
	private void initHMeanOrder() {
		
		// The hmean score for each residue.  Residues with the highest scores should be explored first.
		double hmeanScore [] = new double[numNodesForLevel.length];

		// Compute an hmean score for every residue.
		for(int pos1=0; pos1<numNodesForLevel.length;pos1++){
			// All residues initially have an hmean score of zero.
			hmeanScore[pos1] =0.0;
			for(int pos2=pos1; pos2<hmeanScore.length;pos2++){
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
				
				//Calculate Hmean
				Collections.sort(costs);
				int size = costs.size();
				double minValue = costs.get(0);				
				// The hmean is zero automatically if there is only one element or if the second element is equal in cost to the first.
				double hmean = 0.0;
				if (size > 1 && costs.get(1) != minValue){
					for(int pair = 1; pair < size; pair++){
						hmean += 1.0/costs.get(pair);
					}
					// Invert this for the actual value of the hmean.
					hmean = 1.0/hmean;
				}
				hmeanScore[pos1] += hmean;
			}	
		}
		
		//Get the order
		ArrayIndexComparator comparator = new ArrayIndexComparator(hmeanScore);
		hmean_static_order = comparator.createIndexArray();
		//Sorts in descending order
		Arrays.sort(hmean_static_order,comparator);
			
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

		for (int i=0; i<numTreeLevels; i++){ //initialize for this run
			curConf[i] = -1;
		}

		long numNodesStart = curExpansion.totalNodes;
		
		if (run1) {//if this is the first run of A*, then we need to set-up the empty head node
			//Expand nodes that only have 1 or 2 possibilities so we don't waste our time on them
			//First count how many there are
			int numInitialNodes = 1;
			for(int i=0; i<numNodesForLevel.length;i++)
				if(numNodesForLevel[i] <= preExpandDomainSize ){
					numInitialNodes *= numNodesForLevel[i];
					numPreexpandedLevels++;
				}
			
			int[][] confs = new int[numInitialNodes][];
			makeConfs(0,numNodesForLevel,numTreeLevels,new ArrayList<Integer>(),confs); 
			
			
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

				if(es.useEPIC){
                    //We will only expand terms that have the FS term included
                    //until our lowest lower bound contains an FS term, we will
                    //compute this term for the lowest-bounded nodes we get
                    //and then throw them back in the queue (this is to minimize the number of
                    //times we have to compute the FS term)
                    while(!expNode.FSTermIncluded){
                        expNode.fScore += FSTerm(expNode);
                        expNode.FSTermIncluded = true;
                        curExpansion.insert(expNode);
                        expNode = curExpansion.getMin();
                    }
                    //Now expNode has an FS term so we can expand it
                }
				
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
					double finalScore = fCompute(expNode);

					if(finalScore != retConfE){
						System.out.println("Fscore: "+retConfE +"  Gscore: "+finalScore);
						expNode.fScore = finalScore;
						curExpansion.insert(expNode);
					}
					else{
						done = true;
					}
				}
				else {// Queue not empty, we can continue.

//					long start = System.currentTimeMillis();
//					long validateTime = 0;
					
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
						case HMEANSTATIC: // chose the next residue by the sum of hmeans of all outgoing edges 
							nextLevelNodes = pickNextLevelByHmeanStatic(expNode);
							break;
						//// The following criteria are "dynamic": 
						case MINFSCORE: //choose the variable that has the max min fscore
							nextLevelNodes = pickNextLevelByMinFscore(expNode);
							break;
						case MEDFSCORE: //choose the variable that has the max med fscore
							nextLevelNodes = pickNextLevelByMedFscore(expNode);
							break;
						case HMEANFSCORE: // choose the residue that has the highest harmonic mean (of all rotamers) with respect to the expanded node.
							nextLevelNodes = pickNextLevelByHmeanFscore(expNode);
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
						if(heuristic == HEURISTIC.GUROBI){
//							if(energyTuples.size() > 0)
//								nextLevelNodes[rot].fScore = gurobiTupleFscore(nextLevelNodes[rot], expNode);
//							else
							nextLevelNodes[rot].fScore = gurobiFscore(nextLevelNodes[rot]);
						}else if(heuristic == HEURISTIC.WCSP){
							nextLevelNodes[rot].fScore = wcspFscore(nextLevelNodes[rot],expNode);
						}else if(heuristic == HEURISTIC.MPLP){
							nextLevelNodes[rot].fScore = MPLPfscore(nextLevelNodes[rot]);
						}//Else the fScores already use the A* heuristic

						//Validate the current bound
//						long validStart = System.currentTimeMillis();
//						double actualBound = validate(nextLevelNodes[rot]);
//						KSParser.metrics.updateASMetrics(expNode.nonEmptyLevels.size(), actualBound, nextLevelNodes[rot].fScore);
//						long validEnd = System.currentTimeMillis();
//						validateTime += (validEnd - validStart);
						//End Validate
						
						if(expNode.fScore != 0 && nextLevelNodes[rot].fScore - expNode.fScore < -0.1)
							System.out.println("Something went wrong with fScores");

						curExpansion.insert(nextLevelNodes[rot]);	
						
					}

					
					long stop = System.currentTimeMillis();
//					KSParser.metrics.updateASTimes(expNode.nonEmptyLevels.size(), (stop-start - validateTime), nextLevelNodes.length);
				}	
			}
		}
		
		EMatrixEntryWIndex[] actualConf = getActualConf(expNode.confSoFar, emat);
		expNode.actualConf = actualConf;
		
		System.out.println("Number of A* nodes inserted in the queue: "+(curExpansion.totalNodes -numNodesStart));
		outPS.println("A* returning conformation; lower bound = "+expNode.fScore+" nodes expanded: "+numExpanded+" FS terms evaluated: "+numFS);
		
		return expNode;
	}


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

	private double validate(PGQueueNode node) {
		
		//Solve subproblem with WCSP
		WCSPOptimization optimizer = new WCSPOptimization(node,emat,numNodesForLevel,Double.POSITIVE_INFINITY,twoDTo3D);
		System.out.print(".");
		double E = optimizer.optimize(null);
		optimizer.cleanUp();
		
		return E;
	}

	//KER: Use gurobi LP to get a bound on the score so far
	private double gurobiFscore(PGQueueNode node) {
		int numThreads = 1;
		GurobiOptimization optimizer = new GurobiOptimization(node,emat,numNodesForLevel,twoDTo3D,numTotalNodes,numThreads);
		System.out.print(".");
		return optimizer.optimize();
	}

	//Get MPLP bound for current subproblem
	private double MPLPfscore(PGQueueNode node){
		return mpLP.optimizeEMPLP(node.confSoFar, 100); 
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
			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
			// And for each rotamer at that level
			for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
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
			double[] fscoresAtThisLevel = new double[numNodesForLevel[pos]]; 
			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
			// And for each rotamer at that level
			for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
				newNode.fScore = fCompute(newNode);
				curLevelPGQueueNodes[rot] = newNode;	
				fscoresAtThisLevel[rot] = newNode.fScore;
			}			
			allLevelNodes.add(curLevelPGQueueNodes);

			Arrays.sort(fscoresAtThisLevel);
			int size = numNodesForLevel[pos];
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
			levelToExpand = sortedIndicesNumNodesForLevel[dequeuedNode.emptyLevels.size()-1 ]; 
		}
		else{
			//Sorted in descending order so we want the start from the far left for max if empty
			levelToExpand = sortedIndicesNumNodesForLevel[dequeuedNode.nonEmptyLevels.size() - numPreexpandedLevels];
		}
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
			PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
			newNode.fScore = fCompute(newNode);
			curLevelPGQueueNodes[rot] = newNode;	
		}			
			
		return curLevelPGQueueNodes;
	}

	// Picks the next level by the sum of the harmonic means of residue pair energies.
	private PGQueueNode[] pickNextLevelByHmeanStatic(PGQueueNode dequeuedNode) {

		//Find the level to expand next
		//This assumes that all variables of domain size 1-2 that are preexpanded are to the right of this matrix.
		int levelToExpand=this.hmean_static_order[dequeuedNode.emptyLevels.size()-1];
		
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
			PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
			newNode.fScore = fCompute(newNode);
			curLevelPGQueueNodes[rot] = newNode;	
		}			
			
		return curLevelPGQueueNodes;
	}
	
	// Pick the next level dynamically by the highest harmonic mean fscore.
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelByHmeanFscore(PGQueueNode dequeuedNode) {
		ArrayList<PGQueueNode[]> allLevelNodes = new ArrayList<PGQueueNode[]>();

		// The index and value of the unassigned residue with the highest hmean fscore. 
		int maxHmeanFScoreIndex = -1;
		double maxHmeanFScore = Double.NEGATIVE_INFINITY;

		// Iterate through all residues that have not been assigned in the partial conformation of dequeuedNode.
		int iterationLength = dequeuedNode.emptyLevels.size();
		for(int iter = 0; iter < iterationLength; iter++){
			// pos: the actual residue index.
			int pos = dequeuedNode.emptyLevels.get(iter);
			// Initially  the hmean is zero; it will only be changed for this residue if the minimum fScore for all rotamers is higher than dequeueNode.fScore 
			double hmeanAtThisLevel = 0.0;
			// We store all expanded nodes.
			PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
			// We store the sum of the inverses of each fscore, which is equal to 1/hmean
			double inverseSumHmean = 0.0;
			boolean hmeanIsZero = false; // A flag, in case the minimum fScore is equal to dequeueNode.fScore
			// Go through each rotamer at this residue
			for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
				PGQueueNode newNode = new PGQueueNode(numTreeLevels, dequeuedNode.confSoFar, 0.0, pos, rot);  // No energy yet
				// Compute an fscore for pos 
				newNode.fScore = fCompute(newNode);
				curLevelPGQueueNodes[rot] = newNode;
				// To normalize the hmean, which requires only positive values, we substract the parent's fscore,
				if(newNode.fScore > dequeuedNode.fScore){
					inverseSumHmean += 1.0/(newNode.fScore - dequeuedNode.fScore);
				}
				else{ // If one of the values is equal to the fScore, then the hmean will be zero.
					hmeanIsZero = true;
				}
			}			
			allLevelNodes.add(curLevelPGQueueNodes);
			if(!hmeanIsZero){
				hmeanAtThisLevel = 1.0/inverseSumHmean;
			}

			// Store the highest hmean and its index.
			if(hmeanAtThisLevel > maxHmeanFScore){
				maxHmeanFScore = hmeanAtThisLevel;
				maxHmeanFScoreIndex = iter;
			}			
		}
		// For now return the first level.  This should work the same as the old A*.... 
		if(!allLevelNodes.isEmpty()){
			return allLevelNodes.get(maxHmeanFScoreIndex);			
		}
		else{
			return null;
		}
	}
	// PGC
	// dePGQueueNode is the node at the top of the queue.  This node can be expanded to any of the levels that 
	//	have not been assigned in its own confSofar array.  For each rotamer at each of the unassigned levels,
	// we create a node.  Then, we calculate its fScore.
	private PGQueueNode[] pickNextLevelSequential(PGQueueNode dequeuedNode){		

		//Find the level to expand next
		int levelToExpand=dequeuedNode.emptyLevels.get(0);
		
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
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
		//This assumes that all variables of domain size 1-2 that are preexpanded are to the right of this matrix.
		int levelToExpand=dom_cmed_order[dequeuedNode.emptyLevels.size()-1];
		
		
		int pos = levelToExpand;
		PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
		// And for each rotamer at that level
		for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
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
			bestLevelNodes = new PGQueueNode[numNodesForLevel[maxTuplePos]];
			for(int rot = 0; rot < numNodesForLevel[maxTuplePos]; rot++){				
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
				PGQueueNode curLevelPGQueueNodes[] = new PGQueueNode[numNodesForLevel[pos]];
				// And for each rotamer at that level
				for(int rot = 0; rot < numNodesForLevel[pos]; rot++){				
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


	//Updates and prints the state of the queue
	private void printState(PGQueueNode expNode){

		numExpanded++;

		if (expNode.nonEmptyLevels.size()>topL){
			topL = expNode.nonEmptyLevels.size();
			numTopL = 1;
		}
		else if (expNode.level+1==topL)
			numTopL++;


		if((numExpanded%50)==0){
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




	//////////////////////////////////////////////////////////////////////////
	// THe following is the new versions by PGC
	// First f(X) and g(X)  
	//////////////////////////////////////////////////////////////////////////
	// Having the emptyLevels precomputed can help with some speed.  However, I think emptyLevels should go as a member of the node class.
	private double fCompute(PGQueueNode node){
		ArrayList<LinkedList<EnergyTuple>> tupleOptions; 

		if(energyTuples != null && energyTuples.size() > 0)
			tupleOptions = findTupleOptions(node);
		else
			tupleOptions = new ArrayList<LinkedList<EnergyTuple>>();
		
		double maxScore = Double.NEGATIVE_INFINITY;

		double gScore,hScore,fScore;
		fScore = Double.NEGATIVE_INFINITY;
		LinkedList<EnergyTuple> maxTuples = node.curTuples;
		//if(tupleOptions.size() == 0){
		gScore = gCompute(node);
		hScore = hCompute(node);
		maxScore = gScore + hScore;
		//			fScore = gurobiFscore(node);
		//System.out.println(oldFscore+" "+fScore);
		//}
		//else{

		for(LinkedList<EnergyTuple> curTuples: tupleOptions){
			node.curTuples = curTuples;
			gScore = gCompute(node);
			hScore = hCompute(node);
			fScore = gScore+hScore;
			if(fScore > maxScore){
				maxScore = fScore;
				maxTuples = curTuples;
			}
		}

		//}

		fScore = maxScore;
		node.curTuples = maxTuples;


		return fScore;
	}

	//KER: Use gurobi LP to get a bound on the score so far
//	private double gurobiFscore(PGQueueNode node) {
//		
//		GurobiOptimization optimizer = new GurobiOptimization(node,emat,numNodesForLevel,twoDTo3D,numTotalNodes);
//		System.out.print(".");
//		return optimizer.optimize();
//
//	}
	
	//KER: Use wcsp to get a bound on the score so far
	private double wcspFscore(PGQueueNode node,PGQueueNode parent) {
		WCSPOptimization optimizer = new WCSPOptimization(node,emat,numNodesForLevel,Double.POSITIVE_INFINITY,twoDTo3D);
		System.out.print(".");
		//double E =optimizer.optimize(null);
		double E = optimizer.getBound(null);
		optimizer.cleanUp();
		if(E < parent.fScore)
			E = parent.fScore;
		return E;

	}

	//KER: Use gurobi LP to get a bound on the score so far
//	private double gurobiTupleFscore(PGQueueNode node, PGQueueNode expNode) {
//
//		GurobiOptimization optimizer = new GurobiOptimization(true,node,emat,numNodesForLevel,twoDTo3D,numTotalNodes,0);
//
//		double hScore = optimizer.optimize();
//
//		ArrayList<LinkedList<EnergyTuple>> tupleOptions = findTupleOptions(node);
//		double maxScore = Double.NEGATIVE_INFINITY;
//
//		double gScore,fScore;
//		fScore = Double.NEGATIVE_INFINITY;
//		LinkedList<EnergyTuple> maxTuples = node.curTuples;
//
//		gScore = gCompute(node);
//		maxScore = gScore + hScore;
//
//
//
//		if(tupleOptions.size() > 0){
//
//			for(LinkedList<EnergyTuple> curTuples: tupleOptions){
//				node.curTuples = curTuples;
//				gScore = gCompute(node);
//				fScore = gScore+hScore;
//				if(fScore > maxScore){
//					maxScore = fScore;
//					maxTuples = curTuples;
//				}
//			}
//
//		}
//		node.curTuples = maxTuples;
//		fScore = maxScore;
//
//
//
//		return fScore;
//
//	}

	private double gCompute (PGQueueNode node){
		double gn = 0.0f; 
		Index3 index1;

		double minShellResE;
		double minIndVoxE;		//formula term 1
		double sumMinPairE=0.0;		//formula term 2

		boolean[] excludeLevel = new boolean[numTreeLevels];
		Index3[] curConf = new Index3[numTreeLevels];
		for(int i=0; i<excludeLevel.length;i++){
			excludeLevel[i] = false;
			curConf[i] = null;
		}

		//Setup CurConf
		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			int curLevel = node.nonEmptyLevels.get(curALindex);
			index1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//index of r at i
			//Find Tuples
//			Index3 i1 = rotIndexes[index1];
//			assert curLevel == i1.pos;
			curConf[curLevel] = index1;
		}

		//		LinkedList<EnergyTuple> curTuples = new LinkedList<EnergyTuple>();
		//		int initSize;
		//		do{
		//			initSize = curTuples.size();
		//			EnergyTuplewVal curTuple = null;
		//			ArrayList<Index3> rotPair = new ArrayList<Index3>(2);
		//			rotPair.add(null);
		//			rotPair.add(null);
		//			//Go through all pairs
		//			for(int i=0; i<curConf.length;i++){
		//				if(curConf[i] != null){
		//					rotPair.set(0, curConf[i]);
		//					for(int j=i+1;j<curConf.length;j++){
		//						if(curConf[j] != null){
		//							rotPair.set(1, curConf[j]);
		//							if(energyTuples.containsKey(rotPair)){
		//								EnergyTuplewVal tmpTuple = largestTuple(energyTuples.get(rotPair), curConf, excludeLevel,node);
		//								if( (tmpTuple != null && curTuple == null) || (tmpTuple != null && tmpTuple.et.rots.length > curTuple.et.rots.length) )
		//									curTuple = tmpTuple;
		//								else if(tmpTuple != null && tmpTuple.et.rots.length == curTuple.et.rots.length && tmpTuple.intE > curTuple.intE)
		//									curTuple = tmpTuple;
		//							}
		//						}
		//					}
		//				}
		//			}
		//
		//			if(curTuple != null){
		//				curTuples.add(curTuple.et);
		//				for(int i=0; i<curTuple.et.rots.length;i++){
		//					excludeLevel[curTuple.et.rots[i].pos] = true;
		//				}
		//			}
		//
		//		}while(curTuples.size() > initSize);
		//
		//		//Find most relevant tuple
		//		//(We should be finding all tuples but that requires more complicated energy calculations
		//		/*LinkedList<EnergyTuple> curTuples = new LinkedList<EnergyTuple>();
		//		int initSize;
		//		do{
		//			initSize = curTuples.size();
		//			EnergyTuple curTuple = null;
		//			for (int curALindex=0; curALindex<node.nonEmptyLevels.size(); curALindex++){
		//				int curLevel = node.nonEmptyLevels.get(curALindex);
		//				if(!excludeLevel[curLevel]){
		//					index1 = nodeIndexOffset[curLevel] + node.confSoFar[curLevel];//index of r at i
		//					//Find Tuples
		//					LinkedList<Index3> rots = new LinkedList<Index3>();
		//					Index3 i1 = rotIndexes[index1];
		//					rots.add(i1);
		//					EnergyTuple tmpTuple = null;
		//					//Recurse on the levels until we find the best tuple
		//					tmpTuple = buildTuples(0,node,rots,excludeLevel);
		//
		//					if( (tmpTuple != null && curTuple == null) || (tmpTuple != null && tmpTuple.rots.length > curTuple.rots.length) )
		//						curTuple = tmpTuple;
		//				}
		//			}
		//			if(curTuple != null){
		//				curTuples.add(curTuple);
		//				for(int i=0; i<curTuple.rots.length;i++){
		//					excludeLevel[curTuple.rots[i].pos] = true;
		//				}
		//			}
		//
		//		}while(curTuples.size() > initSize);*/



		//if(curTuples.size() > 0){

		//node.curTuples = curTuples;

		if(node.curTuples != null){
			for(EnergyTuple curTuple:node.curTuples){

				int indexOffset = 0;
				int startPoint = 0;

				//Accumulate positions to exclude so we don't double count anything
				//boolean[] excludeTuple = new boolean[numTreeLevels];
				int[] curTupleOffset = new int[numTreeLevels];
				//for(int i=0; i<excludeTuple.length;i++)
				//	excludeTuple[i] = false;

				//Add Tuple Energies
				for(int i=0; i<curTuple.rots.length;i++){
					//excludeTuple[curTuple.rots[i].pos] = true;
					excludeLevel[curTuple.rots[i].pos] = true;
					for(int j=startPoint;j<=curTuple.rots[i].pos;j++)
						curTupleOffset[j] = indexOffset;
					if(i>0)
						indexOffset--;
					startPoint = curTuple.rots[i].pos+1;
				}
				for(int i=startPoint;i<curTupleOffset.length;i++)
					curTupleOffset[i] = indexOffset;

				//Intra E for tuple
				minIndVoxE = curTuple.intraE;
				//				String tupPrefSingle = "tup_s";
				//				String tupPrefPair = "tup_p";
				//				
				//				for(Index3 rot: curTuple.rots){
				//					tupPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
				//					tupPrefPair += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
				//				}
				//				System.out.println(tupPrefSingle+" "+curTuple.intraE);

				//Pair E for tuple
				int posOffset = 0;
				Index3 index2;
				boolean first = true;
				for(Index3 i3_1: curTuple.rots){
					int level_1 = i3_1.pos;
					Index3 index_1 = twoDTo3D[level_1][node.confSoFar[level_1]];	//s at j
					for (int curALindex=0; curALindex<node.nonEmptyLevels.size(); curALindex++){
						int level = node.nonEmptyLevels.get(curALindex);
	
						if(!excludeLevel[level]){
							index2 = twoDTo3D[level][node.confSoFar[level]];	//s at j
	//						Index3 i3 = rotIndexes[index2];
	//						sumMinPairE += curTuple.E[index2.pos+curTupleOffset[level]][index2.aa][index2.rot];//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
							sumMinPairE += emat.getPairMinE(index_1, index2);
							//System.out.println(tupPrefPair+"_"+(i3.pos+curTupleOffset[level])+"_"+i3.aa+"_"+i3.rot+" "+curTuple.E[i3.pos+curTupleOffset[level]][i3.aa][i3.rot]);
						}
					}
				}
				gn += minIndVoxE + sumMinPairE;
				sumMinPairE = 0;
				minIndVoxE = 0;
			}
		}

		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			minIndVoxE = 0;
			sumMinPairE = 0;

			int curLevel = node.nonEmptyLevels.get(curALindex);
			index1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//index of r at i

			//minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);
			if(!excludeLevel[curLevel]){
				minIndVoxE = emat.getSingleMinE(index1);// pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE(); //the intra-energy is in the last column
				//Index3 i1 = rotIndexes[index1];
				//System.out.println("s_"+i1.pos+"_"+i1.aa+"_"+i1.rot+" "+pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE());
				sumMinPairE = gSumMinPVE(node, curALindex+1, index1,excludeLevel);
			}

			gn += (minIndVoxE + sumMinPairE);
		}

		return gn;
	}


	private double gComputePrint (PGQueueNode node){
		double gn = 0.0f; 
		Index3 index1;

		double minShellResE;
		double minIndVoxE;		//formula term 1
		double sumMinPairE=0.0;		//formula term 2

		boolean[] excludeLevel = new boolean[numTreeLevels];
		Index3[] curConf = new Index3[numTreeLevels];
		for(int i=0; i<excludeLevel.length;i++){
			excludeLevel[i] = false;
			curConf[i] = null;
		}

		//Setup CurConf
		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			int curLevel = node.nonEmptyLevels.get(curALindex);
			index1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//index of r at i
			//Find Tuples
//			Index3 i1 = rotIndexes[index1];
//			assert curLevel == i1.pos;
			curConf[curLevel] = index1;
		}


		if(node.curTuples != null){
			for(EnergyTuple curTuple:node.curTuples){

				int indexOffset = 0;
				int startPoint = 0;

				//Accumulate positions to exclude so we don't double count anything
				//boolean[] excludeTuple = new boolean[numTreeLevels];
				int[] curTupleOffset = new int[numTreeLevels];
				//for(int i=0; i<excludeTuple.length;i++)
				//	excludeTuple[i] = false;

				//Add Tuple Energies
				for(int i=0; i<curTuple.rots.length;i++){
					//excludeTuple[curTuple.rots[i].pos] = true;
					excludeLevel[curTuple.rots[i].pos] = true;
					for(int j=startPoint;j<=curTuple.rots[i].pos;j++)
						curTupleOffset[j] = indexOffset;
					if(i>0)
						indexOffset--;
					startPoint = curTuple.rots[i].pos+1;
				}
				for(int i=startPoint;i<curTupleOffset.length;i++)
					curTupleOffset[i] = indexOffset;

				//Intra E for tuple
				minIndVoxE = curTuple.intraE;
				String tupPrefSingle = "tup_s";
				String tupPrefPair = "tup_p";

				for(Index3 rot: curTuple.rots){
					tupPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
					tupPrefPair += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
				}
				System.out.println(tupPrefSingle+" "+curTuple.intraE);

				//Pair E for tuple
				int posOffset = 0;
				Index3 index2;
				boolean first = true;
				for(Index3 i3_1: curTuple.rots){
					int level_1 = i3_1.pos;
					Index3 index_1 = twoDTo3D[level_1][node.confSoFar[level_1]];	//s at j
					for (int curALindex=0; curALindex<node.nonEmptyLevels.size(); curALindex++){
						int level = node.nonEmptyLevels.get(curALindex);
	
						if(!excludeLevel[level]){
							index2 = twoDTo3D[level][node.confSoFar[level]];	//s at j
	//						Index3 i3 = rotIndexes[index2];
	//						sumMinPairE += curTuple.E[index2.pos+curTupleOffset[level]][index2.aa][index2.rot];//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
							sumMinPairE += emat.getPairMinE(index_1, index2);
							System.out.println(tupPrefPair+"_"+(index2.pos)+"_"+index2.aa+"_"+index2.rot+" "+emat.getPairMinE(index_1, index2));
						}
					}
				}
				gn += minIndVoxE + sumMinPairE;
				sumMinPairE = 0;
				minIndVoxE = 0;
			}
		}

		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			minIndVoxE = 0;
			sumMinPairE = 0;

			int curLevel = node.nonEmptyLevels.get(curALindex);
			index1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//index of r at i

			//minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);
			if(!excludeLevel[curLevel]){
				minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE(); //the intra-energy is in the last column
//				Index3 i1 = rotIndexes[index1];
				System.out.println("s_"+index1.pos+"_"+index1.aa+"_"+index1.rot+" "+emat.getSingleMinE(index1));
				sumMinPairE = gSumMinPVEPrint(node, curALindex+1, index1,excludeLevel);
			}

			gn += (minIndVoxE + sumMinPairE);
		}

		return gn;
	}


	private double gSumMinPVE(PGQueueNode node, int startALindex, Index3 index1,boolean[] excludeLevel){
		Index3 index2;
		double sum = 0.0f;

		for (int levelALIndex = startALindex; levelALIndex < node.nonEmptyLevels.size(); levelALIndex++){
			int level = node.nonEmptyLevels.get(levelALIndex);
			if(!excludeLevel[level]){
				index2 = twoDTo3D[level][node.confSoFar[level]];//nodeIndexOffset[level] + node.confSoFar[level];	//s at j
				if(emat.areNeighbors(index1.pos, index2.pos))//pairwiseMinEnergyMatrix[index1][index2]!=null) //happens if they aren't neighbors
					sum += emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
				//				Index3 i1 = rotIndexes[index1];
				//				Index3 i2 = rotIndexes[index2];
				//				System.out.println("p_"+i1.pos+"_"+i1.aa+"_"+i1.rot+"_"+i2.pos+"_"+i2.aa+"_"+i2.rot+" "+pairwiseMinEnergyMatrix[index1][index2].eme.minE());
			}

		}

		return sum;

	}

	private double gSumMinPVEPrint(PGQueueNode node, int startALindex, Index3 index1,boolean[] excludeLevel){
		Index3 index2;
		double sum = 0.0f;

		for (int levelALIndex = startALindex; levelALIndex < node.nonEmptyLevels.size(); levelALIndex++){
			int level = node.nonEmptyLevels.get(levelALIndex);
			if(!excludeLevel[level]){
				index2 = twoDTo3D[level][node.confSoFar[level]];//nodeIndexOffset[level] + node.confSoFar[level];	//s at j
				if(emat.areNeighbors(index1.pos, index2.pos))//pairwiseMinEnergyMatrix[index1][index2]!=null) //happens if they aren't neighbors
					sum += emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
				
				System.out.println("p_"+index1.pos+"_"+index1.aa+"_"+index1.rot+"_"+
										index2.pos+"_"+index2.aa+"_"+index2.rot+" "+
										emat.getPairMinE(index1, index2));
			}

		}

		return sum;

	}

	/*private EnergyTuple buildTuples(int startLevel, PGQueueNode node,LinkedList<Index3> rots,boolean[] excludeLevel) {
		EnergyTuple retTuple = null;
		EnergyTuple curTuple = null;
		EnergyTuple tmpTuple = null;

		for(int curALindex=startLevel; curALindex<node.nonEmptyLevels.size();curALindex++){

			int curLevel = node.nonEmptyLevels.get(curALindex);
			if(!excludeLevel[curLevel]){
				int index = nodeIndexOffset[curLevel] + node.confSoFar[curLevel];//index of r at i
				Index3 i3 = rotIndexes[index];
				rots.add(i3);
				Collections.sort(rots);

				if(energyTuples.containsKey(rots)){
					curTuple = energyTuples.get(rots);

					if(retTuple == null)
						retTuple = curTuple;

					tmpTuple = buildTuples(startLevel,node,rots,excludeLevel);

					if(tmpTuple != null && tmpTuple.rots.length > retTuple.rots.length)
						retTuple = tmpTuple;
					if(curTuple != null && curTuple.rots.length > retTuple.rots.length)
						retTuple = curTuple;


				}
				rots.remove(i3);
			}

		}

		return retTuple;

	}*/


	//KER: Return largest(longest) tuple, or if there is a tie, the tuple with the largest energy contribution
	private ArrayList<EnergyTuple> allTuples(EnergyTuple parent, ArrayList<EnergyTuple> retTuples, Index3[] curConf, boolean[] excludeLevel,PGQueueNode node) {

		boolean goodParent = true;
		//double tupleEnergyContribution = Double.NEGATIVE_INFINITY;
		for(Index3 rot:parent.rots)
			if(excludeLevel[rot.pos]){
				goodParent = false;
				break;
			}

		//EnergyTuple retTuple = null;
		if(goodParent){
			retTuples.add(parent);
			//tupleEnergyContribution = tupleEnergyContribution(retTuple, node);
		}
		EnergyTuple curTuple = null;
		EnergyTuplewVal tmpTuple = null;

		for(EnergyTuple child:parent.children){
			//Go through each child to find best tuple that doesn't include any level in excludeLevel
			boolean goodChild = true;
			for(Index3 rot:child.rots){
				if(!rot.equals(curConf[rot.pos]))
					goodChild = false;
				if(excludeLevel[rot.pos])
					goodChild = false;

			}
			if(goodChild){
				//KER: If child is good replace parent with child
				retTuples.remove(parent);
				//retTuples.add(child);
				//curTuple = child;

				//				if(retTuple == null){
				//					retTuple = curTuple;
				//					//tupleEnergyContribution = tupleEnergyContribution(retTuple, node);
				//				}

				retTuples = allTuples(child,retTuples,curConf,excludeLevel,node);

				//				if(tmpTuple != null && tmpTuple.et.rots.length > retTuple.rots.length)
				//					retTuple = tmpTuple.et;
				//				else if (tmpTuple != null && tmpTuple.et.rots.length == retTuple.rots.length){
				//					double tmpTEC = tupleEnergyContribution(tmpTuple.et,node);
				//					if(tmpTEC > tupleEnergyContribution){
				//						retTuple = tmpTuple.et;
				//						tupleEnergyContribution = tmpTEC;
				//					}
				//				}
				//				if(curTuple != null && curTuple.rots.length > retTuple.rots.length)
				//					retTuple = curTuple;
				//				else if (curTuple != null && curTuple.rots.length == retTuple.rots.length){
				//					double tmpTEC = tupleEnergyContribution(curTuple,node);
				//					if(tmpTEC > tupleEnergyContribution){
				//						retTuple = curTuple;
				//						tupleEnergyContribution = tmpTEC;
				//					}
				//				}

			}

		}

		//		if(retTuple == null)
		//			return null;
		//		else
		//			return new EnergyTuplewVal(retTuple, tupleEnergyContribution);

		return retTuples;
	}



	//KER: Return largest(longest) tuple, or if there is a tie, the tuple with the largest energy contribution
	private EnergyTuplewVal largestTuple(EnergyTuple parent, Index3[] curConf, boolean[] excludeLevel,PGQueueNode node) {
		boolean goodParent = true;
		double tupleEnergyContribution = Double.NEGATIVE_INFINITY;
		for(Index3 rot:parent.rots)
			if(excludeLevel[rot.pos]){
				goodParent = false;
				break;
			}

		EnergyTuple retTuple = null;
		if(goodParent){
			retTuple = parent;
			tupleEnergyContribution = tupleEnergyContribution(retTuple, node);
		}
		EnergyTuple curTuple = null;
		EnergyTuplewVal tmpTuple = null;

		for(EnergyTuple child:parent.children){
			//Go thorough each child to find best tuple that doesn't include any level in excludeLevel
			boolean goodChild = true;
			for(Index3 rot:child.rots){
				if(!rot.equals(curConf[rot.pos]))
					goodChild = false;
				if(excludeLevel[rot.pos])
					goodChild = false;

			}
			if(goodChild){
				curTuple = child;

				if(retTuple == null){
					retTuple = curTuple;
					tupleEnergyContribution = tupleEnergyContribution(retTuple, node);
				}

				tmpTuple = largestTuple(child,curConf,excludeLevel,node);

				if(tmpTuple != null && tmpTuple.et.rots.length > retTuple.rots.length)
					retTuple = tmpTuple.et;
				else if (tmpTuple != null && tmpTuple.et.rots.length == retTuple.rots.length){
					double tmpTEC = tupleEnergyContribution(tmpTuple.et,node);
					if(tmpTEC > tupleEnergyContribution){
						retTuple = tmpTuple.et;
						tupleEnergyContribution = tmpTEC;
					}
				}
				if(curTuple != null && curTuple.rots.length > retTuple.rots.length)
					retTuple = curTuple;
				else if (curTuple != null && curTuple.rots.length == retTuple.rots.length){
					double tmpTEC = tupleEnergyContribution(curTuple,node);
					if(tmpTEC > tupleEnergyContribution){
						retTuple = curTuple;
						tupleEnergyContribution = tmpTEC;
					}
				}

			}

		}

		if(retTuple == null)
			return null;
		else
			return new EnergyTuplewVal(retTuple, tupleEnergyContribution);

	}

	private double tupleEnergyContribution(EnergyTuple curTuple,PGQueueNode node){

		double gn = 0;
		double minIndVoxE = 0;
		double sumMinPairE = 0;

		boolean[] excludeTuple = new boolean[numTreeLevels];
		int[] curTupleOffset = new int[numTreeLevels];
		for(int i=0; i<excludeTuple.length;i++)
			excludeTuple[i] = false;

		int indexOffset = 0;
		int startPoint = 0;

		//Add Tuple Energies
		for(int i=0; i<curTuple.rots.length;i++){
			excludeTuple[curTuple.rots[i].pos] = true;
			for(int j=startPoint;j<=curTuple.rots[i].pos;j++)
				curTupleOffset[j] = indexOffset;
			if(i>0)
				indexOffset--;
			startPoint = curTuple.rots[i].pos+1;
		}
		for(int i=startPoint;i<curTupleOffset.length;i++)
			curTupleOffset[i] = indexOffset;

		//Intra E for tuple
		minIndVoxE = curTuple.intraE;

		//Pair E for tuple
		int posOffset = 0;
		Index3 index2;
		boolean first = true;
		for(Index3 i3_1: curTuple.rots){
			int level_1 = i3_1.pos;
			Index3 index_1 = twoDTo3D[level_1][node.confSoFar[level_1]];	//s at j
			for (int curALindex=0; curALindex<node.nonEmptyLevels.size(); curALindex++){
				int level = node.nonEmptyLevels.get(curALindex);
	
				if(!excludeTuple[level]){
					index2 = twoDTo3D[level][node.confSoFar[level]];	//s at j
					//sumMinPairE += curTuple.E[index2.pos+curTupleOffset[level]][index2.aa][index2.rot];//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
					sumMinPairE += emat.getPairMinE(index_1, index2);
				}
			}
		}

		gn += minIndVoxE + sumMinPairE;


		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			minIndVoxE = 0;
			sumMinPairE = 0;

			int curLevel = node.nonEmptyLevels.get(curALindex);
			Index3 index1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//index of r at i

			//minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);
			if(!excludeTuple[curLevel]){
				minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE(); //the intra-energy is in the last column
				sumMinPairE = gSumMinPVE(node, curALindex+1, index1,excludeTuple);
			}

			gn += (minIndVoxE + sumMinPairE);
		}


		return gn;
	}


	//////////////////////////////////////////////////////////////////////////
	// Now h(x) 
	//////////////////////////////////////////////////////////////////////////
	//Compute the h(n) score for the new node expanded by expNode
	//		called by doAStar(.)
	// dLevel is the current level that is being evaluated.
	private double hCompute (PGQueueNode node){

		double hn = 0.0f;

		for (int curALindex=0 ;curALindex<node.emptyLevels.size();curALindex++){
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


		int curLevel = node.emptyLevels.get(curALindex);

		for (int i1=0; i1<numNodesForLevel[curLevel];i1++){		//the rotamers at j

			index1 = twoDTo3D[curLevel][i1];	//the index of s at j

			//minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);//pairwiseMinEnergyMatrix[numTotalNodes][index1];//the shell-residue E is in the last row
			minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();//the intra-energy is in the last column

			sumMinPairE = hSumMinPVE (node, index1);
			sumMinMinPairE = sumMinMinPVE(node, curALindex, index1);

			curE = minIndVoxE + sumMinPairE + sumMinMinPairE;
			if (curE<minE)		//compare to the min energy found so far
				minE = curE;
		}

		return minE;

	}

	//Called by EnergyAtLevel(.)
	//  Here we compute the energy with the pairwise energy with respect to things that have already been assigned so far at this queue node (e.g. interaction energies with nodes inside node.nonEmptyNodes)
	private double hSumMinPVE (PGQueueNode node, Index3 index1){

		double sum = 0.0f;
		Index3 index2;


		boolean[] excludeLevel = new boolean[numTreeLevels];
		for(int i=0; i<excludeLevel.length;i++)
			excludeLevel[i] = false;

		//Add energy related to tuple
//		if(node.curTuples != null && node.curTuples.size() > 0){
//
//
//			//Accumulate positions to exclude so we don't double count anything
//
//			int[] curTupleOffset = new int[numTreeLevels];
//
//			for(EnergyTuple curTuple:node.curTuples){
//
//				int indexOffset = 0;
//				int startPoint = 0;
//
//				//Add Tuple Energies
//				for(int i=0; i<curTuple.rots.length;i++){
//					excludeLevel[curTuple.rots[i].pos] = true;
//					for(int j=startPoint;j<=curTuple.rots[i].pos;j++)
//						curTupleOffset[j] = indexOffset;
//					if(i>0)
//						indexOffset--;
//					startPoint = curTuple.rots[i].pos+1;
//				}
//				for(int i=startPoint;i<curTupleOffset.length;i++)
//					curTupleOffset[i] = indexOffset;
//
//				//Pair E for tuple
//				int posOffset = 0;
//				boolean first = true;
////				Index3 i3 = rotIndexes[index1];
//				sum += curTuple.E[index1.pos+curTupleOffset[index1.pos]][index1.aa][index1.rot];//pairwiseMinEnergyMatrix[index1][index2].eme.minE();			
//			}
//		}


		for (int levelALindex=0; levelALindex<node.nonEmptyLevels.size(); levelALindex++){
			int level = node.nonEmptyLevels.get(levelALindex);
			if(!excludeLevel[level]){
				index2 = twoDTo3D[level][node.confSoFar[level]];//nodeIndexOffset[level] + node.confSoFar[level]; //the index of r at i

				if(emat.areNeighbors(index1.pos, index2.pos))//[index2][index1] != null) //This should only happen when the two positions are not neighbors
					sum += emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index2][index1].eme.minE(); //the pairwise energy between the two nodes
			}
		}

		return sum;
	}

	//Called by EnergyAtLevel(.)
	//  Here we compute the energy with the pairwise energy with respect to things that have not been assigned at this queue node (e.g. interaction energies with nodes <= dLevel)
	private double sumMinMinPVE(PGQueueNode node , int jALindex, Index3 index1){

		double sum = 0.0f;


		for (int curALindex=jALindex+1; curALindex<node.emptyLevels.size(); curALindex++){
			sum += indMinMinPVE(node, curALindex, index1);
		}

		return sum;
	}

	//Called by sumMinMinPVE(.)
	private double indMinMinPVE (PGQueueNode node, int kALindex, Index3 index1){

		double minEn = (double)Math.pow(10,30);
		double curEn;
		Index3 secondIndex;
		int kLevel = node.emptyLevels.get(kALindex);

		for (int i2=0; i2<numNodesForLevel[kLevel]; i2++){ //u at k

			secondIndex = twoDTo3D[kLevel][i2];//nodeIndexOffset[kLevel]+i2;
			if(!emat.areNeighbors(index1.pos, secondIndex.pos))//pairwiseMinEnergyMatrix[index1][secondIndex] == null) //These positions aren't neighbors so return 0
				return 0.0;

			curEn = emat.getPairMinE(index1, secondIndex);//pairwiseMinEnergyMatrix[index1][secondIndex].eme.minE();	
			if (curEn<minEn){
				minEn = curEn;
			}
		}

		return minEn;
	}

	@Override
	public void setEnergyTuples(HashMap<ArrayList<Index3>,EnergyTuple> eTups){
		energyTuples = eTups;
		tuplesPerPos = new ArrayList<ArrayList<EnergyTuple>>();

		for(int i=0; i<numTreeLevels;i++){
			tuplesPerPos.add(new ArrayList<EnergyTuple>());
		}

		for(EnergyTuple tuple: energyTuples.values()){
			for(Index3 i3 : tuple.rots){
				tuplesPerPos.get(i3.pos).add(tuple);
			}
		}

	}




	//////////////////////////////////////////////////////////////////////////
	// Old versions
	//////////////////////////////////////////////////////////////////////////
	/**	
	//Compute the h(n) score for the new node expanded by expNode
	//		called by doAStar(.)
	// dLevel is the current level that is being evaluated.
	private double hCompute (int dLevel, int conf[]){

		double hn = 0.0f;

		for (int curLevel=dLevel+1;curLevel<numTreeLevels;curLevel++){
			// For every level after the current one, we calculate the "heuristic" at that level
			hn += EnergyAtLevel(dLevel, curLevel, conf);
		}

		return hn;
	}

	//Called by hCompute(.)
	//  dLevel is called topLevel here.  
	//  At each level we compute the intra energy, the shell energy, the pairwise energy with respect to things that have already been assigned (e.g. interaction energies with nodes <= dLevel)
	//		and the energies with anything further ahead in the tree.
	private double EnergyAtLevel(int topLevel, int curLevel, int conf[]){

		double minE = (double)Math.pow(10,30);
		double curE;		
		int index1;

		double minShellResE;
		double minIndVoxE;			//formula term 1
		double sumMinPairE;			//formula term 2
		double sumMinMinPairE;		//formula term 3

		for (int i1=0; i1<numNodesForLevel[curLevel];i1++){		//the rotamers at j

			index1 = nodeIndexOffset[curLevel]+i1;	//the index of s at j

			minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);//pairwiseMinEnergyMatrix[numTotalNodes][index1];//the shell-residue E is in the last row
			minIndVoxE = pairwiseMinEnergyMatrix[index1][numTotalNodes];//the intra-energy is in the last column
			sumMinPairE = hSumMinPVE (topLevel, index1, conf);
			sumMinMinPairE = sumMinMinPVE(topLevel+1, curLevel, index1);

			curE = minShellResE + minIndVoxE + sumMinPairE + sumMinMinPairE;
			if (curE<minE)		//compare to the min energy found so far
				minE = curE;
		}

		return minE;

	}

	//Called by EnergyAtLevel(.)
	//  Here we compute the energy with the pairwise energy with respect to things that have already been assigned at this node (e.g. interaction energies with nodes <= dLevel)
	private double hSumMinPVE (int topLevel, int index1, int conf[]){

		double sum = 0.0f;
		int index2;

		for (int level=0; level<=topLevel; level++){

			index2 = nodeIndexOffset[level] + conf[level]; //the index of r at i

			sum += pairwiseMinEnergyMatrix[index2][index1]; //the pairwise energy between the two nodes
		}

		return sum;
	}

	//Called by EnergyAtLevel(.)
	private double sumMinMinPVE(int startLevel, int jLevel, int firstIndex){

		double sum = 0.0f;
		for (int level=jLevel+1; level<numTreeLevels; level++){
			sum += indMinMinPVE(level, firstIndex);
		}

		return sum;
	}

	//Called by sumMinMinPVE(.)
	private double indMinMinPVE (int kLevel, int firstIndex){

		double minEn = (double)Math.pow(10,30);
		double curEn;
		int secondIndex;

		for (int i2=0; i2<numNodesForLevel[kLevel]; i2++){ //u at k

			secondIndex = nodeIndexOffset[kLevel]+i2;
			curEn = pairwiseMinEnergyMatrix[firstIndex][secondIndex];			
			if (curEn<minEn){
				minEn = curEn;
			}
		}

		return minEn;
	}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

	//Compute the g(n) score for the new node expanded by expNode;
	//		called by doAStar(.)
	private double gCompute (int dLevel,int conf[]){

		double gn = 0.0f;
		int index1;

		double minShellResE;
		double minIndVoxE;		//formula term 1
		double sumMinPairE;		//formula term 2

		for (int curLevel=0; curLevel<=dLevel; curLevel++){//copute using the formula

			index1 = nodeIndexOffset[curLevel] + conf[curLevel];//index of r at i

			minShellResE = RotamerSearch.getReducedShellRotE(pairwiseMinEnergyMatrix,index1,numTotalNodes);
			minIndVoxE = pairwiseMinEnergyMatrix[index1][numTotalNodes];//the intra-energy is in the last column

			sumMinPairE = gSumMinPVE(dLevel, curLevel+1, index1, conf);

			gn += (minShellResE + minIndVoxE + sumMinPairE);
		}

		return gn;
	}

	//Called by gCompute(.)
	private double gSumMinPVE(int topLevel, int startLevel, int index1, int conf[]){

		int index2;
		double sum = 0.0f;

		for (int level=startLevel; level<=topLevel; level++){

			index2 = nodeIndexOffset[level] + conf[level];	//s at j			
			sum += pairwiseMinEnergyMatrix[index1][index2];

		}

		return sum;
	}
	 **/

	private double FSTerm(PGQueueNode expNode){
        //Get contribution to g-score from fit series
        //include h-score too if useHSer
        //optDOFVals is initially from the parent, to use an initial values; update to be optimal here

    	if(expNode.nonEmptyLevels.size() == 0)
    		return 0;
    	
        //if only want FS terms for full confs
        if( (expNode.level<numTreeLevels-1) && !(es.minPartialConfs) )
            return 0;

        CETObjFunction of = getNodeObjFunc(expNode);
        
        
        CCDMinimizer ccdMin = new CCDMinimizer(of,false);
        
        //TEMPORARILY DE-ACTIVATING OPTFSPOINT
        /*if(expNode.optFSPoint!=null && useHSer ){//we have initial values to use...
            //OUTFIT FOR CHANGING NUMBERS OF DOFS IN !USEHSER THOUGH...
            
            //if we're not getting much improvement by updating the FSTerm
            //even with the old initial values,
            //then the minimized value will not offer significant improvement over the current score
            //so don't spend time on it
            //dof_inh_lazy
            double initE = of.getValue( expNode.optFSPoint );
            if(useHSer && initE < expNode.fScore+0.1)
                return expNode.fScore;
            
            
            ccdMin.singleInitVal = expNode.optFSPoint;
        }*/
            

        DoubleMatrix1D optDOFs = ccdMin.minimize();
        if(optDOFs==null)//for ellipse bounds, if the highest ellipses don't all intersect together,
            //we can exclude the conformations involving them (set bound to inf)
            return Double.POSITIVE_INFINITY;
        
        double LSBE = of.getValue( optDOFs );
        System.out.println("LSBE: "+LSBE+" Fscore: "+expNode.fScore);
        
        if(es.useSVE){
            //cof minimized m...revert to unminimized state
            m.updateCoordinates();
            m.revertPertParamsToCurState();
        }
        
        
        //store initial values for next time
        //TEMPORARILY DE-ACTIVATING OPTFSPOINT
        /*if(useHSer){
            if(expNode.optFSPoint==null)
                expNode.optFSPoint = optDOFs;
            else
                expNode.optFSPoint.assign(optDOFs);
        }*/
        
        /*if(LSBE<-0.001){
            System.err.println("NEGATIVE VALUE ENCOUNTERED FOR LSBE: "+LSBE);
            System.err.println("Outputting LSBObjFunction to LSBOF.dat");
            KSParser.outputObject(of, "LSBOF.dat");
            System.err.println("optDOFs: "+optDOFs);
            System.exit(0);
        }*/
        
        numFS++;

        return LSBE;
    }

    
    
    CETObjFunction getNodeObjFunc(PGQueueNode node) {
        //not for splitBySlack
    	 int[] conf = node.confSoFar;
            
    	 Index3[] indices = new Index3[node.nonEmptyLevels.size()];
         //int rots[] = new int[node.nonEmptyLevels.size()];
    	 for(int i=0; i<node.nonEmptyLevels.size();i++){
    		 int level = node.nonEmptyLevels.get(i);
    		 if(conf[level]>=0){ //Should now always be greater than 0
             	indices[i] = twoDTo3D[level][node.confSoFar[level]];
    		 }
    	 }
    	 
    	 
//         for(int level : node.nonEmptyLevels){
//
//            if(conf[level]>=0){
//            	indices[level] = twoDTo3D[level][node.confSoFar[level]];
//            	Residue r = m.residue[strandMut.allMut[level]];
//            	int str = r.strandNumber;
//            	ResidueConformation rc = ((StrandRCs)strandRot[str]).rcl.getRC(emat.singles.getRot(indices)[0]);
//                //int redRot = pairwiseMinEnergyMatrix.resOffsets[level]+conf[level];
//                //long-form reduced index for rotamer at this level
//
////                AANums[level] = rc.rot.aaType.index;//pairwiseMinEnergyMatrix.indicesEMatrixAA[redRot];
////                rots[level] = rc.rot.rlIndex; //pairwiseMinEnergyMatrix.indicesEMatrixRot[redRot];
//            }
//            else{//special for not-fully-assigned levels
////                AANums[level] = -1;
////                rots[level] = conf[level];
//            }
//        
//        }
        
        //DEBUG!!!!  To use when checking minimization for a particular conformation
        /*
        AANums = new int[] {0, 1, 1, 20, 4, 1, 4};
        rots = new int[] {0, 1, 2, 0, 3, 0, 3};
        */
        
        
        ContSCObjFunction cof = null;//to set DOFs for SVE
            
        
        if(es.useSVE){//will need cof to set DOFs for sve
            //set AA types and rotamers up to now
            //as in RotamerSearch
            applyRotamers(node, indices);
            cof = new ContSCObjFunction(m,m.numberOfStrands,null,strandRot,false,getTransRotStrands(node));
        }
        
        
        
        return CETM.getObjFunc(indices,false,false,cof);//This can include h bounds if they're set up
        
    }
    
    
    boolean[] getTransRotStrands(PGQueueNode node){
        //when minimizing a partial conformation,
        //figure out which strands may be allowed to rotate and translate
        boolean transRotStrands[] = new boolean[m.numberOfStrands];
        
        //if a strand can rotate and translate and has either assigned or template residues,
        //let it rotate and translate
        
        for(int str=0; str<m.numberOfStrands; str++){
            if(m.strand[str].rotTrans){
                
                if( m.strand[str].numberOfResidues > strandMut.numMutPerStrand[str] )//strand has template residues
                    transRotStrands[str] = true;
                
                for(int i: node.nonEmptyLevels){
                    if(str==m.residue[strandMut.allMut[i]].strandNumber)//strand has an assigned residue
                        transRotStrands[str] = true;
                }
                
            }
        }
        
        return transRotStrands;
    }
    
    

    
    
    void applyRotamers(PGQueueNode node, Index3[] indices){
        //apply AA types and rotamers up to the current level
                	
        for (int i=0; i< node.nonEmptyLevels.size();i++){

            RotamerEntry re = emat.singles.getTerm(indices[i]);
            re.applyMutation(m, emat.resByPos, true,true );
			re.applyRC(emat.resByPos, m);
            re.flexible(m, emat.resByPos, true);
           
            
            
            //For Debugging, so delete for increased speed
            Residue r = m.residue[emat.resByPos.get(re.pos).get(0)];
            Rotamer rot = m.strand[r.strandNumber].rcl.getRC(re.r.rotamers[0]).rot;
            int aaInd1 = rot.aaType.index;
            int rotInd1 = rot.aaIndex; 
			System.out.print(re.pos+"_"+aaInd1+"_"+rotInd1+",");
			
            
            //fill in curAANum while we're at it
//            curAANum[m.strand[str].residue[strResNum].moleculeResidueNumber] = AANums[i];
          
            //rotamers
            
//            if(doPerturbations){
                //For DEEPer curRot is the the current RC
            	
//                boolean validRC = ((StrandRCs)strandRot[str]).applyRC(m, strResNum, rc);
//                if(!validRC)
//                    throw new RuntimeException("Error: invalid RC " + rots[i] + " at residue " + strResNum +
//                            " of strand " + str );
//            }
//            else if (strandRot[str].rl.getNumRotForAAtype(AANums[i])!=0)//not GLY or ALA
//                strandRot[str].applyRotamer(m, strResNum, rots[i]);
            
            //for gly or ala don't need to do anything
            
            //also make assigned residues flexible
//            m.strand[str].residue[strResNum].flexible = true;
        }
      
        //make the other residues not flexible
        for(int i: node.emptyLevels){
        	Residue r = m.residue[strandMut.allMut[i]];
            int str = r.strandNumber;
            int strResNum = r.strandResidueNumber;
            m.strand[str].residue[strResNum].flexible = false;
        }
    }
    
}
