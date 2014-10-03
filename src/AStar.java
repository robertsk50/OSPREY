import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;


public abstract class AStar {

	EPICSettings es;
	DegreeOfFreedom[] DOFList;
	CETMatrix CETM;
	
	//number of residues under consideration
	protected int numTreeLevels;
	
	//the offset in the array index for each level
//	protected int nodeIndexOffset[] = null;
	protected Index3 twoDTo3D[][];
	
//	Index3 rotIndexes[];
	protected HashMap<ArrayList<Index3>,EnergyTuple> energyTuples;
	
	double retConfE;
	int numExpanded = 0;
	
	//the leaf nodes visible from the expansions
	protected PGExpansionQueue curExpansion;

	PrintStream outPS = System.out;
	
	abstract PGQueueNode doAStar(boolean run1);

	public void stopSlaves() {
		//Do nothing by default
	}
	
	public void setEnergyTuples(HashMap<ArrayList<Index3>, EnergyTuple> eTups) {
		energyTuples = eTups;
	}
	
	void addNodeBack(PGQueueNode node){
		curExpansion.insert(node);
	}
	
	protected ArrayList<LinkedList<EnergyTuple>> findTupleOptions(PGQueueNode node) {
		boolean[] excludeLevel = new boolean[numTreeLevels];
		Index3[] curConf = new Index3[numTreeLevels];
		for(int i=0; i<excludeLevel.length;i++){
			excludeLevel[i] = false;
			curConf[i] = null;
		}

		int index1;
		//Setup CurConf
		for (int curALindex=0; curALindex< node.nonEmptyLevels.size(); curALindex++){ //compute using the formula
			int curLevel = node.nonEmptyLevels.get(curALindex);
//			index1 = nodeIndexOffset[curLevel] + node.confSoFar[curLevel];//index of r at i
			//Find Tuples
			Index3 i1 = twoDTo3D[curLevel][node.confSoFar[curLevel]];//rotIndexes[index1];
			assert curLevel == i1.pos;
			curConf[curLevel] = i1;
		}

		ArrayList<LinkedList<EnergyTuple>> comboTuples = new ArrayList<LinkedList<EnergyTuple>>();

		if(energyTuples == null || energyTuples.size() ==0)
			return comboTuples;
		//ArrayList<EnergyTuple> tuplesThisRound = new ArrayList<EnergyTuple>();
		//ArrayList<EnergyTuple> tuplesToAddTo = null;

		//EnergyTuplewVal curTuple = null;
		ArrayList<Index3> rotPair = new ArrayList<Index3>(2);
		rotPair.add(null);
		rotPair.add(null);
		//Go through all pairs
		for(int i=0; i<curConf.length;i++){
			if(curConf[i] != null){
				rotPair.set(0, curConf[i]);
				for(int j=i+1;j<curConf.length;j++){
					if(curConf[j] != null){
						rotPair.set(1, curConf[j]);
						if(energyTuples.containsKey(rotPair)){
							ArrayList<EnergyTuple> tmpTuples = new ArrayList<EnergyTuple>();
							tmpTuples = allTuples(energyTuples.get(rotPair), tmpTuples, curConf, excludeLevel, node);
							//Add Every Tuple found as it's own ArrayList
							for(EnergyTuple tup: tmpTuples){
								LinkedList<EnergyTuple> tmpList = new LinkedList<EnergyTuple>();
								tmpList.add(tup);
								comboTuples.add(tmpList);
							}
						}
					}
				}
			}
		}

		//Now we have seed tuples and find if any other tuples fit for the ones I have
		int initSize = 0;
		if(comboTuples.size() > 0){
			do{
				initSize = comboTuples.size();

				ArrayList<LinkedList<EnergyTuple>> toRemove = new ArrayList<LinkedList<EnergyTuple>>();
				ArrayList<LinkedList<EnergyTuple>> toAdd = new ArrayList<LinkedList<EnergyTuple>>();

				for(LinkedList<EnergyTuple> parents: comboTuples){
					boolean[] tmpExclLevel = new boolean[excludeLevel.length];
					System.arraycopy(excludeLevel, 0, tmpExclLevel, 0, excludeLevel.length);
					for(EnergyTuple parent: parents){
						for(int i=0; i<parent.rots.length;i++){
							tmpExclLevel[parent.rots[i].pos] = true;
						}}

					rotPair = new ArrayList<Index3>(2);
					rotPair.add(null);
					rotPair.add(null);
					//Go through all pairs
					for(int i=0; i<curConf.length;i++){
						if(curConf[i] != null){
							rotPair.set(0, curConf[i]);
							for(int j=i+1;j<curConf.length;j++){
								if(curConf[j] != null){
									rotPair.set(1, curConf[j]);
									if(energyTuples.containsKey(rotPair)){
										ArrayList<EnergyTuple> tmpTuples = new ArrayList<EnergyTuple>();
										tmpTuples = allTuples(energyTuples.get(rotPair), tmpTuples, curConf, tmpExclLevel, node);
										//Add Every Tuple found as it's own ArrayList
										if(tmpTuples.size() > 0){
											toRemove.add(parents);
											for(EnergyTuple tup: tmpTuples){
												LinkedList<EnergyTuple> tmpList = new LinkedList<EnergyTuple>();
												for(EnergyTuple parent: parents)
													tmpList.add(parent);
												tmpList.add(tup);
												toAdd.add(tmpList);
											}
										}
									}
								}
							}
						}
					}



				}

				for(LinkedList<EnergyTuple> parent: toRemove){
					comboTuples.remove(parent);
				}
				for(LinkedList<EnergyTuple> child: toAdd){
					comboTuples.add(child);
				}


			}while(comboTuples.size() > initSize);
		}
		return comboTuples;
	}
	
	protected void initialize2Dto3D(Emat emat, int numNodesForLevel[]) {
		twoDTo3D = new Index3[numTreeLevels][];
		SinglesIterator iter = emat.singlesIterator();
		int ctr=0;
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(twoDTo3D[emeWI.pos1()] == null){
				twoDTo3D[emeWI.pos1()] = new Index3[numNodesForLevel[emeWI.pos1()]];
				ctr=0;
			}
			
			if(!emeWI.eme.isPruned()){
				twoDTo3D[emeWI.pos1()][ctr] = new Index3(emeWI.rot1index());
				ctr++;
			}
		}
	}
	
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
	
	//Get the actual numbers for the rotamers of a conformation that is returned by A* by including the information
	//		for the pruned rotamers
	protected EMatrixEntryWIndex [] getActualConf(int curConf[], Emat emat){
		EMatrixEntryWIndex[] retConf = new EMatrixEntryWIndex[curConf.length];

		for (int curLevel=0; curLevel<curConf.length; curLevel++){
			Index3 rot = twoDTo3D[curLevel][curConf[curLevel]];
			retConf[curLevel] = new EMatrixEntryWIndex(emat.singles.getTerm(rot),rot);//arpMatrixRed[posIndex][arpMatrixRed[posIndex].length-1]; //intra entries stored in last column;		
		}
		return retConf;
	}
	
			
}
