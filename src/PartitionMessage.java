import java.io.Serializable;
import java.math.BigDecimal;


/*
 * Author: Kyle Roberts
 * This class is used to give intermediate information about a mutation search
 * to the head node, or for the head node to give intermediate information to
 * a mutation search.
 * 
 * Basically what this is used for is that when doing a drug design the protein
 * won't change so the partition function will be the same. This function is used
 * to relay information about the protein unbound partition function so that it 
 * only needs to be computed once instead of every time a mutation search is done.
 * 
 */

public class PartitionMessage implements Serializable {
	int runNum;
	int mutNum;
	int seqNum;
	BigDecimal q;
	double bestEMin;
	double bestE;
	int q_Time;
	boolean repeatEW;
	boolean allPruned;
	int searchNumConfsTotal;
	int searchNumConfsPrunedByE;
	int searchNumConfsPrunedByS;
	int searchNumConfsEvaluated;
	int searchNumConfsLeft;
	int searchNumPrunedMinDEE;
	double searchBestEnergyFound;
	
	
	PartitionMessage(int runNum, int mutNum, int seqNum,BigDecimal q, double bestEMin,double bestE,int q_Time,
	boolean repeatEW, boolean allPruned,int searchNumConfsTotal,
	int searchNumConfsPrunedByE,int searchNumConfsPrunedByS,int searchNumConfsEvaluated,
	int searchNumConfsLeft,int searchNumPrunedMinDEE,double searchBestEnergyFound){
		this.runNum = runNum;
		this.mutNum = mutNum;
		this.q = q;
		this.bestEMin = bestEMin;
		this.bestE=bestE;
		this.q_Time=q_Time;
		this.repeatEW=repeatEW;
		this.allPruned=allPruned;
		this.searchNumConfsTotal=searchNumConfsTotal;
		this.searchNumConfsPrunedByE=searchNumConfsPrunedByE;
		this.searchNumConfsPrunedByS=searchNumConfsPrunedByS;
		this.searchNumConfsEvaluated=searchNumConfsEvaluated;
		this.searchNumConfsLeft=searchNumConfsLeft;
		this.searchNumPrunedMinDEE=searchNumPrunedMinDEE;
		this.searchBestEnergyFound=searchBestEnergyFound;
	}
	PartitionMessage(PartitionMessage p){
		this.runNum = p.runNum;
		this.mutNum = p.mutNum;
		this.seqNum = p.seqNum;
		this.q = p.q;
		this.bestEMin = p.bestEMin;
		this.bestE=p.bestE;
		this.q_Time=p.q_Time;
		this.repeatEW=p.repeatEW;
		this.allPruned=p.allPruned;
		this.searchNumConfsTotal=p.searchNumConfsTotal;
		this.searchNumConfsPrunedByE=p.searchNumConfsPrunedByE;
		this.searchNumConfsPrunedByS=p.searchNumConfsPrunedByS;
		this.searchNumConfsEvaluated=p.searchNumConfsEvaluated;
		this.searchNumConfsLeft=p.searchNumConfsLeft;
		this.searchNumPrunedMinDEE=p.searchNumPrunedMinDEE;
		this.searchBestEnergyFound=p.searchBestEnergyFound;
		
	}
	PartitionMessage(){
		
	}
	
	
	
}
