import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;


public class GurobiCalcInfo implements Serializable {
	Emat emat;
	int[] numNodesForLevel; 
	int[] nodeIndexOffset;
	int[][] numRotRemainingBySeq;
	int[][] seqIndexOffset;
	int numTotalNodes;
	boolean doInteger;
	int[][] seqIndicesPerlevel;
	double upperE;
	
	
	//BYSubRot Members
	int[] numParentRotPerLvl = null;
	int[][] parentRotIndexPerLvl = null;
	int[][] numSubRotPerParentRot = null;
	//First dim is the pos, second dim is the parent;
	ArrayList<HashMap<Integer,ArrayList<Index3>>> subRotsPerLvlPerParent = null;
	boolean bySubRot = false;
	
	
	public GurobiCalcInfo(Emat emat, int[] numNodesForLevel, 
			int[][] numRotRemainingBySeq, int[][] seqIndexOffset, int numTotalNodes,
			boolean doInteger, int[][] seqIndicesPerlevel, double upperE){
		this.emat = emat;
		this.numNodesForLevel = numNodesForLevel; 
		this.numRotRemainingBySeq = numRotRemainingBySeq;
		this.seqIndexOffset = seqIndexOffset;
		this.numTotalNodes = numTotalNodes;
		this.doInteger = doInteger;
		this.seqIndicesPerlevel = seqIndicesPerlevel;
		this.upperE = upperE;
	}
	
	public GurobiCalcInfo(Emat emat, int[] numNodesForLevel, int numTotalNodes,
			boolean doInteger, double upperE,
			int[] numParentRotPerLvl,
			int[][] parentRotIndexPerLvl, int[][] numSubRotPerParentRot,
			ArrayList<HashMap<Integer,ArrayList<Index3>>> subRotsPerLvlPerParent){
				this.emat = emat;
				this.numNodesForLevel = numNodesForLevel; 
				this.numTotalNodes = numTotalNodes;
				this.doInteger = doInteger;
				this.upperE = upperE;
				this.numParentRotPerLvl = numParentRotPerLvl;
				this.parentRotIndexPerLvl = parentRotIndexPerLvl;
				this.numSubRotPerParentRot = numSubRotPerParentRot;
				this.subRotsPerLvlPerParent = subRotsPerLvlPerParent;
				bySubRot = true;
			}
	
}
