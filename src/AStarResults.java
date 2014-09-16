
public class AStarResults {
	
	static final int NOTDONE = 0;
	static final int DONE = 1;
	static final int ADDTUPLES = 2;
	static final int INCREASEIVAL = 3;
	static final int PARTITIONROT = 4;
	
	double bestE;
	double lowestBound;
	long numConfsEvaluated;
	double lastBound;
	
	int status = NOTDONE;
	
	AStarResults(double bestE, double lowestBound, long numConfsEvaluated, double lastBound){
		this.bestE = bestE;
		this.lowestBound = lowestBound;
		this.numConfsEvaluated = numConfsEvaluated;
		this.lastBound = lastBound;
	}
	
}
