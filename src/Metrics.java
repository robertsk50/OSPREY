import java.util.ArrayList;


public class Metrics {
	
	ArrayList<Double> percPruned;
	
	long numExpanded = 0;
	long totNumNodes = 0;
	long totalNumConfs = 0;
	
	long startTime = -1;
	long EmatTime = 0;
	long DEEtime = 0;
	long EnumTime = 0;
	long AStime = 0;
	long Etime = 0;
	long endTime = -1;
	long totalTime = 0;
	double trueIval;
	long loopStart= -1;
	
	//Measures of A* efficiency
	long[] timePerLevel; 
	double[] boundTightnessPerLevel; //Report how close the bound was to the actual solution (in percentage)
	long[] nodesEvaluatedPerLevel;
	
	public void setStartTime(){
		if(startTime == -1)
			startTime = System.currentTimeMillis();
		
		System.out.println("DEE Start Time: "+startTime);
	}
	
	public void setEndTime(){
		if(endTime == -1)
			endTime = System.currentTimeMillis();
		
		totalTime = endTime - startTime;
	}
	
	public void setLoopStart(){
		if(loopStart == -1)
			loopStart = System.currentTimeMillis();
		
		System.out.println("DEE Start Time: "+loopStart);
	}

	public void print() {
		System.out.println("TRUEIVAL: "+trueIval);
		System.out.println("EmatTime: "+EmatTime);
		System.out.println("DEETime: "+DEEtime);
		System.out.println("ASTime: "+AStime);
		System.out.println("EnumTime: "+EnumTime);
		System.out.println("Etime: "+Etime);
		System.out.println("TotalTime: "+totalTime);
		System.out.println("TotalNumConfs: "+totalNumConfs);
		System.out.println("ASNumExpanded: "+numExpanded);
		System.out.println("ASTotalNodes: "+totNumNodes);
		if(boundTightnessPerLevel != null){
			System.out.print("boundTightnessAvgPerLevel: ");for(int i=0; i<boundTightnessPerLevel.length;i++){System.out.print((boundTightnessPerLevel[i]/nodesEvaluatedPerLevel[i])+" ");}System.out.println("");
		}
	}
	
	public void initASMetrics(int numLevels){
		if(timePerLevel == null){
			timePerLevel = new long[numLevels]; 
			boundTightnessPerLevel = new double[numLevels]; 
			nodesEvaluatedPerLevel = new long[numLevels];
		}
	}

	public void updateASMetrics(int curLevel, double actualBound, double fScore) {
		if(actualBound == Double.POSITIVE_INFINITY){ //Skip these because they mess up the calculation
			nodesEvaluatedPerLevel[curLevel]--;
		}else{
			boundTightnessPerLevel[curLevel] += Math.abs((fScore-actualBound) / actualBound);
		}
	}

	public void updateASTimes(int curLevel, long time, int numNodes) {
		timePerLevel[curLevel] += time;
		nodesEvaluatedPerLevel[curLevel] += numNodes;
		
	}
	
	
	
	
}
