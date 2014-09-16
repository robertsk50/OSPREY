
public class WCSPSearch extends AStar {

	WCSPOptimization wcspOpt;
	Emat emat;
//	EMatrixEntryWIndex[] bestConf; 
	int treeLevels;
	int[] numRotForRes;
	double Ew;
	PGQueueNode dummy;
	
	WCSPSearch (int treeLevels, int numRotForRes[], Emat emat, double Ew){
		this.emat = emat;
		this.treeLevels = treeLevels;
		this.numRotForRes = numRotForRes;
		this.Ew = Ew;
		curExpansion = new PGExpansionQueue(1);
	}
	
	
	
	@Override
	PGQueueNode doAStar(boolean run1) {
		if(run1){
			int[] conf = new int[treeLevels];
			for(int i=0; i<conf.length;i++)
				conf[i] = -1;
			
			dummy = new PGQueueNode (treeLevels, conf, Double.NEGATIVE_INFINITY,0,-1);
			wcspOpt = new WCSPOptimization(dummy, emat, null, null, numRotForRes,Double.POSITIVE_INFINITY);
			wcspOpt.optimize(null);
			
			double bestE = wcspOpt.bestConf.E;
			if(Ew > 0){
				wcspOpt.allConfs.poll(); //Remove top conformation;
				wcspOpt.getAllConfs(Ew,bestE);
			}
		}
		
		RotConf bConf = wcspOpt.allConfs.poll();
		if(bConf == null)
			dummy.actualConf = null;
		else
			dummy.actualConf = bConf.conf;
	
		return dummy;
	}
	
	@Override
	public void stopSlaves() {
		wcspOpt.cleanUp();
	}
	
	

}
