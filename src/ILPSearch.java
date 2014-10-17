

public class ILPSearch extends AStar {

	GurobiOptimization ilpOpt;
	Emat emat;
	int treeLevels;
	int[] numRotForRes;
	double Ew;
	PGQueueNode dummy;
	
	ILPSearch (int treeLevels, int numRotForRes[], Emat emat, double Ew){
		this.emat = emat;
		this.treeLevels = treeLevels;
		this.numRotForRes = numRotForRes;
		this.Ew = Ew;
		curExpansion = new PGExpansionQueue(1);
		ilpOpt = new GurobiOptimization(emat,true);
	}
	
	
	
	@Override
	PGQueueNode doAStar(boolean run1) {
		
		int[] conf = new int[treeLevels];
		for(int i=0; i<conf.length;i++)
			conf[i] = -1;
		
		dummy = new PGQueueNode (treeLevels, conf, Double.NEGATIVE_INFINITY,0,-1);
		
		GurobiOptimization.GurobiConf ilpConf = ilpOpt.optimize(emat,energyTuples);
		dummy.actualConf = ilpConf.conf;
		
		ilpOpt.removeConf(ilpConf.conf);
	
		return dummy;
	}
	

}
