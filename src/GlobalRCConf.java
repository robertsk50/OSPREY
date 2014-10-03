import java.io.Serializable;
import java.util.ArrayList;



class GlobalRCConf implements Serializable{

	int[] globalRCs;
	int[] resNum;
	double[] EforRes;

	
	GlobalRCConf(int[] globalRCs, int[] resNum, double[] EforRes){
		
		this.globalRCs = globalRCs;
		this.resNum = resNum;
		this.EforRes = EforRes;
	}
	
}