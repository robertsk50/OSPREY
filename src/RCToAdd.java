import java.util.ArrayList;


public class RCToAdd {

	ResidueConformation rc;
	int pos;
	ArrayList<double[]> confs;
	
	RCToAdd(ArrayList<double[]> confs, ResidueConformation resConf, int pos ){
		this.confs = confs;
		rc = resConf;
		this.pos = pos;
	}
	
}
