import java.io.Serializable;


public class RotConf implements Comparable,Serializable{

	double E;
	EMatrixEntryWIndex[] conf;
	
	public RotConf(EMatrixEntryWIndex[] conf, double E) {
		this.E = E;
		this.conf = conf;
	}
	
	@Override
	public int compareTo(Object otherObject) {
		if(!(otherObject instanceof RotConf)){
	            throw new ClassCastException("A RotConf object expected.");
	    }
	    RotConf other = (RotConf)otherObject;
	    if(this.E > other.E){
	            return 1;
	    }
	    else if(this.E < other.E){
	            return -1;
	    }
	    else{ // Nodes have the same score
            return checkConf(this, other);
	    }
	}
	
	//Checks if the two given nodes have the same partially assigned conformation
	private int checkConf(RotConf conf1, RotConf conf2){
	
		for (int l=0; l<conf1.conf.length; l++){
			if(conf1.conf[l].index[1] < conf2.conf[l].index[1])
				return -1;
			else if(conf1.conf[l].index[1] > conf2.conf[l].index[1])
				return 1;
			else if(conf1.conf[l].index[2] < conf2.conf[l].index[2])
				return -1;
			else if(conf1.conf[l].index[2] > conf2.conf[l].index[2])
				return 1;
		}
		
		//The partially assigned conformations are the same
		return 0;		
	}
	
	
}
