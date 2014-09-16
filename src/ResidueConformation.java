import java.io.Serializable;


public class ResidueConformation implements Serializable{
	
	int id;
	Rotamer rot;
//	Perturbation[] perts;
	int pertState;
	
	int res;
	
	ResidueConformation(int id, Rotamer rot, int pertState, int res){
		this.id = id;
		this.rot = rot;
//		this.perts = perts;
		this.pertState = pertState;
		this.res = res;
	}
	
	@Override
	public String toString(){
		String s= id+" "+res+" "+rot.rlIndex+" "+pertState;
		
		return s;
		
	}
	
}
