import java.io.Serializable;


public class ResidueConformation implements Serializable{
	
	int id;
	Rotamer rot;
//	Perturbation[] perts;
	int pertState;
	int parent;
	
	int res;
	
	ResidueConformation(int id, Rotamer rot, int pertState, int res){
		this.id = id;
		this.parent = id;
		this.rot = rot;
//		this.perts = perts;
		this.pertState = pertState;
		this.res = res;
	}
	
	public static ResidueConformation newSubRC(int id, Rotamer rot, int pertState, int res, ResidueConformation parent){
		ResidueConformation rc = new ResidueConformation(id, rot, pertState, res);
		rc.parent = parent.parent;
		return rc;
	}
	
	@Override
	public String toString(){
		String s= id+" "+res+" "+rot.rlIndex+" "+pertState;
		
		return s;
		
	}
	
}
