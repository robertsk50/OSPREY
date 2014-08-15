import java.io.Serializable;


public class EMatrixEntrySlim implements Serializable{
	double minE;
	double maxE; 
	boolean pruned;
	//boolean prunedIsSteric;
	int[] index;
	double[] rotDih1;
	double[] rotDih2;
	
	
	EMatrixEntrySlim(EMatrixEntry eme, int[] ind){
		index = ind;
		minE = eme.minE();
		//maxE = eme.maxE(false);
		pruned = eme.isPruned();
		//prunedIsSteric = eme.prunedIsSteric();
	}
	
	EMatrixEntrySlim(RotamerEntry eme, int[] ind){
		index = ind;
		minE = eme.minE();
		maxE = eme.maxE();
		//maxE = eme.maxE(false);
		pruned = eme.isPruned();
		rotDih1 = eme.rotDih;
		//prunedIsSteric = eme.prunedIsSteric();
	}
	
	EMatrixEntrySlim(RotamerPairEntry eme, int[] ind){
		index = ind;
		minE = eme.minE();
		maxE = eme.maxE();
		//maxE = eme.maxE(false);
		pruned = eme.isPruned();
		rotDih1 = eme.rotDih1;
		rotDih2 = eme.rotDih2;
		//prunedIsSteric = eme.prunedIsSteric();
	}
	
	
	EMatrixEntrySlim(double minE, int[] ind){
		index = ind;
		this.minE = minE;
	}
}
