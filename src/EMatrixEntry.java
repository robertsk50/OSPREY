import java.io.Serializable;
import java.util.ArrayList;


public class EMatrixEntry implements Serializable{

	private double minE = Double.POSITIVE_INFINITY;
	private double maxE = Double.POSITIVE_INFINITY;
	//private double maxE = Double.POSITIVE_INFINITY; 
	private boolean pruned = false;
	//private boolean prunedIsSteric = false;
	
	public EMatrixEntry(){
		
	}
	
	public EMatrixEntry(double minE, double maxE, boolean pruned, boolean prunedIsSteric){
		this.minE = minE;
		//this.maxE = maxE;
		this.setPruned(pruned);
		//this.setPrunedIsSteric(prunedIsSteric);
	}
	
	public EMatrixEntry(double minE, boolean pruned){
		this.minE = minE;
		this.setPruned(pruned);
	}
	
	public EMatrixEntry(double minE, double maxE,boolean pruned){
		this.minE = minE;
		this.maxE = maxE;
		this.setPruned(pruned);
	}
	
	public void setMinE(double energy) {
		minE = energy;	
	}

	/*public void setMaxE(double energy) {
		maxE = energy;
	}*/

	/*public double getMaxE(boolean getMin) {
		if(getMin)
			return minE;
		else{
			System.out.println("Code does not support max");
			System.exit(0);
		}
		return maxE;
	}*/
	
	protected void setPruned(boolean p){
		pruned = p;
	}
	
	boolean applyMutation(Molecule m, ArrayList<ArrayList<Integer>> resByPos,boolean addHydrogens, boolean connectResidues){return false;}
	void setEnergyEval(Molecule m, ArrayList<ArrayList<Integer>> resByPos, boolean scEval,boolean bbEval){}
	void setSCEnergyEval(Molecule m,ArrayList<ArrayList<Integer>> resByPos, boolean scEval){}
	void applyRC(ArrayList<ArrayList<Integer>> resByPos, Molecule m){}
	void flexible(Molecule m, ArrayList<ArrayList<Integer>> resByPos, boolean flex){}
	int numNonWT(Molecule m,ArrayList<ArrayList<Integer>> resByPos) {return 0;}
	String printRes(Molecule m,ArrayList<ArrayList<Integer>> resByPos) {return "";}
	public boolean[] transRotStrands(Molecule m,ArrayList<ArrayList<Integer>> resByPos, MutableResParams strandMut) {return null;}
	
	
	public void set(EMatrixEntrySlim eme){
		//maxE = eme.maxE;
		minE = eme.minE;
		setPruned(eme.pruned);
		//setPrunedIsSteric(eme.prunedIsSteric);
	}
	
	/*
	 * Input: if true, returns min instead of max 
	 */
	public double maxE(boolean getMin){
		if(getMin)
			return minE;
		else{
			System.out.println("Code does not support Max energies");
			//System.exit(0);
		}
		return 0.0f;
		//else
			//return maxE;
	}
	
	public double minE(){
		return minE;
	}
	
	public double maxE(){
		return maxE;
	}

	public String getString() {
		String retString = minE+" "+isPruned()+" ";
		return retString;
	}

	public boolean isPruned() {
		return pruned;
	}

	

	/*public void setPrunedIsSteric(boolean prunedIsSteric) {
		this.prunedIsSteric = prunedIsSteric;
	}

	public boolean prunedIsSteric() {
		return prunedIsSteric;
	}*/

	
	
}
