import java.util.ArrayList;

/**
 * 
 * @author kroberts
 * This class is used in the energy matrix to store the 
 * rotamer pair for that position. 
 *
 */
public class RotamerPairEntry extends EMatrixEntry {
	SuperRotamer r1;
	SuperRotamer r2;
	
	double[] rotDih1;
	double[] rotDih2;
	
	int pos1;
	int pos2;
	
	public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2) {
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;
	}
	
	public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2,double minE, double maxE, boolean pruned, boolean prunedIsSteric) {
		super(minE,maxE,pruned,prunedIsSteric);
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;
	}
	
	public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2,double minE, boolean pruned,
			double[] rotDih1, double[] rotDih2) {
		super(minE,pruned);
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;
		
		this.rotDih1 = rotDih1;
		this.rotDih2 = rotDih2;
	}
	
	public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2,double minE, double maxE,boolean pruned,
			double[] rotDih1, double[] rotDih2) {
		super(minE,maxE,pruned);
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;
		
		this.rotDih1 = rotDih1;
		this.rotDih2 = rotDih2;
	}
	
	public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2,double minE, boolean pruned) {
		super(minE,pruned);
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;
	}
	
	/*public RotamerPairEntry(int position1, SuperRotamer rot1, int position2, SuperRotamer rot2,double minE, double maxE, boolean pruned) {
		super(minE,maxE,pruned);
		pos1 = position1;
		pos2 = position2;
		r1 = rot1;
		r2 = rot2;	
	}*/
	
	
	
	public RotamerPairEntry(SuperRotamer rot1, SuperRotamer rot2, double minEnergy, double maxEnergy) {
		r1 = rot1;
		r2 = rot2;
		setMinE(minEnergy);
		//setMaxE(maxEnergy);
	}	
	
	@Override
	public boolean applyMutation(Molecule m, ArrayList<ArrayList<Integer>> resByPos, boolean addHydrogens, boolean connectResidues){
		boolean mut1 = r1.applyMutation(m, resByPos.get(pos1), addHydrogens, connectResidues);
		boolean mut2 = r2.applyMutation(m, resByPos.get(pos2), addHydrogens, connectResidues);
		
		return (mut1 || mut2);
	}

	@Override
	void setEnergyEval(Molecule m,ArrayList<ArrayList<Integer>> resByPos, boolean scEval, boolean bbEval) {
		for(Integer resID:resByPos.get(pos1)){
			m.residue[resID].setEnergyEval(scEval, bbEval);
		}
		for(Integer resID:resByPos.get(pos2)){
			m.residue[resID].setEnergyEval(scEval, bbEval);
		}
	}
	
	@Override
	void setSCEnergyEval(Molecule m,ArrayList<ArrayList<Integer>> resByPos, boolean scEval) {
		for(Integer resID:resByPos.get(pos1)){
			m.residue[resID].setSCEnergyEval(scEval);
		}
		for(Integer resID:resByPos.get(pos2)){
			m.residue[resID].setSCEnergyEval(scEval);
		}
	}

	@Override
	void applyRC(ArrayList<ArrayList<Integer>> resByPos, Molecule m) {
		r1.applyRC(resByPos.get(pos1),m);
		r2.applyRC(resByPos.get(pos2),m);
	}

	@Override
	void flexible(Molecule m,ArrayList<ArrayList<Integer>> resByPos, boolean flex) {
		for(Integer resID:resByPos.get(pos1)){
			m.residue[resID].flexible = flex;
		}
		for(Integer resID:resByPos.get(pos2)){
			m.residue[resID].flexible = flex;
		}
	}
	
	RotamerPairEntry combine(RotamerPairEntry re, int posToMerge){
		RotamerPairEntry retRE = this.copy();
		
		
		SuperRotamer r1ToMerge=null;
		SuperRotamer r2ToMerge=null;
		if(posToMerge == 1){
			r1ToMerge = retRE.r1;
			r2ToMerge = re.r1;
		}
		else if(posToMerge == 2){
			r1ToMerge = retRE.r2;
			r2ToMerge = re.r2;
		}
		
		r1ToMerge.combine(r2ToMerge);
		retRE.setMinE(retRE.minE() + re.minE());
		//retRE.setMaxE(retRE.maxE(false) + re.maxE(false));
		
		retRE.setPruned(isPruned() || re.isPruned());
		//retRE.setPrunedIsSteric(prunedIsSteric() || re.prunedIsSteric());
		
		return retRE;
	}
	
	//This doesn't create a valid energy so we make it infinity to make sure it is reset
	RotamerPairEntry combine(RotamerEntry re, int posToMerge){
		RotamerPairEntry retRE = this.copy();
		
		
		SuperRotamer r1ToMerge=null;
		SuperRotamer r2ToMerge=null;
		if(posToMerge == 1){
			r1ToMerge = retRE.r1;
			r2ToMerge = re.r;
		}
		else if(posToMerge == 2){
			r1ToMerge = retRE.r2;
			r2ToMerge = re.r;
		}
		
		r1ToMerge.combine(r2ToMerge);
		retRE.setMinE(Double.NEGATIVE_INFINITY);
		//retRE.setMaxE(retRE.maxE(false) + re.maxE(false));
		
		retRE.setPruned(isPruned());
		//retRE.setPrunedIsSteric(prunedIsSteric() || re.prunedIsSteric());
		
		return retRE;
	}
	
	RotamerPairEntry copy(){
		RotamerPairEntry re = new RotamerPairEntry(pos1, r1.copy(),pos2, r2.copy());
		re.setMinE( minE());
		//re.setMaxE( maxE(false));
		re.setPruned(isPruned());
		//re.setPrunedIsSteric(prunedIsSteric());
		
		return re;
	}
	
	RotamerPairEntry swappedCopy(){
		RotamerPairEntry re = new RotamerPairEntry(pos2, r2.copy(),pos1, r1.copy());
		re.setMinE(minE());
		//re.setMaxE(maxE(false));
		re.setPruned(isPruned());
		//re.setPrunedIsSteric(prunedIsSteric());
		
		return re;
	}
	
	@Override
	public String getString(){
		String retString = super.getString();
		retString += r1.getString();
		retString += "@ "+ r2.getString();
		return retString;
	}
	
	@Override
	public boolean[] transRotStrands(Molecule m,
		ArrayList<ArrayList<Integer>> resByPos, MutableResParams strandMut) {
		boolean[] transRotStrands = new boolean[m.numberOfStrands];
		
	
		for(int molNum: resByPos.get(pos1)){
			int str = m.residue[molNum].strandNumber;
			if(m.strand[str].rotTrans)
				transRotStrands[str] = true;
		}
		
		for(int molNum: resByPos.get(pos2)){
			int str = m.residue[molNum].strandNumber;
			if(m.strand[str].rotTrans)
				transRotStrands[str] = true;
		}
		

		return transRotStrands;
	}
	
}

