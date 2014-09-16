import java.util.ArrayList;
import java.util.HashMap;


class RotamerEntry extends EMatrixEntry{
		SuperRotamer r;
		double[] rotDih;
		
		int pos;
		
		public RotamerEntry(int position,SuperRotamer rot) {
			pos = position;
			r = rot;
		}
		
		public RotamerEntry(int position,SuperRotamer rot, double minE, double maxE, boolean pruned, boolean prunedIsSteric) {
			super(minE,maxE,pruned,prunedIsSteric);
			pos = position;
			r = rot;
		}
		
		public RotamerEntry(int position,SuperRotamer rot, double minE, boolean pruned, double[] rotDih) {
			super(minE,pruned);
			pos = position;
			r = rot;
			this.rotDih = rotDih;
			
		}
		
		public RotamerEntry(int position,SuperRotamer rot, double minE,double maxE, boolean pruned, double[] rotDih) {
			super(minE,maxE,pruned);
			pos = position;
			r = rot;
			this.rotDih = rotDih;
			
		}
		
		public RotamerEntry(int position,SuperRotamer rot, double minE, boolean pruned) {
			super(minE,pruned);
			pos = position;
			r = rot;
			
			
		}
		public RotamerEntry(int position,SuperRotamer rot, double minE, double maxE,boolean pruned) {
			super(minE,maxE,pruned);
			pos = position;
			r = rot;
		}

		@Override
		public boolean applyMutation(Molecule m, ArrayList<ArrayList<Integer>> resByPos,boolean addHydrogens, boolean connectResidues){
			return r.applyMutation(m, resByPos.get(pos), addHydrogens, connectResidues);
		}

		@Override
		void flexible(Molecule m, ArrayList<ArrayList<Integer>> resByPos,boolean flex) {
			for(int resID: resByPos.get(pos))
				m.residue[resID].flexible = flex;
		}

		@Override
		void setEnergyEval(Molecule m,ArrayList<ArrayList<Integer>> resByPos,boolean scEval, boolean bbEval) {
			for(int resID: resByPos.get(pos))
				m.residue[resID].setEnergyEval(scEval, bbEval);
		}

		@Override
		void applyRC(ArrayList<ArrayList<Integer>> resByPos,Molecule m) {
			r.applyRC(resByPos.get(pos),m);
		}
		
		/*void addEref(HashMap<String,double[]> eRef, Molecule m, Emat emat){//ArrayList<ArrayList<Integer>> resByPos){
			double eRefVal = r.getEref(pos, eRef, m, emat.resByPos.get(pos));
			setMinE(minE() - eRefVal);
			setMaxE(maxE(false) - eRefVal);
		}*/
		
		double getEref(HashMap<String,double[]> eRef, Molecule m, ArrayList<ArrayList<Integer>> resByPos){
			return r.getEref(pos, eRef, m, resByPos.get(pos));
		}
		
		public double getEntropy(Molecule m,
				ArrayList<ArrayList<Integer>> resByPos) {
			return r.getEntropy(pos,m,resByPos.get(pos));
		}

		@Override
		void setSCEnergyEval(Molecule m,ArrayList<ArrayList<Integer>> resByPos,
				boolean scEval) {
			for(int resID: resByPos.get(pos))
				m.residue[resID].setSCEnergyEval(scEval);
			
		}
		
		@Override
		int numNonWT(Molecule m,ArrayList<ArrayList<Integer>> resByPos){
			return r.numNonWT(m, resByPos.get(pos));
		}

		@Override
		String printRes(Molecule m,ArrayList<ArrayList<Integer>> resByPos){
			return r.printRes(m,resByPos.get(pos));
		}
		
		@Override
		public String getString(){
			String retString = super.getString();
			retString += r.getString();
			return retString;
		}
		
		RotamerEntry combine(RotamerEntry re, double extraMinE){
			RotamerEntry retRe = this.copy();
			
			retRe.r.combine(re.r);
			retRe.setMinE(retRe.minE() + re.minE() + extraMinE);
			//retRe.setMaxE(retRe.maxE(false) + re.maxE(false) + extraMaxE);
			
			retRe.setPruned(isPruned() || re.isPruned());
			//retRe.setPrunedIsSteric(prunedIsSteric() || re.prunedIsSteric());
				
			return retRe;
		}
		
		RotamerEntry copy(){
			
			RotamerEntry re = new RotamerEntry(pos, r.copy());
			re.setMinE(minE());
			//re.setMaxE(maxE(false));
			re.setPruned(isPruned());
			//re.setPrunedIsSteric(prunedIsSteric());
			
			return re;
		
		}

		@Override
		public boolean[] transRotStrands(Molecule m,
			ArrayList<ArrayList<Integer>> resByPos, MutableResParams strandMut) {
			boolean[] transRotStrands = new boolean[m.numberOfStrands];
			
			for(int molNum: resByPos.get(pos)){
				int str = m.residue[molNum].strandNumber;
				if(m.strand[str].rotTrans)
					transRotStrands[str] = true;
			}

			//If any strand that can rotate & translate contains part of the template,
			//let it rotate and translate
			for(int str=0; str<m.numberOfStrands; str++){
				if(m.strand[str].rotTrans){
					if( m.strand[str].numberOfResidues > strandMut.numMutPerStrand[str] )//strand has template residues
						transRotStrands[str] = true;
				}
			}
			return transRotStrands;
		}
		
		
	}
