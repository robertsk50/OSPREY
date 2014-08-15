import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;


public class SuperRotamer implements Serializable{
	
	//ArrayList<Residue> residues;
	
	int[] rotamers;
	
	//ArrayList<AARotamerType> aaTypes;
	
	
//	public SuperRotamer() {
//		//rotamers = new int[0];
//	}
	
	public SuperRotamer(int r){
		rotamers = new int[1];
		rotamers[0] = r;
	}
	
	public SuperRotamer(int[] r){
		rotamers = r;
	}
	
	public boolean applyMutation(Molecule m, ArrayList<Integer> residues, boolean addHydrogens, boolean connectResidues){
		boolean neededMut = false;
		int ctr=0;
		if(needsMutation(residues,m)){
			neededMut = true;
			for(Integer resID:residues){
				MutUtils.changeResidueType(m, resID, m.residue[resID].rl.getRot(rotamers[ctr]).aaType.name, addHydrogens,connectResidues);
				ctr++;
			}
		}
		return neededMut;
	}

	public void applyRotamer(ArrayList<Integer> residues, Molecule m) {
		int ctr=0;
		for(int resID:residues){
			MutUtils.applyRotamer(m, m.residue[resID].rl.getRot(rotamers[ctr]),resID);
			ctr++;
		}	
	}
	
	public double getEref(int pos, HashMap<String,double[]> eRef, Molecule m, ArrayList<Integer> residues) {
		int ctr = 0;
		double retEref = 0;
		for(int resID: residues){
			retEref += eRef.get(m.residue[resID].getResNumberString())[m.residue[resID].rl.getRot(rotamers[ctr]).aaType.index];
			ctr++;
		}
		return retEref;
	}
	
	public double getEntropy(int pos, Molecule m, ArrayList<Integer> residues) {
		int ctr = 0;
		double retEntropy = 0;
		for(int resID: residues){
			retEntropy += m.residue[resID].rl.getRot(rotamers[ctr]).aaType.entropy;
			ctr++;
		}
		return retEntropy;
	}
	
	public boolean needsMutation(ArrayList<Integer> residues, Molecule m) {
		int ctr=0;
		boolean mutOnce = false;
		for(Integer resID:residues){
			Residue r = m.residue[resID];
			try{
			if(! r.isSameAA(m.residue[resID].rl.getRot(rotamers[ctr]).aaType.name))
				return true;
			else if(!r.mutatedOnce && m.residue[resID].canMutate){
				mutOnce = true;
			}
			ctr++;
			}catch(Exception E){
				System.out.println(resID + " "+ ctr);
				E.printStackTrace();
			}
		}
		return mutOnce;
	}
	
	public int numNonWT(Molecule m,ArrayList<Integer> residues){
		int numNonWT = 0;
		int ctr=0;
		for(int resID:residues){
			if(!m.residue[resID].rl.getRot(rotamers[ctr]).aaType.name.equalsIgnoreCase(m.residue[resID].defaultAA))
				numNonWT++;
			ctr++;
		}
		return numNonWT;
	}

	public String printRes(Molecule m, ArrayList<Integer> residues) {
		String resStr = "";
		int ctr=0;
		for(int resID: residues){
			Rotamer r = m.residue[resID].rl.getRot(rotamers[ctr]);
			resStr += m.residue[resID].getResNumberString()+" "+r.aaType.name+" "+r.aaIndex+" ";
			ctr++;
		}
		return resStr;
		
	}
	
	public void combine(SuperRotamer sr2){
		int[] combRots = new int[rotamers.length+sr2.rotamers.length];
		
		int ctr=0;
		for(int r:rotamers){
			combRots[ctr] = r;
			ctr++;
		}
		for(int r:sr2.rotamers){
			combRots[ctr] = r;
			ctr++;
		}
		
		rotamers = combRots;
	}
	
	/*public void addRot(int r){
		int[] newRots = new int[rotamers.length+1];
		System.arraycopy(rotamers, 0, newRots, 0, newRots.length);
		newRots[newRots.length-1] = r;
		rotamers = newRots;
	}*/

	public SuperRotamer copy() {
		
		int[] newRotamers = new int[rotamers.length];
		System.arraycopy(rotamers, 0, newRotamers, 0, rotamers.length);
		
		SuperRotamer sr = new SuperRotamer(newRotamers);
		return sr;
	}

	public String getString() {
		String retString = "";
		for(int i:rotamers){
			retString += i+" ";
		}
		return retString;
	}

	public ArrayList<Rotamer> getRotamers(Molecule m, ArrayList<Integer> residues){
		ArrayList<Rotamer> rots = new ArrayList<Rotamer>();
		int ctr=0;
		for(Integer resID:residues){
			rots.add(m.residue[resID].rl.getRot(rotamers[ctr]));
			ctr++;
		}
		
		
		return rots;
	}

	

	
	
}
