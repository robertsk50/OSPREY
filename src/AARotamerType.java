import java.io.Serializable;
import java.util.ArrayList;

/**
 * 
 * @author kroberts
 * Stores the AA type information for each rotamer
 */
public class AARotamerType implements Serializable{
	public String name;                    //AA Name
	public int index;
	public String dihedralAtomNames[][];  // Names of atoms involved in the dihedrals for each amino acid
	
	public ArrayList<Rotamer> rotamers;
	
	double entropy = 0.0;
	String oneLet;
	String threeLet;
	
	public AARotamerType(String AAname, int numRots, String[][] dihedAtomNames, int ind) {
		name = AAname;
		
		dihedralAtomNames = dihedAtomNames;
		rotamers = new ArrayList<Rotamer>();
		
		index = ind;
		this.entropy = EnvironmentVars.getEntropyTerm(name);
	}
	
	public void addRot(Rotamer r){
		rotamers.add(r);
	}
	
	/*public void addRot(int[] dihedVals, int RotRLIndex) {
		Rotamer r = new Rotamer(rotamers.size(),dihedVals,this, RotRLIndex);
		rotamers.add(r);
	}*/
	
	public int numRotamers(){
		return rotamers.size();
		
	}
	
	public int numDihedrals(){
		if(dihedralAtomNames == null)
			return 0;
		else
			return dihedralAtomNames.length;
	}

	/**
	 * Returns only the template rotamers for the aaType
	 * 
	 */
	public ArrayList<Rotamer> getTemplRot() {
		ArrayList<Rotamer> templRot = new ArrayList<Rotamer>();
		
		for(Rotamer r: rotamers)
			if(r.resSpecificRot.equals(Rotamer.TEMPL))
				templRot.add(r);
		
		return templRot;
	}

	public ArrayList<Rotamer> getPosSpecRot(String pdbNumString) {
		ArrayList<Rotamer> templRot = new ArrayList<Rotamer>();
		
		for(Rotamer r: rotamers)
			if(r.resSpecificRot.equals(Rotamer.TEMPL) || r.resSpecificRot.equals(pdbNumString))
				templRot.add(r);
		
		return templRot;
	}
	
	
}
