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
	
	// Clear all rotamers for this amino acid.
	public void clearRotamers(){
		rotamers = new ArrayList<Rotamer>();
	}
	
	
	// Add a new rotamer, provided it is not already in our list.
	// PGC 2014: Optionally add a rotamer only if it one of its dihedrals are DEFAULTMINWIDTH from another rotamer.
	public boolean addRot(Rotamer r){


		if (!EnvironmentVars.CLUSTER_SIMILAR_ROTAMERS){
			rotamers.add(r);
			return true;
		}
		else{
			double clustering_size = r.DEFAULTMINWIDTH;
			// PGC 2014: Do not add a rotamer if it is clustering_size degrees from another rotamer.
			boolean similarRotamerExists = false;
			for( int i = 0; i< rotamers.size(); i++){
				Rotamer existingRotamer = rotamers.get(i);
				// Dihedrals within the minimization range
				int dihedralsInRange = existingRotamer.values.length;
				for(int angleix = 0; angleix < existingRotamer.values.length; angleix++){
					if(Math.abs(r.values[angleix] -existingRotamer.values[angleix]) <= clustering_size){
						dihedralsInRange -= 1;
						
					}				
				}
				// All dihedral angles were within range for a specific rotamer.
				if (dihedralsInRange == 0){
					similarRotamerExists = true;
				}
			}

			// Add this rotamer 
			if (!similarRotamerExists){
				rotamers.add(r);
				return true;
			}
			// Do not add this rotamer because all dihedrals are very close to another dihedral.
			else{
				return false;
			}
		}
	
	}
	
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
