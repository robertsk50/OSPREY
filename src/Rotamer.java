import java.io.Serializable;

/**
 * 
 * @author kroberts
 * This class is meant to represent a single rotamer.
 * The rotamer should store the dihedrals that make it up
 *
 */
public class Rotamer implements Serializable{

	static final String TEMPL = "TEMPL";
	
	double DEFAULTMINWIDTH = 9.0;
	
	//Dihedrals for this specific rotamer
	double[] values; 
	double[] minimizationWidth;
	double volume;
	
	//This string is the a way to tell whether the 
	//rotamer was part of the original LovellRotamer.dat
	//file or whether it was added during the run somewhere.
	//It's value will be the pdbNumber for the residue
	//it was added for.
	String resSpecificRot; 
	
	boolean isWTrot;
	
	AARotamerType aaType;
	
	int aaIndex;
	int rlIndex;
	int parent;
	
	
	public Rotamer(int rotAAIndex,double[] vals, AARotamerType aaType, int rlIndex, String pdbNum, boolean isWTrot) {
		aaIndex = rotAAIndex;
		this.rlIndex = rlIndex;
		this.parent = rlIndex;
		values = vals;
		this.aaType = aaType;
		this.resSpecificRot = pdbNum;
		this.isWTrot = isWTrot;
		
		//KER: for now we set the minimization width to the default
		if(vals != null){
			minimizationWidth = new double[values.length];
			for(int i=0; i<minimizationWidth.length;i++){
				minimizationWidth[i] = DEFAULTMINWIDTH;
			}
		}
		
	}
	
	public Rotamer(int rotAAIndex, double[] vals, AARotamerType aaType, int rlIndex, String pdbNum, double[] minimizationWidth, boolean isWTrot) {
		aaIndex = rotAAIndex;
		this.rlIndex = rlIndex;
		this.parent = rlIndex;
		values = vals;
		this.aaType = aaType;
		this.resSpecificRot = pdbNum;
		this.isWTrot = isWTrot;
		
		if(minimizationWidth == null){
			//KER: for now we set the minimization width to the default
			if(vals != null){
				this.minimizationWidth = new double[values.length];
				for(int i=0; i<this.minimizationWidth.length;i++){
					this.minimizationWidth[i] = DEFAULTMINWIDTH;
				}
			}
		}
		else
			this.minimizationWidth = minimizationWidth;
		
	}
	
	public static Rotamer newSubRot(int rotAAIndex, double[] vals, AARotamerType aaType, int rlIndex, String pdbNum, 
			double[] minimizationWidth, Rotamer parent) {
		Rotamer r = new Rotamer(rotAAIndex, vals,  aaType, rlIndex, pdbNum, minimizationWidth,false);
		r.parent = parent.parent;
		return r;
	}
	
	

	public void setVolume(double vol){
		volume = vol;
	}
	
	@Override
	public String toString(){
		String s="";
		s += rlIndex+" ";
		s += aaType.name+" "+aaIndex+" ";
		
		if(values != null){
			for(int j=0; j<values.length;j++)
				s += values[j]+" ";
			for(int j=0; j<minimizationWidth.length;j++)
				s += minimizationWidth[j]+" ";
		}
		
		s += resSpecificRot;
		
		// PGC 2014: if this is the wildtype rotamer, add "WT" to the end of the aaRots file.
		if(this.isWTrot){
			s += " WT";
		}
		
		return s;
		
	}
	
}
