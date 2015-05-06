import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class to store information about the mutable residues
 * Used to setup the energy matrix
 *  
 *  
 * @author kroberts
 *
 */
public class MutableResParams implements Serializable {

	HashMap<String, Residue> resByPDB;
	//ArrayList<ArrayList<Residue>> resByMut;
	int[] allMut;
	int[] resStrandNum;
	int[] resStrand;
	int[] numMutPerStrand;
	
	public MutableResParams(int numMutPos, int numOfStrands) {
		//resByMut = new ArrayList<ArrayList<Residue>>();
		resByPDB = new HashMap<String,Residue>();
		allMut = new int[numMutPos];
		for(int i=0; i<allMut.length;i++)
			allMut[i] = -1;
		
		numMutPerStrand = new int[numOfStrands];
		
		resStrandNum = new int[numMutPos];
		resStrand = new int[numMutPos];
	}
	
	public void addRes(int mutPos, Residue res, RotamerLibrary rl, boolean addOrigRot, Residue prevRes) {
		res.origMutPos = mutPos;
		
		resByPDB.put(res.getResNumberString(), res);
		for(int i=0; i<allMut.length;i++)
			if(allMut[i]==-1){
				allMut[i] = res.moleculeResidueNumber;
				resStrandNum[i] = res.strandResidueNumber;
				resStrand[i] = res.strandNumber;
				numMutPerStrand[res.strandNumber]++;
				break;
			}
		
		if(res.isMutable == false){ //If we haven't already set up the mutable res, do it now.
			res.isMutable = true;
			res.initializeMutableRes(rl, addOrigRot, prevRes);
		}
		
		
		
	}
	
	public void checkWT(boolean[] strandPresent, ParamSet sParams){
		for(Residue mr: resByPDB.values()){
			mr.updateWT(strandPresent,sParams);
		}
	}

	public int numMutPos() {
		return allMut.length;
	}
	
	
}
