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
public class MutableResParams {

	HashMap<String, Residue> resByPDB;
	//ArrayList<ArrayList<Residue>> resByMut;
	int[] allMut;
	
	public MutableResParams(int numMutPos) {
		//resByMut = new ArrayList<ArrayList<Residue>>();
		resByPDB = new HashMap<String,Residue>();
		allMut = new int[numMutPos];
		for(int i=0; i<allMut.length;i++)
			allMut[i] = -1;
	}
	
	public void addRes(int mutPos, Residue res, RotamerLibrary rl, boolean addOrigRot) {
		res.origMutPos = mutPos;
		
		//ArrayList<Residue> tmpAL = new ArrayList<Residue>();
		//tmpAL.add(res);
		
		//resByMut.add(tmpAL);
		resByPDB.put(res.getResNumberString(), res);
		for(int i=0; i<allMut.length;i++)
			if(allMut[i]==-1){
				allMut[i] = res.moleculeResidueNumber;
				break;
			}
		
		res.isMutable = true;
		res.initializeMutableRes(rl, addOrigRot);
		
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
