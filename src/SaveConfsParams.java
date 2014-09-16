
public class SaveConfsParams {

	int numTopConfs = 0;
	boolean saveTopConfs = false;
	boolean printTopConfs = false;
	boolean savePDBs = false;
	
	public SaveConfsParams(int numSaveConfs, boolean saveConfs,
			boolean printConfs, boolean savePDBs) {
	
		this.numTopConfs = numSaveConfs;
		this.saveTopConfs = saveConfs;
		this.printTopConfs = printConfs;
		this.savePDBs = savePDBs;
	}
	
	
	
}