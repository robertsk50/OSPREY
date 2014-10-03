import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;


public class ResidueConformationLibrary implements Serializable{

	ArrayList<ResidueConformation> allRCs;
		
	ResidueConformationLibrary(){
		allRCs = new ArrayList<ResidueConformation>();
		
	}

	ResidueConformation getRC(int globalID){
		return allRCs.get(globalID);
	}
	
	public ResidueConformation addResidueConformation(Rotamer rot, int pertState, int res){
		ResidueConformation rc = new ResidueConformation(allRCs.size(), rot, pertState, res);
		allRCs.add(rc);
		
		return rc;
	}
	
	public ResidueConformation addSubResidueConformation(Rotamer rot, int pertState, int res, ResidueConformation parent){
		ResidueConformation rc = ResidueConformation.newSubRC(allRCs.size(), rot, pertState, res, parent);
		allRCs.add(rc);
		
		return rc;
	}

	/**
	 * Return all residue conformations for a given pdb number and given amino acid 
	 * 
	 * @param pdbResNum String, the pdb residue number
	 * @param AAname    String, name of the amino acid we want RCs for
	 * @return
	 */
	public ArrayList<ResidueConformation> getRCsPosType(int strResNum, String AAname) {
		//KER: This is innefficient because we loop through allRCs everytime
		//KER: instead of caching the RCs per residue and/or per AA.
		ArrayList<ResidueConformation> rcs = new ArrayList<ResidueConformation>();
		for(ResidueConformation rc : allRCs){
			if(rc.res == strResNum && rc.rot.aaType.name.equalsIgnoreCase(AAname))
				rcs.add(rc);
		}
		return rcs;
	}

	public void renumberRCs() {
		int ctr=0;
		for(ResidueConformation rc: allRCs){
			rc.id = ctr++;
		}
		
	}

	/**
	 * Saves all of the residue conformations in the rc library.
	 * Can be read with loadGlobalRCs
	 * @param fileName
	 * @return
	 */
	public void save(String fileName){
		try {
			
				FileOutputStream fileOutputStream = new FileOutputStream(fileName);
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				PrintStream logPS = new PrintStream( bufferedOutputStream );
				for(ResidueConformation rc: allRCs){
					logPS.println(rc.toString());
				}
				logPS.close();
			
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
	}

	/**
	 * Load in residue conformations from a file
	 * @param rotFilename name of global RC file
	 */
	public void loadGlobalRCs(String rcFilename, RotamerLibrary rl) {
	
		try{
			FileInputStream is = new FileInputStream( rcFilename );
			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;

			

			// Skip over comments (lines starting with !)
			curLine = bufread.readLine();
			while( curLine.charAt(0) == '!' )
				curLine = bufread.readLine();

			
			allRCs = new ArrayList<ResidueConformation>();
			

			while( curLine != null ) {
				if(curLine.charAt(0) == '!'){
					curLine = bufread.readLine();
					continue;
				}


				int id = (new Integer(KSParser.getToken(curLine,1))).intValue();
				int res = (new Integer(KSParser.getToken(curLine, 2))).intValue();
				int rlIndex = (new Integer(KSParser.getToken(curLine, 3))).intValue();
				int pertState = (new Integer(KSParser.getToken(curLine, 4))).intValue();
				
				assert id == allRCs.size();
				
				ResidueConformation rc = addResidueConformation(rl.getRot(rlIndex), pertState, res);

				curLine = bufread.readLine();
			}
			bufread.close();

		}catch(Exception E){
			System.out.println("Couldn't properly load rot library.");
		}

	}

	public int getGlobalIndex(int molResNum, int rlIndex, int state) {
		for(ResidueConformation rc: allRCs){
			if(rc.res == molResNum && rc.rot.rlIndex == rlIndex && rc.pertState == state )
				return rc.id;
		}
		return -1;
	}

	
	

	
}


