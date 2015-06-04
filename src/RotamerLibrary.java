/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.0
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University
	
	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as 
	published by the Free Software Foundation, either version 3 of 
	the License, or (at your option) any later version.
	
	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.
	
	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.
		
	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.
	
	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129 
			USA
			e-mail:   www.cs.duke.edu/brd/
	
	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
*/

////////////////////////////////////////////////////////////////////////////////////////////
// RotamerLibrary.java
//
//  Version:           2.0
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 * 
 */

import java.io.*;
import java.util.*;

/**
 * This class implements a rotamer library reader. It reads from an input file that contains rotamer information
 * for amino acid types or other generic residues.
 */
public class RotamerLibrary implements Serializable {
	

	
	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = false;
	
	
	//First index is the amino acid index, second is the rotamer index for that AA
	private AARotamerType[] aaTypes;
	
	//Holds an array with all rotamers that have been created
	private Rotamer[] allRotamers;
	//String aaNames[];       // Array of amino-acid names
	private int numAAallowed;       // Number of AA types read
	
	//String logFile = "";
	//private int numDihedrals[];     // Number of dihedrals per AA
	//private int numRotamers[];      // Number of rotamers per AA
	
	//private String dihedralAtomNames[][][];  // Names of atoms involved in the dihedrals for each amino acid
	//private int rotamerValues[][][];  // Actual angle values for each rotamer for each amino acid
	private double rotamerVolumes[][];	// Volumes of each rotamer for each amino acid
	
	//private int totalNumRotamers; //AAs with 0 rotamers are counted as 1 rotamer	
	//private int rotamerIndexOffset[] = null; //the rotamer index offset for each amino acid (AAs with 0 rotamers are counted as 1 rotamer)
	
	private String rotFile;
	private String volFilename;
	
	private boolean canMutate = false;
	
	private boolean addedRotamers = false;
	
	public static HashMap<String,String> three2one = null;
	
	public boolean isAddedRotamers() {
		return addedRotamers;
	}

	public void setAddedRotamers(boolean addedRotamers) {
		this.addedRotamers = addedRotamers;
	}
	
	// Generic constructor
	RotamerLibrary(String rotFilename, boolean canMut) {		
		
		canMutate = canMut;
		
		try {
			readRotLibrary(rotFilename);
		}
		catch (Exception e){
			e.printStackTrace();
			System.out.println("ERROR reading rotamer library file: "+e);
			System.exit(1);
		}
		
		if(three2one == null)
			initThree2One();
	}
	
	//Read in all of the rotamers for all amino acids from the rotFilename file
	private void readRotLibrary(String rotFilename) throws Exception {
		
		rotFile = rotFilename;
		volFilename = rotFile + ".vol";
		
		// HANDLE THE NORMAL AAs	
		FileInputStream is = new FileInputStream( rotFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		int curAA = 0;
		numAAallowed = 0;
		
		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
		
		numAAallowed = (new Integer(getToken(curLine,1))).intValue(); //the first non-comment line is the number of AA types
		curLine = bufread.readLine();
		
		ArrayList<Rotamer> allRots = new ArrayList<Rotamer>();
		aaTypes = new AARotamerType[numAAallowed];
		
	  	while( curLine != null ) {
			if(curLine.charAt(0) == '!'){
				curLine = bufread.readLine();
				continue;
			}
			
	  		String aaName = getToken(curLine,1);
			int numDihedrals = (new Integer(getToken(curLine,2))).intValue();
			int numRotamers = (new Integer(getToken(curLine,3))).intValue();
			
			String[][] dihedralAtomNames = new String[numDihedrals][4];
			

			// Read in the actual dihedrals
			for(int q=0;q<numDihedrals;q++) {
				curLine = bufread.readLine();
				dihedralAtomNames[q][0] = getToken(curLine,1);
				dihedralAtomNames[q][1] = getToken(curLine,2);
				dihedralAtomNames[q][2] = getToken(curLine,3);
				dihedralAtomNames[q][3] = getToken(curLine,4);
			}
			
			if(numRotamers == 0) //ALA or GLY
				aaTypes[curAA] = new AARotamerType(aaName, 1, null,curAA);
			else
				aaTypes[curAA] = new AARotamerType(aaName, numRotamers, dihedralAtomNames,curAA);
			
			// Read in the actual rotamers
			for(int q=0;q<numRotamers;q++) {
				curLine = bufread.readLine();
				double[] rotamerValues = new double[numDihedrals];
				for(int w=0;w<numDihedrals;w++) {
					rotamerValues[w] = (new Double(getToken(curLine,(w+1)))).doubleValue();
				}
				Rotamer curRot = new Rotamer(q,rotamerValues,aaTypes[curAA],allRots.size(),Rotamer.TEMPL,false);
				aaTypes[curAA].addRot(curRot);
				allRots.add(curRot);
			}
			if(numRotamers <=0){
				Rotamer curRot = new Rotamer(0,null,aaTypes[curAA],allRots.size(),Rotamer.TEMPL,false);
				aaTypes[curAA].addRot(curRot);
				allRots.add(curRot);
			}
					
			/*totalNumRotamers += numRotamers;
			if (numRotamers<=0) //ALA or GLY
				totalNumRotamers += 1;*/
			
			curAA++;
			curLine = bufread.readLine();
		}
	  	
		bufread.close();
		
		allRotamers = allRots.toArray(new Rotamer[0]);
		
		if (curAA!=numAAallowed){
			System.out.println("ERROR: not all amino acid types read from rotamer library");
			System.exit(1);
		}
	}
	
	//reads in the rotamer volume data from volFilename
	public void loadVolFile (){
		
		try {
			readRotVol(volFilename);
		}
		catch (Exception e){
			System.out.println("Rotamer volumes file not found. Computing it..");
			computeAAVolumes(volFilename);
			try {
				readRotVol(volFilename);
			}
			catch (Exception e1){
				System.out.println("ERROR: "+e1);
				System.exit(1);
			}
		}
	}
	
	//Read in all of the rotamer volumes for all amino acids from the volFilename file
	private void readRotVol(String volFilename) throws Exception {
		
		FileInputStream is = new FileInputStream( volFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		
		is = new FileInputStream( volFilename );
		bufread = new BufferedReader(new InputStreamReader(is));
		
		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
  		
		int curResult = 0;
		while( curLine != null ) {
							
			AARotamerType aaType = getAAType(getToken(curLine,1));
			
				
			//Now we only store one volume per amino acid type
			int numRotVol = 1;//aaType.numRotamers();
			//if (numRotVol==0)
			//	numRotVol++;
			
			for (int j=0; j<numRotVol; j++)
				aaType.volume = (new Double(getToken(curLine,j+2)).doubleValue());
		
			curLine = bufread.readLine();
			curResult++;
		}
		
		bufread.close();
		
		if (curResult!=numAAallowed){
			System.out.println("ERROR: not all amino acid types read from rotamer volumes file");
			System.exit(1);
		}
	}	
	
	// Uses the VolModule class to calculate the volume of each rotamer of each amino acid
	public void computeAAVolumes(String volFileName) {

		Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
		
		Molecule m = new Molecule();
		Residue res = ppr.getResidue("ala"); 
		m.addResidue(0,res);
		VolModule sm = new VolModule(m);
		sm.setResidueTreatment(0,1);		

		//Data Structures for d-amino acids
		Molecule dm = new Molecule();
		Residue dres = ppr.getResidue("dala"); 
		dm.addResidue(0,dres);
		VolModule dsm = new VolModule(dm);
		dsm.setResidueTreatment(0,1);		
		
		
		PrintStream printStream = setupOutputFile(volFileName);

		AARotamerType aaTypes[] = getAAtypesAllowed();
		int numAAs = getNumAAallowed();
		
		Residue r = m.residue[0];
		int molNum = r.moleculeResidueNumber;
		
		Molecule curM;
		VolModule curVM;
		Residue curRes;
		for(int i=0;i<numAAs;i++){
			//Since we don't yet allow mutations from l-d amino acids, we have to switch
			//which molecule we are using for the different amino acids
			if(aaTypes[i].name.length() == 4 && aaTypes[i].name.startsWith("D")) {
				curM = dm;
				curVM = dsm;
			}else{
				curM = m;
				curVM = sm;
			}
			
			if(canMutate)
				//MutUtils.changeResidueType(curM,molNum,aaTypes[i].name,true);
			printStream.print(aaTypes[i].name + " ");
			System.out.println(aaTypes[i].name + " ");
			//for(Rotamer rot: aaTypes[i].rotamers){ 
			//All the volumes are very similar so just use the first rotamer
			
			Rotamer rot = aaTypes[i].rotamers.get(0); 
			MutUtils.applyRotamer(curM,rot,curM.residue[0]);
			double vol = curVM.getMoleculeVolume(0.25,0.0);
			printStream.print(vol + " ");
			System.out.println(vol + " ");
			aaTypes[i].volume = vol;
			//}
			printStream.println();
		}
		printStream.close();		
	}
	
	public double [][] getRotVol(){
		return rotamerVolumes;
	}

	// This function returns the rotamer index for rotamers of
	//  amino acid aaName; returns -1 if name not found
	public int getAARotamerIndex(String aaName) {

		// KER: Allow HID, HIE, and HIP to all be different now
//		if (aaName.equalsIgnoreCase("HID") || aaName.equalsIgnoreCase("HIE") || aaName.equalsIgnoreCase("HIP")) {
//			aaName = "HIS";
//			if(debug)
//				System.out.println("ASSUMING HID/E/P is " + aaName + " for rotamer purposes.");
//		}
		if (aaName.equalsIgnoreCase("CYX")) {
			aaName = "CYS";
			if(debug)
				System.out.println("ASSUMING CYX is " + aaName + " for rotamer purposes.");
		}

		for(int q=0;q<numAAallowed;q++) {
			if (aaTypes[q].name.equalsIgnoreCase(aaName))
				return q;
		}
		return -1;
	}
	
	// This function returns the rotamer index for rotamers of
	//  amino acid aaName; returns -1 if name not found
	public AARotamerType getAAType(String aaName) {

		// KER: Allow HID, HIE, and HIP to all be different now
		/*if (aaName.equalsIgnoreCase("HID") || aaName.equalsIgnoreCase("HIE") || aaName.equalsIgnoreCase("HIP")) {
			aaName = "HIS";
			if(debug)
				System.out.println("ASSUMING HID/E/P is " + aaName + " for rotamer purposes.");
		}*/
		if (aaName.equalsIgnoreCase("CYX")) {
			aaName = "CYS";
			if(debug)
				System.out.println("ASSUMING CYX is " + aaName + " for rotamer purposes.");
		}

		for(int q=0;q<aaTypes.length;q++) {
			if (aaTypes[q].name.equalsIgnoreCase(aaName))
				return aaTypes[q];
		}
		System.out.println("Could not find AA type: "+aaName);
		System.exit(0);
		return null;
	}

	// Returns the name of the amino acid at array index aaIndex
	public String getAAName(int aaIndex){
		return(aaTypes[aaIndex].name);
	}

	public AARotamerType getAAType(int aaIndex){
		return aaTypes[aaIndex];
	}
	
	
	// Returns the number of rotamers for the specified amino acid type
	public int getNumRotForAAtype(int aaTypeInd){
		return(aaTypes[aaTypeInd].numRotamers());
	}

	// This function returns the number of rotamers for a given
	//  amino acid type (by name)
	public int getNumRotamers(String aaName) {
		int aaNum = getAARotamerIndex(aaName);
		return(aaTypes[aaNum].numRotamers());
	}

	// Returns the number of dihedrals for the amino acid
	//  number aaTypeInd. NOTE: aaTypeInd is not an index
	//  into a molecule, it's 0..19
	public int getNumDihedrals(int aaTypeInd){
		if (aaTypeInd != -1)
			return(aaTypes[aaTypeInd].numDihedrals());
		else
			return(0);
	}
	
	// Returns the number of dihedrals for the amino acid
	//  number aaTypeInd. NOTE: aaTypeInd is not an index
	//  into a molecule, it's 0..19
	public int getNumDihedrals(String aaName){
		int aaNum = getAARotamerIndex(aaName);
		return(aaTypes[aaNum].numDihedrals());
	}


	// Returns the residue local atom numbers for dihedral dihedNum of
	//  residue resNum of strand strNum
	public int[] getDihedralInfo(Molecule m, int strNum, int resNum, int dihedNum){

		Residue localResidue = m.strand[strNum].residue[resNum];
		String tmpName = null;
		if(localResidue.name.equalsIgnoreCase("CYX"))
			tmpName = "CYS";
		else
			tmpName = localResidue.name;
				
		AARotamerType aaType = getAAType(tmpName);
		int atNum[] = new int[4];

		// Find atoms involved in the dihedral
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			for(int w=0;w<4;w++) {
				if (localResidue.atom[q].name.equalsIgnoreCase(aaType.dihedralAtomNames[dihedNum][w])) {
					atNum[w] = q;
				}
			}
		}
		return(atNum);
	}
	
	// Returns the residue local atom numbers for dihedral dihedNum of
	//  residue resNum of strand strNum
	public int[] getDihedralInfo(Molecule m, Residue localResidue, int dihedNum){

		//Residue localResidue = m.strand[strNum].residue[resNum];
		String tmpName = null;
		if(localResidue.name.equalsIgnoreCase("CYX"))
			tmpName = "CYS";
		else
			tmpName = localResidue.name;
				
		AARotamerType aaType = getAAType(tmpName);
		int atNum[] = new int[4];

		// Find atoms involved in the dihedral
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			for(int w=0;w<4;w++) {
				if (localResidue.atom[q].name.equalsIgnoreCase(aaType.dihedralAtomNames[dihedNum][w])) {
					atNum[w] = q;
				}
			}
		}
		return(atNum);
	}
	
	// Returns the residue local atom numbers for dihedral dihedNum of
	//  residue resNum of strand strNum
	public int[] getDihedralInfo(Residue res, int dihedNum){

		Residue localResidue = res;
		String tmpName = null;
		if(localResidue.name.equalsIgnoreCase("CYX"))
			tmpName = "CYS";
		else
			tmpName = localResidue.name;
				
		AARotamerType aaType = getAAType(tmpName);
		int atNum[] = new int[4];

		// Find atoms involved in the dihedral
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			for(int w=0;w<4;w++) {
				if (localResidue.atom[q].name.equalsIgnoreCase(aaType.dihedralAtomNames[dihedNum][w])) {
					atNum[w] = q;
				}
			}
		}
		return(atNum);
	}

	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
	
	public int getNumAAallowed(){
		return numAAallowed;
	}
	
	/*public float [][] getRotVol(){
		return rotamerVolumes;
	}*/
	
	/*public String getDihedralAtomNames(int i, int j, int k){
		return dihedralAtomNames[i][j][k];
	}
	
	public int getRotamerValues(int i, int j, int k){
		return rotamerValues[i][j][k];
	}*/
	
	/*public int [] getRotamerIndexOffset(){
		return rotamerIndexOffset;
	}*/
	
	/*public int getTotalNumRotamers(){
		return totalNumRotamers;
	}*/
	
	public AARotamerType [] getAAtypesAllowed(){
		return aaTypes;
	}
	
	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
	}	
	
	//KER: Adding this function so that I can dope the rotamers
	//	with the original rotamer from the input structure
	public Rotamer addRotamer(String AAname, String pdbNum, double[] dihedVals,boolean isWTrot){
		
		return addRotamer(AAname, pdbNum, dihedVals,null, isWTrot);
		
	}
	
	//KER: Adding this function so that I can dope the rotamers
	//	with the original rotamer from the input structure
	public Rotamer addRotamer(String AAname, String pdbNum,double[] dihedVals,double[] minimizationWidth, boolean isWTrot){

		AARotamerType aaType = getAAType(AAname);
		Rotamer newRot = new Rotamer(aaType.rotamers.size(), dihedVals, aaType, allRotamers.length,pdbNum,minimizationWidth,isWTrot);

		// PGC 2014: Optionally, a wildtype rotamer is only added a rotamer if it is beyond the minimization range of another rotamer. 
		if(aaType.addRot(newRot)){
			//aaType.addRot(dihedVals,allRotamers.length);
			Rotamer[] newArray = new Rotamer[allRotamers.length+1];
			System.arraycopy(allRotamers, 0, newArray, 0, allRotamers.length);
			newArray[newArray.length-1] = newRot;
			allRotamers = newArray;
	
			return newRot;
		}else{
			return null;
		}
	}
	
	//KER: Adding this function so that I can dope the rotamers
	//	with the original rotamer from the input structure
	public Rotamer addSubRotamer(String AAname, String pdbNum,double[] dihedVals,double[] minimizationWidth, Rotamer parent){

		AARotamerType aaType = getAAType(AAname);
		Rotamer newRot = Rotamer.newSubRot(aaType.rotamers.size(), dihedVals, aaType, allRotamers.length,pdbNum,minimizationWidth,parent);


		aaType.addRot(newRot);
		//aaType.addRot(dihedVals,allRotamers.length);
		Rotamer[] newArray = new Rotamer[allRotamers.length+1];
		System.arraycopy(allRotamers, 0, newArray, 0, allRotamers.length);
		newArray[newArray.length-1] = newRot;
		allRotamers = newArray;

		return newRot;
	}
	
		// PGC 2015: Load the Dunbrack Rotamers for each residue.  This function clears all rotamers from the 
	//   Backbone independent rotamer library except  the wildtype rotamer.
	// 	The array allflexibleResidues contains all residues.  This is needed to compute phi and psi.  
	public void loadDunbrackRots(Residue allResidues[]){
		
		// We need to remove all previous rotamers EXCEPT Ala and Gly.  It is convenient, however, to keep the wildtype rotamers, so that we do not have to add them again.		
		for (AARotamerType aatype : aaTypes){  // Each amino acid type stores rotamers
			aatype.clearTEMPLRotamers(); // Only backbone independent rotamers, different from Ala and gly are cleared.
		}
		// We will place all new rotamers in this arraylist
		ArrayList<Rotamer> dunbrackAllRots = new ArrayList<Rotamer>();
		// Also store wildtype rotamers in a different list for easier access.
		ArrayList<Rotamer> wtRots = new ArrayList<Rotamer>();
		for (Rotamer r: this.allRotamers){ // And we also store all rotamers in a RotamerLibrary class variable
			if(r.aaType.name.equals("ALA") || r.aaType.name.equals("GLY") || r.isWTrot){// Maintain the wildtype rotamer and GLY and ALA.
				r.rlIndex = dunbrackAllRots.size();
				dunbrackAllRots.add(r);
			}
			if(r.isWTrot){
				wtRots.add(r);
			}
		}
		this.allRotamers = null;
		// Add rotamers for every mutable residue.
		for(int resIx = 0; resIx < allResidues.length; resIx++){
			Residue residue = allResidues[resIx];

			
			if(residue.isMutable){				
				Residue previousRes = allResidues[resIx-1];
				Residue nextRes = allResidues[resIx+1];
				// Now add the dunbrack rotamers.
				ArrayList<Rotamer> rotamersForThisResidue = addDunbrackRotamers(residue, previousRes, nextRes);
				// Now, find the rotamer from the dunbrack library that is closest to the wildtype for this rotamer, and assign its probability to the wildtype rotamer. 
				// Search for the wildtype rotamer for this residue
				Rotamer wtRotForThisRes = null;
				for(Rotamer r: wtRots){
					if(r.isWTrot && (r.resSpecificRot.equals(residue.getResNumberString()))){
						wtRotForThisRes = r;
					}
				}
				// Only do this if the wildtype rotamer was included.
				if(wtRotForThisRes != null){
					int indexOfSmallestRMSD = -1;
					double smallestRMSD = Double.POSITIVE_INFINITY;
					// Compute the RMSD (with respect to dihedrals) of each rotamer with the wildtype rotamer, but ignore those residues of different amino acid identity.
					for(int rIx = 0; rIx < rotamersForThisResidue.size(); rIx++){							
						Rotamer r = rotamersForThisResidue.get(rIx);
						if(r.aaType.name.equals(wtRotForThisRes.aaType.name)){
							double myrmsd = Util.RootMeanSquareDihedrals(wtRotForThisRes.values, r.values);
							if(myrmsd < smallestRMSD){
								smallestRMSD = myrmsd;
								indexOfSmallestRMSD = rIx;
							}							
						}
					}
					// We limit the RMSD so that this is a completely different rotamer, we don't remove any rotamers.
					if(smallestRMSD < EnvironmentVars.REMOVE_CLOSEST_ROTAMER_RMS_CUTOFF){
						// Assign the probability of the closest rotamer to the wildtype rotamer.
						wtRotForThisRes.rotamerProbability = (rotamersForThisResidue.get(indexOfSmallestRMSD)).rotamerProbability;
						if(EnvironmentVars.REMOVE_CLOSEST_ROTAMER_WT_ROTAMER){
							// Remove the rotamer closest to this rotamer.
							rotamersForThisResidue.remove(indexOfSmallestRMSD);
						}

					}
						
					
				}
				// Add all remaining rotamers.
				for(Rotamer r: rotamersForThisResidue){
					AARotamerType myAAtype = getAAType(r.aaType.name);
					r.aaIndex = myAAtype.rotamers.size();
					myAAtype.addRot(r);
					r.rlIndex = dunbrackAllRots.size();
					dunbrackAllRots.add(r);
				}
			}
			
		}
		// Convert the arraylist to an array.
		this.allRotamers = dunbrackAllRots.toArray(new Rotamer[0]);
		
		
	}	
	// PGC 2015: Add all Dunbrack rotamers for residue res.  The previous and next residues are used to compute the backbone angles..  This function computes the phi and psi angles using 
	//   the hydrogen and the oxygen in the residue (i.e. instead of using the atoms of the previous and next residues).
	public ArrayList<Rotamer> addDunbrackRotamers(Residue res, Residue previousRes, Residue nextRes){
		// Get the backbone phi and psi
		Atom Cprev = previousRes.getAtomByName("C");
		Atom N = res.getAtomByName("N");
		Atom CA = res.getAtomByName("CA");
		Atom C = res.getAtomByName("C");
		Atom Nnext = nextRes.getAtomByName("N");
		// We first compute the backbone phi and psi angles because we need them to find the rotamers for this residue in the Dunbrack Library.		
		double Cprev_N_CA_C_torsion = C.torsion(Cprev, N, CA) ;
		double wtphi = Cprev_N_CA_C_torsion;
	
		double N_CA_C_Nnext_torsion = Nnext.torsion(N, CA, C);
		double wtpsi = N_CA_C_Nnext_torsion;
	
		ArrayList<Rotamer> newDunbrackRotamersForThisResidue = new ArrayList<Rotamer>();
		// Now we search the dunbrack rotamer library for the backbone of the wildtype residue, and add all rotamers in it.
		// Read the file name for the Dunbrack Backbone dependent rotamer library. 
		try{
			// Read in the dunbrack library for each residue.
			FileInputStream is = new FileInputStream( EnvironmentVars.dunbrackRotamerLibraryLocation);
			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;
			// Create an arraylist with the allowed amino acids at this residue.
			ArrayList<String> allowedAAatThisResAsString = new ArrayList<String>();
			for(AARotamerType aart: res.AATypesAllowed()){
				allowedAAatThisResAsString.add(aart.name);
			}
			// Round the wildtype phi and wildtype psi to the closest backbone grid (by dividing by 10, rounding, and multiplying by 10 again).
			long wtRoundedPhi = 10* Math.round(wtphi/10.0);
			long wtRoundedPsi = 10* Math.round(wtpsi/10.0);
			// Parse the whole file
			while((curLine = bufread.readLine())!= null){
				// Skip comment lines: those that start with #
				if(!curLine.startsWith("#")){
					// Convert the line into a char array.
					char dbLine [] = curLine.toCharArray(); 
					// First three letters are the amino acid.
					String aaNameNewRot = ""+dbLine[0]+dbLine[1]+dbLine[2];
					// Check that the amino acid is allowed for this position.					
					if(allowedAAatThisResAsString.contains(aaNameNewRot)) {
						// Search for the backbone phi and psi for this rotamer.
						// Phi is in chars 5, 6, 7, and 8 
						String phiForThisEntryAsStr = ""+dbLine[5]+dbLine[6]+dbLine[7]+dbLine[8];
						int phiForThisEntry = Integer.valueOf(phiForThisEntryAsStr.trim());
						// Psi is in chars 10, 11, 12, 13
						String psiForThisEntryAsStr = ""+dbLine[10]+dbLine[11]+dbLine[12]+dbLine[13];
						int psiForThisEntry = Integer.valueOf(psiForThisEntryAsStr.trim());
						if(phiForThisEntry == wtRoundedPhi && psiForThisEntry == wtRoundedPsi){	// We matched a backbone.
							// Fetch the amino acid type AARotamerType object;							
							AARotamerType myAAtype = getAAType(aaNameNewRot);
							// Parse the rotamer probability which is located at position 37 to 44.
							String probabilityStr = ""+dbLine[37]+dbLine[38]+dbLine[39]+dbLine[40]+dbLine[41]+dbLine[42]+dbLine[43]+dbLine[44];
							double myProbability = Double.parseDouble(probabilityStr.trim());
							// Parse the dihedrals in the new rotamer from the dunbrack rotamer library.
							int myNumDihedrals = myAAtype.numDihedrals();
							// The dihedrals of the new rotamer in the rotamer library.
							double myDihedrals[] = new double[myNumDihedrals]; 
							// The first char for each of the dihedrals is, respectively at 47, 55, 63, 71
							int dihedralIndices [] = {47, 55, 63, 71}; 
							for (int chi = 0; chi < myNumDihedrals; chi++){
								String dihedralStr = ""+dbLine[dihedralIndices[chi]]+dbLine[dihedralIndices[chi]+1]+dbLine[dihedralIndices[chi]+2]+
										dbLine[dihedralIndices[chi]+3]+dbLine[dihedralIndices[chi]+4]+dbLine[dihedralIndices[chi]+5];
								myDihedrals[chi] = Double.parseDouble(dihedralStr.trim());
							}							
							// Add a new Rotamer.	
							if(myProbability >= EnvironmentVars.DUNBRACK_PROBABILTY_CUTOFF){
								Rotamer newRot = new Rotamer(myAAtype.rotamers.size(), myDihedrals, myAAtype, 0, res.getResNumberString(), false, myProbability);
								newDunbrackRotamersForThisResidue.add(newRot);
								
								// Add the rotamer diversified by Chi1 in one standard deviation.
								if ( EnvironmentVars.DIVERSIFY_DUNBRACK_ROTAMERS_CHI1){
									// The standard deviation is in fields: 83, 84, 85, 86									
									String standardDeviationChi1Str =  ""+dbLine[83]+dbLine[84]+dbLine[85]+dbLine[86];
									double standardDeviationChi1 = Double.parseDouble(standardDeviationChi1Str);
									// First -Chi1
									// We need a new array for the dihedrals
									double mydihedrals_minusStdDev_Chi1[] = new double[myNumDihedrals];
									for (int dih = 0; dih < myNumDihedrals; dih++){
										mydihedrals_minusStdDev_Chi1[dih] = myDihedrals[dih];
									}
									if(mydihedrals_minusStdDev_Chi1[0]-standardDeviationChi1 < -180 ){ 
										mydihedrals_minusStdDev_Chi1[0] = mydihedrals_minusStdDev_Chi1[0]-standardDeviationChi1+360;
									}
									else{
										mydihedrals_minusStdDev_Chi1[0] -=  standardDeviationChi1;
									}
									Rotamer newRot_MinusStdDev_Chi1 = new Rotamer(myAAtype.rotamers.size(), mydihedrals_minusStdDev_Chi1, myAAtype, 0, res.getResNumberString(), false, myProbability);
									newDunbrackRotamersForThisResidue.add(newRot_MinusStdDev_Chi1);
									// Now +chi1
									double mydihedrals_plusStdDev_Chi1[] = new double[myNumDihedrals];
									for (int dih = 0; dih < myNumDihedrals; dih++){
										mydihedrals_plusStdDev_Chi1[dih] = myDihedrals[dih];
									}									
									if(mydihedrals_plusStdDev_Chi1[0]+standardDeviationChi1 > 180 ){ 
										mydihedrals_plusStdDev_Chi1[0] = mydihedrals_plusStdDev_Chi1[0]+standardDeviationChi1-360;
									}
									else{
										mydihedrals_plusStdDev_Chi1[0] +=  standardDeviationChi1;
									}
									Rotamer newRot_PlusStdDev_Chi1 = new Rotamer(myAAtype.rotamers.size(), mydihedrals_plusStdDev_Chi1, myAAtype, 0, res.getResNumberString(), false, myProbability);
									newDunbrackRotamersForThisResidue.add(newRot_PlusStdDev_Chi1);

								}

									
							}
						}
					}
				}
			}
			is.close();
		}
		catch(Exception e){
			System.out.println("Problem reading the dunbrack rotamer library.");
			System.exit(1);
		}

		return newDunbrackRotamersForThisResidue;
	}	
	// PGC 2015: Add rotamers for the Dunbrack library; rotamers have the probability field.  
	//	 Returns null if rotamer could not be added.
	public Rotamer addRotamer(String AAname, String pdbNum,double[] dihedVals, boolean isWTrot, double rotamerProbability){

		AARotamerType aaType = getAAType(AAname);
		Rotamer newRot = new Rotamer(aaType.rotamers.size(), dihedVals, aaType, allRotamers.length,pdbNum,isWTrot, rotamerProbability);

		if(aaType.addRot(newRot)){
			Rotamer[] newArray = new Rotamer[allRotamers.length+1];
			System.arraycopy(allRotamers, 0, newArray, 0, allRotamers.length);
			newArray[newArray.length-1] = newRot;
			allRotamers = newArray;
	
			return newRot;
		}else{
			return null;
		}
	}
	public Rotamer addOrigRot(Residue res){
		Rotamer wtRot = null;
		if( !res.addedOrigRot){
			//Get Num Dihedrals
			int numDiheds = getNumDihedrals(getAARotamerIndex(res.name));
			if(numDiheds>0){
				double[] diheds = new double[numDiheds]; 
				int atoms[] = new int[4];
				//get all dihedrals
				for(int i=0; i<numDiheds; i++){
					atoms = getDihedralInfo(res, i);
					diheds[i] = res.atom[atoms[3]].torsion(res.atom[atoms[0]], res.atom[atoms[1]], res.atom[atoms[2]]);
				}
				wtRot = addRotamer(res.name, res.getResNumberString(), diheds, true); 
				// PGC 2015: Set to true to diversify X1 or X1,X2				
				boolean X1diversify = false;
				boolean X2diversify = false;
				// PGC2015: Optionally, add the rotamers with +9 and -9 degrees for X1.
				if(X1diversify || (X2diversify && numDiheds == 1)){
					for (double offset = -9; offset <= 9; offset += 18){
						double tmpDiheds []= new double[numDiheds];
						for(int i=0; i<numDiheds; i++){						
							tmpDiheds[i] = diheds[i];
						}
						tmpDiheds[0] += offset;
						addRotamer(res.name, res.getResNumberString(), tmpDiheds, true);
					}		
				}
				// PGC2015: Optionally, add the rotamers with +9 and -9 degrees for X1 and X2.
				if(numDiheds > 1 && X2diversify){
					for (double offset1 = -9; offset1 <= 9; offset1 += 9){
						for (double offset2 = -9; offset2 <= 9; offset2 += 9){
							double tmpDiheds []= new double[numDiheds];
							for(int i=0; i<numDiheds; i++){						
								tmpDiheds[i] = diheds[i];
							}
							tmpDiheds[0] += offset1;
							tmpDiheds[1] += offset2;
							if (offset1 != 0 || offset2 != 0){
								addRotamer(res.name, res.getResNumberString(), tmpDiheds, false);
							}
						}
					}
				}
			}
			res.addedOrigRot = true;
		}
		else{
			wtRot = res.wtRot;
		}
		
		return wtRot;
	}
	
	public static void initThree2One(){
		three2one = new HashMap<String,String>();
		three2one.put("ALA","A");
		three2one.put("CYS","C");
		three2one.put("CYX","C");
		three2one.put("ASP","D");
		three2one.put("GLU","E");
		three2one.put("PHE","F");
		three2one.put("GLY","G");
		three2one.put("HIS","H");
		three2one.put("HIP","H");
		three2one.put("HIE","H");
		three2one.put("HID","H");
		three2one.put("ILE","I");
		three2one.put("LYS","K");
		three2one.put("LEU","L");
		three2one.put("MET","M");
		three2one.put("ASN","N");
		three2one.put("PRO","P");
		three2one.put("GLN","Q");
		three2one.put("ARG","R");
		three2one.put("SER","S");
		three2one.put("THR","T");
		three2one.put("VAL","V");
		three2one.put("TRP","W");
		three2one.put("TYR","Y");
	}
	
	public static String getOneLet(String aa3Name){
		if(aa3Name.length() == 4 && aa3Name.startsWith("D"))
			aa3Name = aa3Name.substring(1);
		
		String res = three2one.get(aa3Name);
		if (res == null)
			res = "X";
		return res;
	}

	public Rotamer getRot(int i) {
		return allRotamers[i];
	}

	public String getRotFile(){
		return rotFile;
	}

	public Rotamer getAARot(String AAname, int i) {
		AARotamerType aaType = getAAType(AAname);
		return aaType.rotamers.get(i);
	}
	
	public Rotamer getAARot(int AAind, int i) {
		AARotamerType aaType = getAAType(getAAName(AAind));
		return aaType.rotamers.get(i);
	}
	
	public void loadGlobalRots(String rotFilename) {

		// PGC 2014: always read in all rotamers.
		//int numRotLoaded = allRotamers.length;
		int numRotLoaded = 0;

		

		// HANDLE THE NORMAL AAs	
		try{
			ArrayList<Rotamer> allRots = new ArrayList<Rotamer>();
			FileInputStream is = new FileInputStream( rotFilename );
			// First, remove all rotamers for all amino acid types. 
			for (AARotamerType aatype : aaTypes){
				aatype.clearRotamers();
			}
			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = bufread.readLine();
			
			
			while( curLine != null ) {
				// Skip over comments (lines starting with !)
				if(curLine.charAt(0) == '!'){
					curLine = bufread.readLine();
					continue;
				}

				// Read the "global"index of the current rotamer, which is the first column.
				int globalRotNum = (new Integer(getToken(curLine,1))).intValue();
				// This is deprecated; in case you don't want to reload rotamers; 
				if(globalRotNum >= numRotLoaded){
					// AA name (e.g. 'ALA')
					String aaName = getToken(curLine,2);
					AARotamerType aa = getAAType(aaName);
					// index of the rotamer for the amino acid
					int aaRotNum = (new Integer(getToken(curLine,3))).intValue();


					// Read in the rotamer
					double[] dihedrals = new double[aa.numDihedrals()];
					double[] minLimits = new double[aa.numDihedrals()];
					// Read dihedrals and minimization boundaries for each dihedral.
					for(int w=0;w<aa.numDihedrals();w++) {
						dihedrals[w] = (new Double(getToken(curLine,(w+4)))).doubleValue();
						minLimits[w] = (new Double(getToken(curLine,(w+4+aa.numDihedrals())))).doubleValue();
					}
					
					// PDB Number for residue specific rotamer, otherwise -1 for global rot;
					String pdbNum = getToken(curLine,4+2*aa.numDihedrals());

					// Is this a wildtype rotamer?					
					boolean isWTrot = false;
					String wtFlag = getToken(curLine,5+2*aa.numDihedrals());
					if(wtFlag.equals("WT")){
						isWTrot = true;
					}
					
					Rotamer curRot = new Rotamer(aaRotNum,dihedrals,aa,globalRotNum,pdbNum,minLimits, isWTrot);
					// Add it first to the amino acid type
					aa.addRot(curRot);
					// And to the global rotamer list.
					allRots.add(curRot);
					
					
				}
				curLine = bufread.readLine();
			}
			bufread.close();
		
			// PGC 2014: replacing the old rotamer library, no matter what it had.		
			allRotamers = allRots.toArray(new Rotamer[0]);

		}catch(Exception E){
			System.out.println("Couldn't properly load rot library.");
		}

	}

	/*
	 * Goes through and finds the rotamer the split rotamer came from
	 */
	public Rotamer remapSplitRot(int globalRotNum){
		
		Rotamer splitRot = getRot(globalRotNum);
		
		for(int i=0; i<=183; i++){
			Rotamer r = getRot(i);
			if(splitRot.aaType.name.equalsIgnoreCase(r.aaType.name)){
				boolean dihMatch = true;
				for(int dih=0; dih<r.aaType.numDihedrals();dih++){
					if(Math.abs(r.values[dih] - splitRot.values[dih]) > 9 ){
						dihMatch = false;
					}
				}
				if(dihMatch)
					return r;
			}
		}
		
		return null;
	}

	/**
	 * Saves all of the rotamers in the rotamer library.
	 * Can be read with loadGlobalRots
	 * @param fileName
	 * @return
	 */
	public void save(String fileName){
		try {
			
				FileOutputStream fileOutputStream = new FileOutputStream(fileName);
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				PrintStream logPS = new PrintStream( bufferedOutputStream );
				for(Rotamer r: allRotamers){
					logPS.println(r.toString());
				}
				logPS.close();
			
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
	}

	public int getTotalNumRotamers() {
		return allRotamers.length;
	}
	
}
