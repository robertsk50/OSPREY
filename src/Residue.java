/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2009 Bruce Donald Lab, Duke University

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

	<signature of Bruce Donald>, 12 Apr, 2009
	Bruce Donald, Professor of Computer Science
 */

////////////////////////////////////////////////////////////////////////////////////////////////
// Residue.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2001-2004);
 * 		other changes made by Ivelin Georgiev (2004-2009)
 *
 * Rewritten by Ryan Lilien based on code by Neill White
 * Many functions have been added, others removed, most have had
 *  at least some parts rewritten. Code rewrites have often been
 *  major to fix bugs or add functionality.
 *
 * Based on software copyrighted, 1999, by Neill White.
 *  The author hereby grants permission to use, copy, modify, and re-distribute
 *  this software and its documentation for any purpose, provided
 *  that existing copyright notices are retained in all copies and that this
 *  notice is included verbatim in any distributions. No written agreement,
 *  license, or royalty fee is required for any of the authorized uses.
 *  Modifications to this software may be distributed provided that
 *  the nature of the modifications are clearly indicated.
 *
 */

import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.io.Serializable;

/**
 * Handles functions and data associated with residues.
 */
public class Residue implements Serializable {

	String name = "";			// The short residue name "GLU"
	String fullName = "";	// The long version of the residue name ie "GLU A  47"
	String defaultName = ""; //The long version of the WT residue name ie "GLU A 47"
	int	numberOfAtoms=0;	// Number of atoms
	int	moleculeResidueNumber=-1;	// The molecule relative residue number
	int	strandResidueNumber=-1;		// The strand relative residue number
	int	strandNumber=-1;	// The number of the strand containing this residue
	Atom	atom[];					// Array of atoms in this residue
	private boolean energyEvalSC = true;	// Should the residue side-chain be included in computing energy
	private boolean energyEvalBB = true;	// Should the residue backbone be included in computing energy
	boolean ffAssigned = false;	// Are forcefield types assigned for atoms in this residue
	boolean flexible = false;	// Is this residue flexible
	boolean cterm =  false; //Is a cterminal residue
	boolean nterm = false; //Is an nterminal residue
	public boolean cofactor = false;
	byte SStype; //Secondary structure type (Set by readDSSP for hbond functionality) Duplicate of the "secondaryStruct" int below that reads PDB entries
	public boolean lAmino = true;
	double[] defaultCB; //Store the WT CB coordinates for when mutating back from Gly
	
	//DEEPer stuff                       
	int perts[]=null;             // Perturbations to which this unit is subject (may be empty; this array consists of indices in m.perts, in ascending order)
	int pertStates[][]=null;      //Defines the perturbation states for this residue:
	//pertStates[perturbation state #][perturbation #] gives which default param. value or interval of the given perturbation is in the given perturbation state of this residue
	//Each residue conformation is defined by a rotameric state and a perturbation state
	//pertStates[0] must contain the unperturbed state

	int affectedPerts[] = null;//Perturbations affected by the conformation of this residue
	//(as opposed to perturbations that can affect by the conformation of this residue, which are in perts,
	//A perturbation is in affectedPerts if it is not in perts but comes after one in perts and overlaps with it at some other residue (or so on recursively)


	boolean pucker = false;//If the residue is a proline this will be ProlineFlip.UP or ProlineFlip.DOWN
	boolean validConf = true;//This is true unless there's an uncloseable proline ring or some other conformational problem



	int secondaryStruct;
	//Types of secondary structure
	final static int HELIX = 0;
	final static int SHEET = 1;
	final static int LOOP = 2;
	
	// PGC 2014: We store the wildtype rotamer coordinates as part of the residue.
	Atom wildtypeAtoms[];

	/*
	 * KER: I use these variables to keep track of what this molecule can mutate to. Since I
	 * want position specific rotamers, each residue has to keep track of the allowed rotamers
	 * and can't just rely on all the rotamers available in the AAtype
	 */
	private ArrayList<AARotamerType> AATypesAllowed = null; //Changed the name so I know it can't be accessed anywhere
	private ArrayList<ResidueConformation> allowedRCs = null;
	public boolean mutatedOnce = false;
	public boolean addedOrigRot = false;
	public ResidueConformation curRC;
	public ResidueConformation wtRC;
	public Rotamer wtRot;
	
	RotamerLibrary rl;
	public boolean canMutate = false; //Is the residue type a type that can mutate (i.e. of type protein).
	public boolean isMutable = false; //Is this residue a mutable residue
	int origMutPos;
	public String defaultAA;
	
	
	Residue(){
		init(); 
	}

	Residue(String resname){
		init();
		name = resname;
	}

	public void init(){
		atom = new Atom[1];
		perts = new int[0];
		affectedPerts = new int[0];
		pertStates = new int[0][];
	}

	// Displays some residue info to System.out for debugging
	public void printResidueInfo(){
		System.out.println("name = *"+name+"*");
		System.out.println("fullName = *"+fullName+"*");
		System.out.println("numberOfAtoms = *"+numberOfAtoms+"*");
		System.out.println("moleculeResidueNumber = *"+moleculeResidueNumber+"*");
		System.out.println("strandResidueNumber = *"+strandResidueNumber+"*");
		System.out.println("strandNumber = *"+strandNumber+"*");
	}

	// Adds specified atom to this residue
	// Appropriately updates all Residue and Atom fields except for
	//  Atom.moleculeAtomNumber which should be updated in the molecule class
	// *This function should rarely be called. Atoms should be added via the
	//  molecule containing this residue. This function is important however
	//  as it is called by the molecule and strand addAtom functions
	public int addAtom(Atom newAtom){
		int newAtomNumber = numberOfAtoms + 1;
		int newAtomNumberx3 = newAtomNumber * 3;

		Atom largerAtomArray[] = new Atom[newAtomNumber];
		System.arraycopy(atom, 0, largerAtomArray, 0, atom.length);
		atom = largerAtomArray;
		atom[numberOfAtoms] = newAtom;

		newAtom.residueAtomNumber = numberOfAtoms;
		newAtom.moleculeResidueNumber = moleculeResidueNumber;
		newAtom.strandResidueNumber = strandResidueNumber;
		newAtom.strandNumber = strandNumber;

		return(numberOfAtoms++);
	}

	// Deletes the specified atom
	// *This function should rarely be called. Atoms should be deleted via the
	//  molecule containing this residue. This function is important however
	//  as it is called by the molecule and strand deleteAtom functions
	public int deleteAtom(int atomNumber){
		int newAtomNumber = numberOfAtoms - 1;
		int newAtomNumberx3 = newAtomNumber * 3;
		int atomNumberx3 = atomNumber * 3;
		int moleculeAtomNumber = atom[atomNumber].moleculeAtomNumber;

		// If there is no molecule atom number then this is a lone
		//  residue (ie. it's not in a molecule, so use the simple
		//  numbering as the moleculeAtomNumber
		if (moleculeAtomNumber == -1)
			moleculeAtomNumber = atomNumber;
		// update residueAtomNumbers of higher numbered atoms
		for(int i=atomNumber; i<numberOfAtoms; i++)
			atom[i].residueAtomNumber -= 1;

		Atom smallerAtomArray[] = new Atom[numberOfAtoms-1];
		System.arraycopy(atom, 0, smallerAtomArray, 0, atomNumber);
		if (atomNumber<newAtomNumber)
			System.arraycopy(atom, atomNumber+1, smallerAtomArray,
					atomNumber, atom.length-atomNumber-1);
		atom = smallerAtomArray;

		numberOfAtoms--;

		// update bond indices
		for(int i=0;i<numberOfAtoms;i++) {
			for(int m=0; m<atom[i].numberOfBonds; m++){
				if (atom[i].bond[m] == moleculeAtomNumber){
					atom[i].deleteBond(m);
					m--;
				}
				else if (atom[i].bond[m] > moleculeAtomNumber){
					atom[i].bond[m] -= 1;
				}
			}
		}

		return(numberOfAtoms);
	}

	// This function takes a residue and renumbers its
	//  internal indices such that the first atom is atom 0;
	// It adjusts atom members moleculeAtomNumber and residueAtomNumber;
	// It resets atom bond[] and numberOfBonds, so the bond information must be updated
	//  using the determineBonds() method from Molecule.java;
	// This function is used when we pull a residue out
	//  from one molecule and we want to add it back
	//  into another one.
	public void renumberResidue(){
		for(int i=0;i<numberOfAtoms;i++){
			atom[i].moleculeAtomNumber = i;
			atom[i].residueAtomNumber = i;
			atom[i].bond = null;
			atom[i].numberOfBonds = 0;
		}
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
				System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())
			return(st.nextToken());
		return(new String(""));
	} // end getToken

	// This function returns the residue number (not the sequential one,
	//  but the 'true' one from the pdb file).  If the pdb doesn't have
	//  one (unlikely) then use the sequential numbering
	//  PGC 2013: Residue number from the BD is now a String
	public String getResNumber() {
		if (fullName.length() > 5)
			return( (getToken(fullName.substring(5),1)) );
		return Integer.toString(moleculeResidueNumber+1);
	}

	// This function rotates the specified atoms in the residue
	//  by thetaDeg degrees around axis dx, dy, dz
	// at1 is the pivot atom (the 3rd atom if we were doing dihedrals)
	//  where the 4th atom+ would be the moving atoms
	// atomList is a list of atom numbers to rotate
	// numAtoms is the number of atoms to rotate (size of atomList)
	public void rotateResidue(Atom at1, double dx, double dy,
			double dz, double thetaDeg, int atomList[], int numAtoms) {

		double fx,fy,fz, tx,ty,tz;
		fx = (new Double(dx)).doubleValue();
		fy = (new Double(dy)).doubleValue();
		fz = (new Double(dz)).doubleValue();
		Atom pseudoAtom = new Atom("a", fx, fy, fz);

		int atomNumber;
		if (numAtoms == 0)
			return;

		int numberOfCoordinatesx3 = numAtoms*3;
		double temporaryCoordinates[] = new double[numberOfCoordinatesx3];
		int qx3;
		for(int q=0;q<numAtoms;q++) {
			qx3 = q*3;
			temporaryCoordinates[qx3+0]=atom[atomList[q]].coord[0] - at1.coord[0];
			temporaryCoordinates[qx3+1]=atom[atomList[q]].coord[1] - at1.coord[1];
			temporaryCoordinates[qx3+2]=atom[atomList[q]].coord[2] - at1.coord[2];
		}

		double[][] rot_mtx = new double[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(fx,fy,fz,(double) thetaDeg,rot_mtx);

		for(int q=0;q<numAtoms;q++) {
			qx3 = q*3;
			tx = temporaryCoordinates[qx3+0];
			ty = temporaryCoordinates[qx3+1];
			tz = temporaryCoordinates[qx3+2];

			atom[atomList[q]].coord[0] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + at1.coord[0];
			atom[atomList[q]].coord[1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + at1.coord[1];
			atom[atomList[q]].coord[2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + at1.coord[2];
		}

	}

	//Returns the distance (the minimum distance between a pair of atoms) between the side-chains of this residue and res2;
	//If bbAt is false, then only side-chain atoms are considered
	public double getDistCountH(Residue res2, boolean bbAt){
		
		double minDist = Double.POSITIVE_INFINITY;
		
		for (int a1=0; a1<numberOfAtoms; a1++){
			
			//String a1type = atom[a1].elementType;
			
			if ( /*(!a1type.equalsIgnoreCase("H")) &&*/ ( bbAt || (!atom[a1].getIsBBatom()) ) ){
				
				for (int a2=0; a2<res2.numberOfAtoms; a2++){
					
					//String a2type = res2.atom[a2].elementType;
					
					if ( /*(!a2type.equalsIgnoreCase("H")) &&*/ (bbAt || (!res2.atom[a2].getIsBBatom()) ) ){
						
						double curDist = getDist(atom[a1],res2.atom[a2]);
						minDist = Math.min(minDist, curDist);
					}
				}
			}
		}
		
		return minDist;
	}
	
	
	//Returns the distance (the minimum distance between a pair of non-hydrogen side-chain atoms) between the side-chains of this residue and res2;
	//If bbAt is false, then only side-chain atoms are considered
	public double getDist(Residue res2, boolean bbAt){

		double minDist = (double)Math.pow(10, 10);

		for (int a1=0; a1<numberOfAtoms; a1++){

			String a1type = atom[a1].elementType;

			if ( (!a1type.equalsIgnoreCase("H")) && ( bbAt || (!atom[a1].getIsBBatom()) ) ){

				for (int a2=0; a2<res2.numberOfAtoms; a2++){

					String a2type = res2.atom[a2].elementType;

					if ( (!a2type.equalsIgnoreCase("H")) && (bbAt || (!res2.atom[a2].getIsBBatom()) ) ){

						double curDist = getDist(atom[a1],res2.atom[a2]);
						minDist = Math.min(minDist, curDist);
					}
				}
			}
		}

		return minDist;
	}

	//Returns the distance between the two atoms
	private double getDist(Atom a1, Atom a2){

		double rijx, rijy, rijz, rij, rij2;

		rijx = a1.coord[0] - a2.coord[0];
		rijy = a1.coord[1] - a2.coord[1];
		rijz = a1.coord[2] - a2.coord[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		rij = (double)Math.sqrt(rij2);

		return rij;
	}

	public void setEnergyEval(boolean scEval, boolean bbEval){
		energyEvalSC = scEval;
		energyEvalBB = bbEval;
	}

	public void setSCEnergyEval(boolean scEval){
		energyEvalSC = scEval;
	}

	public boolean getEnergyEvalSC(){
		return energyEvalSC;
	}

	public boolean getEnergyEvalBB(){
		return energyEvalBB;
	}

	//Returns the moleculeAtomNumber for the atom of this residue that has the same name as the parameter n
	public int getAtomNameToMolnum(String n){

		for (int i=0; i<numberOfAtoms; i++){
			if (atom[i].name.equalsIgnoreCase(n))
				return atom[i].moleculeAtomNumber;
		}
		return -1;
	}

	// This function returns the residue number (not the sequential one,
	//  but the 'true' one from the pdb file).  If the pdb doesn't have
	//  one (unlikely) then use the sequential numbering
	public String getResNumberString() {
		if (fullName.length() > 5)
			return( getToken(fullName.substring(5),1) ) ;
		return new Integer(moleculeResidueNumber+1).toString();
	}


	//DEEPer functions

	//Returns the first (and thus presumably the only) atom in the unit named n
	//Return null if there is no atom with this name
	public Atom getAtomByName(String n){

		for (int i=0; i<numberOfAtoms; i++){
			if (atom[i].name.equalsIgnoreCase(n))
				return atom[i];
		}
		return null;
	}


	//Get an array of molecule atom numbers for the specified parts of the residue
	//SC (sidechain) classified as in Atom.setIsBBAtom() (so alpha hydrogens are counted as sidechain)
	//CA is just the alpha carbon
	//Proline is different: we treat the amide group and the CD and its hydrogens together as a rigid
	//body during perturbations, so the CD and its hydrogens are included in the amide atom list
	public int[] getAtomList(boolean includeAmide, boolean includeSC,
			boolean includeCA, boolean includeCarbonyl){

		int listSize = 0;

		if( getAtomByName("H3") != null ){//Checking if N-terminal 
			if(includeAmide)
				listSize += 4;
			if(includeSC)
				listSize += atom.length - 7;
			if(includeCarbonyl)
				listSize += 2;
		}
		else if( getAtomByName("OXT") != null ){//Checking if C-terminal
			if(includeAmide)
				listSize += 2;
			if(includeSC)
				listSize += atom.length - 6;
			if(includeCarbonyl)
				listSize += 3;
		}
		else{
			if(includeAmide)
				listSize += 2;
			if(includeSC)
				listSize += atom.length - 5;
			if(includeCarbonyl)
				listSize += 2;
		}
		if(includeCA)
			listSize += 1;//Same for any kind of residue

		if( name.equalsIgnoreCase("PRO") ){
			//Two more amide atoms; thus also two less sidechain atoms since the
			//number of SC atoms was calculated using atom.length
			if(includeAmide)
				listSize+=2;
			if(includeSC)
				listSize-=2;
		}

		int atomList[] = new int[listSize];
		int a=0;//Index in atomList

		if(includeAmide){
			for(int b=0;b<atom.length;b++){
				if( ( atom[b].name.equalsIgnoreCase("N") ) || ( atom[b].name.equalsIgnoreCase("H") )
						|| ( atom[b].name.equalsIgnoreCase("H1") ) || ( atom[b].name.equalsIgnoreCase("H2") )
						|| ( atom[b].name.equalsIgnoreCase("H3") ) ){
					atomList[a] = atom[b].moleculeAtomNumber;
					a++;
				}

				else if( name.equalsIgnoreCase("PRO") ){
					if( ( atom[b].name.equalsIgnoreCase("CD") ) || ( atom[b].name.contains("HD") ) ){
						atomList[a] = atom[b].moleculeAtomNumber;
						a++;
					}
				}

			}
		}
		if(includeSC){
			for(int b=0;b<atom.length;b++){

				if( !atom[b].isBBatom ){

					if( name.equalsIgnoreCase("PRO") ){
						if( ( atom[b].name.equalsIgnoreCase("CD") ) || ( atom[b].name.contains("HD") ) )
							continue;//Don't count CD or HD's in the sidechain for proline
					}

					atomList[a] = atom[b].moleculeAtomNumber;
					a++;
				}
			}
		}
		if(includeCA){
			atomList[a] = getAtomByName("CA").moleculeAtomNumber;
			a++;
		}
		if(includeCarbonyl){
			for(int b=0;b<atom.length;b++){
				if( ( atom[b].name.equalsIgnoreCase("C") ) || ( atom[b].name.equalsIgnoreCase("O") )
						|| ( atom[b].name.equalsIgnoreCase("OXT") ) ){
					atomList[a] = atom[b].moleculeAtomNumber;
					a++;
				}
			}
		}

		return atomList;

	}

	// PGC 2013: The following methods will allow the comparison of two residues by their
	//      PDB number even if they have Kabat numbering.
	// Returns true if "r2_pdbResNum" precedes r1 anywhere in the PDB chain. Used for helix purposes.
	public static boolean lessThanInPDBChain(String r2_pdbResNum, Residue r1){
		String r1_pdbResNum = r1.getResNumber();
		int r1_integerPart = getResNumIntegerPart(r1_pdbResNum);//Integer.valueOf(r1_pdbResNum.split("[0-9]")[0]);
		String r1_kabatCharacter = getResNumKabatLetter(r1_pdbResNum);//r1_pdbResNum.split("[A-Z]*")[1];
		int r2_integerPart = getResNumIntegerPart(r2_pdbResNum);//Integer.valueOf(r2_pdbResNum.split("[0-9]")[0]);
		String r2_kabatCharacter = getResNumKabatLetter(r2_pdbResNum);//r1_pdbResNum.split("[A-Z]*")[1];
		if(r2_integerPart < r1_integerPart){
			return true;
		}
		else if(r1_integerPart < r2_integerPart){
			return false;
		}
		else if( r2_kabatCharacter.compareTo(r1_kabatCharacter) == -1){
			return true;
		}
		else{
			return false;
		}

	}
	// Returns true if r2_pdbResNum precedes r1 in the PDB chain. (e.g. residue 19 precedes residue 20g)
	public static boolean lessThanOrEqualInPDBChain(String r2_pdbResNum, Residue r1){
		String r1_pdbResNum = r1.getResNumber();
		if (r2_pdbResNum.compareTo(r1_pdbResNum)  == 0){
			return true;
		}
		else{
			return lessThanInPDBChain(r2_pdbResNum, r1);
		}
	}

	//splitting a residue number in string form into integer and character parts
	public static int getResNumIntegerPart(String resNumString){
		String resNumTrimmed = resNumString.trim();
		if( Character.isDigit( resNumTrimmed.charAt(resNumTrimmed.length()-1) ) )//resNumTrimmed ends with a digit
			return Integer.valueOf(resNumTrimmed);//so it's all the integer part
		else//remove last character, take the rest as integer part
			return Integer.valueOf( resNumTrimmed.substring(0, resNumTrimmed.length()-1) );
	}

	public static String getResNumKabatLetter(String resNumString){
		String resNumTrimmed = resNumString.trim();
		if( Character.isDigit( resNumTrimmed.charAt(resNumTrimmed.length()-1) ) )//resNumTrimmed ends with a digit
			return "";
		else//remove last character, take the rest as integer part
			return resNumTrimmed.substring(resNumTrimmed.length()-1, resNumTrimmed.length());
	}

	public void assignHandedness() {

		//Find the CA,N,C,CB atoms and calculate the dihedral between them
		Atom CA = null;
		Atom N = null;
		Atom C = null;
		Atom CB = null;
		for(Atom a: atom){
			if(a.name.equalsIgnoreCase("CA")){
				CA = a;
			}
			else if(a.name.equalsIgnoreCase("N")){
				N = a;
			}
			else if(a.name.equalsIgnoreCase("C")){
				C = a;
			}
			else if(a.name.equalsIgnoreCase("CB")){
				CB = a;
			}else if( (name.equalsIgnoreCase("GLY") || name.equalsIgnoreCase("DGLY")) && (a.name.equalsIgnoreCase("HA3") || a.name.equalsIgnoreCase("3HA"))){
				CB = a;
			}
		}

		double torsion = 1;
		if(CA != null && N != null && C != null)
			torsion = CB.torsion(CA, N, C);

		if(torsion > 0)
			lAmino  = true;
		else{
			lAmino = false; //dAmino acid
			//    				System.out.println("Residue: "+fullName+" is a d-amino acid!!");

			if(name.length() ==3){
				name = "D" + name;
			}
		}
	}


	public void reflect(){
		Atom C = null;
		Atom CA = null;
		Atom N = null;

		//if(!name.equalsIgnoreCase("GLY")){
		for(Atom a: atom){
			if(a.name.equalsIgnoreCase("CA")){
				CA = a;
			}
			else if(a.name.equalsIgnoreCase("N")){
				N = a;
			}
			else if(a.name.equalsIgnoreCase("C")){
				C = a;
			}
		}

		double[] C_CA = new double[3];
		double[] N_CA = new double[3];

		C_CA[0] = C.coord[0] - CA.coord[0];
		C_CA[1] = C.coord[1] - CA.coord[1];
		C_CA[2] = C.coord[2] - CA.coord[2];

		N_CA[0] = N.coord[0] - CA.coord[0];
		N_CA[1] = N.coord[1] - CA.coord[1];
		N_CA[2] = N.coord[2] - CA.coord[2];

		double[] normal = MathUtils.cross(C_CA, N_CA);

		double denom = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];
		//denom = (float) Math.sqrt(denom);

		double d = -1 * (normal[0]*CA.coord[0]+normal[1]*CA.coord[1]+normal[2]*CA.coord[2]);


		double[] tmpCoord = new double[3];
		for(Atom a: atom){
			if(!a.isBBatom || a.name.equals("HA")){

				tmpCoord[0] = a.coord[0] - normal[0]*2*(normal[0]*a.coord[0]+normal[1]*a.coord[1]+normal[2]*a.coord[2]+d)/denom;
				tmpCoord[1] = a.coord[1] - normal[1]*2*(normal[0]*a.coord[0]+normal[1]*a.coord[1]+normal[2]*a.coord[2]+d)/denom;
				tmpCoord[2] = a.coord[2] - normal[2]*2*(normal[0]*a.coord[0]+normal[1]*a.coord[1]+normal[2]*a.coord[2]+d)/denom;
				a.coord[0] = tmpCoord[0];
				a.coord[1] = tmpCoord[1];
				a.coord[2] = tmpCoord[2];
			}
		}

		lAmino = !lAmino; //We reflected the amino so it changes its chirality

	}

	public String threeLet(){
		if(!lAmino){
			return name.substring(1);
		}
		return name;
	}

	public int numAllowedAATypes() {
		return AATypesAllowed.size();
	}

	/**
	 * Given the amino acid type index (0 .. AAtypesAllowed.length)
	 * return the rotamers for that amino acid type that are allowed
	 * 
	 * @param i AA allowed index
	 * @return ArrayList of Rotamers that are allowed
	 */
//	public ArrayList<Rotamer> getRotsForType(int i) {
//		
//		ArrayList<Rotamer> rots = new ArrayList<Rotamer>();
//		AARotamerType aaType = AATypesAllowed.get(i);
//		for(Rotamer r: allowedRotamers){
//			if(r.aaType.index == aaType.index){
//				rots.add(r);
//			}
//		}
//		
//		return rots;
//	}

	/**
	 * Given the amino acid type index (0 .. AAtypesAllowed.length)
	 * return the <code>ResidueConformation</code> for that amino acid type that are allowed
	 * 
	 * @param i AA allowed index
	 * @return ArrayList of Residue Conformations that are allowed
	 */
	public ArrayList<ResidueConformation> getRCsForType(int i) {
		
		ArrayList<ResidueConformation> RCs = new ArrayList<ResidueConformation>();
		AARotamerType aaType = AATypesAllowed.get(i);
		for(ResidueConformation r: allowedRCs){
			if(r.rot.aaType.index == aaType.index){
				RCs.add(r);
			}
		}
		
		return RCs;
	}
	
	/**
	 * Given the amino acid type 
	 * return the <code>ResidueConformation</code> for that amino acid type that are allowed
	 * 
	 * @param Type of amino acid to find the allowed RCs for
	 * @return ArrayList of Residue Conformations that are allowed
	 */
	public ArrayList<ResidueConformation> getRCsForType(AARotamerType aaType) {
		
		ArrayList<ResidueConformation> RCs = new ArrayList<ResidueConformation>();
		for(ResidueConformation r: allowedRCs){
			if(r.rot.aaType.index == aaType.index){
				RCs.add(r);
			}
		}
		
		return RCs;
	}

	/**
	 * Given the amino acid type 
	 * return the <code>ResidueConformation</code> for that amino acid type that are allowed
	 * 
	 * @param Type of amino acid to find the allowed RCs for
	 * @return ArrayList of Residue Conformations that are allowed
	 */
	public ArrayList<ResidueConformation> getRCsForType(String AA) {
		
		ArrayList<ResidueConformation> RCs = new ArrayList<ResidueConformation>();
		for(ResidueConformation r: allowedRCs){
			if(r.rot.aaType.name.equals(AA)){
				RCs.add(r);
			}
		}
		
		return RCs;
	}
	
	/**
	 * Given a Residue Conformation index, check whether the residue 
	 * is allowed to mutate to that rotamer's AA type
	 * @param i Global residue conformation library index for rotamer
	 * @return
	 */
	public boolean isResAllowed(ResidueConformationLibrary rcl, int i) {
		ResidueConformation curRC = rcl.getRC(i);
		for(AARotamerType aa:AATypesAllowed){
			if(aa.index == curRC.rot.aaType.index)
				return true;
		}
		return false;
	}
	
	/** PGC 2014: After we load the individual rotamer libraries from the file for 
	 * each unbound/bound entity, we reset the rotamer library for each residue
	 */
	public void setRotamerLibrary(RotamerLibrary aRL){
		this.rl = aRL;
	}

	public boolean isSameAA(String name){
		boolean isLamino = true;
		if(name.length() > 3){
			if(name.startsWith("D"))
				isLamino = false;
			name = name.substring(1);
		}
		if(lAmino != isLamino)
			return false;
		
		if(threeLet().equalsIgnoreCase(name)){
			return true;
		}
		
		return false;
	}

	
	public int getNumRCs() {
		return allowedRCs.size();
	}
	
	/**
	 * Turns this residue into a mutable residue. A mutable residue
	 * knows that it can mutate and stores the rotamers that it's
	 * allowed to mutate to.
	 */
	public void initializeMutableRes(RotamerLibrary rl, boolean addOrigRot){
		if(defaultAA == null)
			defaultAA = name;
		this.rl = rl;
		if(addOrigRot){
			wtRot = rl.addOrigRot(this);
		}
		AATypesAllowed = new ArrayList<AARotamerType>();
		allowedRCs = new ArrayList<ResidueConformation>();
		
		//Save initial coordinates for WT rotamer
		saveWTcoords();
		
		//resRots = new ResRotamers(rl, this, addOrigRot);
	}
	
	/**
	 * Allows the residue to mutate to all of the template rotamers 
	 * for the amino acid given 
	 * 
	 * @param AAname The name of the amino acid that should be added
	 */
	void setAllowable(String AAname){
		AARotamerType aaType = rl.getAAType(AAname);
		
		boolean alreadyAdded = false;
		for(AARotamerType curType: AATypesAllowed){
			if(curType.index == aaType.index)
				alreadyAdded = true;
		}
		
		if(!alreadyAdded)
			AATypesAllowed.add(aaType);
	}
	
	/**
	 * 
	 * Allows the residue to mutate to the given residue conformation
	 * (Checks if the residue conformation has already been added)
	 * 
	 * @param rot The {@link ResidueConformation} that the residue is allowed to mutate to 
	 */
	void setAllowable(ResidueConformation rc){
		
		boolean alreadyAdded = false;
		
		for(ResidueConformation curRC: allowedRCs){
			if(curRC.id == rc.id){
				alreadyAdded = true;
			}
		}
		
		if(!alreadyAdded)
			allowedRCs.add(rc);

	}
	
	void updateWT(boolean[] strandPresent, ParamSet sParams) {
		if(strandPresent[strandNumber]){
			String tempResAllow = (String)sParams.getValue("RESALLOWED"+getResNumberString(),"");
			if(KSParser.numTokens(tempResAllow)==1)
				defaultAA = KSParser.getToken(tempResAllow, 1);
		}
			
	}
	
	/**
	 * Clears all the allowed rotamers that the residue is
	 * allowed to mutate to
	 */
	public void clearAllowable(){
		AATypesAllowed = new ArrayList<AARotamerType>();
		allowedRCs = new ArrayList<ResidueConformation>();
	}

	/**
	 * Be very careful about using this function. Should only be used
	 * if you really only need the AA types allowed. Otherwise rely on
	 * the isRotAllowed functions
	 * @return
	 */
	public ArrayList<AARotamerType> AATypesAllowed() {
		return AATypesAllowed;
	}
	
//	public ArrayList<Rotamer> allowedRotamers(){
//		return allowedRotamers;
//	}
	
	public ArrayList<ResidueConformation> allowedRCs(){
		return allowedRCs;
	}

	public void removeRC(ResidueConformation rc) {
		Iterator<ResidueConformation> iter = allowedRCs.iterator();
		while(iter.hasNext()){ //TODO: Should just implement rc equals function to do this
			ResidueConformation curRC = iter.next();
			if(curRC.id == rc.id)
				iter.remove();
		}
	}

	/**
	 * Copies over the mutable information to the input residue. 
	 * Note the mutable information is not deep copied, so this should
	 * only be used when you are "throwing" away the input residue.
	 * @param res Source Residue to copy info from.
	 */
	public void copyMutInfo(Residue res) {
		AATypesAllowed = res.AATypesAllowed();
		allowedRCs = res.allowedRCs();
		addedOrigRot = res.addedOrigRot;
		wtRC = res.wtRC;
		wtRot = res.wtRot;
		rl = res.rl;
		
		canMutate = res.canMutate; //Is the residue type a type that can mutate (i.e. of type protein).
		isMutable = res.isMutable; //Is this residue a mutable residue
		origMutPos = res.origMutPos;
		defaultAA = res.defaultAA;
		
	}

	/**
	 * PGC 2014: Save the coordinates of all the atoms in the 
	 * input structure to use as a wildtype rotamer;
	 * This prevents problems in recovery because of bond angle and bond length non-ideal values.
	 * KER: Store the coordinates relative to the CA atom, so we can just put them on the current
	 * backbone conformations.
	 * @param 
	 * @return
	 */
	public void saveWTcoords(){
		Atom CA = null;
		
		this.wildtypeAtoms = new Atom[this.atom.length];
		for(int aIx = 0; aIx < this.atom.length; aIx++){
			this.wildtypeAtoms[aIx] = this.atom[aIx].copy();
			if(this.atom[aIx].name.equalsIgnoreCase("CA"))
				CA = this.wildtypeAtoms[aIx];
		}
		
		double CAx = CA.coord[0];
		double CAy = CA.coord[1];
		double CAz = CA.coord[2];
		if(CA == null){
			System.out.println("Wasn't able to find CA atom in residue: "+name+" when saving WT coords.");
			//System.exit(0);
		}else{
			for(Atom a: this.wildtypeAtoms){
				a.coord[0] -= CAx;
				a.coord[1] -= CAy;
				a.coord[2] -= CAz;
			}
		}
		
		
	}
	
	/**
	 * PGC 2014: Set the atom coordinates to the wildtype coordinates. 
	 * KER: Atom coordinates are stored relative to the CA atom, so we need
	 * to find the CA atom and add the coordinates to it.
	 * @param skipBB Don't change the BB atoms (if BB has moved we don't want to move it)
	 * @return
	 */
	public void setResidueToWTcoords(boolean skipBB){

		Atom CA = null;
		for(Atom a: atom){
			if(a.name.equalsIgnoreCase("CA")){
				CA = a;
				break;
			}
		}
		if(CA == null){
			System.out.println("Could not find CA for WT rotamer. Using (0,0,0)");
			double[] zeros = {0.0,0.0,0.0};
			CA = new Atom(zeros);
		}
		
		int nonBBatoms = 0;
		int atomsChanged = 0;
//		System.out.println("Restoring wildtype rotamer coordinates for residue: "+this.fullName);
		try{
		for(int aIx1 = 0; aIx1 < this.atom.length; aIx1++){
			if(!this.atom[aIx1].isBBatom || !skipBB ){
				nonBBatoms++;
				boolean atomMatched = false;
				for(int aIx2 = 0; aIx2 < this.wildtypeAtoms.length; aIx2++){
					if(this.atom[aIx1].isNameEqualTo(this.wildtypeAtoms[aIx2].name)){
	//					System.out.println(this.atom[aIx1]+" changed to :"+this.wildtypeAtoms[aIx2]);
						this.atom[aIx1].coord[0] = CA.coord[0] + this.wildtypeAtoms[aIx2].coord[0];
						this.atom[aIx1].coord[1] = CA.coord[1] + this.wildtypeAtoms[aIx2].coord[1];
						this.atom[aIx1].coord[2] = CA.coord[2] + this.wildtypeAtoms[aIx2].coord[2];
						atomsChanged ++;
						atomMatched = true;				
					}				
				}
				if(!atomMatched){
					System.out.println("Couldn't find a match for "+this.atom[aIx1].name);
				}
			}
		}
		}catch(Exception E){
			E.printStackTrace();
		}
		if (atomsChanged != nonBBatoms || (!skipBB && (atomsChanged != this.wildtypeAtoms.length || atomsChanged != this.atom.length) ) ){
			System.out.println("Something went wrong when trying to copy atom coordinates");
			System.exit(1);
			
		}
	}	

	
}
