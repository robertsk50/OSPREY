/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
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

////////////////////////////////////////////////////////////////////////////////////////////////
// SaveMolecule.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * This Function Rewritten by Ryan Lilien (2001-2004) as the function by Neill White
 *  had several errors and did not conform to the standard pdb format
 *  
 *  Minor modifications by Ivelin Georgiev (2004-2009)
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

import java.io.PrintStream;
import java.io.BufferedInputStream;
import java.util.Hashtable;

/**
 * Saves a molecule to a pdb format file.
 */
class SaveMolecule {

	PrintStream pw;


	// Saves molecule m (atom[][] coordinates, not actualCoordinates[]) to the printstream pw. Params include
	//  comment: (string) a one line comment for header
	//  printSegID: (boolean) print atom seg ids
	//  showConnect: (boolean) add connect terms after atom entries
	SaveMolecule (Molecule inM, PrintStream pw, Hashtable params) throws Exception {
		Molecule m = inM;
		inM = null;
		m.resolveCoordinates(); //copy from actualCoordinates to atom coordinates
		DoSaveMolecule(m,pw,params);
	}

	private void DoSaveMolecule (Molecule m, PrintStream pw, Hashtable params) throws Exception {
		this.pw = pw;
		int atomCounter = 1;
		int residueCounter = 0;
		int prevResNum = -1;
		if (m.numberOfAtoms == 0)
			return;

		if (m.connectivity12Valid == false)
			m.establishConnectivity(true);

		// Pull out parameters
		String comment = (String) params.get("comment");
		boolean printSegID = ((Boolean) params.get("printSegID")).booleanValue();
		boolean showConnect = ((Boolean) params.get("showConnect")).booleanValue();
		double energy = ((Double) params.get("energy")).doubleValue();

		m.resolveCoordinates();
		pw.println("AUTHOR generated by OSPREY");
		pw.println("REMARK   6 " + comment);
		pw.println("REMARK   7 " + energy);

		char[] tmpChr = new char[80];
		String tmpStg;
		Integer tmpInt;
		int tmpLen;
		for(int q=0; q<80; q++)  // clear the array
			tmpChr[q]=' ';

		for(int i=0; i<m.numberOfStrands; i++){
			Strand strand = m.strand[i];
			if (strand.numberOfAtoms == 0)
				continue;

			if (i>0)
				pw.println("TER");

			// copy over the one character strand identifier
			tmpStg = strand.name.toUpperCase();
			tmpLen = tmpStg.length();
			if (tmpLen != 0)
				tmpStg.getChars(0,1,tmpChr,21);
			else
				tmpChr[21] = ' ';

			residueCounter = 0;

			for(int j=0; j<strand.numberOfResidues; j++){
				Residue residue = strand.residue[j];
				if (residue.numberOfAtoms == 0)
					continue;
				residueCounter++;

				tmpStg = residue.name.toUpperCase();
				tmpLen = 3-tmpStg.length();
				if (tmpLen <= 0)
					if(tmpLen < 0 && tmpStg.startsWith("D"))
						tmpStg.getChars(1, 4, tmpChr, 17);
					else
						tmpStg.getChars(0,3,tmpChr,17);
				else {
					// from 17-19 (zero based)
					tmpStg.getChars(0,tmpStg.length(),tmpChr,17+tmpLen);
					for(int q=tmpStg.length();q<3;q++)  // add whitespace on left
						tmpChr[19-q]=' ';
				}
				tmpStg = residue.getResNumber();
				tmpLen = 4-tmpStg.length();
				String kabatL = " ";
				// PGC 2013: support for Kabat formatting: if tmpStg ends in a letter, then this letter will go into field 26.
				if (tmpStg.matches(".*[A-Z]$")){
					kabatL = tmpStg.substring(tmpStg.length()-1, tmpStg.length() );
					tmpStg = tmpStg.substring(0, tmpStg.length()-1);
					tmpLen = 4-tmpStg.length();
				}

				if (tmpLen <= 0)
					tmpStg.getChars(0,4,tmpChr,22);
				else {
					// from 22-25 (zero based)
					tmpStg.getChars(0,tmpStg.length(),tmpChr,22+tmpLen);
					for(int q=tmpStg.length();q<4;q++)  // add whitespace on left
						tmpChr[25-q]=' ';
					kabatL.getChars(0, 1, tmpChr, 26);
				}

				// If the residue number decreases in the middle
				//  of a strand change the strand name to "R"
				//				if(tmpInt.intValue() < prevResNum){
				//					tmpChr[21] = 'R';
				//				}
				//prevResNum = tmpInt.intValue();

				for(int k=0; k<residue.numberOfAtoms; k++){
					Atom atom = residue.atom[k];
					tmpStg = "ATOM  ";
					tmpStg.getChars(0,6,tmpChr,0);  // Copy over the ATOM term
					tmpStg = "1.00";
					tmpStg.getChars(0,4,tmpChr,56); // Make the occupancy 1.00
					tmpStg = "0.00";
					tmpStg.getChars(0,4,tmpChr,62); // and the temp factor 0.0

					tmpInt = new Integer(atomCounter);
					tmpStg = tmpInt.toString();
					tmpLen = 5-tmpStg.length();
					if (tmpLen <= 0)
						tmpStg.getChars(0,5,tmpChr,6);
					else {
						// from 6-10 (zero based)
						tmpStg.getChars(0,tmpStg.length(),tmpChr,6+tmpLen);
						for(int q=tmpStg.length();q<5;q++)  // add whitespace on left
							tmpChr[10-q]=' ';
					}
					atom.modelAtomNumber = atomCounter++;

					// writing the atom name is a little fuzzy, although the
					//  atom name is allocated columns 12-15(zero based), rasmol
					//  likes and other people essentially start with column 13
					//  leaving column 12 blank. So we only allow 3 characters for
					//  the atom name and it should be left justified
					//  unless the first character is a number in which case we
					//  start with column 12
					// there are also exceptions when the atom has a two letter
					//  element code

					tmpStg = getAtomField(atom);
					tmpStg.getChars(0,4,tmpChr,12);

					// Write the x coordinate
					tmpStg = coordinate(atom.coord[0]);
					tmpLen = 8-tmpStg.length();
					if (tmpLen < 0) {
						System.out.println("ERROR: coordinate exceeds pdb format size");
						tmpStg.getChars(0,8,tmpChr,30);
					}
					else {
						// from 30-37 (zero based)
						tmpStg.getChars(0,tmpStg.length(),tmpChr,30+tmpLen);
						for(int q=tmpStg.length();q<8;q++)  // add whitespace on left
							tmpChr[37-q]=' ';
					}

					// Write the y coordinate
					tmpStg = coordinate(atom.coord[1]);
					tmpLen = 8-tmpStg.length();
					if (tmpLen < 0) {
						System.out.println("ERROR: coordinate exceeds pdb format size");
						tmpStg.getChars(0,8,tmpChr,38);
					}
					else {
						// from 38-45 (zero based)
						tmpStg.getChars(0,tmpStg.length(),tmpChr,38+tmpLen);
						for(int q=tmpStg.length();q<8;q++)  // add whitespace on left
							tmpChr[45-q]=' ';
					}

					// Write the z coordinate
					tmpStg = coordinate(atom.coord[2]);
					tmpLen = 8-tmpStg.length();
					if (tmpLen < 0) {
						System.out.println("ERROR: coordinate exceeds pdb format size");
						tmpStg.getChars(0,8,tmpChr,46);
					}
					else {
						// from 46-53 (zero based)
						tmpStg.getChars(0,tmpStg.length(),tmpChr,46+tmpLen);
						for(int q=tmpStg.length();q<8;q++)  // add whitespace on left
							tmpChr[53-q]=' ';
					}
					if (printSegID) {
						tmpStg = atom.segID;
						tmpLen = 4-tmpStg.length();
						tmpStg.getChars(0,1,tmpChr,21); // Make the chain identifier the first segID char
						if (tmpLen < 0) {
							System.out.println("ERROR: segment id exceeds pdb format size");
							tmpStg.getChars(0,4,tmpChr,72);
						}
						else {
							// from 72-75 (zero based)
							tmpStg.getChars(0,tmpStg.length(),tmpChr,72);
							for(int q=(72+tmpStg.length());q<76;q++) {
								tmpChr[q]=' ';
							}
						}
					}

					pw.println(new String(tmpChr));
				}
			}
		}

		if (showConnect) {
			for (int i=0; i<m.numberOfAtoms; i++){
				int numberOfConnections = m.connected[i][0];
				if (numberOfConnections > 0)
					pw.print("CONECT");
				else 
					continue;
				printAtom(m.atom[ i ].modelAtomNumber);	
				for(int j=1; j<=numberOfConnections; j++){
					int bondedToAtomNumber = m.connected[i][j];
					Atom bondedToAtom = m.atom[bondedToAtomNumber];
					printAtom(bondedToAtom.modelAtomNumber);
				}
				pw.println();
			}
		}
		pw.println("END");
	}

	private void printAtom(int atomNumber){
		if (atomNumber < 10)
			pw.print("    " + atomNumber );
		else if (atomNumber < 100)
			pw.print("   " + atomNumber );
		else if (atomNumber < 1000)
			pw.print("  " + atomNumber);
		else if (atomNumber < 10000)
			pw.print(" " + atomNumber);
		else if (atomNumber < 100000)
			pw.print(atomNumber);
	}

	private String coordinate(double coord){
		if((coord<0.001) && (coord>-0.001))
			coord = 0.0f;
		String coordString = String.valueOf(coord);
		String intPart = " ";
		String doublePart = " ";
		String returnString = null;
		int radix = coordString.indexOf(".");
		if (radix != -1){
			intPart = coordString.substring(0, radix);
			if (radix != coordString.length())
				doublePart = coordString.substring(radix+1);
		}
		else 
			return(coordString);
		if (intPart.length() == 1)
			returnString = "   " + intPart;
		else if (intPart.length() == 2)
			returnString = "  " + intPart;
		else if (intPart.length() == 3)
			returnString = " " + intPart;
		else if (intPart.length() == 4)
			returnString = intPart;
		else if (intPart.length() > 4)
			returnString = intPart.substring(0, 3);
		returnString += ".";
		if (doublePart.length() == 1)
			returnString += doublePart + "00";
		else if (doublePart.length() == 2)
			returnString += doublePart + "0";
		else if (doublePart.length() == 3)
			returnString += doublePart;
		else if (doublePart.length() > 3)
			returnString += doublePart.substring(0, 3);
		return(returnString);
	}

	//Returns the atom name field
	private String getAtomField(Atom at){
		if(at.elementType.length()==1) {
			if (at.name.length() == 1)
				return(" " + at.name + "   ");
			else if (at.name.length() == 2){
				if ( (at.name.charAt(0)>='0') && (at.name.charAt(0)<='9') )
					return(at.name + "   ");
				else
					return(" " + at.name + "  ");
			}
			else if (at.name.length() == 3){
				if ( (at.name.charAt(0)>='0') && (at.name.charAt(0)<='9') )
					return(at.name + "  ");
				else
					return(" " + at.name + " ");
			}
			else // (at.name.length() > 3)
				return(at.name.substring(0, 4) + " ");
		}
		else
		{
			if (at.name.length() == 2)
				return(at.name + "   ");
			else if (at.name.length() == 3)
				return(at.name + "  ");
			else // (at.name.length() > 3)
				return(at.name.substring(0, 4) + " ");
		}
	}
}
