

public class MutUtils {

	static final int N = 0;
	static final int CA = 1;
	static final int C = 2;
	static final int O = 3;
	static final int H = 4;
	static final int CB = 5;
	static final int HA = 6;
	
	private static boolean debug = false;
	private static int ctr=0;
	
	// This function converts residue (resNum) to the conformation specified by rotamer rotNum
	public static boolean applyRotamer(Molecule m, Rotamer rot, Residue localResidue) {

		
		
		//If the rotamer is already set we don't need to set it again
		if(localResidue.curRC != null && localResidue.curRC.rot != null && localResidue.curRC.rot.rlIndex == rot.rlIndex){
			return true;
		}
		
		if(EnvironmentVars.STORE_FULL_WT_ROT && rot.isWTrot){
			localResidue.setResidueToWTcoords(true);
			
			for(int a=0;a<localResidue.numberOfAtoms;a++){
                if( !localResidue.atom[a].isBBatom )
                    m.updateCoordinates(localResidue.atom[a].moleculeAtomNumber);
            }
			
			//BB conformational changes might make this orientation inappropriate though
	        //So idealize the sidechain to fix it
	        //This won't return false because it's not a proline
	        m.idealizeResSidechain(localResidue);
	
	        //Finally make sure the (generalized) chi1 is correct
	        Perturbation.setGenChi1(m, localResidue.moleculeResidueNumber, localResidue.WTGenChi1 );

			return true;
		}
		
		AARotamerType aaType = rot.aaType;
		
		if(rot.values == null)
			return true;
		
		int aaDihedrals = aaType.numDihedrals();
		
		int atomList[] = new int[localResidue.numberOfAtoms];
		int alSize = 0;

		for(int i=0;i<aaDihedrals;i++) {
	
			Atom at[] = new Atom[4];
			at[0] = at[1] = at[2] = at[3] = null;
			int atNum[] = new int[4];

			// Find atoms involved in the dihedral
			for(int q=0;q<localResidue.numberOfAtoms;q++) {
				for(int w=0;w<4;w++) {
					if (localResidue.atom[q].name.equalsIgnoreCase(aaType.dihedralAtomNames[i][w])) {
						at[w] = localResidue.atom[q];
						atNum[w] = q;
					}
				}		
			}

			// If we didn't find all the atoms return false
			if ((at[0] == null) || (at[1] == null) || (at[2] == null) ||
				(at[3] == null))
				return false;

//			//If the rotamer is already set we don't need to set it again
//			if(localResidue.curRot != null && localResidue.curRot.rlIndex == rot.rlIndex){
//				
//				//If the rotamer values are still valid then we can skip to next dihedral
//				//(Need to check because they might have been minimized)
//				double angleDiff = Math.abs(at[3].torsion(at[0], at[1], at[2])-rot.values[i]);
//				if(angleDiff > 180)
//					angleDiff = 360-angleDiff;
//							
//				if(angleDiff < 0.0001){
//					continue;
//				}
//			}
		
			
			for(int q=0;q<localResidue.numberOfAtoms;q++) {
				atomList[q]=1;
			}
			atomList[atNum[1]]=0;
			atomList[atNum[2]]=0;
				
			// Now find all atoms in the rotamer that are more distal than at[2]
			getAtomsMoreDistal(m,localResidue.moleculeResidueNumber,at[2],atomList);

			// Copy atoms over in atomList, ie. if current atomList[q]==2 then
			//  the new atomList[] should include q (ie. atom q counts)
			// Skip the 4th atom of the torsion as it's automatically rotated
			//  by m.setTorsion
			alSize=0;
			for(int q=0;q<localResidue.numberOfAtoms;q++){
				if ((atomList[q]==2) && (q != atNum[3]))
					atomList[alSize++]=localResidue.atom[q].moleculeAtomNumber;
			}

			// Perform the rotation
			// at[0], at[1], and at[2] don't move, at[3] and the atoms
			//  in the atomList rotate
			try{
			m.setTorsion(at[0].moleculeAtomNumber,at[1].moleculeAtomNumber,
				at[2].moleculeAtomNumber,at[3].moleculeAtomNumber,
				rot.values[i],atomList,alSize);
			}catch(Exception E){
				System.out.println("DELETE ME");
			}
		}
		
		return true;
	}

	
	
	// This function changes the residue type (it performs a mutation)
	// The residue resNum of strand standNumber of molecule m is changed
	//  to a residue with name newResType
	// If addHydrogens is true then hydrogens are added to the new
	//  residue
	public static void changeResidueType(Molecule m, int resID, String newResType, boolean addHydrogens) {

		// call with connectResidue = true and useOldBBatoms = false
		changeResidueType(m,resID,newResType,addHydrogens,true,true);
	}
	
	public static void changeResidueType(Molecule m, int resID, String newResType, boolean addHydrogens, boolean connectResidue) {

		// call with useOldBBatoms = false
		changeResidueType(m,resID,newResType,addHydrogens,connectResidue,true);
	}

	// This function changes the residue type (it performs a mutation)
	// The residue resNum of strand strandNumber of molecule m is changed
	//  to a residue with name newResType
	// If addHydrogens is true then hydrogens are added to the new
	//  residue, otherwise hydrogens are stripped
	// If connectResidue is true then AA's are connected if their residue
	//  numbers are sequential.
	// If useOldBBatoms is true, then the *exact* coordinates for the old backbone
	//	CA, C, and CB atoms are copied for the new residue
	//
	// This function makes use of the CB carbon, unfortunately if we are changing
	//  to or from Gly we don't have a CB, in that case we match with HA3
	// So the overall scheme is:
	//  1) Translate N's to overlap
	//  2) Rotate CA's to overlap
	//  3) Rotate around N-CA axis to get CB's to overlap
	//
	// In step 3), we could align the backbone C's, but due to differences between 
	//		the PPR amino acid template and actual residue conformations, the CB may 
	//		differ significantly in such a case. To preserve the CA-CB orientation, 
	//		the CB's are aligned instead.
	//
	// This function is rather complicated due to the bookkeeping
	//  required to maintain our messy molecule datastructure.
	public static void changeResidueType(Molecule m, int resID, String newResType, boolean addHydrogens, boolean connectResidue, boolean useOldBBatoms) {

		Residue localResidue = m.residue[resID];
		localResidue.curRC = null;
		
		RotMatrix rm = new RotMatrix();
		
		int resNum = localResidue.strandResidueNumber;
		int strandNumber = localResidue.strandNumber;
		
		boolean newDaa    = false;
		boolean oldDaa    = false;
		boolean newResGly = false;
		boolean oldResGly = false;
		
		boolean newResPro = false;
		boolean oldResPro = false;
		boolean proMutation = false;
		
		if (newResType.length() == 4 && newResType.startsWith("D")) {
			newDaa = true;
		}
		if (!localResidue.lAmino){
			oldDaa = true;
		}
		
		// If the old or new residues are glycine then we must do special things as we treat the H as CB
		if(newResType.equalsIgnoreCase("gly") ||
				newResType.equalsIgnoreCase("dgly")){
				newResGly = true;
				if(oldDaa == true)
					newDaa = true;
		}
		if(localResidue.name.equalsIgnoreCase("gly") ||
				localResidue.name.equalsIgnoreCase("dgly"))
				oldResGly = true;

		//Proline is also a special case:
		if(localResidue.name.equalsIgnoreCase("pro") || localResidue.name.equalsIgnoreCase("dpro"))
			oldResPro = true;
		else if(newResType.equalsIgnoreCase("pro") || newResType.equalsIgnoreCase("dpro"))
			newResPro = true;
		
		//glyMutation = newResGly || oldResGly;	
		proMutation = newResPro || oldResPro;
		
		// We assume a standard backbone ordering of: N,CA,C,O

		int savedMoleculeResidueNumber = localResidue.atom[0].moleculeResidueNumber;
		int savedStrandResidueNumber = localResidue.atom[0].strandResidueNumber;
		int savedStrandNumber = localResidue.atom[0].strandNumber;
		String savedSegID = localResidue.atom[0].segID;

		// Get the new residue from the templates
 		Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
		Residue r = ppr.getResidue(newResType);
		
				
		// Residue r = ppr.getResidue("Lala");
		Molecule m2 = new Molecule();
		m2.addResidue(0,r);
		m2.establishConnectivity(false);
		
		// First get a handle on the backbone N, CA, C, O, H, and CB atoms
		Atom at[] = null;
		at = getBBatoms(r); //for the new residue
		Atom NNew = at[0]; Atom CANew = at[1]; Atom CNew = at[2]; Atom ONew = at[3]; Atom HNew = at[4]; Atom CBNew = at[5];Atom HANew = at[6];
		
		at = getBBatoms(localResidue); //for the old residue
		Atom NOld = at[0]; Atom CAOld = at[1]; Atom COld = at[2]; Atom OOld = at[3]; Atom HOld = at[4]; Atom CBOld = at[5];Atom HAOld = at[6];
		
		//KER: DAA switching residues
		//There are several scenarios which is complicated by the fact that we use L-aa templates
		//Now we just use the right template and only have to switch residues when swapping between the types
		//Case1: Mutate L to L -- Everything valid no switching
		//Case2: Mutate L to D -- Template not ok so we just switch oldCB and oldHA
		//Case3: Mutate D to L -- Template is  ok so we just switch oldCB and oldHA 
		//Case4: Mutate D to D -- Template not ok so we switch newCB and newHA
		
		if(oldDaa!=newDaa){ //if we are switching a L-aa with a D-aa
			System.out.println("Switching from D-L or L-D has not been fully tested...Exiting");
			System.exit(0);
			// if we are replacing a L-aa with a D-aa we need to align the new CB with the old HA
			// CBOld = localResidue.atom[HA];				
		}		

		//Reset CB atom for mutation from gly or pro
		if (oldResGly || oldResPro){ 
			setOrigAtom(CBOld, OOld, COld, CAOld, localResidue.CBplace);
		}
		
		double newNHLength = 0.0;
		if(proMutation) //mutation to or from Pro
			newNHLength = rm.norm( rm.subtract( HNew.coord, NNew.coord ) );//New amide NH or N-CD bond length
		
		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}
		
		// START ALIGNMENT
		// Translate N's to overlap
		double Ntrans[] = new double[3];
		Ntrans[0] = NNew.coord[0] - NOld.coord[0];
		Ntrans[1] = NNew.coord[1] - NOld.coord[1];
		Ntrans[2] = NNew.coord[2] - NOld.coord[2];
		for(int q=0;q<r.numberOfAtoms;q++) {
			r.atom[q].coord[0] -= Ntrans[0];
			r.atom[q].coord[1] -= Ntrans[1];
			r.atom[q].coord[2] -= Ntrans[2];
		}

		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}
		
		int atomList[] = null;
		double thetaDeg[] = new double[1];
		double rotAxis[] = new double[3];
		int numAtoms = r.numberOfAtoms;
		atomList = new int[numAtoms];

		// Rotate CAs to overlap
		if (CAOld.distance(CANew) > 0.0001) {
			getRotationInfoA(CAOld, NOld, CANew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(NNew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}
		
		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}
		
		// Rotate CBs to overlap - now the residue backbones should be aligned		
		if ( (!newResGly) && (CBOld.distance(CBNew) > 0.0001) ) { //not a to- or from- Gly mutation
			getRotationInfoB(CBOld, CAOld, NOld, CBNew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(CANew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}
		
		if (newResGly) //mutation to Gly
			alignCBNewGly(CBOld,CBNew,CAOld,CANew,NOld,HAOld, r);
	
		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}
		
		// Remove hydrogens if we don't want them
		if (!addHydrogens) {
			for(int q=0;q<r.numberOfAtoms;q++) {
				if (r.atom[q].elementType.equalsIgnoreCase("H"))
					m2.deleteAtom(r.atom[q].moleculeAtomNumber);
			}
		}
		else {
			// ELSE KEEP ALL HYDROGENS
		}
		
		//Make the positions of the new backbone atoms coincide *exactly* with the old ones;
		// 	NNew already has the same coordinates as NOld; 
		// 	the CAs, CBs, and CBs are already aligned, but may differ due to differences between the template PPR and the residue backbone in the PDB
		if (useOldBBatoms){
			setAtomCoord(CANew,CAOld);
			//Don't count CB as part of the backbone if we're idealizing side chains.
			//if ( !newResGly && !Perturbation.idealizeSC ){ //not a from/to Gly mutation (Gly does not have a CB)
			//	moveCB(m2.residue[0],CBNew,CBOld);
			//}
		}
		setAtomCoord(CNew,COld);
		setAtomCoord(ONew,OOld);
		
		//Set the H's for the residue
		if(addHydrogens){
			if( newResPro ){ //if mutating to proline rescale the N-CD vector
				double NVec[] = RotMatrix.subtract(HOld.coord, NOld.coord);//N- to H or CD bond vector
				NVec = RotMatrix.scale(NVec, newNHLength/rm.norm(NVec) );
				HNew.coord = RotMatrix.add(NVec, NNew.coord);
			}
			else if( oldResPro  ){ //if mutating from proline, set the N-H vector back to the original position
				setOrigAtom(HNew, COld, CAOld, NOld, localResidue.Hplace);
			}
			else { //otherwise, just set the H position to the previous H position
				if(HOld!=null)
					setAtomCoord(HNew,HOld);
			}
			if(HAOld!=null && HANew!=null) //Set the HA atom to the old value
				setAtomCoord(HANew,HAOld);
		}
		
		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}
		//////////////////////////////////////////////////////////////////////////////////
		if (localResidue.nterm){ //we are changing the nterm residue
			//we must add the H1, H2, H3 atoms (and delete the HN atom) to the new residue, 
			//	since the PPR templates only handle polypeptide residues			
			
			Residue r1 = getNtermRes(r,localResidue);
			
			m2 = new Molecule();
			m2.addResidue(0, r1);
			m2.establishConnectivity(false);
		}
		if (localResidue.cterm){ //we are changing the cterm residue
			//we must add the OXT or H1, H2, H3 atoms (and delete the HN atom) to the new residue, 
			//	since the PPR templates only handle polypeptide residues			
			
			Residue r1 = getCtermRes(r,localResidue);
			
			m2 = new Molecule();
			m2.addResidue(0, r1);
			m2.establishConnectivity(false);
		}
		////////////////////////////////////////////////////////////////////////////////
		
		if(debug){
			m2.saveMolecule("mutPDB_"+ ctr++ + ".pdb", 0.0,true);
		}

		// Copy the new residue information into the old residue
		int changeInAtoms = r.numberOfAtoms - localResidue.numberOfAtoms;
		int baseChangedAtom = localResidue.atom[0].moleculeAtomNumber;
		// first atomnum in next residue
		int nextResidueBase = -1;
		// the first atomnum in this residue
		int thisResidueBase = -1;
		// atomnumber of the C in the last residue
		int lastResidueC = -1;
		// Determine if the residue before and the one after are have sequential
		//  residue numbers. If they do and if we're interested in connecting
		//  sequential residues then gather information so we can make the
		//  appropriate bonds.
		boolean connectedResidue = false;
		boolean connectNextResidue = false;
		boolean connectLastResidue = false;
		thisResidueBase = m.strand[strandNumber].residue[resNum].atom[0].moleculeAtomNumber;
		nextResidueBase = thisResidueBase + localResidue.numberOfAtoms;

		if ((resNum+1) < m.strand[strandNumber].numberOfResidues) {
			if (connectResidue)
				connectNextResidue = m.residuesAreBBbonded(strandNumber, resNum, strandNumber, resNum+1);
		}
		if (resNum > 0) {
			if (connectResidue)
				connectLastResidue = m.residuesAreBBbonded(strandNumber, resNum-1, strandNumber, resNum);	
			if (connectLastResidue) {
				for(int q=0;q<m.strand[strandNumber].residue[resNum-1].numberOfAtoms;q++) {
					if (m.strand[strandNumber].residue[resNum-1].atom[q].name.equalsIgnoreCase("C"))
						lastResidueC = m.strand[strandNumber].residue[resNum-1].atom[q].moleculeAtomNumber;
				}
			}
			if (lastResidueC == -1)
				connectLastResidue = false;
		}
		

		localResidue.name = r.name; //Will include "D" if a D amino acid
		if (localResidue.fullName.length()  > 4)
			localResidue.fullName = r.name.substring(0,3) + " " + localResidue.fullName.substring(4);
		else
			localResidue.fullName = r.name.substring(0,3);
		localResidue.numberOfAtoms = r.numberOfAtoms;
		localResidue.atom = r.atom;
		localResidue.lAmino = r.lAmino;
		
		for(int j=0;j<localResidue.numberOfAtoms;j++){
			localResidue.atom[j].moleculeResidueNumber = savedMoleculeResidueNumber;
			localResidue.atom[j].strandResidueNumber = savedStrandResidueNumber;
			localResidue.atom[j].strandNumber = savedStrandNumber;
			localResidue.atom[j].segID = new String(savedSegID);
		}
		
		// Update atoms in residues of this strand after this residue
		//  as well as the bookkeeping in the molecule itself
		int curAtom = 0;
		int linkfrom=-1, linkto=-1;
		m.numberOfAtoms += changeInAtoms;
		m.numberOfAtomsx3 = m.numberOfAtoms * 3;
		m.atom = new Atom[m.numberOfAtoms];
		m.actualCoordinates = new double[m.numberOfAtomsx3];
		for(int j=0;j<m.numberOfStrands;j++) {
			for(int q=0;q<m.strand[j].numberOfResidues;q++) {
				for(int w=0;w<m.strand[j].residue[q].numberOfAtoms;w++) {
					m.strand[j].residue[q].atom[w].moleculeAtomNumber = curAtom;
					try{
						m.atom[curAtom++] = m.strand[j].residue[q].atom[w];
					}catch(Exception E){
						System.out.println("DELETE ME.");
						E.printStackTrace();
					}
					int tmpIntAry[];
					tmpIntAry = m.strand[j].residue[q].atom[w].bond;
					if ((q==resNum) && (j==strandNumber)) {
						if (tmpIntAry != null) {
							for(int i=0;i<m.strand[j].residue[q].atom[w].bond.length;i++) {
								tmpIntAry[i] += baseChangedAtom;
							}
						}
						if (m.strand[j].residue[q].atom[w].name.equalsIgnoreCase("C") && (nextResidueBase != -1) && connectNextResidue) {
							// mark the atoms to bond so we can join when we're done
							//  we can't join them now as the second atom doesn't
							//  yet exist in the atom[] array
							connectedResidue = true;
							linkfrom = nextResidueBase+changeInAtoms;
							linkto = m.strand[j].residue[q].atom[w].moleculeAtomNumber;
						}
						if (m.strand[j].residue[q].atom[w].name.equalsIgnoreCase("N") && (lastResidueC != -1) && connectLastResidue) {
							m.addBondBetween(lastResidueC,m.strand[j].residue[q].atom[w].moleculeAtomNumber);
						}
					}
					else {
						if ((tmpIntAry != null) && (nextResidueBase != -1)) {
							for(int i=0;i<m.strand[j].residue[q].atom[w].bond.length;i++) {
								// Because the position of this residue's backbone C atom can move
								//  we have to be careful about simply updating the bond numbers.
								// The backbone N of numRes+1 if it bonds back to this residue needs
								//  special attention
								// Basically if there is a bond to the old residue we'll delete that
								//  bond. The only types of bonds into the new residue should be
								//  the previous and next peptide bond which we handle ourselves.
								if (tmpIntAry[i] >= nextResidueBase)
									tmpIntAry[i] += changeInAtoms;
								else if (tmpIntAry[i] >= thisResidueBase){
									m.strand[j].residue[q].atom[w].deleteBond(i);
									tmpIntAry = m.strand[j].residue[q].atom[w].bond;
									i--;
								}
							}
						}
					}
					
				}
			}
		}
		if (connectedResidue)
			m.addBondBetween(linkfrom,linkto);
		
		
		
		// Establish all connectivity including non-bonded interactions
		//KER: Remember that if you try to get rid of this establishConnectivity, you have to check for connectivity in Amber96ext.calculateTypesWithTemplates
		m.connectivityValid = false;
		//m.establishConnectivity(false);
		m.updateNumAtoms();
		
		localResidue.ffAssigned = false;

		// Copy atom coordinates back into actualCoordinates
		for(int q=0;q<m.numberOfAtoms;q++)
			m.updateCoordinates(q);
		
		//Idealize the sidechain if needed.
		if(newResPro || Perturbation.idealizeSC){//Idealize the sidechain since proline has an unusual geometry
			//This is especially important if we are mutating to proline (some bonds are probably way off length then:
			//the idealization reconstructs the ideal ring given the backbone, CB, and CD coordinates)

			m.idealizeResSidechain(localResidue);//this will also enforce the specified pucker (thus matching the original if we mutated away from and then back to proline)
			int firstAtom = localResidue.atom[0].moleculeAtomNumber;
			for(int a=0; a<localResidue.numberOfAtoms; a++)
				m.resolveCoordinates(firstAtom+a);//Copy the idealized coordinates back into the Atom.coord arrays

			if(!newResPro)//No longer proline, so ring closure is no longer an issue
				localResidue.validConf = true;
			
			if(debug){
				m.saveMolecule("idealize_"+ctr+".pdb", 0.0);
			}
		}
		
	
	}
	
	
	/**
	 * Move the CB atom and all atoms that are downstream of the CB atom
	 * 
	 * @param cBNew
	 * @param cBOld
	 */
	private static void moveCB(Residue r, Atom CBNew, Atom CBOld) {
		
		double[] delta = new double[3];
		
		delta[0] = CBOld.coord[0] - CBNew.coord[0];
		delta[1] = CBOld.coord[1] - CBNew.coord[1];
		delta[2] = CBOld.coord[2] - CBNew.coord[2];
		
		for(Atom a:r.atom){
			if(!a.isBBatom){
				translateAtomCoord(a, delta);
			}
		}
		
	}




	//Returns a handle on the backbone N, CA, C, O, H, and CB atoms for residue res
	public static Atom [] getBBatoms(Residue res){
		
		Atom at[] = new Atom[7];
		
		for(int q=0;q<res.numberOfAtoms;q++) {
			if (res.atom[q].name.equalsIgnoreCase("N"))
				at[N] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("CA"))
				at[CA] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("C"))
				at[C] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("O")||res.atom[q].name.equalsIgnoreCase("O1"))
				at[O] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("H") || ( res.atom[q].name.equalsIgnoreCase("CD") && res.name.equalsIgnoreCase("PRO") ))
				at[H] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("CB")||res.atom[q].name.equalsIgnoreCase("HA3")||res.atom[q].name.equalsIgnoreCase("3HA")) //HA3 for gly
				at[CB] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("HA")||res.atom[q].name.equalsIgnoreCase("HA2")||res.atom[q].name.equalsIgnoreCase("2HA")) //HA2 for gly
				at[HA] = res.atom[q];
		}
		
		if(at[CB] == null){
			System.out.println("HA3 on "+res.fullName+" couldn't be found. Please check labelling.");
			System.out.println("Going to assume we want HA2");
			for(int q=0;q<res.numberOfAtoms;q++) {
				if ( (res.atom[q].name.equalsIgnoreCase("HA2")) || (res.atom[q].name.equalsIgnoreCase("2HA")) )
					at[CB] = res.atom[q];
				if ( (res.atom[q].name.equalsIgnoreCase("HA")) || (res.atom[q].name.equalsIgnoreCase("HA1")) )
					at[HA] = res.atom[q];
			}
		}
		
		return at;
	}
	
	// This function is called when NNew is aligned with NOld. We want
	//  the the axis and angle that will rotate CANew into CAOld
	// In the function call at1 = CAOld, at2 = NOld = NNew,
	//  and at3 = CANew
	// The axis is returned in axis and the angle is in thetaAngle
	// This is one of the more obfuscated functions
	private static void getRotationInfoA(Atom at1, Atom at2, Atom at3,
		double thetaAngle[], double axis[]){
	
		// Theta is the angle between at1-at2-at3
		double	R2D = 57.29577951308232090712;
		double x12,y12,z12,x32,y32,z32,l12,l32,dp,dx,dy,dz;
		double fx,fy,fz,fp,thetaDeg2;
		double thetaDeg;

		x12 = at1.coord[0] - at2.coord[0];
		y12 = at1.coord[1] - at2.coord[1];
		z12 = at1.coord[2] - at2.coord[2];
		x32 = at3.coord[0] - at2.coord[0];
		y32 = at3.coord[1] - at2.coord[1];
		z32 = at3.coord[2] - at2.coord[2];		
		l12 = Math.sqrt(x12*x12 + y12*y12 + z12*z12);
		l32 = Math.sqrt(x32*x32 + y32*y32 + z32*z32);
		if ((l12 == 0.0) || (l32 == 0.0)) {
			thetaDeg = 0;
		}
		else {
			dp = (x12*x32 + y12*y32 + z12*z32) / (l12*l32);
			if (dp < -1.0)
				dp = -1.0;
			else if (dp > 1.0)
				dp = 1.0;
			thetaDeg = R2D * Math.acos(dp);
		}
		// To exactly pin down the angle we need to take the
		//  cross product of 12 and 32
		fx = y12*z32 - z12*y32;
		fy = z12*x32 - x12*z32;
		fz = x12*y32 - y12*x32;
		fp = (Math.sqrt(fx*fx + fy*fy + fz*fz))/(l12*l32);
		if (fp < -1.0)
		  fp = -1.0;
		else if (fp > 1.0)
		  fp = 1.0;
		thetaDeg2 = R2D * Math.asin(fp);
		// If the sign of the angle from the asin and acos do not agree
		//  then make the sign of the cos the same as that of the sin
		if (((thetaDeg2 > 0) && (thetaDeg < 0)) || 
				((thetaDeg2 < 0) && (thetaDeg > 0)))
			thetaDeg = -thetaDeg;

		// The rotation axis is the cross product of 12 cross 32
		// this is fx,fy,fz
			
		thetaAngle[0]=thetaDeg;
		axis[0]=fx;
		axis[1]=fy;
		axis[2]=fz;
	}

	// This function is called when CANew is already aligned with CAOld and NNew is
	//  aligned with NOld. We want the axis and angle that will rotate
	//  CBOld into CBNew
	// In the function call at1 = CBOld, at2 = CAOld = CANew,
	//  at3 = NOld = NNew, and at4 = CBNew
	// The axis is returned in axis and the angle is in thetaAngle
	// This is one of the more obfuscated functions
	private static void getRotationInfoB(Atom at1, Atom at2, Atom at3, Atom at4,
		double thetaAngle[], double axis[]){

		double R2D = 57.29577951308232090712;
		double x12,y12,z12,x32,y32,z32,l12,l32,dp,dx,dy,dz;
		double x12b,y12b,z12b,x32b,y32b,z32b,ex,ey,ez;
		double fx,fy,fz,fp,thetaDeg2;
		double thetaDeg;

		x12 = at1.coord[0] - at2.coord[0];
		y12 = at1.coord[1] - at2.coord[1];
		z12 = at1.coord[2] - at2.coord[2];
		x32 = at3.coord[0] - at2.coord[0];
		y32 = at3.coord[1] - at2.coord[1];
		z32 = at3.coord[2] - at2.coord[2];		
		// d is cross product of vectors 12 and 32
		dx = y12*z32 - z12*y32;
		dy = z12*x32 - x12*z32;
		dz = x12*y32 - y12*x32;
		x12b = at4.coord[0] - at2.coord[0];
		y12b = at4.coord[1] - at2.coord[1];
		z12b = at4.coord[2] - at2.coord[2];
		x32b = at3.coord[0] - at2.coord[0];
		y32b = at3.coord[1] - at2.coord[1];
		z32b = at3.coord[2] - at2.coord[2];		
		// e is cross product of vectors 12b and 32b
		ex = y12b*z32b - z12b*y32b;
		ey = z12b*x32b - x12b*z32b;
		ez = x12b*y32b - y12b*x32b;

		// Now determine dot product between vectors d and e
		//  to get cos of angle
		l12 = Math.sqrt(dx*dx + dy*dy + dz*dz);
		l32 = Math.sqrt(ex*ex + ey*ey + ez*ez);
		if ((l12 == 0.0) || (l32 == 0.0)) {
			thetaDeg = 0;
		}
		else {
			dp = (dx*ex + dy*ey + dz*ez) / (l12*l32);
			if (dp < -1.0)
				dp = -1.0;
			else if (dp > 1.0)
				dp = 1.0;
			thetaDeg = R2D * Math.acos(dp);
		}
		// To exactly pin down the angle we need to take the
		//  cross product of d and e
		fx = dy*ez - dz*ey;
		fy = dz*ex - dx*ez;
		fz = dx*ey - dy*ex;
		fp = (Math.sqrt(fx*fx + fy*fy + fz*fz))/(l12*l32);

		if (fp < -1.0)
		  fp = -1.0;
		else if (fp > 1.0)
		  fp = 1.0;
		thetaDeg2 = R2D * Math.asin(fp);

		// If the sign of the angle from the asin and acos do not agree
		//  then make the sign of the cos the same as that of the sin
		if (((thetaDeg2 > 0) && (thetaDeg < 0)) ||
				((thetaDeg2 < 0) && (thetaDeg > 0)))
			thetaDeg = -thetaDeg;
		thetaAngle[0]=thetaDeg;
		
		if(debug)
			if ((Math.abs(fp) < 0.00001) && (Math.abs(thetaDeg-180) > 0.1))
				System.out.println("Warning: LovellRotamers: fp is small but rotation is not 180!");
				
		// if we're doing a 180 degree rotation then we can't
		//  get the axis from the cross product above as the
		//  cross product is zero, let the axis be the 32
		//  vector
		if (Math.abs(fp) < 0.00001) {
			System.out.println("Warning: GetRotationInfoB, crossproduct = 0");
			axis[0]=x32;
			axis[1]=y32;
			axis[2]=z32;
		}
		else {		 
			axis[0]=fx;
			axis[1]=fy;
			axis[2]=fz;
		}
	}

	/**
	 * Takes in the HA3 location from a glycine and transforms it into the correct 
	 * location for a CB for the corresponding amino acid.
	 * 
	 * @param CBOld HA3 atom that will be converted into a CB location
	 * @param CBNew The CB atom from the new amino acid template
	 * @param CAOld The CA atom from the glycine
	 * @param CANew The CA atom from the new amino acid template
	 */
	private static void moveCBOldGly(Atom CBOld, Atom CBNew, Atom CAOld, Atom CANew){
		double magOld = CBOld.distance(CAOld);
		double magNew = CBNew.distance(CANew);
		
		double distToMove = magNew - magOld;
		
		double[] dirToMove = new double[3];
		
		dirToMove[0] = CBOld.coord[0] - CAOld.coord[0];
		dirToMove[1] = CBOld.coord[1] - CAOld.coord[1];
		dirToMove[2] = CBOld.coord[2] - CAOld.coord[2];
		
		double scale = distToMove / magOld;
		
		dirToMove[0] *= scale;
		dirToMove[1] *= scale;
		dirToMove[2] *= scale;
		
		CBOld.coord[0] += dirToMove[0];
		CBOld.coord[1] += dirToMove[1];
		CBOld.coord[2] += dirToMove[2];
		
	}
	
	/**
	 * Takes in the a1,a2,a3 location from a residue and the original
	 * parameters for the placement of an atom and sets the the correct 
	 * original position of the destAtom
	 *  
	 * @param destAtom Atom where the coordinates will be updated atom coords will be set
	 * @param a1 First  atom in dihedral
	 * @param a2 Second atom in dihedral
	 * @param a3 Third  atom in dihedral
	 */
	private static void setOrigAtom(Atom destAtom, Atom a1, Atom a2, Atom a3, AtomPlacement atomPlace){
		double[] coords = RotMatrix.get4thPoint(a1.coord, a2.coord, a3.coord, atomPlace.len, atomPlace.ang, atomPlace.dihe);
		for(int i=0; i<destAtom.coord.length;i++){
			destAtom.coord[i] = coords[i];
		}
	}
	
	//Performs the CB alignment (changes the coordinates of the new residue r) for the new residue when the old residue is Gly
	private static void alignCBOldGly(Atom CBOld, Atom CBNew, Atom CAOld, Atom CANew, Atom NOld, Residue r){
		
	  try{
		int numAtoms = -1;
		int atomList[] = null;
		double thetaDeg[] = new double[1];
		double rotAxis[] = new double[3];
		numAtoms = r.numberOfAtoms;
		atomList = new int[numAtoms];
		
		// Compute the angle between CBOld-CA and CBNew-CA; here, CBOld is actually HA3
		double magOld = CBOld.distance(CAOld);
		double magNew = CBNew.distance(CAOld);
		double dotProd = ((CBNew.coord[0] - CAOld.coord[0]) * (CBOld.coord[0] - CAOld.coord[0])) +
			((CBNew.coord[1] - CAOld.coord[1]) * (CBOld.coord[1] - CAOld.coord[1])) +
			((CBNew.coord[2] - CAOld.coord[2]) * (CBOld.coord[2] - CAOld.coord[2]));
		double angle = 180.0 / 3.1415 * Math.acos(dotProd / magOld / magNew);
		if (Math.abs(angle) > 0.1){
			getRotationInfoB(CBOld, CAOld, NOld, CBNew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(CANew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}
	  }
	  catch(Exception e){
		  System.out.println("Why are you failing?");
	  }
	}
	
	//Performs the CB alignment (changes the coordinates of the new residue r) when the new residue is Gly
	private static void alignCBNewGly(Atom CBOld, Atom CBNew, Atom CAOld, Atom CANew, Atom NOld, 
			Atom HAOld, Residue r){
		
		//If the new residue is a Gly, we snap the HA2 and HA3 atoms
		
		// Snap the HA2
		for(int q=0;q<r.numberOfAtoms;q++){
			if( r.atom[q].name.equalsIgnoreCase("HA2") || r.atom[q].name.equalsIgnoreCase("2HA") ){
				r.atom[q].coord[0] = HAOld.coord[0];
				r.atom[q].coord[1] = HAOld.coord[1];
				r.atom[q].coord[2] = HAOld.coord[2];
			}
		}
		
		// Snap the HA3 (but note that we can't just copy coordinates
		//  because the C-H bond length is different than the C-C
		//  bond length, so we have to scale
		int localH = -1;
		for(int q=0;q<r.numberOfAtoms;q++){
			if( r.atom[q].name.equalsIgnoreCase("HA3") || r.atom[q].name.equalsIgnoreCase("3HA") )
				localH = q;
		}
		Atom tempH = r.atom[localH];
		double magNew = Math.sqrt(((tempH.coord[0] - CANew.coord[0]) * (tempH.coord[0] - CANew.coord[0])) +
			((tempH.coord[1] - CANew.coord[1]) * (tempH.coord[1] - CANew.coord[1])) +
			((tempH.coord[2] - CANew.coord[2]) * (tempH.coord[2] - CANew.coord[2])));
		double magOld = Math.sqrt(((CBOld.coord[0] - CANew.coord[0]) * (CBOld.coord[0] - CANew.coord[0])) +
			((CBOld.coord[1] - CANew.coord[1]) * (CBOld.coord[1] - CANew.coord[1])) +
			((CBOld.coord[2] - CANew.coord[2]) * (CBOld.coord[2] - CANew.coord[2])));
		r.atom[localH].coord[0] = (float)(CANew.coord[0]+(CBOld.coord[0]-CANew.coord[0])*magNew/magOld);
		r.atom[localH].coord[1] = (float)(CANew.coord[1]+(CBOld.coord[1]-CANew.coord[1])*magNew/magOld);
		r.atom[localH].coord[2] = (float)(CANew.coord[2]+(CBOld.coord[2]-CANew.coord[2])*magNew/magOld);
	}
	
	//Sets the coordinates of atom toAt to be the same as the coordinates of atom fromAt
	private static void setAtomCoord(Atom toAt, Atom fromAt){
		toAt.coord[0] = fromAt.coord[0];
		toAt.coord[1] = fromAt.coord[1];
		toAt.coord[2] = fromAt.coord[2];
	}
	
	//Translates the atom by delta
	private static void translateAtomCoord(Atom a, double[] delta){
		a.coord[0] += delta[0];
		a.coord[1] += delta[1];
		a.coord[2] += delta[2];
	}
	
	
	//This function returns a modified version of residue r by replacing the single backbone H atom
	//		with the H1, H2, and H3 atoms from localResidue, and updaing the bond information
	private static Residue getNtermRes(final Residue r, final Residue localResidue){
		
		Residue r1 = r;
		
		//Delete the HN atom from the new residue
		for(int q=0;q<r1.numberOfAtoms;q++) {
			if (r1.atom[q].name.equalsIgnoreCase("H")){
				r1.deleteAtom(q);
			}
		}
		
		//Find the old H1, H2, H3 atoms and add them to the new residue
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			if ((localResidue.atom[q].name.equalsIgnoreCase("H1"))||
					(localResidue.atom[q].name.equalsIgnoreCase("H2"))||
					(localResidue.atom[q].name.equalsIgnoreCase("H3"))){
				
				Atom lAtom = localResidue.atom[q];					
				Atom at1 = new Atom(lAtom.name,lAtom.coord[0],lAtom.coord[1],lAtom.coord[2],lAtom.charge,lAtom.forceFieldType);					
				r1.addAtom(at1);
			}
		}
		
		//Add the new bonds between the H1, H2, and H3 atoms and the N atom
		for(int q=0;q<r1.numberOfAtoms;q++) {
			if ((r1.atom[q].name.equalsIgnoreCase("H1"))||(r.atom[q].name.equalsIgnoreCase("H2"))||
					(r1.atom[q].name.equalsIgnoreCase("H3"))){
				
				r1.atom[q].addBond(0); //add the bond for this H atom to the N atom
				r1.atom[0].addBond(q); //add the bond for the N atom to this H atom
			}
		}		
		
		return r1;
	}
	
	private static Residue getCtermRes(final Residue r, final Residue localResidue){
		
		Residue r1 = r;
		
		
		//Find the old H1, H2, H3 atoms and add them to the new residue
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			if (localResidue.atom[q].name.equalsIgnoreCase("OXT")||
					localResidue.atom[q].name.equalsIgnoreCase("O2")){
				
				Atom lAtom = localResidue.atom[q];					
				Atom at1 = new Atom(lAtom.name,lAtom.coord[0],lAtom.coord[1],lAtom.coord[2],lAtom.charge,lAtom.forceFieldType);					
				r1.addAtom(at1);
			}
		}
		
		int cAtom = -1;
		//find which atom is the C to connect the OXT to
		for(int q=0; q<r1.numberOfAtoms;q++){
			if(r1.atom[q].name.equalsIgnoreCase("C")){
				cAtom = q;
				break;
			}
		}
		
		//update the bond information for the new atom
		for(int q=0; q<r1.numberOfAtoms;q++){
			if (r1.atom[q].name.equalsIgnoreCase("OXT")||r1.atom[q].name.equalsIgnoreCase("O2")){
				r1.atom[q].addBond(cAtom); //add the bond to the new C atom
				r1.atom[cAtom].addBond(q);
			}
		}		
		
		return r1;
	}
	
	// This function returns a list of atoms from residue resNum that are
	//  further distal than at2, atomList for a2 == 0 so we ignore it
	// atomList elements are 0=ignore, 1=unchecked, 2=include
	public static void getAtomsMoreDistal(Molecule m, int resNum, Atom at2, int atomList[]) {
		
		if (at2.bond != null) {
			for(int q=0;q<at2.bond.length;q++) {
					Atom at = m.atom[at2.bond[q]];
					if (at.moleculeResidueNumber==resNum){
					if (atomList[at.residueAtomNumber] == 1) {
						atomList[at.residueAtomNumber] = 2;
						getAtomsMoreDistal(m,resNum,at,atomList);
					}
				}
			}
		}
	}
	
	// This function converts residue resNum to the residue conformation specified by RCNum
	//It will start with the actualCoordinates but will change both the molecule actualCoordinates and the atom coord arrays
	//(like applyRotamer)
	public static boolean applyRC(Molecule m, Residue res, ResidueConformation rc) {


		AARotamerType AAType = rc.rot.aaType; //rl.getAARotamerIndex(m.strand[strandNumber].residue[resNum].name);
		Rotamer rot = rc.rot;//RCRots[resNum][AANum][RCNum];
		int pertState = rc.pertState;//RCPertStates[resNum][AANum][RCNum];

		boolean outcome = true;

		//                if(rotNum == -2){//WT rotamer
		//
		//                    Residue localRes=m.strand[strandNumber].residue[resNum];
		//
		//                    if( WTRes[resNum] == null ){
		//                        System.err.println("Error: WT rotamer being called but AddWTRotamers option is off");
		//                        System.exit(1);
		//                    }
		//
		//                    if( !isProtein ){
		//                        System.err.println("Error: WT rotamer not supported for non-protein: " + localRes.fullName);
		//                        System.exit(1);
		//                    }
		//
		//                    //MH 2012: It's an error to be assigning a WT rotamer to a mutant residue
		//                    //But sometimes CYX gets called
		//                    if( !localRes.name.equalsIgnoreCase( WTRes[resNum].name ) ){
		//                        System.err.println("Error: WT rotamer assigned to mutant residue "+localRes.fullName+", wild type " + WTRes[resNum].name );
		//                        System.exit(1);
		//                    }
		//
		//                    RotMatrix r = new RotMatrix();
		//                    double CACoord[] = m.getActualCoord( localRes.getAtomNameToMolnum("CA") );
		//
		//                    //First put the atoms in in the original PDB orientation
		//                    for(int a=0;a<localRes.numberOfAtoms;a++){
		//                        Atom localAt = localRes.atom[a];
		//                        if( !localAt.isBBatom ){//BB atoms shouldn't be restored
		//                            Atom WTAt = WTRes[resNum].atom[a];
		//                            if( ! WTAt.name.equalsIgnoreCase( localAt.name ) )//Different ordering of atoms
		//                                WTAt  = WTRes[resNum].getAtomByName(localAt.name);
		//                            localAt.coord = r.add(WTAt.coord, CACoord);
		//                        }
		//                    }
		//
		//
		//                    for(int a=0;a<localRes.numberOfAtoms;a++){
		//                        if( !localRes.atom[a].isBBatom )
		//                            m.updateCoordinates(localRes.atom[a].moleculeAtomNumber);
		//                    }
		//                    //Update the actualCoordinates
		//
		//
		//                    //BB conformational changes might make this orientation inappropriate though
		//                    //So idealize the sidechain to fix it
		//                    //This won't return false because it's not a proline
		//                    outcome = outcome && m.idealizeResSidechain(localRes);
		//
		//                    //Finally make sure the (generalized) chi1 is correct
		//                    Perturbation.setGenChi1(m, localRes.moleculeResidueNumber, WTGenChi1[resNum] );
		//
		//                    curRotNum[resNum] = -2;
		//                }
		//                else
		applyRotamer(m, rot, res);
		//applyRotamer returns false for conditions like the AA having no rotamers,
		//e.g. applying "rotamer 0" of glycine, which do not actually mean RC application failed
		//So we do not store its return value

		outcome = outcome && applyPertState(m, res, pertState);

		m.resolveCoordinates();//Get the perturbed coordinates into the Atom.coord arrays


		res.curRC = rc;

		return outcome;
	}
	
	//Apply a residue perturbation state to the m.actualCoordinates
	//The false return value here is used to detect bad perturbation states
	public static boolean applyPertState(Molecule m, Residue localRes, int pertState){

//		Residue localRes=m.strand[strandNumber].residue[resNum];

		boolean outcome = true;//Indicates success or failure

		if( localRes.perts.length>0 ){//Apply perturbation state

			int pertInd = localRes.perts.length - 1;//Index in localRes.perts
			int affectedInd = localRes.affectedPerts.length - 1;//Index in localRes.affectedPerts

			//First see what perturbations we have to change
			int firstPert = m.perts.length;

			for(int a=0;a<localRes.perts.length;a++){
				Perturbation pert = m.perts[localRes.perts[a]];
				int thisState = localRes.pertStates[pertState][a];
				if( pert.curParam != ( pert.minParams[thisState] + pert.maxParams[thisState] ) / 2 ){//The perturbation is not at the default parameter value for the desired state
					firstPert = localRes.perts[a];
					break;
				}
			}

			for(int pertNum=m.perts.length-1; pertNum>=firstPert; pertNum--){//Undoing perturbations, in reverse order
				if(pertNum == localRes.perts[pertInd]){
					m.perts[pertNum].undo();
					pertInd--;
				}
				else if(affectedInd > -1){//Still have affected perturbations to look for
					if( pertNum == localRes.affectedPerts[affectedInd] ){
						m.perts[pertNum].undo();
						affectedInd--;
					}
				}
			}

			if(pertInd < localRes.perts.length - 1)
				pertInd++;//It's now the first perturbation to reapply
			if(affectedInd < localRes.affectedPerts.length - 1 )
				affectedInd++;//The first affected perturbation to reapply


			for(int pertNum=firstPert; pertNum<m.perts.length;pertNum++){//Redoing them in the correct state
				if(pertNum==localRes.perts[pertInd]){
					Perturbation pert = m.perts[pertNum];

					if(outcome){
						int newState = localRes.pertStates[pertState][pertInd];
						outcome = pert.applyPerturbation( ( pert.minParams[newState] + pert.maxParams[newState] ) / 2 );
						if(outcome)
							pert.curState = newState;//The actual state of the perturbation (not the residue perturbation state)
							else
								pert.curState = 0;
					}
					else{
						pert.applyPerturbation(0);//Don't waste time applying this invalid perturbation state
						pert.curState = 0;//Put it into an unperturbed state
					}

					if(pertInd < localRes.perts.length - 1)
						pertInd++;
				}
				else if(affectedInd > -1){
					if(pertNum==localRes.affectedPerts[affectedInd]){

						if(outcome)
							outcome = m.perts[pertNum].applyPerturbation(m.perts[pertNum].curParam);
						else{
							m.perts[pertNum].applyPerturbation(0);
							m.perts[pertNum].curState = 0;
						}


						if( affectedInd < localRes.affectedPerts.length - 1 )
							affectedInd++;
					}
				}
			}

		}

//		curPertState[resNum] = pertState;

		return outcome;
	}
	
	
}
