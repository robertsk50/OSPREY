import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.StringTokenizer;


public class CompareEmats {
	final static int COMPLEX = -1;
	static Molecule[] mols = new Molecule[3];
	public static void main(String[] args) {
		
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(args[0]); //read system parameters
		sParams.addParamsFromFile(args[1]); //read mutation search parameters
		
		
//		String runName = (String)sParams.getValue("RUNNAME");
//		boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS"))).booleanValue();//Triggers DEEPer
//		String pertFile = (String)sParams.getValue("PERTURBATIONFILE");//Input file giving perturbation information
//		
//		MolParameters mp = loadMolecule(sParams, -1, false, -1,true);
//		Molecule m = mp.m;
//		int numberMutable = mp.strandMut.numMutPos();
//		MutableResParams strandMut = mp.strandMut;
//		int numOfStrands = mp.m.numberOfStrands;
//
//		StrandRotamers[] strandRot = new StrandRotamers[mp.m.numberOfStrands];
//		for(int i=0; i<mp.m.numberOfStrands;i++){
//			if(doPerturbations)
//				strandRot[i] = new StrandRCs(m.rotLibForStrand(i),m.strand[i]);
//			else
//				strandRot[i] = new StrandRotamers(m.rotLibForStrand(i),m.strand[i]);
//		}
//		
//		//Load rotamer and residue conformation libraries
//		String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
//		m.aaRotLib.loadGlobalRots(runNameEMatrixMin+"_COM.dat.aaRots");
//		if(m.genRotLib != null)
//			m.genRotLib.loadGlobalRots(runNameEMatrixMin+"_COM.dat.genRots");
//		PertFileHandler.readPertFile(pertFile, m, strandRot,true);
//		for(Strand strand : m.strand){
//			if(strand.isProtein)
//				strand.rcl.loadGlobalRCs(runNameEMatrixMin+"_COM.dat.rcl_"+strand.number, m.aaRotLib);
//			else
//				strand.rcl.loadGlobalRCs(runNameEMatrixMin+"_COM.dat.rcl_"+strand.number, m.genRotLib);
//		}
//		
		RotamerLibrary aaRotLib = new RotamerLibrary("/usr/project/dlab/Users/kroberts/Troubleshooting/Kyle/082714/DEEPERepic/dataFiles/LovellRotamer.dat", false);
//		PertFileHandler.readPertFile(pertFile, m, strandRot,true);
		ResidueConformationLibrary rcl = new ResidueConformationLibrary();
		rcl.loadGlobalRCs("/usr/project/dlab/Users/kroberts/Troubleshooting/Kyle/082714/DEEPERepic/merged/1CC8_pertminM_COM.dat.rcl_0", aaRotLib);
		
		
		Emat emat = new Emat("1CC8_pertminM_COM.dat",false,null);
		PairwiseEnergyMatrix pairEmat =  new PairwiseEnergyMatrix();
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream("/usr/project/dlab/Users/kroberts/Troubleshooting/Kyle/082714/DEEPERepic/github/1CC8_pertminM_COM.dat"));
			pairEmat.eMatrix = (double [][][][][][])in.readObject();
			in.close();
		}catch(Exception e){}
				
				
		
		for(int p1=0; p1<pairEmat.eMatrix.length-1;p1++){
			for(int a1=0; a1<pairEmat.eMatrix[p1].length;a1++){
				if(pairEmat.eMatrix[p1][a1] == null)
					continue;
				for(int r1=0; r1<pairEmat.eMatrix[p1][a1].length;r1++){
					double e1 = pairEmat.getIntraAndShellE(p1, a1, r1);
					AARotamerType aaType = aaRotLib.getAAType(a1);
					int state = 0;
					int rot = r1;
					while(rot > aaType.numRotamers()-1){
						rot -= aaType.numRotamers();
						state++;
					}
					
					int globalRC = -1; 
					while(globalRC == -1){
						globalRC = rcl.getGlobalIndex(emat.resByPos.get(p1).get(0),aaType.rotamers.get(rot).rlIndex, state );
						//If globalRC doesn't exist in rcl then it was not a valid state (the previous while loop breaks when more than one pert is invalid)
						state++;
					}
					
					
					double e2=0.0;
					Index3 RCind;
					
					RCind = emat.getGlobalRotForPos(p1, globalRC);
					if(RCind == null){
						if(e1 < 100)
							System.out.println("Energy: "+e1 +" was pruned");
					}else{
						e2 = emat.getSingleMinE(RCind);
						
						double diff = e1 - e2;
						if(Math.abs(diff) > 0.01){
							System.out.println(p1+"_"+a1+"_"+r1+" Diff: "+diff);
						}
						
						for(int p2=p1+1; p2<pairEmat.eMatrix.length-1;p2++){
							for(int a2=0; a2<pairEmat.eMatrix[p1][a1][r1][p2].length;a2++){
								if(pairEmat.eMatrix[p1][a1][r1][p2][a2] == null)
									continue;
								for(int r2=0; r2<pairEmat.eMatrix[p1][a1][r1][p2][a2].length;r2++){
									double pair_e1 = pairEmat.getPairwiseE(p1, a1, r1,p2,a2,r2);
									AARotamerType aaType2 = aaRotLib.getAAType(a2);
									int state2 = 0;
									int rot2 = r2;
									while(rot2 > aaType2.numRotamers()-1){
										rot2 -= aaType2.numRotamers();
										state2++;
									}
									
									int globalRC2 = -1; 
									while(globalRC2 == -1){
										globalRC2 = rcl.getGlobalIndex(emat.resByPos.get(p2).get(0),aaType2.rotamers.get(rot2).rlIndex, state2 );
										//If globalRC doesn't exist in rcl then it was not a valid state (the previous while loop breaks when more than one pert is invalid)
										state2++;
									}
									
									
									
									Index3 RCind2 = emat.getGlobalRotForPos(p2, globalRC2);
									if(RCind2 != null){
										double pair_e2 = emat.getPairMinE(RCind, RCind2);
										double pair_diff = pair_e1 - pair_e2;
										
										if(Math.abs(pair_diff) > 0.01)
											System.out.println(p1+"_"+a1+"_"+r1+"_"+p2+"_"+a2+"_"+r2+" Diff: "+pair_diff);
									}
								}
							}
						}
					}
					
				}
			}
		}
				
				
	}
	
	public static MolParameters loadMolecule(ParamSet sParams, int curStrForMatrix, 
			boolean neighborList, double distCutoff, boolean useCache){	

		MolParameters mp = new MolParameters();

		loadStrandParams(sParams, mp, curStrForMatrix);

		//Setup the molecule system

		//Setup the molecule system
		if(mols[curStrForMatrix+1] == null || !useCache){
			mp.m = new Molecule();
			mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits, true);
			mols[curStrForMatrix+1] = mp.m;
		}else{
			mp.m = mols[curStrForMatrix+1];
			mp.loadedFromCache = true;
		}

		loadMutationParams(sParams, mp);

		if(neighborList){
			mp.m.genDistNeighborList(distCutoff);
			mp.m.dumpNeighborList();
		}

		return mp;

	}


	private int getNumberMutable(int[][] strandMut){
		int numberMutable = 0;
		for (int i=0; i<strandMut.length;i++)
			numberMutable += strandMut[i].length;

		return numberMutable;
	}

	private static void loadStrandParams(ParamSet sParams, MolParameters mp, int curStrForMatrix){
		mp.numOfStrands = (new Integer((String)sParams.getValue("NUMOFSTRANDS"))).intValue();	
		mp.strandLimits = new String[mp.numOfStrands][2];
		mp.strandPresent = new boolean[mp.numOfStrands];
		mp.strandsPresent = 0;
		for (int i=0; i<mp.numOfStrands; i++){
			String strandLimit = (String)sParams.getValue("STRAND"+i);
			String limit1 = getToken(strandLimit,1);
			String limit2 = getToken(strandLimit,2);
			mp.strandLimits[i][0] = limit1;
			mp.strandLimits[i][1] = limit2;
			if(curStrForMatrix == COMPLEX)
				mp.strandPresent[i] = true;
			else if(curStrForMatrix == i)
				mp.strandPresent[i] = true;
			else
				mp.strandPresent[i] = false;
			//strandPresent[i] = (new Boolean((String)sParams.getValue("STRANDPRESENT" + i))).booleanValue();
			if(mp.strandPresent[i]){
				mp.strandsPresent++;
			}
		}
	}



	//KER: mp must already have the following terms set:
	//mp.numOfStrands	
	//mp.strandLimits
	//mp.strandPresent
	//mp.strandsPresent
	private static void loadMutationParams(ParamSet sParams, MolParameters mp) {

		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","true"))).booleanValue();
		/**********Get the regions of each strand that are mutable****/
		String strandMutNums = (String)sParams.getValue("STRANDMUTNUMS");
		boolean addOrigRots = (new Boolean((String)sParams.getValue("ADDWTROTS", "false"))).booleanValue();
		int totalNumMut = 0;
		for(int i=0; i<mp.strandPresent.length; i++){
			if(mp.strandPresent[i])
				totalNumMut += (new Integer(KSParser.getToken(strandMutNums,i+1))).intValue();
		}
		mp.strandMut = new MutableResParams(totalNumMut,mp.m.numberOfStrands);  //taking the place of resMap and ligMap

		int flatCtr = 0;
		int strCtr = 0;
		try{
			for(int i=0; i<mp.strandPresent.length; i++){
				if(mp.strandPresent[i]){
					int numberOfMutables = (new Integer(getToken(strandMutNums,i+1))).intValue();
					//mp.strandMut[ctr] = new int[numberOfMutables];
					String strandMutResNum = (String)sParams.getValue("STRANDMUT"+i);
					for(int j=0; j<numberOfMutables; j++){
						String strandMutRes = getToken(strandMutResNum,j+1);
						Residue r = mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)];
						Residue prevRes = null;
						if(mp.m.residuesAreBBbonded(r.moleculeResidueNumber-1,r.moleculeResidueNumber))
							prevRes = mp.m.residue[r.moleculeResidueNumber-1];
						mp.strandMut.addRes(flatCtr,r,mp.m.rotLibForStrand(strCtr),addOrigRots,prevRes);
						if(mp.m.strand[r.strandNumber].isProtein)
							mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)].canMutate = true;
						flatCtr++;
					}
					strCtr++;
				}
			}
			strCtr = 0;
			//Say rotamers have been set
			for(int i=0; i<mp.strandPresent.length;i++){
				if(mp.strandPresent[i]){
					mp.m.rotLibForStrand(strCtr).setAddedRotamers(true);
					strCtr++;
				}
			}
			if(!addWT)
				mp.strandMut.checkWT(mp.strandPresent, sParams);
		}

		catch(Exception E){
			System.out.println("PROBLEM Loading the strand mut numbers (Check System.cfg) ");
			E.printStackTrace();
			return;
		}
	}
	
	static // This function returns the xth token in string s
	String getToken(String s, int x) {

		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
				st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	
	//Sets up the molecule system and returns the number of ligand rotamers
		private static Molecule setupMolSystem(Molecule m, ParamSet sParams, boolean[] strandPresent, String[][] strandLimits, boolean keepCofactor){

			int ligStrNum = -1;
			int cofStrNum = -1;

			try{
				FileInputStream is = new FileInputStream((String)sParams.getValue("PDBNAME"));
				new PDBChemModel(m, is);
			}
			catch (Exception e){
				System.out.println("WARNING: An error occurred while reading file");
				System.out.println(e);
				e.printStackTrace();
				System.exit(1);
			}

			/*m.strand[sysStrNum].isProtein = strandPresent[sysStrNum];		// main protein
			m.strand[sysStrNum].rotTrans = (new Boolean((String)sParams.getValue("STRANDROTTRANS"+sysStrNum))).booleanValue();
			strandLength[sysStrNum] = m.mapPDBresNumToMolResNum(strandLimits[0][1])-m.mapPDBresNumToMolResNum(strandLimits[0][0])+1;*/
			int strNum = 0; //the current strand number; 0 is reserved for the protein strand, the ligand strand is 1 (if present)

			//get the number of strands that are present
			int numPresent = 0;
			for(int str=0;str<strandPresent.length;str++)
				if(strandPresent[str])
					numPresent++;

			Molecule newMol = new Molecule();

			//int curPresStr = 0;
			for(int i=0; i<strandLimits.length;i++){
				String pdbStart = strandLimits[i][0]; 
				String pdbEnd   = strandLimits[i][1];
				//if (pdbEnd>=0){ //with ligand in PDB
				int molStartNum = m.mapPDBresNumToMolResNum(pdbStart);
				int molEndNum   = m.mapPDBresNumToMolResNum(pdbEnd);
				if(molStartNum <0 || molEndNum < 0){
					System.out.println("Please make sure strand "+i+"'s begin and end are set properly.");
					System.exit(0);
				}
				int numInStrand = molEndNum-molStartNum+1;
				//strandLength[i] = numInStrand;
				Residue myStrand[] = new Residue[numInStrand];
				for (int j=(numInStrand-1); j>=0; j--){
					myStrand[j] = m.residue[molStartNum+j]; // pull out the ligand
					m.deleteResidue(molStartNum+j);
					myStrand[j].renumberResidue();
				}
				if (strandPresent[i]) { //ligand will be used in design

					newMol.addStrand(""+i);
					strNum = newMol.strand.length-1;
					newMol.strand[strNum].rotTrans = (new Boolean((String)sParams.getValue("STRANDROTTRANS"+i))).booleanValue();

					for(int j=0; j<numInStrand;j++)
						newMol.addResidue(strNum,myStrand[j],false);

					newMol.strand[strNum].isProtein = (new Boolean((String)sParams.getValue("STRANDAA"+i))).booleanValue();
					if (newMol.strand[strNum].isProtein) //use the AA rotamer library for the ligand
						newMol.setAARotLib(EnvironmentVars.aaRotLibFile);
					else //use the non-AA rotamer library for the ligand
						newMol.setGenRotLib(sParams.getValue("GROTFILE","GenericRotamers.dat"));

					//change ligand to the specified residue type
					//if ( m.strand[ligStrNum].isProtein && !m.strand[ligStrNum].residue[0].name.equalsIgnoreCase(ligType) ) //not the same ligand type
					//	(new StrandRotamers(grl,m.strand[ligStrNum])).changeResidueType(m,0,ligType,true);

					strNum++;
					//curPresStr++;

				}
				//}
				/*else if (strandPresent[i]){
				System.out.println("ERROR: Attempting to use a ligand, but ligand not found in system config file");
				System.exit(1);
			}*/
			}


			//Get the cofactors (if present)
			if(keepCofactor){
				int str=-1; 
				for(int i=0; i<strandPresent.length; i++){
					if(strandPresent[i]){
						str++;
						String cofMapString = sParams.getValue("COFMAP"+i, "-1");
						//int numCofactorRes = (new Integer((String)sParams.getValue("NUMCOFRES"))).intValue();

						Residue cof;
						for(int res=0;res<numTokens(cofMapString);res++){
							int cofactorRes = m.mapPDBresNumToMolResNum(getToken(cofMapString, res));
							if(cofactorRes >= 0){
								cof = m.residue[cofactorRes];
								cof.cofactor = true;
								m.deleteResidue(cofactorRes);
								newMol.addResidue(str,cof,false);
							}
						}
					}
				}
			}

			newMol.determineBonds();
			newMol.establishConnectivity(false);

			return newMol;

			//Determine the number of rotamers for the ligand (if used)
			/*int numLigRotamers = 0;
			if (useLig) {
				numLigRotamers = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
				if (numLigRotamers == 0)
					numLigRotamers = 1;
			}
			return numLigRotamers;*/
		}
		
		static // This function returns the number of tokens in string s
		int numTokens(String s) {

			int curNum = 0;	
			StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

			while (st.hasMoreTokens()) {
				curNum++;
				st.nextToken();
			}
			return(curNum);
		}

	
}

class MolParameters{
	String[][] strandLimits;
	int numOfStrands;
	boolean[] strandPresent;
	int strandsPresent;
	MutableResParams strandMut;
	String[][] strandDefault;
	Molecule m;
//	int numberMutable;
	boolean loadedFromCache = false;
	
	MolParameters(){
		
	}
	
	public MolParameters copy(){
		MolParameters retMP = new MolParameters();
		retMP.strandLimits = strandLimits;
		retMP.numOfStrands = numOfStrands;
		retMP.strandPresent = strandPresent;
		retMP.strandsPresent = strandsPresent;
		retMP.strandMut = strandMut;
		retMP.strandDefault = strandDefault;
		retMP.m = m;
		retMP.loadedFromCache = loadedFromCache;
		return retMP;
	}
	
}