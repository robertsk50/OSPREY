import static org.junit.Assert.*;

import java.io.IOException;
import java.net.URL;

import org.junit.Before;
import org.junit.Test;


public class BackrubMutationTest {

	KSParser ksParser;
	KSParser.MolParameters mp;
	ParamSet sParams;
	ParamSet rParams;
	Settings settings;
	double tol = 0.001;

	@Before
	public void loadParams(){
		ksParser = new KSParser();
		ksParser.cfgName = "bin/staticFiles/KStar.cfg";

		ksParser.setConfigPars(); //Setup the KStar.cfg parameters

		sParams = new ParamSet();
		try{
			sParams.addParamsFromFile(ClassLoader.class.getResourceAsStream("/staticFiles/brMut/System.cfg"));
			sParams.addParamsFromFile(ClassLoader.class.getResourceAsStream("/staticFiles/brMut/MutSearch.cfg"));
		}catch(IOException e){
			System.out.println("Trouble reading param files");
		}

		Settings settings = new Settings();
		String runName = settings.getRunName(sParams);
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		Settings.Emat ematSettings = settings.new Emat(sParams, runName, true);
		mp = ksParser.loadMolecule(sParams, KSParser.COMPLEX, graphSettings.neighborList,graphSettings.distCutoff,true);

		ksParser.selectPerturbations(mp, true, "bin/staticFiles/brMut/pertFile.pert", false, ematSettings.addWTRot, sParams);

	}


	/*
	 * Test whether the native backbone coordinates can be recovered after mutation all amino acids
	 */
	@Test
	public void simpleBackrub() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("38"); //pdb number 38

		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		ResidueConformation[] brRC = new ResidueConformation[mp.strandMut.allMut.length];
		//Find BR Residue Conformation
		//The goal is to turn pert 0 (i.e. 38-39-40 5deg backrub on)
		//and have everything else turned off

		int ctr = 0;
		for(ResidueConformation rc: mp.m.strand[0].rcl.allRCs){
			if(ctr < mp.strandMut.allMut.length && rc.res == mp.strandMut.allMut[ctr] && rc.rot.isWTrot){
				Residue.PertPair[] pertPair = mp.m.residue[mp.strandMut.allMut[ctr]].pertsForPertState(rc.pertState);
				boolean correctBR = true;
				for(int i=0; i<pertPair.length;i++){
					if( (pertPair[i].pertID == 0 && pertPair[i].status != 2) ||
							pertPair[i].pertID != 0 && pertPair[i].status != 0)
						correctBR = false;
				}

				if(correctBR)
					brRC[ctr++] = rc;
			}
		}

		mp.m.saveMolecule("beforeRC.pdb", 0.0);
		for(int i=0; i<brRC.length; i++){
			MutUtils.applyRC(mp.m, mp.m.residue[mp.strandMut.allMut[i]], brRC[i]);
			mp.m.saveMolecule("simpleBackrub_"+i+".pdb", 0.0);
		}
		mp.m.saveMolecule("afterRC.pdb", 0.0);

		//		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "PRO", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "GLY", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "VAL", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "LEU", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "ILE", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "PHE", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "TYR", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "TRP", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "CYS", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "MET", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "SER", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "THR", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "LYS", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "ARG", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "HIP", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "HID", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "HIE", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "ASP", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "GLU", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "ASN", true);
		//		MutUtils.changeResidueType(mp.m, mutRes, "GLN", true);
		//		
		//		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
		//
		//		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		//		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		//		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		//		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		//		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		//		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		//		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		//		
		//		assertEquals("H is recovered.",0,Hdist,tol);
		//		assertEquals("N is recovered.",0,Ndist,tol);
		//		assertEquals("CA is recovered.",0,CAdist,tol);
		//		assertEquals("HA is recovered.",0,HAdist,tol);
		//		assertEquals("C is recovered.",0,Cdist,tol);
		//		assertEquals("O is recovered.",0,Odist,tol);
		//		assertEquals("CB is recovered.",0,CBdist,tol);

		fail("Test not implemented yet");

	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation all amino acids
	 */
	@Test
	public void undoPerturbations() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("38"); //pdb number 38

		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		ResidueConformation[] brRC = new ResidueConformation[mp.strandMut.allMut.length];
		//Find BR Residue Conformation
		//The goal is to turn pert 0 (i.e. 38-39-40 5deg backrub on)
		//and have everything else turned off

		int ctr = 0;
		for(ResidueConformation rc: mp.m.strand[0].rcl.allRCs){
			if(ctr < mp.strandMut.allMut.length && rc.res == mp.strandMut.allMut[ctr] && rc.rot.isWTrot){
				Residue.PertPair[] pertPair = mp.m.residue[mp.strandMut.allMut[ctr]].pertsForPertState(rc.pertState);
				boolean correctBR = true;
				for(int i=0; i<pertPair.length;i++){
					if( (pertPair[i].pertID == 0 && pertPair[i].status != 0) ||
							pertPair[i].pertID != 0 && pertPair[i].status != 0)
						correctBR = false;
				}

				if(correctBR)
					brRC[ctr++] = rc;
			}
		}

		mp.m.saveMolecule("undoPerturbations_before.pdb", 0.0);
		for(int i=0; i<brRC.length; i++){
			MutUtils.applyRC(mp.m, mp.m.residue[mp.strandMut.allMut[i]], brRC[i]);
			mp.m.saveMolecule("undoPerturbations_"+i+".pdb", 0.0);
		}
		mp.m.saveMolecule("undoPerturbations_after.pdb", 0.0);

		fail("Test not implemented yet");

	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation all amino acids
	 */
	@Test
	public void undoThenBackrub() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("38"); //pdb number 38

		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		ResidueConformation[] brRC = new ResidueConformation[mp.strandMut.allMut.length];
		//Find BR Residue Conformation
		//The goal is to turn pert 0 (i.e. 38-39-40 5deg backrub on)
		//and have everything else turned off

		int ctr = 0;
		for(ResidueConformation rc: mp.m.strand[0].rcl.allRCs){
			if(ctr < mp.strandMut.allMut.length && rc.res == mp.strandMut.allMut[ctr] && rc.rot.isWTrot){
				Residue.PertPair[] pertPair = mp.m.residue[mp.strandMut.allMut[ctr]].pertsForPertState(rc.pertState);
				boolean correctBR = true;
				for(int i=0; i<pertPair.length;i++){
					if( (pertPair[i].pertID == 0 && pertPair[i].status != 0) ||
							pertPair[i].pertID != 0 && pertPair[i].status != 0)
						correctBR = false;
				}

				if(correctBR)
					brRC[ctr++] = rc;
			}
		}

		for(int i=0; i<brRC.length; i++){
			MutUtils.applyRC(mp.m, mp.m.residue[mp.strandMut.allMut[i]], brRC[i]);
		}

		brRC = new ResidueConformation[mp.strandMut.allMut.length];
		//Find BR Residue Conformation
		//The goal is to turn pert 0 (i.e. 38-39-40 5deg backrub on)
		//and have everything else turned off

		ctr = 0;
		for(ResidueConformation rc: mp.m.strand[0].rcl.allRCs){
			if(ctr < mp.strandMut.allMut.length && rc.res == mp.strandMut.allMut[ctr] && rc.rot.isWTrot){
				Residue.PertPair[] pertPair = mp.m.residue[mp.strandMut.allMut[ctr]].pertsForPertState(rc.pertState);
				boolean correctBR = true;
				for(int i=0; i<pertPair.length;i++){
					if( (pertPair[i].pertID == 0 && pertPair[i].status != 2) ||
							(pertPair[i].pertID == 1 && pertPair[i].status != 0) ||
							(pertPair[i].pertID == 2 && pertPair[i].status != 2) ||
							(pertPair[i].pertID >= 3 && pertPair[i].status != 0))
						correctBR = false;
				}

				if(correctBR)
					brRC[ctr++] = rc;
			}
		}

		for(int i=0; i<brRC.length; i++){
			MutUtils.applyRC(mp.m, mp.m.residue[mp.strandMut.allMut[i]], brRC[i]);
		}
		mp.m.saveMolecule("redoBackrub_after.pdb",0.0);
		
		fail("Test not implemented yet");

	}

}
