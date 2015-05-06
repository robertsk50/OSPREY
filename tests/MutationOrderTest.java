import static org.junit.Assert.*;

import java.io.IOException;
import java.net.URL;
import java.util.HashMap;

import org.junit.Before;
import org.junit.Test;


public class MutationOrderTest {

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
		sParams.addParamsFromFile(ClassLoader.class.getResourceAsStream("/staticFiles/System_mutOrder.cfg"));
		sParams.addParamsFromFile(ClassLoader.class.getResourceAsStream("/staticFiles/MutSearch_mutOrder.cfg"));
		}catch(IOException e){
			System.out.println("Trouble reading param files");
		}
		
		Settings settings = new Settings();
		
		//InteractionGraph Settings
		Settings.InteractionGraph graphSettings = settings.new InteractionGraph(sParams);
		//InteractionGraph Settings
		Settings.Minimization minSettings = settings.new Minimization(sParams);
		
		mp = ksParser.loadMolecule(sParams, KSParser.COMPLEX, graphSettings.neighborList,graphSettings.distCutoff,true);
		
		//If molecule was reloaded we need to update the mutable information
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
		if(!addWT)
			mp.strandMut.checkWT(mp.strandPresent, sParams);
		int molStrand = 0;
		for(int resID:mp.strandMut.allMut){
			ksParser.setAllowablesHelper(sParams, addWT, mp.m.residue[resID]);
		}
		
		RotamerSearch.setupRCs(minSettings.doPerturbations, mp.m, minSettings.pertFile);
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation to a glycine
	 */
	@Test
	public void glyMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41
		
		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
//		mp.m.saveMolecule("beforeMut.pdb", 0.0);
		
		MutUtils.changeResidueType(mp.m, mutRes, "GLY", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		
//		mp.m.saveMolecule("afterMut.pdb", 0.0);
		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		
		assertEquals("H is recovered.",0,Hdist,tol);
		assertEquals("N is recovered.",0,Ndist,tol);
		assertEquals("CA is recovered.",0,CAdist,tol);
		assertEquals("HA is recovered.",0,HAdist,tol);
		assertEquals("C is recovered.",0,Cdist,tol);
		assertEquals("O is recovered.",0,Odist,tol);
		assertEquals("CB is recovered.",0,CBdist,tol);
		
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation to a proline
	 */
	@Test
	public void proMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41
		
		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
		
		MutUtils.changeResidueType(mp.m, mutRes, "PRO", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		
		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		
		assertEquals("H is recovered.",0,Hdist,tol);
		assertEquals("N is recovered.",0,Ndist,tol);
		assertEquals("CA is recovered.",0,CAdist,tol);
		assertEquals("HA is recovered.",0,HAdist,tol);
		assertEquals("C is recovered.",0,Cdist,tol);
		assertEquals("O is recovered.",0,Odist,tol);
		assertEquals("CB is recovered.",0,CBdist,tol);
		
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation to a 
	 * glycine followed by a proline.
	 */
	@Test
	public void glyProMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41
		
		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
		
		MutUtils.changeResidueType(mp.m, mutRes, "GLY", true);
		MutUtils.changeResidueType(mp.m, mutRes, "PRO", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		
		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		
		assertEquals("H is recovered.",0,Hdist,tol);
		assertEquals("N is recovered.",0,Ndist,tol);
		assertEquals("CA is recovered.",0,CAdist,tol);
		assertEquals("HA is recovered.",0,HAdist,tol);
		assertEquals("C is recovered.",0,Cdist,tol);
		assertEquals("O is recovered.",0,Odist,tol);
		assertEquals("CB is recovered.",0,CBdist,tol);
		
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation to a proline
	 * followed by a glycine
	 */
	@Test
	public void proGlyMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41
		
		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
		
		MutUtils.changeResidueType(mp.m, mutRes, "PRO", true);
		MutUtils.changeResidueType(mp.m, mutRes, "GLY", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		
		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		
		assertEquals("H is recovered.",0,Hdist,tol);
		assertEquals("N is recovered.",0,Ndist,tol);
		assertEquals("CA is recovered.",0,CAdist,tol);
		assertEquals("HA is recovered.",0,HAdist,tol);
		assertEquals("C is recovered.",0,Cdist,tol);
		assertEquals("O is recovered.",0,Odist,tol);
		assertEquals("CB is recovered.",0,CBdist,tol);
		
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after mutation all amino acids
	 */
	@Test
	public void everyMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41
		
		Atom[] prevBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);
		
		MutUtils.changeResidueType(mp.m, mutRes, "ALA", true);
		MutUtils.changeResidueType(mp.m, mutRes, "PRO", true);
		MutUtils.changeResidueType(mp.m, mutRes, "GLY", true);
		MutUtils.changeResidueType(mp.m, mutRes, "VAL", true);
		MutUtils.changeResidueType(mp.m, mutRes, "LEU", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ILE", true);
		MutUtils.changeResidueType(mp.m, mutRes, "PHE", true);
		MutUtils.changeResidueType(mp.m, mutRes, "TYR", true);
		MutUtils.changeResidueType(mp.m, mutRes, "TRP", true);
		MutUtils.changeResidueType(mp.m, mutRes, "CYS", true);
		MutUtils.changeResidueType(mp.m, mutRes, "MET", true);
		MutUtils.changeResidueType(mp.m, mutRes, "SER", true);
		MutUtils.changeResidueType(mp.m, mutRes, "THR", true);
		MutUtils.changeResidueType(mp.m, mutRes, "LYS", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ARG", true);
		MutUtils.changeResidueType(mp.m, mutRes, "HIP", true);
		MutUtils.changeResidueType(mp.m, mutRes, "HID", true);
		MutUtils.changeResidueType(mp.m, mutRes, "HIE", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ASP", true);
		MutUtils.changeResidueType(mp.m, mutRes, "GLU", true);
		MutUtils.changeResidueType(mp.m, mutRes, "ASN", true);
		MutUtils.changeResidueType(mp.m, mutRes, "GLN", true);
		
		Atom[] curBBAtoms = MutUtils.getBBatoms(mp.m.residue[mutRes]);

		double Hdist = prevBBAtoms[MutUtils.H].distance(curBBAtoms[MutUtils.H]);
		double Ndist = prevBBAtoms[MutUtils.N].distance(curBBAtoms[MutUtils.N]);
		double CAdist = prevBBAtoms[MutUtils.CA].distance(curBBAtoms[MutUtils.CA]);
		double HAdist = prevBBAtoms[MutUtils.HA].distance(curBBAtoms[MutUtils.HA]);
		double Cdist = prevBBAtoms[MutUtils.C].distance(curBBAtoms[MutUtils.C]);
		double Odist = prevBBAtoms[MutUtils.O].distance(curBBAtoms[MutUtils.O]);
		double CBdist = prevBBAtoms[MutUtils.CB].distance(curBBAtoms[MutUtils.CB]);
		
		assertEquals("H is recovered.",0,Hdist,tol);
		assertEquals("N is recovered.",0,Ndist,tol);
		assertEquals("CA is recovered.",0,CAdist,tol);
		assertEquals("HA is recovered.",0,HAdist,tol);
		assertEquals("C is recovered.",0,Cdist,tol);
		assertEquals("O is recovered.",0,Odist,tol);
		assertEquals("CB is recovered.",0,CBdist,tol);
		
		
	}
	
	/*
	 * Test whether the native backbone coordinates can be recovered after 
	 * rotamer is changed to the WT rotamer and back to the original rotamer
	 */
	@Test
	public void wtRotamerMutation() {
		int mutRes = mp.m.mapPDBresNumToMolResNum("41"); //pdb number 41 (WT is Glu)
		
		
		
		AARotamerType glu = mp.m.aaRotLib.getAAType("GLU");
		
		ResidueConformation wtRC = null;
		ResidueConformation testRC = null;
		for(ResidueConformation rc: mp.m.strand[0].rcl.getRCsPosType(mutRes, "GLU")){
			if(rc.rot.isWTrot){
				wtRC = rc;
			}else if(rc.rot.aaIndex == 0)
				testRC = rc;
		}
		
		MutUtils.changeResidueType(mp.m, mutRes, "GLU", true);
		MutUtils.applyRC(mp.m, mp.m.residue[mutRes], testRC);
		
		HashMap<String,Atom> prevAtoms = MutUtils.copyAllAtomPos(mp.m.residue[mutRes]);
		
		MutUtils.applyRC(mp.m, mp.m.residue[mutRes],wtRC); //Apply WT Rotamer
		MutUtils.applyRC(mp.m, mp.m.residue[mutRes],testRC);
		
		HashMap<String,Atom> curAtoms = MutUtils.copyAllAtomPos(mp.m.residue[mutRes]);
		
		
		for(Atom a:mp.m.residue[mutRes].atom){
			assertEquals(a.name+" is recovered.",0,prevAtoms.get(a.name).distance(curAtoms.get(a.name)),tol);
		}
		
		
	}

}
