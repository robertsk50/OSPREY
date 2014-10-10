import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.math.BigInteger;
import java.security.AllPermission;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.TreeSet;


public class Emat implements Serializable {


	//KER: Since objects take up lots of space I'm just going
	//KER: to have lots of matrices that all have primitives
	PairMats pairs;
	SingleMats singles;

	HashMap<String,double[]> eRef = null;
	
	double templ_E = Double.POSITIVE_INFINITY;

	ArrayList<ArrayList<Integer>> resByPos;
	/**
	 * Have the reference energies been added to the energy matrix
	 */
	boolean hasEref = false;
	boolean hasEntropy = false;

	public Emat(){}
	
	public Emat(Molecule m,MutableResParams strandMut, boolean intraOnly,boolean doDih) {
		//templTemplE = new EMatrixEntry();
		resByPos = new ArrayList<ArrayList<Integer>>();
		initializePairEMatrix(m,strandMut,-1,-1,intraOnly,doDih);
	}

	public Emat(Molecule m,MutableResParams strandMut, int pos1, int pos2, boolean doDih) {
		//templTemplE = new EMatrixEntry();
		resByPos = new ArrayList<ArrayList<Integer>>();
		initializePairEMatrix(m,strandMut,pos1,pos2,false,doDih);
	}

	public Emat(Emat emat, int mutPos) {
		templ_E = emat.templ_E;
		resByPos = emat.resByPos;

		//intraE = emat.intraE;//new EMatrixEntry[emat.intraE.length][][];
		singles = emat.singles;

		if(emat.pairs != null && emat.pairs.E!=null && mutPos >= 0){
			pairs = emat.pairs.getSinglePosMat(mutPos);
		}

	}

	public Emat(Emat emat, boolean intraOnly, EmatCalcParams runParams) {
		templ_E = emat.templ_E;
		resByPos = emat.resByPos;

		//intraE = emat.intraE;//new EMatrixEntry[emat.intraE.length][][];
		singles = emat.singles;

		if(!intraOnly && emat.pairs != null && emat.pairs.E!=null){
			pairs = emat.pairs.getOnlyDualPosMat(runParams.pos1,runParams.pos2);
		}

	}
	
	public static Emat dualPosMat(Emat emat, boolean intraOnly, int pos1, int pos2){
		Emat newEmat = new Emat();
		newEmat.templ_E = emat.templ_E;
		newEmat.resByPos = emat.resByPos;

		//intraE = emat.intraE;//new EMatrixEntry[emat.intraE.length][][];
		newEmat.singles = emat.singles;

		if(!intraOnly && emat.pairs != null && emat.pairs.E!=null){
			newEmat.pairs = emat.pairs.getDualPosMat(pos1,pos2);
		}
		return newEmat;
	}
	
	public Emat(Emat emat, boolean intraOnly, int[] mutPos) {
		templ_E = emat.templ_E;
		resByPos = emat.resByPos;

		//intraE = emat.intraE;//new EMatrixEntry[emat.intraE.length][][];
		singles = emat.singles;

		if(!intraOnly && emat.pairs != null && emat.pairs.E!=null){
			pairs = emat.pairs.getOnlyDualPosMat(mutPos[0],mutPos[1]);
		}

	}
	


	public Emat(Emat emat,Molecule m,boolean doDih ){
		templ_E = emat.templ_E;
		fullEmat(emat,m,doDih);
	}

	public Emat(String fileName, boolean doDih, Molecule m){
		open(fileName,doDih, m);
	}

	public Emat(double templ_E, ArrayList<ArrayList<Integer>> resByPos,SingleMats singles, PairMats pairs){
		this.templ_E = templ_E;
		try {
			this.resByPos = (ArrayList<ArrayList<Integer>>) MPItoThread.deepCopy(resByPos);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.singles = singles;
		this.pairs = pairs;
	}

	public ArrayList<Integer> allMutRes(){

		ArrayList<Integer> mutRes = new ArrayList<Integer>();
		Iterator<ArrayList<Integer>> iter = resByPos.iterator();
		while(iter.hasNext()){
			ArrayList<Integer> resAtPos = iter.next();
			Iterator<Integer> iter2 = resAtPos.iterator();
			while(iter2.hasNext()){
				mutRes.add(iter2.next());
			}
		}

		return mutRes;
	}
	
	public static HashMap<String,double[]> loadErefMatrix(String fileName) {
		HashMap<String,double[]> curEref = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
			curEref = (HashMap<String,double[]>) in.readObject();
			in.close();
		}
		catch(Exception e){
			System.out.println("Couldn't find Eref");
		}
		return curEref;
	}

	public void convertToRigid(){
		singles.doDih = false;
		singles.maxE = null;
		singles.rotDih = null;

		pairs.doDih = false;
		pairs.maxE = null;
		pairs.rotDih1 = null;
		pairs.rotDih2 = null;
	}

	/*public void setPairE(EMatrixEntrySlim eme){
		pairs.setE(eme);
	}*/

	public void setSingleE(EMatrixEntrySlim eme){
		singles.setE(eme);
	}

	/*public void setPairE(int[] i, double E){
		pairs.setE(i,E);
	}*/

	public void setE(int[] i, double E){
		if(i.length == 3)
			singles.setE(i,E);
		else if(i.length == 6)
			pairs.setE(i,E);
		else
			System.out.println("Index length not valid");
	}

	public void addE(int[] i, double E){
		if(i.length == 3)
			singles.addE(i,E);
		else if(i.length == 6)
			pairs.addE(i,E);
		else
			System.out.println("Index length not valid");
	}

	public void setSingleE(int[] i, double E){
		singles.setE(i,E);
	}

	/*public void setPairMaxE(EMatrixEntrySlim eme){
		//pairs.setE(eme);
		System.out.println("Code does not support max E");
	}

	public void setSingleMaxE(EMatrixEntrySlim eme){
		//singles.setE(eme);
		System.out.println("Code does not support max E");
	}

	public void setPairMaxE(int[] i, double E){
		//pairs.setE(eme);
		System.out.println("Code does not support max E");
	}

	public void setMaxE(int[] i, double E){
		//pairs.setE(eme);
		System.out.println("Code does not support max E");
	}

	public void setSingleMaxE(int[] i, double E){
		//singles.setE(eme);
		System.out.println("Code does not support max E");
	}*/

	public void setSinglePruned(int[] i, boolean pruned){
		singles.pruned[i[0]][i[1]][i[2]] = pruned;
	}

	public void setSinglePruned(Index3 i, boolean pruned){
		singles.pruned[i.pos][i.aa][i.rot] = pruned;
	}
	
	public void setSinglePruned(int pos, int aa, int rot, boolean pruned){
		singles.pruned[pos][aa][rot] = pruned;
	}

	public double getSingleMinE(int[] i){
		return singles.E[i[0]][i[1]][i[2]];
	}

	public double getSingleMinE(Index3 i){
		return singles.E[i.pos][i.aa][i.rot];
	}
	
	public double getSingleMinE(int pos, int aa, int rot){
		return singles.E[pos][aa][rot];
	}

	/*public void setPairE(EMatrixEntrySlim eme, int[] i){
		pairs.pairE[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]].set(eme);
		//Symmetric so set the other symmetry

		pairE[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = ((RotamerPairEntry)pairE[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]).swappedCopy();
	}*/

	/*public void setSingle(EMatrixEntrySlim eme, int[] i){
		intraE[i[0]][i[1]][i[2]].set(eme);
	}*/



	/*public void initializePairEMatrix(MutableResParams strandMut){
		intraE = new ArrayList<ArrayList<EMatrixEntry>>();
		pairE = new ArrayList<ArrayList<ArrayList<EMatrixEntry>>>();
		//Initialize Array
		for(int i=0; i<strandMut.numMutPos();i++){
			pairE.add(new ArrayList<ArrayList<EMatrixEntry>>());
			intraE.add(new ArrayList<EMatrixEntry>());
			for(int j=0; j<strandMut.numMutPos();j++){
				pairE.get(i).add(new ArrayList<EMatrixEntry>());
			}
		}


		for(int i=0; i<strandMut.allMut.length;i++){
			ArrayList<Residue> tmpRes = new ArrayList<Residue>();
			tmpRes.add(strandMut.allMut[i]);
			resByPos.add(tmpRes);
			for(Rotamer r1: strandMut.allMut[i].getRotamers()){
				SuperRotamer sr1 = new SuperRotamer(r1);
				intraE.get(i).add(new RotamerEntry(i,sr1));
				for(int j=i; j<strandMut.numMutPos();j++){
					for(Rotamer r2: strandMut.allMut[i].getRotamers()){
						RotamerPairEntry rp = new RotamerPairEntry(i,sr1,j,new SuperRotamer(r2));
						pairE.get(i).get(j).add(rp);
						pairE.get(j).get(i).add(rp);
					}

				}

			}
		}

	}*/

	/**
	 * Unprunes all the rotamers that are currently pruned from the matrix
	 */
	public void unPrune(){

		for(int p=0;p<singles.pruned.length;p++){
			for(int a=0;a<singles.pruned[p].length;a++){
				for(int r=0;r<singles.pruned[p][a].length;r++){
					singles.pruned[p][a][r] = false;
					for(int p2=0;p2<pairs.pruned[p][a][r].length;p2++){
						if(p != p2 && areNeighbors(p, p2)){
							for(int a2=0;a2<pairs.pruned[p][a][r][p2].length;a2++){
								for(int r2=0;r2<pairs.pruned[p][a][r][p2][a2].length;r2++){
									pairs.pruned[p][a][r][p2][a2][r2] = false;
								}
							}
						}
					}
				}
			}
		}

	}


	//KER: Wrapper functions to change way Emat is written and read
	public void open(String fileName, boolean doDih, Molecule m){
		//readFromFile(fileName);
		read(fileName,doDih,m);
	}

	public void save(String fileName, Molecule m){
		//writeToFile(fileName);
		write(fileName, m);
	}


	public void write(String fileName,Molecule m){

		//KER: make pdbs directory if it doesn't exist
		File tmpFile = new File(fileName);
		File ematDir = tmpFile.getParentFile();
		if(ematDir != null && !ematDir.exists()){
			ematDir.mkdir();
		}

		pairs.write(fileName);
		singles.write(fileName);

		KSParser.outputObject(templ_E, fileName+".templE");

		KSParser.outputObject(resByPos,fileName+".resByPos");
		/**
		 * Have the reference energies been added to the energy matrix
		 */
		KSParser.outputObject(hasEref,fileName+".hasEref");
		KSParser.outputObject(hasEntropy,fileName+".hasEntropy");

		//The energy matrix is specific rotamer library so save that rotamer library with the emat
		m.aaRotLib.save(fileName+".aaRots");
		if(m.genRotLib != null)
		m.genRotLib.save(fileName+".genRots");
		
		//Save the residue conformation libraries;
		for(Strand s: m.strand)
			s.rcl.save(fileName+".rcl_"+s.number);

		//Save the reference energies
		KSParser.outputObject(eRef, fileName+".eref");
		
	}

	public void read(String file,boolean doDih, Molecule m){


		pairs = PairMats.read(file,doDih);
		singles = SingleMats.read(file,doDih);

		try{
			//			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file+".templE"));
			templ_E = (Double)KSParser.loadObject(file+".templE");//in.readObject();
			//			in.close();
			//			in = new ObjectInputStream(new FileInputStream(file+".resByPos"));
			resByPos = (ArrayList<ArrayList<Integer>>)KSParser.loadObject(file+".resByPos");//in.readObject();
			//			in.close();
			//			in = new ObjectInputStream(new FileInputStream(file+".hasEref"));
			hasEref = (Boolean)KSParser.loadObject(file+".hasEref");//in.readObject();
			//			in.close();
		}
		catch (Exception e){
			//e.printStackTrace();
			System.out.println("Couldn't find resByPos for Ematrix");
			//			return aaRotLib;
		}	

		try{
			//			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file+".hasEntropy"));
			hasEntropy = (Boolean)KSParser.loadObject(file+".hasEntropy");//in.readObject();
			//			in.close();
		}
		catch (Exception e){
			hasEntropy = false;
		}

		//Load the amino acid rotamer library
		try{
			m.aaRotLib.loadGlobalRots(file+".aaRots");
		}
		catch (Exception e){
			System.out.println("Couldn't properly load amino acid rot library for Emat.");
		}
		
		//Load the general rotamer library 
		try{
			m.genRotLib.loadGlobalRots(file+".genRots");
		}
		catch (Exception e){
			System.out.println("Couldn't properly load rot library for Emat.");
		}

		//Load the strand dependent residue conformation libraries
		for(Strand s : m.strand){
			s.rcl.loadGlobalRCs(file+".rcl_"+s.number,m.rotLibForStrand(s.number));
		}
		
		//Read Eref File
		eRef = loadErefMatrix(file+".eref");

	}

	public void writeToFile(String fileName){
		try {
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			PrintStream logPS = new PrintStream( bufferedOutputStream );

			logPS.println("resByPos "+resByPos.size());
			Iterator<ArrayList<Integer>> iter = resByPos.iterator();
			while(iter.hasNext()){
				ArrayList<Integer> positions = iter.next();
				for(int pos:positions){
					logPS.print(pos+" ");
				}
				logPS.println("");
			}
			logPS.flush();
			logPS.println("intraE");
			//Print out all the dimensions
			logPS.println(singles.E.length);
			for(int i=0; i<singles.E.length;i++){
				logPS.print(singles.E[i].length+" ");
				for(int j=0; j<singles.E[i].length;j++){
					logPS.print(singles.E[i][j].length+" ");
				}
				logPS.println("");
			}
			SinglesIterator rotiter = singlesIterator();
			while(rotiter.hasNext()){
				EMatrixEntryWIndex emeWI = rotiter.next();
				for(int i:emeWI.index){
					logPS.print(i+" ");
				}
				logPS.println(emeWI.eme.getString());
			}
			logPS.flush();
			logPS.println("pairE");
			PairsIterator pairiter = pairsIterator();
			while(pairiter.hasNext()){
				EMatrixEntryWIndex emeWI = pairiter.next();
				//if(emeWI.eme.minE() < Double.POSITIVE_INFINITY){
				for(int i:emeWI.index){
					logPS.print(i+" ");
				}
				logPS.println(emeWI.eme.getString());
				//}
			}
			logPS.flush();
			logPS.println("templTemplE");
			logPS.println(templ_E);
			logPS.flush();
			logPS.close();
		}
		catch (Exception ex) {
			ex.printStackTrace();
			System.out.println("ERROR: An exception occured while writing emat file");
		}
	}

	public void readFromFile(String fileName,boolean doDih){
		try{	
			FileInputStream is = new FileInputStream( fileName );
			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;

			//ResByPos
			resByPos = new ArrayList<ArrayList<Integer>>();
			curLine = bufread.readLine(); //First line is 
			KSParser.getToken(curLine, 2);
			curLine = bufread.readLine();
			while(!curLine.startsWith("intraE")){
				StringTokenizer st = new StringTokenizer(curLine," ,;\t\n\r\f");
				ArrayList<Integer> residues = new ArrayList<Integer>();
				while(st.hasMoreTokens()){
					int res = Integer.parseInt(st.nextToken());
					residues.add(res);
				}
				resByPos.add(residues);
				curLine = bufread.readLine();
			}

			//Read in dimensions for matrices
			singles = new SingleMats(doDih);
			pairs = new PairMats(doDih);
			curLine = bufread.readLine();
			int firstDim = Integer.parseInt(curLine);
			int[] p1arr = new int[0];
			singles.addDim(p1arr, firstDim);

			for(int pos1 = 0; pos1 < firstDim;pos1++){
				curLine = bufread.readLine();
				StringTokenizer st = new StringTokenizer(curLine," ,;\t\n\r\f");
				int[] p1a1 = {pos1};
				int secondDim = Integer.parseInt(st.nextToken());
				singles.addDim(p1a1, secondDim);
				for(int aa1 = 0; aa1 < secondDim;aa1++){
					int[] p1a1r1 = {pos1,aa1};
					int thirdDim = Integer.parseInt(st.nextToken());
					singles.addDim(p1a1r1, thirdDim);
				}
			}

			//Set the rest of the dimensions for the pairs matrices
			pairs.addDim(p1arr,singles.E.length);
			for(int p1=0; p1<singles.E.length;p1++){
				int[] p1a1 = {p1};
				pairs.addDim(p1a1,singles.E[p1].length);
				for(int a1=0; a1<singles.E[p1].length;a1++){
					int[] p1a1r1 = {p1,a1};
					pairs.addDim(p1a1r1,singles.E[p1][a1].length);
					for(int r1=0; r1<singles.E[p1][a1].length;r1++){
						int[] p1a1r1p2 = {p1,a1,r1};
						pairs.addDim(p1a1r1p2,singles.E.length);
						for(int p2=0; p2<singles.E.length;p2++){
							if(p1 == p2)
								continue;
							int[] p1a1r1p2a2 = {p1,a1,r1,p2};
							pairs.addDim(p1a1r1p2a2,singles.E[p2].length);
							for(int a2=0; a2<singles.E[p2].length;a2++){
								int[] p1a1r1p2a2r2 = {p1,a1,r1,p2,a2};
								pairs.addDim(p1a1r1p2a2r2,singles.E[p2][a2].length);	
							}	
						}
					}
				}
			}


			curLine = bufread.readLine();
			while(!curLine.startsWith("pairE")){
				StringTokenizer st = new StringTokenizer(curLine," ,;\t\n\r\f");
				int i1 = Integer.parseInt(st.nextToken());
				int i2 = Integer.parseInt(st.nextToken());
				int i3 = Integer.parseInt(st.nextToken());
				double minE = Double.parseDouble(st.nextToken());
				//double maxE = Double.parseDouble(st.nextToken());
				boolean pruned = Boolean.parseBoolean(st.nextToken());
				//boolean prunedIsSteric = Boolean.parseBoolean(st.nextToken());

				int[] rotamers = new int[resByPos.get(i1).size()];
				for(int i=0; i<rotamers.length;i++){
					rotamers[i] = new Integer(st.nextToken());
				}

				SuperRotamer sr = new SuperRotamer(rotamers);
				//sr.rotamers = rotamers;

				//intraE[i1][i2][i3] = new RotamerEntry(i1,sr,minE,maxE,pruned,prunedIsSteric);
				singles.setAll(i1, i2, i3, minE, pruned, rotamers);
				curLine = bufread.readLine();
			}
			curLine = bufread.readLine();
			while(!curLine.startsWith("templTemplE")){
				StringTokenizer st = new StringTokenizer(curLine," ,;\t\n\r\f");
				int i1 = Integer.parseInt(st.nextToken());
				int i2 = Integer.parseInt(st.nextToken());
				int i3 = Integer.parseInt(st.nextToken());
				int i4 = Integer.parseInt(st.nextToken());
				int i5 = Integer.parseInt(st.nextToken());
				int i6 = Integer.parseInt(st.nextToken());
				double minE = Double.parseDouble(st.nextToken());
				//double maxE = Double.parseDouble(st.nextToken());
				boolean pruned = Boolean.parseBoolean(st.nextToken());
				//boolean prunedIsSteric = Boolean.parseBoolean(st.nextToken());

				/*int[] rots1 = new int[resByPos.get(i1).size()];
				for(int i=0; i<rots1.length;i++){
					rots1[i] = Integer.parseInt(st.nextToken());
				}
				String atSymb = st.nextToken();
				int[] rots2 = new int[resByPos.get(i4).size()];
				for(int i=0; i<rots2.length;i++){
					rots2[i] = Integer.parseInt(st.nextToken());
				}*/

				pairs.setAll(i1,i2,i3,i4,i5,i6,minE,pruned);

				curLine = bufread.readLine();
			}
			curLine = bufread.readLine();
			//StringTokenizer st = new StringTokenizer(curLine," ,;\t\n\r\f");
			double minE = Double.parseDouble(curLine);
			//double maxE = Double.parseDouble(st.nextToken());
			//boolean pruned = Boolean.parseBoolean(st.nextToken());
			//boolean prunedIsSteric = Boolean.parseBoolean(st.nextToken());
			templ_E = minE; //new EMatrixEntry(minE,maxE,pruned,prunedIsSteric);			
		}
		catch(Exception e){
			e.printStackTrace();
			System.out.println("ERROR: An exception occured while reading emat file.");
		}
	}

	public void fullEmat(Emat emat, Molecule m,boolean doDih){

		this.singles = new SingleMats(doDih);
		this.pairs = new PairMats(doDih);

		this.resByPos = emat.resByPos;
		this.singles.addDim(new int[0], emat.singles.E.length);
		this.pairs.addDim(new int[0],emat.singles.E.length);

		//Runtime runtime = Runtime.getRuntime();
		for(int p1=0; p1<emat.singles.E.length;p1++){
			int[] p1ind = {p1};
			this.singles.addDim(p1ind,emat.singles.E[p1].length);
			this.pairs.addDim(p1ind,emat.singles.E[p1].length);

			for(int a1=0;a1<emat.singles.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				this.singles.addDim(p1a1ind,emat.singles.E[p1][a1].length);
				this.pairs.addDim(p1a1ind,emat.singles.E[p1][a1].length);
				for(int r1=0; r1<emat.singles.E[p1][a1].length;r1++){
					int[] p1a1r1ind = {p1,a1,r1};
					this.singles.copy(p1a1r1ind, emat.singles,p1a1r1ind);
					this.pairs.addDim(p1a1r1ind,emat.singles.E.length);
					for(int p2=0; p2<emat.singles.E.length;p2++){
						if(p1!=p2 && areNeighbors(p1,p2,m)  ){
							int[] p1a1r1p2ind = {p1,a1,r1,p2};
							this.pairs.addDim(p1a1r1p2ind,emat.singles.E[p2].length);
							for(int a2=0;a2<emat.singles.E[p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
								this.pairs.addDim(p1a1r1p2a2ind,emat.singles.E[p2][a2].length);
								/*for(int r2=0; r2<emat.singles.E[p2][a2].length;r2++){
									//int[] p1a1r1p2a2r2ind = {p1,a1,r1,p2,a2,r2};
									//RotamerPairEntry rp = new RotamerPairEntry(p1,((RotamerEntry)emat.singles[p1][a1][r1]).r,p2,((RotamerEntry)emat.singles[p2][a2][r2]).r);
									this.pairs.setSupRot(p1, a1, r1, p2, a2, r2, emat.singles.supRot[p1][a1][r1], emat.singles.supRot[p2][a2][r2]);
								}*/
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Check if two mutable positions are neighbors. They
	 * are considered neighbors if any residues between the two positions
	 * are neighbors.
	 */
	public boolean areNeighbors(int p1, int p2, Molecule m) {

		for(int res1 : resByPos.get(p1))
			for(int res2 : resByPos.get(p2))
				if(m.areNeighbors(res1, res2))
					return true;

		return false;
	}

	public void initializePairEMatrix(Molecule m,MutableResParams strandMut, int pos1, int pos2, boolean intraOnly,boolean doDih){

		this.singles = new SingleMats(doDih);
		this.pairs = new PairMats(doDih);

		this.singles.addDim(new int[0], strandMut.numMutPos());//intraE = new EMatrixEntry[strandMut.numMutPos()][][];
		if(!intraOnly)
			this.pairs.addDim(new int[0],strandMut.numMutPos());//pairE = new EMatrixEntry[strandMut.numMutPos()][][][][][];

		//Runtime runtime = Runtime.getRuntime();

		for(int p1=0; p1<strandMut.allMut.length;p1++){
			int[] p1ind = {p1};
			///long maxMemory = runtime.maxMemory();  
			//long allocatedMemory = runtime.totalMemory();  
			//long freeMemory = runtime.freeMemory();
			//long totalFreeMemory = (freeMemory + (maxMemory - allocatedMemory));
			//System.out.println("Pos: "+p1+" FreeMem: " + totalFreeMemory);
			Residue res1 = m.residue[strandMut.allMut[p1]];
			ArrayList<Integer> tmpRes = new ArrayList<Integer>();
			tmpRes.add(res1.moleculeResidueNumber);
			resByPos.add(tmpRes);
			if(pos1 == -1 || pos1==p1 || pos2 == p1){
				this.singles.addDim(p1ind,res1.numAllowedAATypes());//intraE[p1] = new EMatrixEntry[res1.allowedAATypes.size()][];
				if(!intraOnly)
					this.pairs.addDim(p1ind,res1.numAllowedAATypes());//pairE[p1] = new EMatrixEntry[res1.allowedAATypes.size()][][][][];

				for(int a1=0;a1<res1.numAllowedAATypes();a1++){
					int[] p1a1ind = {p1,a1};
					//AARotamerType aaType1 = res1.getAllowedAAType(a1);
					ArrayList<ResidueConformation> rotsForAAType1 = res1.getRCsForType(a1); 
					this.singles.addDim(p1a1ind,rotsForAAType1.size());//intraE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()];
					if(!intraOnly)
						this.pairs.addDim(p1a1ind,rotsForAAType1.size());//pairE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()][][][];
					for(int r1=0; r1<rotsForAAType1.size();r1++){
						int[] p1a1r1ind = {p1,a1,r1};
						ResidueConformation rc1 = rotsForAAType1.get(r1);
						//SuperRotamer sr1 = new SuperRotamer(rot1.rlIndex);
						int[] rot1Array = {rc1.id};
						this.singles.setSupRot(p1, a1, r1, rot1Array);//intraE[p1][a1][r1] = new RotamerEntry(p1,sr1);
						if(!intraOnly){
							this.pairs.addDim(p1a1r1ind,strandMut.numMutPos());//pairE[p1][a1][r1] = new EMatrixEntry[strandMut.numMutPos()][][];
							for(int p2=0; p2<strandMut.allMut.length;p2++){
								if(p1!=p2){
									int[] p1a1r1p2ind = {p1,a1,r1,p2};
									if( pos2 == -1 || pos2==p2 ){
										Residue res2 = m.residue[strandMut.allMut[p2]];
										if(m.areNeighbors(res1.moleculeResidueNumber,res2.moleculeResidueNumber)){
											this.pairs.addDim(p1a1r1p2ind,res2.numAllowedAATypes());//pairE[p1][a1][r1][p2] = new EMatrixEntry[res2.allowedAATypes.size()][];
											for(int a2=0;a2<res2.numAllowedAATypes();a2++){
												int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
												//AARotamerType aaType2 = res2.allowedAATypes.get(a2);
												ArrayList<ResidueConformation> rotsForAAType2 = res2.getRCsForType(a2);
												this.pairs.addDim(p1a1r1p2a2ind,rotsForAAType2.size());//pairE[p1][a1][r1][p2][a2] = new EMatrixEntry[aaType2.numRotamers()];
												
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

	}


	public int numMutPos() {
		return resByPos.size();
	}

	public SinglesIterator singlesIterator() {
		return new SinglesIterator(this);
	}

	public SinglesIterator singlesIterator(int p1) {
		return new SinglesIterator(this,p1);
	}

	public PairsIterator pairsIterator() {
		return new PairsIterator(this);
	}

	public PairsIterator pairsIterator(int p1, int p2) {
		return new PairsIterator(this,p1,p2);
	}

	public PairsIterator pairsIterator(int p1, int p2, TreeSet<Integer> AAs1, TreeSet<Integer> AAs2){
		return new PairsIterator(this,p1,p2,AAs1,AAs2);
	}

	/*public EMatrixEntry getPairE(int[] index){
		assert index.length == 6;
		return pairs.getE(index);
	}

	public EMatrixEntry getSymmetricPairE(int[] index){
		assert index.length == 6;
		return pairE[index[3]][index[4]][index[5]][index[0]][index[1]][index[2]];
	}

	public EMatrixEntry getPairE(int[] i1, int[] i2){
		assert i1.length == 3;
		return pairE[i1[0]][i1[1]][i1[2]][i2[0]][i2[1]][i2[2]];
	}

	public EMatrixEntry getPairE(Index3 i1, Index3 i2){
		return pairE[i1.pos][i1.aa][i1.rot][i2.pos][i2.aa][i2.rot];
	}*/

	public double getTemplMinE(){
		return templ_E;
	}

	public void setTemplMinE(double E){
		templ_E = E;
	}

	public void setTemplMaxE(double E){
		System.out.println("Code does not support Max energies");
		//System.exit(0);
	}

	public double getTemplMaxE(boolean doIMinDEE){
		if(doIMinDEE){
			return templ_E;
		}
		else{
			System.out.println("Code does not support Max energies");
			//System.exit(0);
		}
		return 0.0f;
	}

	public double getIntraMinE(int[] index){
		assert index.length ==3;
		return singles.E[index[0]][index[1]][index[2]];
	}

	public double getIntraMaxE(int[] index, boolean doIMinDEE){
		if(doIMinDEE){
			assert index.length ==3;
			return singles.E[index[0]][index[1]][index[2]];
		}
		else{
			System.out.println("Code does not support Max energies");
			System.exit(0);
		}
		return 0.0f;
	}

	public boolean getPairPruned(int[] i1, int[] i2){
		assert i1.length == 3;
		return pairs.pruned[i1[0]][i1[1]][i1[2]][i2[0]][i2[1]][i2[2]];
	}

	public boolean getPairPruned(Index3 i1, Index3 i2){
		return pairs.pruned[i1.pos][i1.aa][i1.rot][i2.pos][i2.aa][i2.rot];
	}

	public boolean getPairPruned(int p1, int a1, int r1, int p2, int a2, int r2){
		return pairs.pruned[p1][a1][r1][p2][a2][r2];
	}

	public double getPairMinE(int[] i1, int[] i2){
		assert i1.length == 3;
		return pairs.E[i1[0]][i1[1]][i1[2]][i2[0]][i2[1]][i2[2]];
	}
	
	public double getPairMinE(int[] i1){
		assert i1.length == 6;
		return pairs.E[i1[0]][i1[1]][i1[2]][i1[3]][i1[4]][i1[5]];
	}
	
	public double getPairwiseE(int[] i1){
		assert i1.length == 6;
		return pairs.E[i1[0]][i1[1]][i1[2]][i1[3]][i1[4]][i1[5]];
	}
	
	public double getPairwiseE(int p1, int a1, int r1, int p2, int a2, int r2){
		return pairs.E[p1][a1][r1][p2][a2][r2];
	}

	public double getPairMinE(Index3 i1, Index3 i2){
		return pairs.E[i1.pos][i1.aa][i1.rot][i2.pos][i2.aa][i2.rot];
	}

	public EMatrixEntry getPairTerm(int[] i) {
		return pairs.getTerm(i, singles.supRot);
	}

	public boolean getSinglePruned(int[] i){
		return singles.pruned[i[0]][i[1]][i[2]];
	}

	public boolean getSinglePruned(Index3 i){
		return singles.pruned[i.pos][i.aa][i.rot];
	}

	public boolean getSinglePruned(int pos,int aa, int rot){
		return singles.pruned[pos][aa][rot];
	}

	public double getPairMaxE(int[] i1, int[] i2, boolean doIMinDEE){
		if(doIMinDEE){
			assert i1.length == 3;
			return pairs.E[i1[0]][i1[1]][i1[2]][i2[0]][i2[1]][i2[2]];
		}
		else{
			System.out.println("Code does not support Max energies");
			System.exit(0);
		}
		return 0.0f;
	}

	/*public void setPairPruned(int[] i, boolean pruned) {
		pairs.pruned[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = pruned;
	}*/

	public void setPruned(int[] i, boolean pruned) {
		if(i.length == 3)
			singles.pruned[i[0]][i[1]][i[2]] = pruned;
		else if(i.length == 6)
			pairs.pruned[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = pruned;
		else
			System.out.println("Index not valid!");
	}

	public void setPairPruned(int[] i1, int[] i2, boolean pruned) {
		pairs.pruned[i1[0]][i1[1]][i1[2]][i2[0]][i2[1]][i2[2]] = pruned;
	}
	
	public void setPairPruned(int p1, int a1, int r1, int p2, int a2, int r2, boolean pruned) {
		pairs.pruned[p1][a1][r1][p2][a2][r2] = pruned;
	}
	
	public void setPairPruned(int[] i1, boolean pruned) {
		pairs.pruned[i1[0]][i1[1]][i1[2]][i1[3]][i1[4]][i1[5]] = pruned;
	}
	
	public void setPairwiseE(int[] i1, double E) {
		pairs.E[i1[0]][i1[1]][i1[2]][i1[3]][i1[4]][i1[5]] = E;
	}

	public void setSymmetricPairPruned(int[] i, boolean pruned) {
		pairs.pruned[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = pruned;
	}


	/*public EMatrixEntry getIntraE(Index3 index){
		return intraE[index.pos][index.aa][index.rot];
	}*/
	
	//Marks all rotamer pairs for which the min energy matrix entry is greater than cutoff as having a steric clash;
	//		If the max energy matrix exists, the corresponding entries are marked in the same way
	public void preprocessPairs(double cutoff, double stericE){
			
		PairsIterator iter = pairsIterator();
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(getPairwiseE( emeWI.index ) > cutoff){
				setPairwiseE( emeWI.index, stericE );
			}
		}
			
	}
	
	//Prune rotamer or residue-conformation pairs that are definitely impossible due to a steric clash or to parametric incompatibility
	//as measured by their interaction energy being greater than the cutoff
	//Particularly useful for DEEPer
	public void pruneRidiculousPairs(double cutoff){

		int numPruned = 0;
		PairsIterator iter = pairsIterator();
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(getPairwiseE( emeWI.index ) >= cutoff){
				setPairPruned(emeWI.index, true);
				numPruned++;
			}
		}

		System.out.println("Number of pairs pruned due to steric clashes or parametric incompatibility: "+numPruned);
	}

	public int numRot(int pos){
		int numRot = 0;
		for(int aa = 0; aa<singles.pruned[pos].length;aa++){
			numRot += singles.pruned[pos][aa].length;
		}
		return numRot;
	}

	public int numRot(int pos,int aa){

		return singles.pruned[pos][aa].length;

	}

	/*
	 * Given the flat rotamer val for a pos, return back array with
	 * AAnum at [0] and RotNum at [1]
	 */
	public int[] rot2AArot(int pos, int rotNum){
		int[] retVal = new int[2];
		int tmpRot = 0;
		for(int aa = 0; aa<singles.pruned[pos].length;aa++){
			if(rotNum >= tmpRot+singles.pruned[pos][aa].length)
				tmpRot += singles.pruned[pos][aa].length;
			else{
				retVal[0] = aa;
				retVal[1] = rotNum-tmpRot;
				return retVal;
			}
		}
		retVal[0] = -1;retVal[1] = -1;
		return retVal;
	}


	public int[] remainingRot(){
		int[] rotsRemaining = new int[numMutPos()];
		for(int i=0; i<rotsRemaining.length;i++)
			rotsRemaining[i] = 0;

		for(int pos = 0; pos<singles.pruned.length;pos++){
			for(int aa=0; aa<singles.pruned[pos].length;aa++)
				for(int rot=0; rot<singles.pruned[pos][aa].length;rot++)
					if(!singles.pruned[pos][aa][rot])
						rotsRemaining[pos]++;
		}


		return rotsRemaining;
	}

	//This returns how many rotamers there are per position (both pruned and unpruned)
	public int[] numRotPerPos(){
		int[] rotsRemaining = new int[numMutPos()];
		for(int i=0; i<rotsRemaining.length;i++)
			rotsRemaining[i] = 0;

		for(int pos = 0; pos<singles.pruned.length;pos++){
			for(int aa=0; aa<singles.pruned[pos].length;aa++)
				rotsRemaining[pos] += singles.pruned[pos][aa].length;
		}


		return rotsRemaining;
	}

	//This returns how many rotamers there are per position (both pruned and unpruned)
	public int[] numAANotPruned(){
		int[] AARemaining = new int[numMutPos()];
		for(int i=0; i<AARemaining.length;i++)
			AARemaining[i] = 0;

		for(int pos = 0; pos<singles.pruned.length;pos++){
			for(int aa=0; aa<singles.pruned[pos].length;aa++)
				for(int rot=0; rot<singles.pruned[pos][aa].length;rot++){
					if(!singles.pruned[pos][aa][rot]){
						AARemaining[pos]++;
						break;
					}
				}
		}
		return AARemaining;
	}
	
	public int[][] seqIndicesPerLevel(int[] numSeqForLevel){
		int[][] seqIndicesPerLevel = new int[numMutPos()][];
		for(int i=0; i<seqIndicesPerLevel.length;i++)
			seqIndicesPerLevel[i] = new int[numSeqForLevel[i]];

		for(int pos = 0; pos<singles.pruned.length;pos++){
			int ctr=0;
			for(int aa=0; aa<singles.pruned[pos].length;aa++)
				for(int rot=0; rot<singles.pruned[pos][aa].length;rot++){
					if(!singles.pruned[pos][aa][rot]){
						seqIndicesPerLevel[pos][ctr] = aa;
						ctr++;
						break;
					}
				}
			
		}
		return seqIndicesPerLevel;
	}
	
	//This returns how many rotamers there are per position (only unpruned)
	public int[][] numRotRemainingBySeq(){
		int[][] rotRemaining = new int[numMutPos()][];
		
		for(int pos = 0; pos<singles.pruned.length;pos++){
			ArrayList<Integer> rotPerAA = new ArrayList<Integer>();
			for(int aa=0; aa<singles.pruned[pos].length;aa++){
				int numRot = 0;
				for(int rot=0; rot<singles.pruned[pos][aa].length;rot++){
					if(!singles.pruned[pos][aa][rot]){
						numRot++;
					}
				}
				if(numRot > 0)
					rotPerAA.add(numRot);
			}
			rotRemaining[pos] = new int[rotPerAA.size()];
			for(int i=0; i<rotRemaining[pos].length;i++)
				rotRemaining[pos][i] = rotPerAA.get(i);
		}
		return rotRemaining;
	}
	
	/**
	 * Since the energy matrix now contains objects instead of indexes
	 * when it is loaded it needs to be reconnected to the current molecule
	 * and the current rotamerLibrary  
	 */
	/*public void reconnect(Molecule m) {
		for(int i=0; i<resByPos.size();i++){
			for(int j=0; j<resByPos.get(j).size();j++){
				Residue curRes = resByPos.get(i).get(j);
				resByPos.get(i).set(j, m.residue[curRes.moleculeResidueNumber]);
			}
		}

	}*/

	public void removePrunedRotMoreMem(boolean intraOnly, Emat emat){

		SingleMats newSingles = new SingleMats(singles.doDih);
		PairMats newPair = new PairMats(pairs.doDih);

		newSingles.addDim(new int[0], singles.E.length);//EMatrixEntry[][][] newIntra = new EMatrixEntry[intraE.length][][];
		if(!intraOnly)
			newPair.addDim(new int[0], pairs.E.length);//newPair = new EMatrixEntry[pairE.length][][][][][];

		for(int p1=0;p1<singles.E.length;p1++){
			int[] p1ind = {p1};
			//else we just copy over the array from the Emat array
			newSingles.addDim(p1ind,singles.E[p1].length);//newIntra[p1] = new EMatrixEntry[intraE[p1].length][];
			if(!intraOnly)
				newPair.addDim(p1ind,pairs.E[p1].length);//newPair[p1] = new EMatrixEntry[pairE[p1].length][][][][];

			for(int a1=0;a1<singles.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				//Else just copy over the array
				newSingles.addDim(p1a1ind, getNumRotUnpruned(p1,a1));//newIntra[p1][a1] = new EMatrixEntry[getNumRotUnpruned(p1,a1)];
				if(!intraOnly)
					newPair.addDim(p1a1ind, getNumRotUnpruned(p1,a1));//newPair[p1][a1] = new EMatrixEntry[getNumRotUnpruned(p1,a1)][][][];

				int offset1 = 0;
				for(int r1ctr = 0; r1ctr<singles.E[p1][a1].length;r1ctr++){
					int r1 = r1ctr+offset1;
					int[] p1a1r1ind = {p1,a1,r1};
					int[] p1a1r1ctrind = {p1,a1,r1ctr};
					//We are getting rid of position mutPos2 so make the length 1 less
					if(singles.pruned[p1][a1][r1ctr]){
						offset1--;
					}
					else{
						newSingles.copy(p1a1r1ind, singles,p1a1r1ctrind); //newIntra[p1][a1][r1] = intraE[p1][a1][r1ctr];
						if(!intraOnly){
							newPair.addDim(p1a1r1ind, pairs.E[p1][a1][r1ctr].length);//newPair[p1][a1][r1] = new EMatrixEntry[pairE[p1][a1][r1ctr].length][][];

							for(int p2 = 0; p2<pairs.E[p1][a1][r1].length;p2++){
								int[] p1a1r1p2ind = {p1,a1,r1,p2};
								//int[] p1a1r1ctrp2ind = {p1,a1,r1ctr,p2};

								if(p1 == p2 || !areNeighbors(p1, p2)){
									//do nothing here;
									continue;
								}
								else{
									newPair.addDim(p1a1r1p2ind, pairs.E[p1][a1][r1ctr][p2].length);//newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1][a1][r1ctr][p2].length][];
								}

								for(int a2 = 0; a2<pairs.E[p1][a1][r1][p2].length;a2++){
									int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
									newPair.addDim(p1a1r1p2a2ind, getNumRotUnpruned(p2,a2));//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[getNumRotUnpruned(p2,a2)];

									int offset2=0;
									for(int r2ctr=0;r2ctr<pairs.E[p1][a1][r1][p2][a2].length;r2ctr++){

										int r2 = r2ctr + offset2;
										int[] p1a1r1p2a2r2ind = {p1,a1,r1,p2,a2,r2};
										int[] p1a1r1ctrp2a2r2ctrind = {p1,a1,r1ctr,p2,a2,r2ctr};
										if(singles.pruned[p2][a2][r2ctr]){
											offset2--;
										}
										else{
											newPair.copy(p1a1r1p2a2r2ind, pairs, p1a1r1ctrp2a2r2ctrind);//newPair[p1][a1][r1][p2][a2][r2] = pairE[p1][a1][r1ctr][p2][a2][r2ctr];
										}
									}

								}								

							}
						}
					}

				}
			}

		}

		emat.pairs = newPair;
		emat.singles = newSingles;

		//Contract eliminatedRotAtRes
		//rotamerIndex.contract(eliminatedRotAtRes.getPrunedArray());
		//eliminatedRotAtRes.contract(eliminatedRotAtRes.getPrunedArray());
	}


	public Emat unprunedMatrix(){
		Emat unprunedMat = new Emat();
		unprunedMat.resByPos = resByPos;
		unprunedMat.templ_E = templ_E;
		
		removePrunedRotMoreMem(false,unprunedMat);
		
		return unprunedMat;
	}


	/*
	 * This function tries to remove excess rotamers while not "remaking" the entire
	 * pairs matrix but rather "shrinking" the existing matrix
	 */
	public void removePrunedRotReducedMem(boolean intraOnly){

		SingleMats newSingles = new SingleMats(singles.doDih);

		//First make the new singles matrix 
		newSingles.addDim(new int[0], singles.E.length);//EMatrixEntry[][][] newIntra = new EMatrixEntry[intraE.length][][];

		for(int p1=0;p1<singles.E.length;p1++){
			int[] p1ind = {p1};
			//else we just copy over the array from the Emat array
			newSingles.addDim(p1ind,singles.E[p1].length);//newIntra[p1] = new EMatrixEntry[intraE[p1].length][];
			for(int a1=0;a1<singles.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				//Else just copy over the array
				newSingles.addDim(p1a1ind, getNumRotUnpruned(p1,a1));//newIntra[p1][a1] = new EMatrixEntry[getNumRotUnpruned(p1,a1)];

				int offset1 = 0;
				for(int r1ctr = 0; r1ctr<singles.E[p1][a1].length;r1ctr++){
					int r1 = r1ctr+offset1;
					int[] p1a1r1ind = {p1,a1,r1};
					int[] p1a1r1ctrind = {p1,a1,r1ctr};
					//We are getting rid of position mutPos2 so make the length 1 less
					if(singles.pruned[p1][a1][r1ctr]){
						offset1--;
					}
					else{
						newSingles.copy(p1a1r1ind, singles,p1a1r1ctrind); //newIntra[p1][a1][r1] = intraE[p1][a1][r1ctr];

					}
				}

			}
		}


		//Now shrink all of the dim=6 matrices
		if(!intraOnly){

			for(int p1=0;p1<singles.E.length;p1++){
				for(int a1=0;a1<singles.E[p1].length;a1++){
					for(int r1 = 0; r1<singles.E[p1][a1].length;r1++){
						for(int p2 = 0; p2<pairs.E[p1][a1][r1].length;p2++){
							if(p1 == p2 || !areNeighbors(p1, p2)){
								//do nothing here;
								continue;
							}

							for(int a2 = 0; a2<pairs.E[p1][a1][r1][p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};

								boolean[] tmpPruned = new boolean[pairs.E[p1][a1][r1][p2][a2].length];
								for(int r2ctr=0;r2ctr<pairs.E[p1][a1][r1][p2][a2].length;r2ctr++){

									try{
									if(singles.pruned[p2][a2][r2ctr]){
										tmpPruned[r2ctr] = true;
									}
									else{
										tmpPruned[r2ctr] = false;
									}}catch(Exception E){
										E.printStackTrace();
										System.out.println("DELETE ME");
									}

								}
								pairs.shrinkDim6(p1a1r1p2a2ind,tmpPruned);

							}
						}

					}								

				}
			}
		}

		//Now shrink all of the dim=3 matrices
		if(!intraOnly){

			for(int p1=0;p1<singles.E.length;p1++){
				for(int a1=0;a1<singles.E[p1].length;a1++){
					int[] p1a1ind = {p1,a1};
					boolean[] tmpPruned = new boolean[pairs.E[p1][a1].length];
					for(int r1 = 0; r1<singles.E[p1][a1].length;r1++){
						if(singles.pruned[p1][a1][r1])
							tmpPruned[r1] = true;
						else
							tmpPruned[r1] = false;
					}
					pairs.shrinkDim3(p1a1ind,tmpPruned);

				}
			}
		}


		singles = newSingles;

		//Contract eliminatedRotAtRes
		//rotamerIndex.contract(eliminatedRotAtRes.getPrunedArray());
		//eliminatedRotAtRes.contract(eliminatedRotAtRes.getPrunedArray());
	}


	private int getNumRotUnpruned(int p1, int a1) {
		int numUnpruned = 0;
		//RotInfo<Boolean> ri = new RotInfo<Boolean>(p1,a1,0,true);
		for(int r1=0; r1<singles.pruned[p1][a1].length;r1++){
			if (!singles.pruned[p1][a1][r1]){ //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
				numUnpruned++;
			}
		}
		return numUnpruned;
	}


	public void combinePos(int mutPos1, int mutPos2){

		SingleMats newSingles = new SingleMats(singles.doDih);
		PairMats newPair = new PairMats(pairs.doDih);
		System.out.println("Combining pos: "+mutPos1+" and "+mutPos2);

		newSingles.addDim(new int[0], singles.E.length-1);//EMatrixEntry[][][] newIntra = new EMatrixEntry[intraE.length-1][][];
		newPair.addDim(new int[0], pairs.E.length-1);//EMatrixEntry[][][][][][] newPair = new EMatrixEntry[pairE.length-1][][][][][];

		int offSet1 = 0;
		for(int p1old=0;p1old<pairs.E.length;p1old++){
			int p1 = p1old+offSet1;
			int[] p1ind = {p1};
			int[] maxAAs = new int[2];maxAAs[0] = singles.E[mutPos1].length; maxAAs[1] = singles.E[mutPos2].length;
			if(p1old == mutPos1){
				//If we are at the position we're merging we need to multiply the numAAs from mutPos1 and mutPos2
				int numAAs = getNumCombinedAA(mutPos1, mutPos2, maxAAs);
				newSingles.addDim(p1ind, numAAs);//newIntra[p1] = new EMatrixEntry[numAAs][];
				newPair.addDim(p1ind, numAAs);//newPair[p1] = new EMatrixEntry[numAAs][][][][];
			}
			else if(p1old == mutPos2){
				//If we are at position mutPos2 we need to skip this position
				offSet1--;
				continue;
			}
			else{
				//else we just copy over the array from the Emat array
				newSingles.addDim(p1ind, singles.E[p1old].length);//newIntra[p1] = new EMatrixEntry[intraE[p1old].length][];
				newPair.addDim(p1ind, pairs.E[p1old].length);//newPair[p1] = new EMatrixEntry[pairE[p1old].length][][][][];
			}

			int[] curAAs = new int[2];curAAs[0] = 0; curAAs[1] = -1;
			for(int a1=0;a1<newPair.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				int[] maxRots = {0,0};
				int[] curRots = {0,0};
				if(p1 == mutPos1){
					//If we are at the position we're merging multiple numRot for the aas at this index
					//KER: We only want to copy over the non-pruned ones
					do{
						Emat.incrementCtr(curAAs, maxAAs);
					}while(getCombinedRotUnpruned(mutPos1,mutPos2,curAAs)==0);
					int numRot = getCombinedRotUnpruned(mutPos1,mutPos2,curAAs);
					newSingles.addDim(p1a1ind, numRot);//newIntra[p1][a1] = new EMatrixEntry[numRot];//new EMatrixEntry[intraE[mutPos1][curAAs[0]].length * intraE[mutPos2][curAAs[1]].length];
					newPair.addDim(p1a1ind, numRot);//newPair[p1][a1] = new EMatrixEntry[numRot][][][];//new EMatrixEntry[pairE[mutPos1][curAAs[0]].length * pairE[mutPos2][curAAs[1]].length][][][];
					maxRots = new int[2];maxRots[0] = singles.E[mutPos1][curAAs[0]].length;maxRots[1] = singles.E[mutPos2][curAAs[1]].length;
					curRots = new int[2];curRots[0] = 0; curRots[1] = -1;

				}
				else{
					//Else just copy over the array
					newSingles.addDim(p1a1ind, singles.E[p1old][a1].length);//newIntra[p1][a1] = new EMatrixEntry[intraE[p1old][a1].length];
					newPair.addDim(p1a1ind, pairs.E[p1old][a1].length);//newPair[p1][a1] = new EMatrixEntry[pairE[p1old][a1].length][][][];
				}

				for(int r1 = 0; r1<newPair.E[p1][a1].length;r1++){
					int[] p1a1r1ind = {p1,a1,r1};
					if(p1old == mutPos1){ //Combine the intra terms
						do{
							Emat.incrementCtr(curRots, maxRots);
						} while(singles.pruned[mutPos1][curAAs[0]][curRots[0]] || singles.pruned[mutPos2][curAAs[1]][curRots[1]] || pairs.pruned[mutPos1][curAAs[0]][curRots[0]][mutPos2][curAAs[1]][curRots[1]]);
						RotamerPairEntry pairTerm = pairs.getTerm(mutPos1,curAAs[0],curRots[0],mutPos2,curAAs[1],curRots[1], singles.supRot);
						RotamerEntry re1 = singles.getTerm(mutPos1,curAAs[0],curRots[0]);
						RotamerEntry re2 = singles.getTerm(mutPos2,curAAs[1],curRots[1]);
						newSingles.setAll(p1,a1,r1,re1.combine(re2,pairTerm.minE()));//newIntra[p1][a1][r1] = re1.combine(re2,pairTerm.minE(),pairTerm.maxE(false));
					}
					else{ //just copy over the value
						RotamerEntry re1 = singles.getTerm(p1old,a1,r1);
						re1.pos = p1;
						newSingles.setAll(p1,a1,r1,re1);//newIntra[p1][a1][r1] = re1;
					}


					//We are getting rid of position mutPos2 so make the length 1 less
					newPair.addDim(p1a1r1ind, pairs.E.length-1);//newPair[p1][a1][r1] = new EMatrixEntry[pairE.length-1][][];
					int offSet2 = 0;
					for(int p2old = 0; p2old<pairs.E.length;p2old++){
						int p2 = p2old+offSet2;
						int[] p1a1r1p2ind = {p1,a1,r1,p2};
						if(p2old == mutPos2){
							offSet2--;
							continue;
						}
						else if((p1old == mutPos1 && p2old == mutPos1) 
								|| (p2old==mutPos1 && !areNeighbors(p1old, mutPos1) && !areNeighbors(p1old, mutPos2))
								|| (p1old == mutPos1 && !areNeighbors(p2old, mutPos1) && !areNeighbors(p2old, mutPos2))){
							//Do Nothing
							continue;
						}
						else if(p2old == mutPos1){
							//If we are at the position we're merging, mult number of AA
							//Figure out how many AA combinations have viable rotamers (i.e. non-zero number of rotamers)
							int numAAs = getNumCombinedAA(mutPos1, mutPos2, maxAAs);
							newPair.addDim(p1a1r1p2ind, numAAs);//newPair[p1][a1][r1][p2] = new EMatrixEntry[numAAs][];
						}
						else if(p1old == mutPos1){
							if(areNeighbors(mutPos1, p2old))
								newPair.addDim(p1a1r1p2ind, pairs.E[mutPos1][curAAs[0]][curRots[0]][p2old].length);//newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old].length][];
							else if(areNeighbors(mutPos2,p2old))
								newPair.addDim(p1a1r1p2ind, pairs.E[mutPos2][curAAs[1]][curRots[1]][p2old].length);//newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old].length][];
						}
						else{
							if(pairs.E[p1old][a1][r1][p2old] != null)
								newPair.addDim(p1a1r1p2ind, pairs.E[p1old][a1][r1][p2old].length);  //newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1old][a1][r1][p2old].length][];
							else
								continue;
						}

						int[] curAAs2 = new int[2];curAAs2[0] = 0; curAAs2[1] = -1;
						for(int a2 = 0; a2<newPair.E[p1][a1][r1][p2].length;a2++){
							int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
							int[] maxRots2 = {0,0};
							int[] curRots2 = {0,0};

							if(p2old == mutPos1){
								do{
									Emat.incrementCtr(curAAs2, maxAAs);
								}while(getCombinedRotUnpruned(mutPos1,mutPos2,curAAs2)==0);

								int numRot = getCombinedRotUnpruned(mutPos1,mutPos2,curAAs2);
								newPair.addDim(p1a1r1p2a2ind, numRot);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[numRot];

								maxRots2 = new int[2];maxRots2[0] = singles.E[mutPos1][curAAs2[0]].length;maxRots2[1] = singles.E[mutPos2][curAAs2[1]].length;
								curRots2 = new int[2];curRots2[0] = 0; curRots2[1] = -1;
							}
							else if(p1old == mutPos1){
								if(areNeighbors(mutPos1, p2old))
									newPair.addDim(p1a1r1p2a2ind, pairs.E[mutPos1][curAAs[0]][curRots[0]][p2old][a2].length);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old][a2].length];
								else if(areNeighbors(mutPos2,p2old))
									newPair.addDim(p1a1r1p2a2ind, pairs.E[mutPos2][curAAs[1]][curRots[1]][p2old][a2].length);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old][a2].length];
							}
							else{
								//Just copy over the array
								if(p1old!=p2old){
									newPair.addDim(p1a1r1p2a2ind, pairs.E[p1old][a1][r1][p2old][a2].length);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[pairE[p1old][a1][r1][p2old][a2].length];
									//System.arraycopy(pairE[p1old][a1][r1][p2old][a2], 0, newPair[p1][a1][r1][p2][a2], 0, pairE[p1old][a1][r1][p2old][a2].length);
									int[] p1olda1r1p2olda2ind = {p1old,a1,r1,p2old,a2};
									newPair.copyRow(p1a1r1p2a2ind,pairs,p1olda1r1p2olda2ind);
								}

							}	


							for(int r2=0;r2<newPair.E[p1][a1][r1][p2][a2].length;r2++){									
								//If we are here p1old = mutPos1 || p2old = mutPos1
								//We want to only merge the pair once and then set it to both places
								if(p1old == mutPos1){
									//The new pairwise energy is (mut1-r2) + (mut2-r2)
									RotamerPairEntry re1,re2;
									if(areNeighbors(mutPos1, p2old))
										re1 = pairs.getTerm(mutPos1,curAAs[0],curRots[0],p2old,a2,r2, singles.supRot);
									else //If they are not neighbors there is no interaction (make a dummy entry that is not pruned and has 0 energy)
										re1 = new RotamerPairEntry(mutPos1,new SuperRotamer(singles.supRot[mutPos1][curAAs[0]][curRots[0]]),p2old,new SuperRotamer(singles.supRot[p2old][a2][r2]),0.0,false);
									if(areNeighbors(mutPos2,p2old))
										re2 = pairs.getTerm(mutPos2,curAAs[1],curRots[1],p2old,a2,r2, singles.supRot);
									else //If they are not neighbors there is no interaction (make a dummy entry that is not pruned and has 0 energy)
										re2 = new RotamerPairEntry(mutPos2,new SuperRotamer(singles.supRot[mutPos2][curAAs[1]][curRots[1]]),p2old,new SuperRotamer(singles.supRot[p2old][a2][r2]),0.0,false);

									newPair.setAll(p1,a1,r1,p2,a2,r2,re1.combine(re2, 1));
								}
								else if(p2old==mutPos1){
									do{
										Emat.incrementCtr(curRots2, maxRots2);
									}while(singles.pruned[mutPos1][curAAs2[0]][curRots2[0]] || singles.pruned[mutPos2][curAAs2[1]][curRots2[1]] || pairs.pruned[mutPos1][curAAs2[0]][curRots2[0]][mutPos2][curAAs2[1]][curRots2[1]]);
									RotamerPairEntry re1,re2;
									if(areNeighbors(p1old,mutPos1))
										re1 = pairs.getTerm(p1old,a1,r1,mutPos1,curAAs2[0],curRots2[0], singles.supRot);
									else
										re1 = new RotamerPairEntry(p1old,new SuperRotamer(singles.supRot[p1old][a1][r1]),mutPos1,new SuperRotamer(singles.supRot[mutPos1][curAAs2[0]][curRots2[0]]),0.0,false);
									if(areNeighbors(p1old,mutPos2))
										re2 = pairs.getTerm(p1old,a1,r1,mutPos2,curAAs2[1],curRots2[1], singles.supRot);
									else
										re2 = new RotamerPairEntry(p1old,new SuperRotamer(singles.supRot[p1old][a1][r1]),mutPos2,new SuperRotamer(singles.supRot[mutPos2][curAAs2[1]][curRots2[1]]),0.0,false);
									newPair.setAll(p1, a1, r1, p2, a2, r2, re1.combine(re2,2));
								}

								//we are just copying over so need to fix the pos
								//((RotamerPairEntry)newPair[p1][a1][r1][p2][a2][r2]).pos1 = p1;
								//((RotamerPairEntry)newPair[p1][a1][r1][p2][a2][r2]).pos2 = p2;


							}


						}
					}

				}
			}

		}	

		mergeSpots(mutPos1, mutPos2);

		singles = newSingles;
		pairs = newPair;
	}

	public Emat combineRots(Index3 mutRot1, Index3 mutRot2, Molecule m){

		SingleMats newSingles = new SingleMats(singles.doDih);
		PairMats newPair = new PairMats(pairs.doDih);


		newSingles.addDim(new int[0], singles.E.length-1);//EMatrixEntry[][][] newIntra = new EMatrixEntry[intraE.length-1][][];
		newPair.addDim(new int[0], pairs.E.length-1);//EMatrixEntry[][][][][][] newPair = new EMatrixEntry[pairE.length-1][][][][][];

		int offSet1 = 0;
		for(int p1old=0;p1old<pairs.E.length;p1old++){
			int p1 = p1old+offSet1;
			int[] p1ind = {p1};
			int[] maxAAs = new int[2];maxAAs[0] = singles.E[mutRot1.pos].length; maxAAs[1] = singles.E[mutRot2.pos].length;
			if(p1old == mutRot2.pos){
				//If we are at position mutPos2 we need to skip this position
				offSet1--;
				continue;
			}
			else if(p1old != mutRot1.pos){
				//Copy over singles
				newSingles.addDim(p1ind, singles.E[p1old].length);//newIntra[p1] = new EMatrixEntry[intraE[p1ctr].length][];
			
				for(int a1=0;a1<newSingles.E[p1].length;a1++){
					int[] p1a1ind = {p1,a1};
					newSingles.addDim(p1a1ind, singles.E[p1old][a1].length);//newIntra[p1][a1] = new EMatrixEntry[intraE[p1ctr][a1].length];
					for(int r1=0;r1<newSingles.E[p1][a1].length;r1++){
						RotamerEntry re1 = singles.getTerm(p1old,a1,r1);
						re1.pos = p1;
						newSingles.setAll(p1,a1,r1,re1);//newIntra[p1][a1][r1] = re1;
					}
				}
			}
			else if(p1old == mutRot1.pos){
				//If we are at the position we're merging we add the one new rotamer
				int numAAs = 1;
				newSingles.addDim(p1ind, numAAs);//newIntra[p1] = new EMatrixEntry[numAAs][];
				newPair.addDim(p1ind, numAAs);//newPair[p1] = new EMatrixEntry[numAAs][][][][];

				int a1 = 0;
				int[] p1a1ind = {p1,a1};
					
					if(p1 == mutRot1.pos){
						int numRot = 1; //Only have 1 rotamer;
						newSingles.addDim(p1a1ind, numRot);//newIntra[p1][a1] = new EMatrixEntry[numRot];//new EMatrixEntry[intraE[mutPos1][curAAs[0]].length * intraE[mutPos2][curAAs[1]].length];
						newPair.addDim(p1a1ind, numRot);//newPair[p1][a1] = new EMatrixEntry[numRot][][][];//new EMatrixEntry[pairE[mutPos1][curAAs[0]].length * pairE[mutPos2][curAAs[1]].length][][][];
						
					}

					int r1=0;
					int[] p1a1r1ind = {p1,a1,r1};
						if(p1old == mutRot1.pos){ //Combine the intra terms
							//do{
							//	Emat.incrementold(curRots, maxRots);
							//} while(singles.pruned[mutPos1][curAAs[0]][curRots[0]] || singles.pruned[mutPos2][curAAs[1]][curRots[1]] || pairs.pruned[mutPos1][curAAs[0]][curRots[0]][mutPos2][curAAs[1]][curRots[1]]);
							RotamerPairEntry pairTerm = pairs.getTerm(mutRot1,mutRot2, singles.supRot);
							RotamerEntry re1 = singles.getTerm(mutRot1);
							RotamerEntry re2 = singles.getTerm(mutRot2);
							newSingles.setAll(p1,a1,r1,re1.combine(re2,pairTerm.minE()));//newIntra[p1][a1][r1] = re1.combine(re2,pairTerm.minE(),pairTerm.maxE(false));
//						}



						//We are getting rid of position mutPos2 so make the length 1 less
						newPair.addDim(p1a1r1ind, pairs.E.length-1);//newPair[p1][a1][r1] = new EMatrixEntry[pairE.length-1][][];
						int offSet2 = 0;
						for(int p2old = 0; p2old<pairs.E.length;p2old++){
							int p2 = p2old+offSet2;
							int[] p1a1r1p2ind = {p1,a1,r1,p2};
							if(p2old == mutRot2.pos){ //skip the position we are removing
								offSet2--;
								continue;
							}
							else if((p1old == mutRot1.pos && p2old == mutRot1.pos) //Same position
|| (p2old==mutRot1.pos && !areNeighbors(p1old, mutRot1.pos,m) && !areNeighbors(p1old, mutRot2.pos,m)) //The mutPos isn't neighbors with our current position
								|| (p1old == mutRot1.pos && !areNeighbors(p2old, mutRot1.pos,m) && !areNeighbors(p2old, mutRot2.pos,m))){
								//Do Nothing
								continue;
							}
							else if(p1old == mutRot1.pos){
								if(areNeighbors(mutRot1.pos,p2old,m))
									newPair.addDim(p1a1r1p2ind, pairs.E[mutRot1.pos][mutRot1.aa][mutRot1.rot][p2old].length);//newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old].length][];
								else if(areNeighbors(mutRot2.pos,p2old,m))
									newPair.addDim(p1a1r1p2ind, pairs.E[mutRot2.pos][mutRot2.aa][mutRot2.rot][p2old].length);
							}
							else{
								
//								if(pairs.E[p1old][a1][r1][p2old] != null)
//									newPair.addDim(p1a1r1p2ind, pairs.E[p1old][a1][r1][p2old].length);  //newPair[p1][a1][r1][p2] = new EMatrixEntry[pairE[p1old][a1][r1][p2old].length][];
//								else
									continue;
							}

							int[] curAAs2 = new int[2];curAAs2[0] = 0; curAAs2[1] = -1;
							for(int a2 = 0; a2<newPair.E[p1][a1][r1][p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
								int[] maxRots2 = {0,0};
								int[] curRots2 = {0,0};

								if(p2old == mutRot1.pos){
									//KER: Shouldn't happen
									continue;
									//do{
									//	Emat.incrementold(curAAs2, maxAAs);
									//}while(getCombinedRotUnpruned(mutPos1,mutPos2,curAAs2)==0);

									//int numRot = getCombinedRotUnpruned(mutPos1,mutPos2,curAAs2);
//									int numRot = 1;
//									newPair.addDim(p1a1r1p2a2ind, numRot);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[numRot];
//
//									maxRots2 = new int[2];maxRots2[0] = singles.E[mutRot1.pos][curAAs2[0]].length;maxRots2[1] = singles.E[mutRot2.pos][curAAs2[1]].length;
//									curRots2 = new int[2];curRots2[0] = 0; curRots2[1] = -1;
								}
								else if(p1old == mutRot1.pos){
									if(areNeighbors(mutRot1.pos, p2old,m))
										newPair.addDim(p1a1r1p2a2ind, pairs.E[mutRot1.pos][mutRot1.aa][mutRot1.rot][p2old][a2].length);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[pairE[p1old][curAAs[0]][curRots[0]][p2old][a2].length];
									else if(areNeighbors(mutRot2.pos,p2old,m))
										newPair.addDim(p1a1r1p2a2ind, pairs.E[mutRot2.pos][mutRot2.aa][mutRot2.rot][p2old][a2].length);
								}
								else{
									
									//Just copy over the array
									if(p1old!=p2old){
										newPair.addDim(p1a1r1p2a2ind, pairs.E[p1old][a1][r1][p2old][a2].length);//newPair[p1][a1][r1][p2][a2] = new EMatrixEntry[pairE[p1old][a1][r1][p2old][a2].length];
										//System.arraycopy(pairE[p1old][a1][r1][p2old][a2], 0, newPair[p1][a1][r1][p2][a2], 0, pairE[p1old][a1][r1][p2old][a2].length);
										int[] p1olda1r1p2olda2ind = {p1old,a1,r1,p2old,a2};
										newPair.copyRow(p1a1r1p2a2ind,pairs,p1olda1r1p2olda2ind);
									}

								}	


//								for(int r2=0;r2<newPair.E[p1][a1][r1][p2][a2].length;r2++){									
//									//If we are here p1old = mutPos1 || p2old = mutPos1
//									//We want to only merge the pair once and then set it to both places
//									if(p1old == mutRot1.pos){
//										//The new pairwise energy is (mut1-r2) + (mut2-r2)
//										
//										
////										RotamerPairEntry re1 = pairs.getTerm(mutRot1.pos,mutRot1.aa,mutRot1.rot,p2old,a2,r2,singles.supRot);
////										re2 = new RotamerPairEntry(mutRot2.pos,new SuperRotamer(singles.supRot[mutRot2.pos]))
////									if(areNeighbors(mutRot2.pos,p2old,m))
////										re2 = pairs.getTerm(mutRot2.pos,0,0,p2old,a2,r2, singles.supRot);
////									else //If they are not neighbors there is no interaction (make a dummy entry that is not pruned and has 0 energy)
////										re2 = new RotamerPairEntry(mutRot2.pos,new SuperRotamer(singles.supRot[mutRot2.pos][0][0]),p2old,new SuperRotamer(singles.supRot[p2old][a2][r2]),0.0,false);
//										//RotamerPairEntry re1 = pairs.getTerm(mutRot1.pos,mutRot1.aa,mutRot1.rot,p2old,a2,r2,singles.supRot);
//
//
//										//IMPORTANT: Energies will no longer be correct but that's ok
//										//since we are recalculating energies after this function
////										RotamerEntry re2 = singles.getTerm(mutRot2.pos,mutRot2.aa,mutRot2.rot);
//
////										newPair.setAll(p1,a1,r1,p2,a2,r2,re1.combine(re2, 1));
//									}
//									/*else if(p2old==mutRot1.pos){
//									do{
//										Emat.incrementCtr(curRots2, maxRots2);
//									}while(singles.pruned[mutPos1][curAAs2[0]][curRots2[0]] || singles.pruned[mutPos2][curAAs2[1]][curRots2[1]] || pairs.pruned[mutPos1][curAAs2[0]][curRots2[0]][mutPos2][curAAs2[1]][curRots2[1]]);
//RotamerPairEntry re1,re2;
//									if(areNeighbors(p1old,mutPos1))
//										re1 = pairs.getTerm(p1old,a1,r1,mutPos1,curAAs2[0],curRots2[0], singles.supRot);
//									else
//										re1 = new RotamerPairEntry(p1old,new SuperRotamer(singles.supRot[p1old][a1][r1]),mutPos1,new SuperRotamer(singles.supRot[mutPos1][curAAs2[0]][curRots2[0]]),0.0,false);
//									if(areNeighbors(p1old,mutPos2))
//										re2 = pairs.getTerm(p1old,a1,r1,mutPos2,curAAs2[1],curRots2[1], singles.supRot);
//									else
//										re2 = new RotamerPairEntry(p1old,new SuperRotamer(singles.supRot[p1old][a1][r1]),mutPos2,new SuperRotamer(singles.supRot[mutPos2][curAAs2[1]][curRots2[1]]),0.0,false);
//									newPair.setAll(p1, a1, r1, p2, a2, r2, re1.combine(re2,2));
////									RotamerPairEntry re1 = pairs.getTerm(p1old,a1,r1,mutPos1,curAAs2[0],curRots2[0]);
////									RotamerPairEntry re2 = pairs.getTerm(p1old,a1,r1,mutPos2,curAAs2[1],curRots2[1]);
////									newPair.setAll(p1, a1, r1, p2, a2, r2, re1.combine(re2,2));
//								}*/
//
//									//we are just copying over so need to fix the pos
//									//((RotamerPairEntry)newPair[p1][a1][r1][p2][a2][r2]).pos1 = p1;
//									//((RotamerPairEntry)newPair[p1][a1][r1][p2][a2][r2]).pos2 = p2;
//
//
//								}


							}
						}
					}
				}
			}



		Emat newMat = new Emat(templ_E,resByPos,newSingles,newPair);
		newMat.mergeSpots(mutRot1.pos, mutRot2.pos);

		return newMat;
	}

	private int getNumCombinedAA(int mutPos1, int mutPos2, int[] maxAAs) {
		int numAA = 0;
		int[] curAAs = new int[2];curAAs[0] = 0; curAAs[1] = -1;
		int totalAAposs = pairs.E[mutPos1].length*pairs.E[mutPos2].length; 
		for(int i=0; i<totalAAposs;i++){
			Emat.incrementCtr(curAAs, maxAAs);
			int numRot = getCombinedRotUnpruned(mutPos1, mutPos2, curAAs);
			if(numRot>0)
				numAA++;
		}
		return numAA;
	}

	private int getCombinedRotUnpruned(int mutPos1, int mutPos2, int[] curAAs) {
		int numUnpruned = 0;
		for(int r1=0; r1<pairs.E[mutPos1][curAAs[0]].length;r1++){
			for(int r2=0; r2<pairs.E[mutPos1][curAAs[0]][r1][mutPos2][curAAs[1]].length;r2++){
				if(!pairs.pruned[mutPos1][curAAs[0]][r1][mutPos2][curAAs[1]][r2])
					numUnpruned++;
			}
		}
		return numUnpruned;
	}

	void mergeSpots(int mutPos1,int mutPos2){
		for(Integer r:resByPos.get(mutPos2)){
			resByPos.get(mutPos1).add(r);
		}
		resByPos.remove(mutPos2);
	}

	//This isn't quite correct because the max should be updated if the curAA changes,
	//but this works due to the matrix structure that if the curAA is defined rotamer 0 will be defined.
	static boolean incrementCtr(int[] ctr, int[] max) {
		boolean updateMax = false;

		ctr[ctr.length-1]++;
		for(int i=ctr.length-1;i>0; i--){
			if(ctr[i] >= max[i]){
				ctr[i] = 0;
				ctr[i-1]++;
				updateMax = true; //Changed ctr in array need to update max;
			}
		}
		return updateMax;
	}

	//Checks if all the pairs or all the singles for a given position are pruned
	public void checkIfAllPruned(Molecule m) {
		//Singles
		for(int i=0; i<resByPos.size();i++){
			Iterator<EMatrixEntryWIndex> iter1 = singlesIterator(i);
			boolean singlesPruned = true;
			while(iter1.hasNext()){
				EMatrixEntryWIndex pos1entry = iter1.next();
				if(!pos1entry.eme.isPruned()){
					singlesPruned = false;
					break;
				}
			}
			if(singlesPruned){
				System.out.println("All rotamers at position: "+i+" are pruned");
				System.exit(0);
			}
		}

		//Pairs
		for(int i=0; i<resByPos.size();i++){
			for(int j=i+1;j<resByPos.size();j++){
				Iterator<EMatrixEntryWIndex> iter1 = pairsIterator(i, j);
				boolean pairsPruned = true;
				while(iter1.hasNext()){
					EMatrixEntryWIndex pos1entry = iter1.next();
					if(!pos1entry.eme.isPruned() && !getSinglePruned(pos1entry.rot1index()) && !getSinglePruned(pos1entry.rot2index())){
						pairsPruned = false;
						break;
					}
				}
				if(pairsPruned && areNeighbors(i, j)){
					System.out.println("All pairs for positions: "+i+" "+j+" are pruned");
					if(m != null){
						System.out.println("The following residues are involved: ");
						int posToCheck[] = {i,j};
						for(int pos : posToCheck)
							for(int resNum: resByPos.get(pos)){
								System.out.println(m.residue[resNum].fullName);
							}
					}
					System.exit(0);
				}
			}
		}

	}

	public int numMutRes() {
		int numMutRes = 0;
		for(ArrayList<Integer> resAtPos: resByPos){
			numMutRes += resAtPos.size();
		}
		return numMutRes;
	}

	/**
	 *	IMPORTANT: This function assumes no superrotamers (i.e. nothing has been combined) 
	 * @param resNum
	 * @param globalRot
	 * @return
	 */
	public Index3 getGlobalRot(int resNum, int globalRot) {

		//KER: First find position
		int pos = 0;
		for(ArrayList<Integer> residues: resByPos){
			if(residues.contains(resNum))
				break;
			pos++;
		}

		SinglesIterator iter = singlesIterator(pos);
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(singles.getRot(emeWI.index)[0] == globalRot){
				return new Index3(emeWI.index);
			}
		}

		return null;
	}

	
	public Index3 getGlobalRotForPos(int pos, int globalRot) {

		SinglesIterator iter = singlesIterator(pos);
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(singles.getRot(emeWI.index)[0] == globalRot){
				return new Index3(emeWI.index);
			}
		}

		return null;
	}
	
//	ArrayList<Rotamer> getRotamers(Molecule m, Index3 rot){
//		return singles.getTerm(rot).r.getRotamers(m, resByPos.get(rot.pos));
//	}
	
	ArrayList<ResidueConformation> getRCs(Molecule m, Index3 rot){
		return singles.getTerm(rot).r.getRCs(m, resByPos.get(rot.pos));
	}

	//KER: go through the Emat and for every rotamer that's not pruned
	//KER: split that rotamer do that it minimize half the distance it did before
	//IMPORTANT: This function assumes no superRotamer actually exist
	public void splitRotamers(Molecule m){


		//First just going to add the rotamer to the 6th dimension
		for(int p1=0; p1<pairs.E.length;p1++){
			for(int a1=0; a1<pairs.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				int resID = resByPos.get(p1).get(0);
				String pdbNum = m.residue[resID].getResNumberString();
				if(singles.E[p1][a1].length>0){
					AARotamerType aaT = m.residue[resID].rl.getRot(singles.supRot[p1][a1][0][0]).aaType;

					//We are splitting rotamers so we double the length of the array
					//KER: Except we shouldn't split ALA or GLY res
					if(aaT.numDihedrals() > 0){

						int prevLength = singles.E[p1][a1].length;
						singles.extendDim(p1a1ind, prevLength);
						for(int r1=0; r1<prevLength;r1++){
							int dihedDepth = 0;

							Rotamer curRot = m.residue[resID].rl.getRot(singles.supRot[p1][a1][r1][0]);

							for(int i=0; i<curRot.minimizationWidth.length-1;i++)
								if(curRot.minimizationWidth[i+1] > curRot.minimizationWidth[i]){
									dihedDepth = i+1;
								}

							double rotInterval = (curRot.minimizationWidth[dihedDepth])/(2);
							//For each rotamer in the original singles matrix we need to split it
							for(int i=0; i<2;i++){
								double[] dihedVals = new double[curRot.values.length];
								double[] minimizationWidth = new double[curRot.minimizationWidth.length];
								for(int q=0; q<dihedVals.length;q++){
									dihedVals[q] = curRot.values[q];
									minimizationWidth[q] = curRot.minimizationWidth[q];
								}
								int dir = i*2-1;
								dihedVals[dihedDepth] += dir*rotInterval;
								minimizationWidth[dihedDepth] = rotInterval;

								Rotamer subRot = m.residue[resID].rl.addRotamer(curRot.aaType.name, pdbNum,dihedVals,minimizationWidth,false);
								int[] subRotArray = new int[1];
								subRotArray[0] = subRot.rlIndex;
								singles.setSupRot(p1, a1, r1+prevLength*i, subRotArray);
							}


						}}}
			}
		}

		//Now we build the entire pairs matrix based on the new singles array setup
		this.pairs = new PairMats(pairs.doDih);
		this.pairs.addDim(new int[0],singles.E.length);//pairE = new EMatrixEntry[strandMut.numMutPos()][][][][][];

		//Runtime runtime = Runtime.getRuntime();

		for(int p1=0; p1<singles.E.length;p1++){
			int[] p1ind = {p1};

			this.pairs.addDim(p1ind,singles.E[p1].length);//pairE[p1] = new EMatrixEntry[res1.allowedAATypes.size()][][][][];

			for(int a1=0;a1<singles.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				this.pairs.addDim(p1a1ind,singles.E[p1][a1].length);//pairE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()][][][];
				for(int r1=0; r1<singles.E[p1][a1].length;r1++){
					int[] p1a1r1ind = {p1,a1,r1};

					this.pairs.addDim(p1a1r1ind,singles.E.length);//pairE[p1][a1][r1] = new EMatrixEntry[strandMut.numMutPos()][][];
					for(int p2=0; p2<singles.E.length;p2++){
						if(p1!=p2){
							int[] p1a1r1p2ind = {p1,a1,r1,p2};
							this.pairs.addDim(p1a1r1p2ind,singles.E[p2].length);//pairE[p1][a1][r1][p2] = new EMatrixEntry[res2.allowedAATypes.size()][];
							for(int a2=0;a2<singles.E[p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
								this.pairs.addDim(p1a1r1p2a2ind,singles.E[p2][a2].length);//pairE[p1][a1][r1][p2][a2] = new EMatrixEntry[aaType2.numRotamers()];
								//for(int r2=0; r2<singles.E[p2][a2].length;r2++){
								//int[] rot2Array = {rot2.rlIndex};
								//this.pairs.setSupRot(p1, a1, r1, p2, a2, r2, singles.supRot[p1][a1][r1], singles.supRot[p2][a2][r2]);//pairE[p1][a1][r1][p2][a2][r2] = rp;
								//}
							}
						}
					}
				}
			}
		}


	}

	//KER: Similar to addRotamers except that we are updating the rotamers to a known position
	public void updateRotamer(Molecule m,Index3 mutRotamer, ArrayList<ResidueConformation> rotToUpdate, boolean intraOnly){

		ArrayList<Residue> resForPos1 = new ArrayList<Residue>();
		int ctr=0;
		for(ResidueConformation r : rotToUpdate){
			resForPos1.add(m.residue[resByPos.get(mutRotamer.pos).get(ctr)]);
			ctr++;
		}

		int[] rotArrayToAdd = new int[rotToUpdate.size()];
		for(int i=0; i<rotArrayToAdd.length;i++)
			rotArrayToAdd[i] = rotToUpdate.get(i).id;

		//Run through every pair and then update rotamer

		for(int p1 = 0; p1<pairs.E.length;p1++){
			for(int a1=0;a1<pairs.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				//Extend Dimension by one for both pair and single
				for(int r1 = 0;r1<pairs.E[p1][a1].length;r1++){
					int[] p1a1r1ind = {p1,a1,r1};
					Index3 p1a1r1ind3 = new Index3(p1,a1,r1);
					//Add rot to singles
					if(p1a1r1ind3.equals(mutRotamer)){
						System.arraycopy(rotArrayToAdd, 0, singles.supRot[p1][a1][r1], 0, rotArrayToAdd.length);
						/*for(int p2=0; p2<pairs.E[p1][a1][r1].length;p2++){
							if(p2 != p1){
								for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){
									for(int r2=0;r2<pairs.E[p1][a1][r1][p2][a2].length;r2++){
										System.arraycopy(rotArrayToAdd, 0, pairs.supRot1[p1][a1][r1][p2][a2][r2], 0, rotArrayToAdd.length);
									}
								}
							}
						}*/
					}
					/*else if(p1 != mutRotamer.pos){
						//KER: set the rot 6 dimension
						System.arraycopy(rotArrayToAdd, 0, pairs.supRot2[p1][a1][r1][mutRotamer.pos][mutRotamer.aa][mutRotamer.rot], 0, rotArrayToAdd.length);
					}*/
				}



			}
		}


	}


	public ArrayList<Index3> addRotamers(Molecule m,int mutPos, ArrayList<ArrayList<ResidueConformation>> rotsToAdd, boolean intraOnly){

		ArrayList<Index3> returnTuples = new ArrayList<Index3>();


		//KER: Loop through each rotamer adding the correct dimensions
		for(ArrayList<ResidueConformation> rotsForPos:rotsToAdd){

			ArrayList<Residue> resForPos1 = new ArrayList<Residue>();
			int ctr=0;
			for(ResidueConformation r : rotsForPos){
				resForPos1.add(m.residue[resByPos.get(mutPos).get(ctr)]);
				ctr++;
			}

			int[] rot2Array = new int[rotsForPos.size()];
			for(int i=0; i<rot2Array.length;i++)
				rot2Array[i] = rotsForPos.get(i).id;

			//First just going to add the rotamer to the 6th dimension
			for(int p1=0; p1<pairs.E.length;p1++){
				if(p1 != mutPos && areNeighbors(p1, mutPos, m)){
					for(int a1=0; a1<pairs.E[p1].length;a1++){
						for(int r1=0; r1<pairs.E[p1][a1].length;r1++){

							int p2 = mutPos;
							boolean addedRot = false;
							for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){
								boolean correctAA = true;
								if(singles.supRot[p2][a2].length == 0)
									continue;

								int[] rotInd2 = singles.supRot[p2][a2][0]; //Just want to extract the AA type (and all rotamers have same type so just grab first one)
								int a2ctr=0;
								for(Residue res1 : resForPos1){

									AARotamerType aaType2 = m.strand[res1.strandNumber].rcl.getRC(rotInd2[a2ctr]).rot.aaType; 

									if(!aaType2.name.equalsIgnoreCase(rotsForPos.get(a2ctr).rot.aaType.name)){
										correctAA = false;
									}	
									a2ctr++;
								}

								if(correctAA){
									int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
									pairs.extendDim(p1a1r1p2a2ind, 1);
									addedRot = true;



									//KER: Add the rotamer
									//pairs.setSupRot(p1, a1, r1, p2, a2, pairs.E[p1][a1][r1][p2][a2].length-1, rot1Array, rot2Array);

								}

							}
							if(!addedRot)
								System.out.println("Error adding rot");
						}
					}
				}
			}

			//Now we add the rotamer to the singles array and matrix to the 3 dim

			int[] rot1Array = new int[rotsForPos.size()];
			for(int i=0; i<rot1Array.length;i++)
				rot1Array[i] = rotsForPos.get(i).id;

			int p1 = mutPos;
			for(int a1=0;a1<pairs.E[p1].length;a1++){
				boolean correctAA = true;
				if(singles.supRot[p1][a1].length==0)
					continue;
				int[] rotInd1 = singles.supRot[p1][a1][0]; //Just want to extract the AA type (and all rotamers have same type so just grab first one)
				int a1ctr=0;
				for(Residue res1 : resForPos1){
					AARotamerType aaType1 = m.strand[res1.strandNumber].rcl.getRC(rotInd1[a1ctr]).rot.aaType;
					//AARotamerType aaType1 = res1.allowedAATypes.get(a1);
					if(!aaType1.name.equalsIgnoreCase(rotsForPos.get(a1ctr).rot.aaType.name)){
						correctAA = false;
					}	
					a1ctr++;
				}

				if(correctAA){	
					int[] p1a1ind = {p1,a1};
					//Extend Dimension by one for both pair and single
					this.singles.extendDim(p1a1ind,1);//intraE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()];
					if(!intraOnly)
						this.pairs.extendDim(p1a1ind,1);//pairE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()][][][];

					int r1 = singles.E[p1][a1].length-1;
					int[] p1a1r1ind = {p1,a1,r1};
					returnTuples.add(new Index3(p1,a1,r1));
					//Add rot to singles
					this.singles.setSupRot(p1, a1, r1, rot1Array);
					if(!intraOnly){
						this.pairs.addDim(p1a1r1ind,resByPos.size());//pairE[p1][a1][r1] = new EMatrixEntry[strandMut.numMutPos()][][];
					}
					for(int p2=0;p2<pairs.E[p1][a1][r1].length;p2++){
						if(p2 != p1 && areNeighbors(p1, p2, m)){
							int[] p1a1r1p2ind = {p1,a1,r1,p2};
							this.pairs.addDim(p1a1r1p2ind,singles.E[p2].length);
							for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
								this.pairs.addDim(p1a1r1p2a2ind,singles.E[p2][a2].length);
								for(int r2=0;r2<pairs.E[p1][a1][r1][p2][a2].length;r2++){
									rot2Array = new int[singles.supRot[p2][a2][r2].length];
									for(int i=0; i<rot2Array.length;i++)
										rot2Array[i] = singles.supRot[p2][a2][r2][i]; 
									//pairs.setSupRot(p1,a1,r1,p2,a2,r2,rot1Array,rot2Array);
								}
							}
						}
					}



				}
			}
		}

		return returnTuples;

	}



	/**
	 * Takes in a list of rotamers to add at every position.
	 * This function makes sure that you only need to extend each array once instead of
	 * calling the above function (addRotamers).
	 * @param m
	 * @param rotsToAdd - Rotamers to add by position,AA
	 * @param intraOnly
	 * @return
	 */
	public void addAllRotamers(Molecule m, ArrayList<ArrayList<ArrayList<ArrayList<Rotamer>>>> rotsToAddByPosAA, boolean intraOnly){


		for(int mutPos=0; mutPos<rotsToAddByPosAA.size();mutPos++){
			ArrayList<ArrayList<ArrayList<Rotamer>>> rotsToAddByAA = rotsToAddByPosAA.get(mutPos);
			if(rotsToAddByAA.size() == 0)
				continue;
			//KER: Loop through each rotamer adding the correct dimensions
			//		for(ArrayList<Rotamer> rotsForPos:rotsToAdd){

			//			ArrayList<Residue> resForPos1 = new ArrayList<Residue>();
			//			int ctr=0;
			//			for(Rotamer r : rotsForPos){
			//				resForPos1.add(m.residue[resByPos.get(mutPos).get(ctr)]);
			//				ctr++;
			//			}
			//
			//			int[] rot2Array = new int[rotsForPos.size()];
			//			for(int i=0; i<rot2Array.length;i++)
			//				rot2Array[i] = rotsForPos.get(i).rlIndex;
			if(!intraOnly){
				//First just going to add the rotamers to the 6th dimension
				for(int p1=0; p1<pairs.E.length;p1++){
					if(p1 != mutPos && areNeighbors(p1, mutPos, m)){
						for(int a1=0; a1<pairs.E[p1].length;a1++){
							for(int r1=0; r1<pairs.E[p1][a1].length;r1++){

								//							int[] rot1Array = new int[singles.supRot[p1][a1][r1].length];
								//							for(int i=0; i<rot1Array.length;i++)
								//								rot1Array[i] = singles.supRot[p1][a1][r1][i];

								int p2 = mutPos;
								boolean addedRot = false;
								for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){

									if(singles.supRot[p2][a2].length == 0)
										continue;

									int numRotsToAdd = rotsToAddByAA.get(a2).size();


									if(numRotsToAdd > 0){
										int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
										pairs.extendDim(p1a1r1p2a2ind, numRotsToAdd);
										addedRot = true;
									}

								}
								if(!addedRot)
									System.out.println("Error adding rot");
							}
						}
					}
				}
			}

			//Now we add the rotamer to the singles array and matrix to the 3 dim

			//			int[] rot1Array = new int[rotsForPos.size()];
			//			for(int i=0; i<rot1Array.length;i++)
			//				rot1Array[i] = rotsForPos.get(i).rlIndex;


			int p1 = mutPos;
			for(int a1=0;a1<singles.E[p1].length;a1++){

				if(singles.supRot[p1][a1].length==0)
					continue;

				ArrayList<ArrayList<Rotamer>> curRotsToAdd = rotsToAddByAA.get(a1);


				if(curRotsToAdd.size() > 0){
					int[] p1a1ind = {p1,a1};
					int r1 = singles.E[p1][a1].length;

					//Extend Dimension by one for both pair and single
					this.singles.extendDim(p1a1ind,curRotsToAdd.size());//intraE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()];
					if(!intraOnly)
						this.pairs.extendDim(p1a1ind,curRotsToAdd.size());//pairE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()][][][];

					for(ArrayList<Rotamer> supRot: curRotsToAdd){
						//Add rot to singles
						int[] rot1Array = new int[supRot.size()];
						for(int i=0; i<supRot.size();i++)
							rot1Array[i] = supRot.get(i).rlIndex;

						int[] p1a1r1ind = {p1,a1,r1};
						this.singles.setSupRot(p1, a1, r1, rot1Array);
						if(!intraOnly){
							this.pairs.addDim(p1a1r1ind,resByPos.size());//pairE[p1][a1][r1] = new EMatrixEntry[strandMut.numMutPos()][][];
							for(int p2=0;p2<pairs.E[p1][a1][r1].length;p2++){
								if(p2 != p1 && areNeighbors(p1, p2, m)){
									int[] p1a1r1p2ind = {p1,a1,r1,p2};
									this.pairs.addDim(p1a1r1p2ind,singles.E[p2].length);
									for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){
										int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
										this.pairs.addDim(p1a1r1p2a2ind,singles.E[p2][a2].length);
										//									for(int r2=0;r2<pairs.E[p1][a1][r1][p2][a2].length;r2++){
										//										int[] rot2Array = new int[singles.supRot[p2][a2][r2].length];
										//										for(int i=0; i<rot2Array.length;i++)
										//											rot2Array[i] = singles.supRot[p2][a2][r2][i]; 
										//										//pairs.setSupRot(p1,a1,r1,p2,a2,r2,rot1Array,rot2Array);
										//									}
									}
								}
							}
						}

						r1++;
					}
				}
			}
			//		}
		}

		//		return returnTuples;

	}
	
	/**
	 * Takes in a list of rotamers to add at every position.
	 * This function makes sure that you only need to extend each array once instead of
	 * calling the above function (addRotamers).
	 * @param m
	 * @param rotsToAdd - Rotamers to add by position,AA
	 * @param intraOnly
	 * @return
	 */
	public void addAllRCs(Molecule m, ArrayList<ArrayList<ArrayList<ArrayList<ResidueConformation>>>> rcsToAddByPosAA, boolean intraOnly){


		for(int mutPos=0; mutPos<rcsToAddByPosAA.size();mutPos++){
			ArrayList<ArrayList<ArrayList<ResidueConformation>>> rcsToAddByAA = rcsToAddByPosAA.get(mutPos);
			if(rcsToAddByAA.size() == 0)
				continue;
			
			if(!intraOnly){
				//First just going to add the rotamers to the 6th dimension
				for(int p1=0; p1<pairs.E.length;p1++){
					if(p1 != mutPos && areNeighbors(p1, mutPos, m)){
						for(int a1=0; a1<pairs.E[p1].length;a1++){
							for(int r1=0; r1<pairs.E[p1][a1].length;r1++){
								int p2 = mutPos;
								boolean addedRot = false;
								for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){

									if(singles.supRot[p2][a2].length == 0)
										continue;

									int numRotsToAdd = rcsToAddByAA.get(a2).size();


									if(numRotsToAdd > 0){
										int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
										pairs.extendDim(p1a1r1p2a2ind, numRotsToAdd);
										addedRot = true;
									}

								}
								if(!addedRot)
									System.out.println("Error adding rot");
							}
						}
					}
				}
			}

			int p1 = mutPos;
			for(int a1=0;a1<singles.E[p1].length;a1++){

				if(singles.supRot[p1][a1].length==0)
					continue;

				ArrayList<ArrayList<ResidueConformation>> curRCsToAdd = rcsToAddByAA.get(a1);


				if(curRCsToAdd.size() > 0){
					int[] p1a1ind = {p1,a1};
					int r1 = singles.E[p1][a1].length;

					//Extend Dimension by one for both pair and single
					this.singles.extendDim(p1a1ind,curRCsToAdd.size());//intraE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()];
					if(!intraOnly)
						this.pairs.extendDim(p1a1ind,curRCsToAdd.size());//pairE[p1][a1] = new EMatrixEntry[aaType1.numRotamers()][][][];

					for(ArrayList<ResidueConformation> supRot: curRCsToAdd){
						//Add rot to singles
						int[] rot1Array = new int[supRot.size()];
						for(int i=0; i<supRot.size();i++)
							rot1Array[i] = supRot.get(i).id;

						int[] p1a1r1ind = {p1,a1,r1};
						this.singles.setSupRot(p1, a1, r1, rot1Array);
						if(!intraOnly){
							this.pairs.addDim(p1a1r1ind,resByPos.size());//pairE[p1][a1][r1] = new EMatrixEntry[strandMut.numMutPos()][][];
							for(int p2=0;p2<pairs.E[p1][a1][r1].length;p2++){
								if(p2 != p1 && areNeighbors(p1, p2, m)){
									int[] p1a1r1p2ind = {p1,a1,r1,p2};
									this.pairs.addDim(p1a1r1p2ind,singles.E[p2].length);
									for(int a2=0;a2<pairs.E[p1][a1][r1][p2].length;a2++){
										int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
										this.pairs.addDim(p1a1r1p2a2ind,singles.E[p2][a2].length);
									}
								}
							}
						}

						r1++;
					}
				}
			}

		}

	}


	/**
	 * Prune all rotamers that are not listed as allowed in the residue
	 * @param m
	 */
	public void pruneNotAllowed(Molecule m) {

		for(int pos = 0; pos<resByPos.size();pos++){
			SinglesIterator iter = singlesIterator(pos);
			while(iter.hasNext()){
				EMatrixEntryWIndex rot = iter.next();

				int[] rcs = singles.getRot(rot.index);
				for(int i=0; i<rcs.length;i++){
					Residue r = m.residue[resByPos.get(pos).get(i)];

					if(!r.isResAllowed(m.strand[r.strandNumber].rcl, rcs[i])){
						setPruned(rot.index,true);
					}
					
				}
			}
		}
	}

	/**
	 * Determines if two positions in the Emat are neighbors
	 * For now I just check if the pair matrix is defined for
	 * these two positions
	 * @param pos1
	 * @param pos2
	 * @return areNeighbors
	 */
	public boolean areNeighbors(int pos1, int pos2) {
		//Find a valid rotamer for pos1 and then check
		//if pos2 is defined for the rotamer
		for(int a1=0; a1<pairs.pruned[pos1].length;a1++){
			for(int r1=0; r1<pairs.pruned[pos1][a1].length;r1++){
				if(pairs.pruned[pos1][a1][r1][pos2] == null)
					return false;
				else
					return true;
			}
		}

		return true;
	}

	public boolean doDih(){
		return singles.doDih;
	}

	public void setMaxE(int[] i, double maxE){
		if(i.length == 3)
			singles.setMaxE(i,maxE);
		else if(i.length == 6)
			pairs.setMaxE(i,maxE);
		else
			System.out.println("Index length not valid");
	}

	public void setDihedrals(int[] i, double[][] diheds){
		if(i.length == 3)
			singles.setDihedrals(i,diheds[0]);
		else if(i.length == 6)
			pairs.setDihedrals(i,diheds);
		else
			System.out.println("Index length not valid");
	}

	/**
	 * Retuns an ArrayList<ArrayList<?> that has dimensions equal to the first two dimensions of 
	 * the singles matrix
	 * @return
	 */
	public ArrayList<ArrayList<ArrayList<ArrayList<Rotamer>>>> getEmptyRotsToAdd() {

		ArrayList<ArrayList<ArrayList<ArrayList<Rotamer>>>> retArray = new ArrayList<ArrayList<ArrayList<ArrayList<Rotamer>>>>();

		for(int p = 0; p<singles.E.length;p++){
			ArrayList<ArrayList<ArrayList<Rotamer>>> tmp = new ArrayList<ArrayList<ArrayList<Rotamer>>>();
			for(int a=0; a<singles.E[p].length;a++){
				tmp.add(new ArrayList<ArrayList<Rotamer>>());
			}
			retArray.add(tmp);
		}

		return retArray;
	}
	
	/**
	 * Retuns an ArrayList<ArrayList<?> that has dimensions equal to the first two dimensions of 
	 * the singles matrix
	 * @return
	 */
	public ArrayList<ArrayList<ArrayList<ArrayList<ResidueConformation>>>> getEmptyRCsToAdd() {

		ArrayList<ArrayList<ArrayList<ArrayList<ResidueConformation>>>> retArray = new ArrayList<ArrayList<ArrayList<ArrayList<ResidueConformation>>>>();

		for(int p = 0; p<singles.E.length;p++){
			ArrayList<ArrayList<ArrayList<ResidueConformation>>> tmp = new ArrayList<ArrayList<ArrayList<ResidueConformation>>>();
			for(int a=0; a<singles.E[p].length;a++){
				tmp.add(new ArrayList<ArrayList<ResidueConformation>>());
			}
			retArray.add(tmp);
		}

		return retArray;
	}

	public int getRotAALoc(Molecule m, int pos, ArrayList<Rotamer> supRot) {



		int p1 = pos;
		for(int a1=0;a1<singles.E[p1].length;a1++){

			if(singles.supRot[p1][a1].length==0)
				continue;
			int[] rotInd1 = singles.supRot[p1][a1][0]; //Just want to extract the AA type (and all rotamers have same type so just grab first one)
			String[] aaName1 = new String[rotInd1.length];
			//AA names of the current emat rotamer
			for(int a1ctr=0; a1ctr<aaName1.length;a1ctr++)
				aaName1[a1ctr] = m.residue[resByPos.get(p1).get(a1ctr)].rl.getRot(rotInd1[a1ctr]).aaType.name;

			//Figure out how many Rotamers are of this AA type
			ArrayList<ArrayList<Rotamer>> curRotsToAdd = new ArrayList<ArrayList<Rotamer>>();
			boolean correctAA = true;
			int a1ctr = 0;

			for(Rotamer rot : supRot){

				//Is this AAtype the same as the rotamer we want to add
				if(!aaName1[a1ctr].equalsIgnoreCase(rot.aaType.name)){
					correctAA = false;
				}	
				a1ctr++;
			}

			if(correctAA){	
				return a1;
			}
		}



		return -1;

	}

	public int getRCAALoc(Molecule m, int pos, ArrayList<ResidueConformation> supRot) {



		int p1 = pos;
		for(int a1=0;a1<singles.E[p1].length;a1++){

			if(singles.supRot[p1][a1].length==0)
				continue;
			int[] rotInd1 = singles.supRot[p1][a1][0]; //Just want to extract the AA type (and all rotamers have same type so just grab first one)
			String[] aaName1 = new String[rotInd1.length];
			//AA names of the current emat rotamer
			for(int a1ctr=0; a1ctr<aaName1.length;a1ctr++)
				aaName1[a1ctr] = m.strand[m.residue[resByPos.get(p1).get(a1ctr)].strandNumber].rcl.getRC(rotInd1[a1ctr]).rot.aaType.name;

			//Figure out how many Rotamers are of this AA type
			ArrayList<ArrayList<Rotamer>> curRotsToAdd = new ArrayList<ArrayList<Rotamer>>();
			boolean correctAA = true;
			int a1ctr = 0;

			for(ResidueConformation rc : supRot){

				//Is this AAtype the same as the rotamer we want to add
				if(!aaName1[a1ctr].equalsIgnoreCase(rc.rot.aaType.name)){
					correctAA = false;
				}	
				a1ctr++;
			}

			if(correctAA){	
				return a1;
			}
		}



		return -1;

	}
	
	public Index3[] getConf(int[] globalRots) {
		Index3[] conf = new Index3[globalRots.length];
		for(int i=0; i<resByPos.size();i++){
			int resNum = resByPos.get(i).get(0);
			conf[i] = getGlobalRot(resNum, globalRots[i]);
		}
		return conf;
	}
	
	public ArrayList<ArrayList<Index3>> getRotsPerRes(int[] globalRots,Molecule m) {
		ArrayList<ArrayList<Index3>> conf = new ArrayList<ArrayList<Index3>>();
		for(int i=0; i<resByPos.size();i++){
			int resNum = resByPos.get(i).get(0);
			conf.add(getGlobalParents(resNum, globalRots[i],m));
		}
		return conf;
	}

	/**
	 *	IMPORTANT: This function assumes no superrotamers (i.e. nothing has been combined) 
	 * @param resNum
	 * @param globalRot
	 * @return
	 */
	public ArrayList<Index3> getGlobalParents(int resNum, int globalRot,Molecule m) {
		ArrayList<Index3> rots = new ArrayList<Index3>();
		//KER: First find position
		int pos = 0;
		for(ArrayList<Integer> residues: resByPos){
			if(residues.contains(resNum))
				break;
			pos++;
		}

		SinglesIterator iter = singlesIterator(pos);
		while(iter.hasNext()){
			EMatrixEntryWIndex emeWI = iter.next();
			if(m.residue[resNum].rl.getRot(singles.getRot(emeWI.index)[0]).parent == globalRot){
				rots.add(new Index3(emeWI.index));
			}
		}

		return rots;
	}

	public boolean[][][][][][] copyPairPruned() {
		return pairs.copyPairPruned();
	}

	public void setPairPruned(boolean[][][][][][] savedSplitFlags) {
		pairs.setPairPruned(savedSplitFlags);
	}

	public void setShellShellE(double minE) {
		templ_E = minE;
		
	}

	public void setShellRotE(int res1, int res1aaNum, int res1RotNum,
			double minEnergy) {
		singles.E[res1][res1aaNum][res1RotNum] = minEnergy;
	}
	
	public void setShellRotE(int[] i, double minEnergy) {
		singles.E[i[0]][i[1]][i[2]] = minEnergy;
	}

	public void setPairwiseE(int res1, int res1aaNum, int res1RotNum, int res2,
			int res2aaNum, int res2RotNum, double minEnergy) {
		pairs.E[res1][res1aaNum][res1RotNum][res2][res2aaNum][res2RotNum] = minEnergy;
		
	}
	
	public void setIntraE(int pos, int aa, int rot, double E){
		singles.E[pos][aa][rot] = E;
	}
	
	public void addToIntraE(int pos, int aa, int rot, double E){
		singles.E[pos][aa][rot] += E;
	}

	public void addToShellRotE(int pos, int aa, int rot, double E) {
		singles.E[pos][aa][rot] += E;
		
	}
	
	public void addToShellRotE(int[] i, double E) {
		singles.E[i[0]][i[1]][i[2]] += E;
		
	}

	public double getIntraAndShellE(int pos, int aa, int rot) {
		return singles.E[pos][aa][rot];
	}

	public double getShellShellE() {
		return templ_E;
	}

	public ArrayList<ArrayList<AARotamerType>> getAATypes(Molecule m, int[] i) {
		ArrayList<ArrayList<AARotamerType>> aaTypes = new ArrayList<ArrayList<AARotamerType>>();
		
		int[] rcs = singles.supRot[i[0]][i[1]][i[2]];
		ArrayList<AARotamerType> aaTypes1 = new ArrayList<AARotamerType>();
		int ctr = 0;
		for(int resID: resByPos.get(i[0]) ){
			Residue r = m.residue[resID];
			ResidueConformation rcl = m.strand[r.strandNumber].rcl.getRC(rcs[ctr]);
			aaTypes1.add(rcl.rot.aaType);
			ctr++;
		}
		
		aaTypes.add(aaTypes1);
		
		if(i.length > 3){
			rcs = singles.supRot[i[3]][i[4]][i[5]];
			ArrayList<AARotamerType> aaTypes2 = new ArrayList<AARotamerType>();
			ctr = 0;
			for(int resID: resByPos.get(i[3]) ){
				Residue r = m.residue[resID];
				ResidueConformation rcl = m.strand[r.strandNumber].rcl.getRC(rcs[ctr]);
				aaTypes2.add(rcl.rot.aaType);
				ctr++;
			}
			aaTypes.add(aaTypes2);
		}
		
		
		return aaTypes;
	}
	
	/**
	 * Goes through all pairs of mutable spots and finds the ones with the least number of pairs
	 */
	public int[] findPositionsToContractLeastPairs() {

		int bestPos1 = -1;
		int bestPos2 = -1;
		double minPairs = Double.MAX_VALUE;

		for(int p1 = 0; p1<numMutPos();p1++){
			for(int p2 = p1+1; p2<numMutPos();p2++){
				//loop through splitFlags to find all pairs for two positions that have been pruned

				//KER: since we now remove rotamers from the energy matrix the base pruning is
				//the number of rotamers still left in the eliminatedRotAtRes array
				int pairs = 0;
				//int totalNumPairs = 0;
				Iterator<EMatrixEntryWIndex> iter = pairsIterator(p1, p2);
				while(iter.hasNext()){
					EMatrixEntryWIndex eme = iter.next();
					int[] index1 = {eme.index[0],eme.index[1],eme.index[2]};
					int[] index2 = {eme.index[3],eme.index[4],eme.index[5]};
					if(!getSinglePruned(index1) && !getSinglePruned(index2) && !eme.eme.isPruned()){
						pairs++;
					}
					
				}
				if(pairs < minPairs){
					bestPos1 = p1;
					bestPos2 = p2;
					minPairs = pairs;
				}
			}
		}

		int[] pos = {bestPos1,bestPos2};
		return pos;
	}

	/**
	 * Given the residue find the residue with the closest CB atom to its CB atom
	 * This approach is hoping that things close should minimize together.
	 * @param res
	 * @return
	 */
	public int[] findPositionsToContractClosestRes(int pos, Molecule m) {

		Atom[] CBs = new Atom[resByPos.get(pos).size()];
		for(int i=0; i<CBs.length;i++)
			for(Atom a: m.residue[resByPos.get(pos).get(i)].atom)
				if(a.name.equalsIgnoreCase("CB"))
					CBs[i] = a;

		int closestPos = -1;
		double minDist = Double.POSITIVE_INFINITY;
		int ctr = 0;
		Iterator<ArrayList<Integer>> iter1 = resByPos.iterator();
		while(iter1.hasNext()){
			Iterator<Integer> iter2 = iter1.next().iterator();
			if(ctr != pos){
				while(iter2.hasNext()){
					Integer curRes = iter2.next();
					Atom curCB = null;
					for(Atom a: m.residue[curRes].atom)
						if(a.name.equalsIgnoreCase("CB"))
							curCB = a;

					for(int i=0; i<CBs.length;i++){
						if(CBs[i] != null && curCB != null && CBs[i].distance(curCB) < minDist){
							minDist = CBs[i].distance(curCB);
							closestPos = ctr;
						}
					}
				}
			}
			ctr++;
		}

		int[] posToMut = new int[2];
		if(closestPos < pos){
			posToMut[0] = closestPos;
			posToMut[1] = pos;
		}
		else{
			posToMut[0] = pos;
			posToMut[1] = closestPos;
		}

		return posToMut;
	}

	/**
	 * Goes through all pairs of mutable spots and finds the ones with the largest fraction
	 * of pairs pruned
	 */
	public int[] findPositionsToContract() {

		int bestPos1 = -1;
		int bestPos2 = -1;
		double maxPercent = 0;

		for(int p1 = 0; p1<numMutPos();p1++){
			for(int p2 = p1+1; p2<numMutPos();p2++){
				//loop through splitFlags to find all pairs for two positions that have been pruned

				//KER: since we now remove rotamers from the energy matrix the base pruning is
				//the number of rotamers still left in the eliminatedRotAtRes array
				int pairsPruned = 0;
				int totalNumPairs = 0;
				Iterator<EMatrixEntryWIndex> iter = pairsIterator(p1, p2);
				while(iter.hasNext()){
					EMatrixEntryWIndex eme = iter.next();
					int[] index1 = {eme.index[0],eme.index[1],eme.index[2]};
					int[] index2 = {eme.index[3],eme.index[4],eme.index[5]};
					if(getSinglePruned(index1) || getSinglePruned(index2) || eme.eme.isPruned()){
						pairsPruned++;
					}
					totalNumPairs++;
				}
				double fractionPruned = (double)pairsPruned/(double)totalNumPairs;
				if(fractionPruned > maxPercent){
					bestPos1 = p1;
					bestPos2 = p2;
					maxPercent = fractionPruned;
				}
				//KER: BYPASS PERCENT SELECTION IF totalNumPairs == 1
				if(totalNumPairs == 1){
					bestPos1 = p1;
					bestPos2 = p2;
					maxPercent = 100;
				}
				//System.out.println(p1+" "+p2+" "+totalNumPairs +" "+pairsPruned+" "+fractionPruned);
			}
		}

		//contractMutSpots(bestPos1, bestPos2);
		int[] pos = {bestPos1,bestPos2};
		return pos;
	}

	/**
	 * Goes through all pairs of mutable spots and finds the ones with the largest fraction
	 * of pairs pruned
	 */
	public int[] findPositionsToContractPercentLeast() {


		int[] bestPos1 = new int[5];
		int[] bestPos2 = new int[5];
		int[] pairs = new int[5];
		double[] maxPercent = new double[5];

		for(int i=0; i<maxPercent.length;i++){
			maxPercent[i] = 0.0;
			bestPos1[i] = -1;
			bestPos2[i] = -1;
		}

		System.out.println("Begin Find Contract Positions");
		for(int p1 = 0; p1<numMutPos();p1++){
			for(int p2 = p1+1; p2<numMutPos();p2++){
				if(areNeighbors(p1, p2)){
					//loop through splitFlags to find all pairs for two positions that have been pruned

					//KER: since we now remove rotamers from the energy matrix the base pruning is
					//the number of rotamers still left in the eliminatedRotAtRes array
					int pairsPruned = 0;
					int totalNumPairs = 0;
					Iterator<EMatrixEntryWIndex> iter = pairsIterator(p1, p2);
					while(iter.hasNext()){
						EMatrixEntryWIndex eme = iter.next();
						int[] index1 = {eme.index[0],eme.index[1],eme.index[2]};
						int[] index2 = {eme.index[3],eme.index[4],eme.index[5]};
						if(getSinglePruned(index1) || getSinglePruned(index2) || eme.eme.isPruned()){
							pairsPruned++;
						}
						totalNumPairs++;
					}
					double fractionPruned = (double)pairsPruned/(double)totalNumPairs;
					String fp = String.format("%.2f",fractionPruned );
					System.out.print(fp+" ");
					double[] minAndIndex = Util.indexOfMin(maxPercent);
					int index = (int) minAndIndex[1];
					if(fractionPruned > minAndIndex[0]){
						bestPos1[index] = p1;
						bestPos2[index] = p2;
						maxPercent[index] = fractionPruned;
						pairs[index] = totalNumPairs - pairsPruned;
					}
					//KER: BYPASS PERCENT SELECTION IF totalNumPairs == 1
					if(totalNumPairs == 1){
						int[] pos = {p1,p2};
						return pos;	
					}
					//System.out.println(p1+" "+p2+" "+totalNumPairs +" "+pairsPruned+" "+fractionPruned);
				}
			}System.out.println("");
		}
		int[] minAndIndex = Util.indexOfMin(pairs);
		int index = (int) minAndIndex[1];

		//contractMutSpots(bestPos1, bestPos2);
		int[] pos = {bestPos1[index],bestPos2[index]};
		return pos;
	}


}

class SinglesIterator implements Iterator<EMatrixEntryWIndex>{
	private Emat emat;
	private boolean hasNextItem = false;
	private EMatrixEntry nextItem = null;
	private int[] curI;
	private int pos = -1;

	public SinglesIterator(Emat energyMat){
		emat = energyMat;
		curI = new int[3];
		for(int i=0; i<curI.length;i++)
			curI[i] = 0;

		if(emat.singles != null){
			hasNextItem = true;
			while(emat.singles.E[curI[0]][curI[1]].length == 0){
				curI[1]++;
			}
			nextItem = emat.singles.getTerm(curI);

		}
		else{
			hasNextItem = false;
		}

	}

	// This iterator returns an iterator of all rotamers at position p.   
	public SinglesIterator(Emat energyMat, int p){
		pos = p;
		emat = energyMat;
		// 3D index for Position, AA, and rotamer.
		curI = new int[3];
		for(int i=0; i<curI.length;i++)
			curI[i] = 0;

		curI[0] = p;

		if(emat.singles != null){
			hasNextItem = true;
			try{
				// Find the first valid rotamer of a position.
				while(emat.singles.E[curI[0]][curI[1]].length == 0){
					curI[1]++;
				}
			}catch(Exception E){
				System.out.println("An error occurred while generating rotamer iterator.  Possibly, all rotamers were pruned at residue position "+pos+".  Try increasing the pruning energy so that not all rotamers are pruned.");
				E.printStackTrace();
				System.exit(1);

			}
			nextItem = emat.singles.getTerm(curI);

		}
		else{
			hasNextItem = false;
		}

	}

	public boolean hasNext() {
		return hasNextItem;
	}

	// Returns the next rotamer entry either pruned or unpruned. 
	public EMatrixEntryWIndex next() {
		int[] index = new int[curI.length];
		for(int i=0; i<index.length;i++)
			index[i] = curI[i];
		EMatrixEntryWIndex ret = new EMatrixEntryWIndex(nextItem, index);
		calcNext();
		return ret;
	}

	// Find the next rotamer; called by next function
	private void calcNext() {
		hasNextItem = false;

		int[] ctr = curI; 
		int[] max = {emat.singles.E.length,emat.singles.E[curI[0]].length,emat.singles.E[curI[0]][curI[1]].length};
		do{
			if(Emat.incrementCtr(ctr,max)){
				//Gone past the end of the matrix
				if(ctr[0] == emat.singles.E.length)
					return;
				else if(pos != -1 && ctr[0] != pos)
					return;
				max = getMax(ctr);
			}
		}while(emat.singles.E[curI[0]][curI[1]].length==0);

		hasNextItem = true;
		nextItem = emat.singles.getTerm(ctr);
	}

	// Returns the max number of rotamers, amino acids, and positions, for the current "rotamer entry"
	// For example, if their are 4 positions, 6 amino acids at the current position and 9 rotamers at the current "amino acid",
	// then it returns [4, 6, 9].
	private int[] getMax(int[] ctr){
		int[] maxes = new int[ctr.length];
		maxes[0] = emat.singles.E.length;
		maxes[1] = emat.singles.E[ctr[0]].length;
		maxes[2] = emat.singles.E[ctr[0]][ctr[1]].length;

		return maxes;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

}

class PairsIterator implements Iterator<EMatrixEntryWIndex>{
	static final int POS1 = 0;
	static final int AA1 = 1;
	static final int ROT1 = 2;
	static final int POS2 = 3;
	static final int AA2 = 4;
	static final int ROT2 = 5;

	private final Emat emat;
	private boolean hasNextItem = false;
	private EMatrixEntry nextItem = null;
	private int[] curI;
	private int pos1 = -1;
	private int pos2 = -1;
	private TreeSet<Integer> AAs1 = null; //Allowed AAs for pos 1
	private TreeSet<Integer> AAs2 = null; //Allowed AAs for pos 2
	private TreeSet<Integer> Rots1 = null; //Allowed Rots for pos 1 //This only works if 1 AA1 is allowed
	//private Iterator<Integer> AAs1_iter=null;
	//private Iterator<Integer> AAs2_iter=null;

	public PairsIterator(Emat pairEnergyMat){
		emat = pairEnergyMat;
		curI = new int[6];

		for(int i=0; i<curI.length;i++)
			curI[i] = 0;

		//KER: set pos2 to 1 since 0-0 doesn't exist in pairs
		//curI[POS2] = 0;

		curI[ROT2] = -1;
		calcNext();

		/*try{
		if(emat.pairs.E != null && emat.pairs.E.length>1){
			hasNextItem = true;
			while(emat.pairs.E[curI[POS1]][curI[AA1]].length == 0){
				curI[AA1]++;
			}
			if(emat.pairs.E[curI[POS1]][curI[AA1]][curI[ROT1]] == null){
				hasNextItem = false;
				return;
			}
			while(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]][curI[4]].length == 0){
				curI[AA2]++;
			}
			nextItem = emat.pairs.getTerm(curI);
		}
		else{
			hasNextItem = false;
		}
}
catch(Exception E){
	System.out.println("DELETE ME");
}*/
	}

	/**
	 * Initializes the allowed AA treesets to everything
	 */
	/*private void initializeAAs() {
		AAs1 = new TreeSet<Integer>();
		for(int i=0; i<emat.singles.supRot[curI[POS1]].length;i++)
			AAs1.add(i);

		AAs1_iter = AAs1.iterator();

		AAs2 = new TreeSet<Integer>();
		for(int i=0; i<emat.singles.supRot[curI[POS2]].length;i++)
			AAs2.add(i);
		AAs2_iter = AAs2.iterator();

	}*/

	public PairsIterator(Emat pairEnergyMat,int p1, int p2){
		emat = pairEnergyMat;
		pos1 = p1;
		pos2 = p2;
		curI = new int[6];
		for(int i=0; i<curI.length;i++)
			curI[i] = 0;

		curI[POS1] = p1;
		curI[POS2] = p2;
		curI[ROT2] = -1;

		calcNext();
		
//		if(emat.pairs.E != null){
//			hasNextItem = true;
//			while(emat.pairs.E[curI[0]][curI[1]].length == 0){
//				curI[AA1]++;
//			}
//			if(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]] == null){
//				hasNextItem = false;
//			}
//			else{
//				while(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]][curI[4]].length == 0){
//					curI[AA2]++;
//				}
//				nextItem = emat.pairs.getTerm(curI, emat.singles.supRot);
//			}
//		}
//		else{
//			hasNextItem = false;
//		}
		
		

	}

	public PairsIterator(Emat pairEnergyMat,int p1, int p2,TreeSet<Integer> AAs1,TreeSet<Integer> AAs2){
		emat = pairEnergyMat;
		pos1 = p1;
		pos2 = p2;
		this.AAs1 = AAs1;
		this.AAs2 = AAs2;
		Iterator<Integer> AAs1_iter = AAs1.iterator();
		Iterator<Integer> AAs2_iter = AAs2.iterator();

		curI = new int[6];
		for(int i=0; i<curI.length;i++)
			curI[i] = 0;

		//Set Positions
		curI[0] = p1;
		curI[3] = p2;
		//Set AAs
		curI[1] = AAs1_iter.next();
		curI[4] = AAs2_iter.next();

		if(emat.pairs.E != null){
			hasNextItem = true;
			while(emat.pairs.E[curI[0]][curI[1]].length == 0){
				curI[AA1] = AAs1_iter.next();
			}
			if(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]] == null){
				hasNextItem = false;
			}
			else{
				while(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]][curI[4]].length == 0){
					curI[AA2] = AAs2_iter.next();
				}
				nextItem = emat.pairs.getTerm(curI, emat.singles.supRot);
			}
		}
		else{
			hasNextItem = false;
		}

	}

	//One rot vs. All pairs
	public PairsIterator(Emat pairEnergyMat,Index3 rot1){
		emat = pairEnergyMat;
		pos1 = rot1.pos;
		pos2 = 0;
		this.AAs1 = new TreeSet<Integer>();
		this.AAs1.add(rot1.aa);
		this.Rots1 = new TreeSet<Integer>();
		this.Rots1.add(rot1.rot);
		//this.AAs2 = null; 
		Iterator<Integer> AAs1_iter = AAs1.iterator();
		//Iterator<Integer> AAs2_iter = AAs2.iterator();
		Iterator<Integer> Rots1_iter = Rots1.iterator();

		curI = new int[6];
		/*for(int i=0; i<curI.length;i++)
			curI[i] = 0;*/

		//Set Positions
		curI[POS1] = rot1.pos;
		curI[POS2] = 0;
		//Set AAs
		curI[AA1] = AAs1_iter.next();
		curI[AA2] = 0;
		//Set Rots
		curI[ROT1] = Rots1_iter.next();
		curI[ROT2] = 0;

		if(emat.pairs.E != null){
			hasNextItem = true;
			while(emat.pairs.E[curI[0]][curI[1]].length == 0){
				curI[AA1] = AAs1_iter.next();
			}
			if(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]] == null){
				hasNextItem = false;
			}
			else{
				while(emat.pairs.E[curI[0]][curI[1]][curI[2]][curI[3]][curI[4]].length == 0){
					curI[AA2]++;
				}
				nextItem = emat.pairs.getTerm(curI, emat.singles.supRot);
			}
		}
		else{
			hasNextItem = false;
		}

	}

	public boolean hasNext() {
		return hasNextItem;
	}

	public EMatrixEntryWIndex next() {
		int[] index = new int[curI.length];
		for(int i=0; i<index.length;i++)
			index[i] = curI[i];
		EMatrixEntryWIndex ret = new EMatrixEntryWIndex(nextItem, index);
		calcNext();
		return ret;
	}

	/*public void calcNext() {
	    hasNextItem = false;

	    int[] ctr = {curI[0],curI[1],curI[2],curI[3],curI[4],curI[5]}; 
	    int[] max = getMax(curI);


	    do{
	    	if(Emat.incrementCtr(ctr,max)){
		    	//incrementCtr(ctr,max);
		    	//Gone past the end of the matrix
		    	if(ctr[0] == emat.pairE.length)
		    		return;
		    	else if(pos2 != -1 && ctr[3] != pos2){
		    		int ctr1[] = {ctr[0],ctr[1],ctr[2]};
		    		int max1[] = {max[0],max[1],max[2]};
		    		Emat.incrementCtr(ctr1, max1);

		    		if(ctr1[0] != pos1)
		    			return;

		    		ctr[0] = ctr1[0];ctr[1] = ctr1[1];ctr[2] = ctr1[2];
		    		ctr[3] = pos2; ctr[4] = 0; ctr[5] = 0;
		    	}
		    	max = getMax(ctr);
		    }
    	}while(ctr[0] == ctr[3] || emat.pairE[ctr[0]][ctr[1]].length==0 || emat.pairE[ctr[0]][ctr[1]][ctr[2]][ctr[3]][ctr[4]].length==0);


	    hasNextItem = true;
    	nextItem = emat.getPairE(ctr);
    	curI = ctr;

	  }*/

	//This function isn't efficient for mutation because it doesn't go through all rotamers before going through all aas
	/*public void calcNext() {
		    hasNextItem = false;

		    int[] ctr1 = {curI[0],curI[1],curI[2],curI[3],curI[4],curI[5]+1}; 
		    //int[] max = getMax(curI);


		    for(int p1=ctr1[0]; p1<emat.pairE.length;p1++){
		    	if(emat.pairE[p1] != null && (pos1 == -1 || p1==pos1)){
		    	for(int a1=ctr1[1];a1<emat.pairE[p1].length;a1++){
		    		if(emat.pairE[p1][a1] != null){
		    		for(int r1=ctr1[2];r1<emat.pairE[p1][a1].length;r1++){
		    			if(emat.pairE[p1][a1][r1] != null){
		    			for(int p2=ctr1[3]; p2<emat.pairE[p1][a1][r1].length;p2++){
		    				if(emat.pairE[p1][a1][r1][p2] != null && (pos2 == -1 || p2==pos2) && (p1 != p2)){
		    		    	for(int a2=ctr1[4];a2<emat.pairE[p1][a1][r1][p2].length;a2++){
		    		    		if(emat.pairE[p1][a1][r1][p2][a2] != null){
		    		    		for(int r2=ctr1[5];r2<emat.pairE[p1][a1][r1][p2][a2].length;r2++){
		    		    			hasNextItem = true;
		    		    			int[] ctr = {p1,a1,r1,p2,a2,r2};
		    		    	    	nextItem = emat.getPairE(ctr);
		    		    	    	curI = ctr;
		    		    	    	return;
		    		    		} ctr1[5] = 0;
		    		    		}
    		    			} ctr1[4] = 0;
		    				}
		    			}ctr1[3] = 0;
		    			}
		    		}ctr1[2] = 0;
		    		}
		    	}ctr1[1] = 0;
		    	}
		    }ctr1[0] = 0;





		  }*/

	//Goes through array in order of amino acid to minimize number of mutations  
	public void calcNext() {
		hasNextItem = false;

		int[] ctr1 = {curI[0],curI[1],curI[2],curI[3],curI[4],curI[5]+1}; 
		//int[] max = getMax(curI);


		for(int p1=ctr1[0]; p1<emat.pairs.E.length;p1++){
			if(emat.pairs.E[p1] != null && (pos1 == -1 || p1==pos1)){
				for(int a1=ctr1[1];a1<emat.singles.E[p1].length;a1++){
					if(emat.pairs.E[p1][a1]!=null && (AAs1 == null || AAs1.contains(a1))){
						for(int p2=ctr1[3]; p2<emat.singles.E.length;p2++){
							if((pos2 == -1 || p2==pos2) && (p1 != p2)){
								for(int a2=ctr1[4];a2<emat.singles.E[p2].length;a2++){
									if(AAs2 == null || AAs2.contains(a2)){
										for(int r1=ctr1[2];r1<emat.pairs.E[p1][a1].length;r1++){
											if(Rots1 == null || Rots1.contains(r1)){
												if(emat.pairs.E[p1][a1][r1] != null && emat.pairs.E[p1][a1][r1][p2] != null && emat.pairs.E[p1][a1][r1][p2][a2] != null){
													for(int r2=ctr1[5];r2<emat.pairs.E[p1][a1][r1][p2][a2].length;r2++){
														hasNextItem = true;
														int[] ctr = {p1,a1,r1,p2,a2,r2};
														nextItem = emat.pairs.getTerm(ctr, emat.singles.supRot);
														curI = ctr;
														return;
													} ctr1[5] = 0;
												}}
										} ctr1[2] = 0;
									}
								} ctr1[4] = 0;
							}
						}ctr1[3] = 0;
					}
				}ctr1[1] = 0;
			}
		}ctr1[0] = 0;
	}

	/*public void calcNextExp() {
		hasNextItem = false;

		Iterator<Integer> Pos1_iter = AAs1.iterator();
		Iterator<Integer> Pos2_iter = AAs1.iterator();
		Iterator<Integer> AA1_iter = AAs1.iterator();
		Iterator<Integer> AA2_iter = AAs1.iterator();
		Iterator<Integer> Rot1_iter = AAs1.iterator();
		Iterator<Integer> Rot2_iter = AAs1.iterator();

		TreeSet<Integer> Pos1 = AAs1;
		TreeSet<Integer> Pos2 = AAs1;
		TreeSet<Integer> AA1 = AAs1;
		TreeSet<Integer> AA2 = AAs1;
		TreeSet<Integer> Rot1 = AAs1;
		TreeSet<Integer> Rot2 = AAs1;


		//int[] ctr1 = {curI[0],curI[1],curI[2],curI[3],curI[4],curI[5]};

		while(Pos1_iter.hasNext()){
			while(Pos2_iter.hasNext()){
				while(AA1_iter.hasNext()){
					while(AA2_iter.hasNext()){
						while(Rot1_iter.hasNext()){
							while(Rot2_iter.hasNext()){
								curI[ROT2] = Rot2_iter.next();
								hasNextItem = true;
								nextItem = emat.pairs.getTerm(curI);
							}
							//Reset Rot2 iterator
							Rot2_iter = Rot2.iterator();
							//Increment Rot1
							curI[ROT1] = Rot1_iter.next();
						}
						//Reset Rot1 iterator
						Rot1_iter = Rot1.iterator();
						//Increment AA2
					}
				}
			}
		}



		for(int p1=ctr1[0]; p1<emat.pairs.E.length;p1++){
			if(emat.pairs.E[p1] != null && (pos1 == -1 || p1==pos1)){
				for(int a1=ctr1[1];a1<emat.singles.E[p1].length;a1++){
					if(emat.pairs.E[p1][a1]!=null){
						for(int p2=ctr1[3]; p2<emat.singles.E.length;p2++){
							if((pos2 == -1 || p2==pos2) && (p1 != p2)){
								for(int a2=ctr1[4];a2<emat.singles.E[p2].length;a2++){
									for(int r1=ctr1[2];r1<emat.pairs.E[p1][a1].length;r1++){
										if(emat.pairs.E[p1][a1][r1] != null && emat.pairs.E[p1][a1][r1][p2] != null && emat.pairs.E[p1][a1][r1][p2][a2] != null){
											for(int r2=ctr1[5];r2<emat.pairs.E[p1][a1][r1][p2][a2].length;r2++){
												hasNextItem = true;
												int[] ctr = {p1,a1,r1,p2,a2,r2};
												nextItem = emat.pairs.getTerm(ctr);
												curI = ctr;
												return;
											} ctr1[5] = 0;
										}
									} ctr1[2] = 0;
								} ctr1[4] = 0;
							}
						}ctr1[3] = 0;
					}
				}ctr1[1] = 0;
			}
		}ctr1[0] = 0;
	}*/

	private int[] getMax(int[] ctr){
		int[] maxes = new int[ctr.length];
		if(emat.pairs.E == null)
			return maxes;
		maxes[0] = emat.pairs.E.length;
		if(emat.pairs.E.length == 0 || emat.pairs.E[ctr[0]]==null)
			return maxes;
		maxes[1] = emat.pairs.E[ctr[0]].length;
		if(emat.pairs.E[ctr[0]].length == 0 || emat.pairs.E[ctr[0]][ctr[1]]==null)
			return maxes;
		maxes[2] = emat.pairs.E[ctr[0]][ctr[1]].length;
		if(emat.pairs.E[ctr[0]][ctr[1]].length==0 || emat.pairs.E[ctr[0]][ctr[1]][ctr[2]]==null)
			return maxes;
		maxes[3] = emat.pairs.E[ctr[0]][ctr[1]][ctr[2]].length;
		if(emat.pairs.E[ctr[0]][ctr[1]][ctr[2]].length==0 || emat.pairs.E[ctr[0]][ctr[1]][ctr[2]][ctr[3]]==null)
			return maxes;
		maxes[4] = emat.pairs.E[ctr[0]][ctr[1]][ctr[2]][ctr[3]].length;
		if(emat.pairs.E[ctr[0]][ctr[1]][ctr[2]][ctr[3]].length==0 || emat.pairs.E[ctr[0]][ctr[1]][ctr[2]][ctr[3]][ctr[4]]==null)
			return maxes;
		maxes[5] = emat.pairs.E[ctr[0]][ctr[1]][ctr[2]][ctr[3]][ctr[4]].length; 

		return maxes;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

	//This isn't quite correct because the max should be updated if the curAA changes,
	//but this works due to the matrix structure that if the curAA is defined rotamer 0 will be defined.
	void incrementCtr(int[] ctr, int[] max) {
		boolean[] changedPos = {false,false,false,false,false,false};
		ctr[ctr.length-1]++;
		for(int i=ctr.length-1;i>0; i--){
			if(ctr[i] >= max[i]){
				ctr[i] = 0;
				ctr[i-1]++;
			}
		}
	}

}
