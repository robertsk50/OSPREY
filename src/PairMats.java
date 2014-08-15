import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.Serializable;


public class PairMats implements Serializable {

	
	double[][][][][][] E;
	boolean[][][][][][] pruned;
	
	//These last three are only present if doDih = true
	double[][][][][][] maxE;
	//First 6 dim are indices, last dim are rotamers
	double[][][][][][][] rotDih1;
	double[][][][][][][] rotDih2;
	
	
	
	boolean doDih;
	
	
	public PairMats(boolean doDih){
		this.doDih = doDih;
	}
	
	public PairMats getSinglePosMat(int pos){
		PairMats newPairs = new PairMats(doDih);
		
		newPairs.E = new double[E.length][][][][][];
		
		newPairs.pruned = new boolean[E.length][][][][][];
		if(doDih){
			newPairs.maxE = new double[E.length][][][][][];
			newPairs.rotDih1 = new double[E.length][][][][][][];
			newPairs.rotDih2 = new double[E.length][][][][][][];
		}
		//newPairs.supRot2 = new int[E.length][][][][][][];
		
		for(int i=0;i<E.length;i++){
			if(i == pos){
				//intraE[i] = emat.intraE[i];
				newPairs.E[i] = E[i];
				newPairs.pruned[i] = pruned[i];
				if(doDih){
					newPairs.maxE[i] = maxE[i];
					newPairs.rotDih1[i] = rotDih1[i];
					newPairs.rotDih2[i] = rotDih2[i];
				}
				//newPairs.supRot2[i] = supRot2[i];
			}
		}
		
		return newPairs;
	}
	
	public PairMats getDualPosMat(int pos1, int pos2){
		PairMats newPairs = new PairMats(doDih);
		
		newPairs.E = new double[E.length][][][][][];
		
		newPairs.pruned = new boolean[E.length][][][][][];
		if(doDih){
			newPairs.maxE = new double[E.length][][][][][];
			newPairs.rotDih1 = new double[E.length][][][][][][];
			newPairs.rotDih2 = new double[E.length][][][][][][];
		}
		//newPairs.supRot2 = new int[E.length][][][][][][];
		
		for(int i=0;i<E.length;i++){
			if(i == pos1 || i == pos2 ){
				//intraE[i] = emat.intraE[i];
				newPairs.E[i] = E[i];
				newPairs.pruned[i] = pruned[i];
				if(doDih){
					newPairs.maxE[i] = maxE[i];
					newPairs.rotDih1[i] = rotDih1[i];
					newPairs.rotDih2[i] = rotDih2[i];
				}
				//newPairs.supRot2[i] = supRot2[i];
			}
		}
		
		return newPairs;
	}
	
	
	
	public PairMats getOnlyDualPosMat(int pos1,int pos2){
		PairMats newPairs = new PairMats(doDih);
		
		newPairs.E = new double[E.length][][][][][];
		newPairs.pruned = new boolean[E.length][][][][][];
		if(doDih){
			newPairs.maxE = new double[E.length][][][][][];
			newPairs.rotDih1 = new double[E.length][][][][][][];
			newPairs.rotDih2 = new double[E.length][][][][][][];
		}
		//newPairs.supRot2 = new int[E.length][][][][][][];
		
		for(int p1=0;p1<E.length;p1++){
			if(p1 == pos1){
				
				newPairs.E[p1] = new double[E[p1].length][][][][];
				newPairs.pruned[p1] = new boolean[E[p1].length][][][][];
				if(doDih){
					newPairs.maxE[p1] = new double[E[p1].length][][][][];
					newPairs.rotDih1[p1] = new double[E[p1].length][][][][][];
					newPairs.rotDih2[p1] = new double[E[p1].length][][][][][];
				}
				for(int a1=0; a1<E[p1].length;a1++){
					newPairs.E[p1][a1] = new double[E[p1][a1].length][][][];
					newPairs.pruned[p1][a1] = new boolean[E[p1][a1].length][][][];
					if(doDih){
						newPairs.maxE[p1][a1] = new double[E[p1][a1].length][][][];
						newPairs.rotDih1[p1][a1] = new double[E[p1][a1].length][][][][];
						newPairs.rotDih2[p1][a1] = new double[E[p1][a1].length][][][][];
					}
					for(int r1=0; r1<E[p1][a1].length;r1++){
						newPairs.E[p1][a1][r1] = new double[E[p1][a1][r1].length][][];
						newPairs.pruned[p1][a1][r1] = new boolean[E[p1][a1][r1].length][][];
						if(doDih){
							newPairs.maxE[p1][a1][r1] = new double[E[p1][a1][r1].length][][];
							newPairs.rotDih1[p1][a1][r1] = new double[E[p1][a1][r1].length][][][];
							newPairs.rotDih2[p1][a1][r1] = new double[E[p1][a1][r1].length][][][];
						}
						
						for(int p2=0; p2<E[p1][a1][r1].length;p2++){
							if(p2 == pos2){
								newPairs.E[p1][a1][r1][p2] = E[p1][a1][r1][p2];
								newPairs.pruned[p1][a1][r1][p2] = pruned[p1][a1][r1][p2];
								if(doDih){
									newPairs.maxE[p1][a1][r1][p2] = maxE[p1][a1][r1][p2];
									newPairs.rotDih1[p1][a1][r1][p2] = rotDih1[p1][a1][r1][p2];
									newPairs.rotDih2[p1][a1][r1][p2] = rotDih2[p1][a1][r1][p2];
								}
							}
						}
					}
				}
			}
		}
		
		return newPairs;
	}
	
	
	public void setDihedrals(int[] i, double[][] diheds){
		this.rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = diheds[0];
		this.rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = diheds[1];
		//Also set the symmetric E
		//this.E[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = E;
	}
	
	public void setE(int[] i, double E){
		this.E[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = E;
		//Also set the symmetric E
		//this.E[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = E;
	}
	
	public void setMaxE(int[] i, double E){
		this.maxE[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = E;
		//Also set the symmetric E
		//this.E[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = E;
	}
	
	public void addE(int[] i, double E){
		this.E[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] += E;
		//Also set the symmetric E
		//this.E[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = E;
	}
	
	public void setE(EMatrixEntrySlim eme){
		this.E[eme.index[0]][eme.index[1]][eme.index[2]][eme.index[3]][eme.index[4]][eme.index[5]] = eme.minE;
		//Also set the symmetric E
		this.E[eme.index[3]][eme.index[4]][eme.index[5]][eme.index[0]][eme.index[1]][eme.index[2]] = eme.minE;
	}
	public void setMaxE(EMatrixEntrySlim eme){
		this.maxE[eme.index[0]][eme.index[1]][eme.index[2]][eme.index[3]][eme.index[4]][eme.index[5]] = eme.maxE;
		//Also set the symmetric E
		this.maxE[eme.index[3]][eme.index[4]][eme.index[5]][eme.index[0]][eme.index[1]][eme.index[2]] = eme.maxE;
	}
	
	public void setDihed(EMatrixEntrySlim eme){
		this.rotDih1[eme.index[0]][eme.index[1]][eme.index[2]][eme.index[3]][eme.index[4]][eme.index[5]] = eme.rotDih1;
		//Also set the symmetric E
		this.rotDih1[eme.index[3]][eme.index[4]][eme.index[5]][eme.index[0]][eme.index[1]][eme.index[2]] = eme.rotDih2;
		
		this.rotDih2[eme.index[0]][eme.index[1]][eme.index[2]][eme.index[3]][eme.index[4]][eme.index[5]] = eme.rotDih2;
		//Also set the symmetric E
		this.rotDih2[eme.index[3]][eme.index[4]][eme.index[5]][eme.index[0]][eme.index[1]][eme.index[2]] = eme.rotDih1;
	}
	
	public void setENoSym(EMatrixEntrySlim eme){
		this.E[eme.index[0]][eme.index[1]][eme.index[2]][eme.index[3]][eme.index[4]][eme.index[5]] = eme.minE;
			
	}
	
	public void setAll(int p1, int a1, int r1, int p2, int a2, int r2, double E, boolean pruned){
		this.E[p1][a1][r1][p2][a2][r2] = E;
		this.pruned[p1][a1][r1][p2][a2][r2] = pruned;
	}
	
	public void setAll(int p1, int a1, int r1, int p2, int a2, int r2, RotamerPairEntry rpe){
		this.E[p1][a1][r1][p2][a2][r2] = rpe.minE();
		this.pruned[p1][a1][r1][p2][a2][r2] = rpe.isPruned();
	}
	
	public void addDim(int[] i, int length){
		int dim = i.length;
		switch(dim){
			case 0:
				E = new double[length][][][][][];
				pruned = new boolean[length][][][][][];
				if(doDih){
					maxE = new double[length][][][][][];
					rotDih1 = new double[length][][][][][][];
					rotDih2 = new double[length][][][][][][];
				}
				//supRot2 = new int[length][][][][][][];
				break;
			case 1:
				E[i[0]] = new double[length][][][][];
				
				pruned[i[0]] = new boolean[length][][][][];
				if(doDih){
					maxE[i[0]] = new double[length][][][][];
				rotDih1[i[0]] = new double[length][][][][][];
				rotDih2[i[0]] = new double[length][][][][][];
				}
				//supRot2[i[0]] = new int[length][][][][][];
				break;
			case 2:
				E[i[0]][i[1]] = new double[length][][][];
				
				pruned[i[0]][i[1]] = new boolean[length][][][];
				if(doDih){
					maxE[i[0]][i[1]] = new double[length][][][];
				rotDih1[i[0]][i[1]] = new double[length][][][][];
				rotDih2[i[0]][i[1]] = new double[length][][][][];
				}
//				supRot2[i[0]][i[1]] = new int[length][][][][];
				break;
			case 3:
				E[i[0]][i[1]][i[2]] = new double[length][][];
				pruned[i[0]][i[1]][i[2]] = new boolean[length][][];
				if(doDih){
					maxE[i[0]][i[1]][i[2]] = new double[length][][];
				rotDih1[i[0]][i[1]][i[2]] = new double[length][][][];
				rotDih2[i[0]][i[1]][i[2]] = new double[length][][][];
				}
//				supRot2[i[0]][i[1]][i[2]] = new int[length][][][];
				break;
			case 4:
				E[i[0]][i[1]][i[2]][i[3]] = new double[length][];
				pruned[i[0]][i[1]][i[2]][i[3]] = new boolean[length][];
				if(doDih){
					maxE[i[0]][i[1]][i[2]][i[3]] = new double[length][];
				rotDih1[i[0]][i[1]][i[2]][i[3]] = new double[length][][];
				rotDih2[i[0]][i[1]][i[2]][i[3]] = new double[length][][];
				}
//				supRot2[i[0]][i[1]][i[2]][i[3]] = new int[length][][];
				break;
			case 5:
				E[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[length];
				pruned[i[0]][i[1]][i[2]][i[3]][i[4]] = new boolean[length];
				if(doDih){
					maxE[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[length];
				rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[length][];
				rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[length][];
				}
//				supRot2[i[0]][i[1]][i[2]][i[3]][i[4]] = new int[length][];
				break;
			default:
				System.out.println("Index too high for pair mats.");
				break;
		}
	}
		
	
	
	public void copy(int[] newI, PairMats pairs, int[] oldI){
		this.E[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = pairs.E[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]];
		this.pruned[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = pairs.pruned[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]];
		
		if(doDih){
			this.maxE[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = pairs.maxE[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]];
		this.rotDih1[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = new double[pairs.rotDih1[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length];
		System.arraycopy(pairs.rotDih1[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]], 0, this.rotDih1[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]], 0, pairs.rotDih1[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length);
		
		this.rotDih2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = new double[pairs.rotDih2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length];
		System.arraycopy(pairs.rotDih2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]], 0, this.rotDih2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]], 0, pairs.rotDih2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length);
		}
//		
//		this.supRot2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]] = new int[pairs.supRot2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length];
//		System.arraycopy(pairs.supRot2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]], 0, this.supRot2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]][newI[5]], 0, pairs.supRot2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]][oldI[5]].length);
	}
	
	public void copyRow(int[] newI, PairMats pairs,
			int[] oldI) {
		System.arraycopy(pairs.E[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.E[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.E[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
		System.arraycopy(pairs.pruned[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.pruned[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.pruned[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
		if(doDih){
			System.arraycopy(pairs.maxE[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.maxE[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.maxE[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
		System.arraycopy(pairs.rotDih1[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.rotDih1[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.rotDih1[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
		System.arraycopy(pairs.rotDih2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.rotDih2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.rotDih2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
		}
//		System.arraycopy(pairs.supRot2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]], 0, this.supRot2[newI[0]][newI[1]][newI[2]][newI[3]][newI[4]], 0, pairs.supRot2[oldI[0]][oldI[1]][oldI[2]][oldI[3]][oldI[4]].length);
	}


//	public void setSupRot(int p1, int a1, int r1, int p2, int a2, int r2,
//			int[] supRot1, int[] supRot2) {
//		this.supRot1[p1][a1][r1][p2][a2][r2] = new int[supRot1.length];
//		System.arraycopy(supRot1, 0, this.supRot1[p1][a1][r1][p2][a2][r2], 0, supRot1.length);
//		
//		this.supRot2[p1][a1][r1][p2][a2][r2] = new int[supRot2.length];
//		System.arraycopy(supRot2, 0, this.supRot2[p1][a1][r1][p2][a2][r2], 0, supRot2.length);
//		
//	}


	public EMatrixEntry getE(int[] i, int[][][][] globalRots) {
		
//		SuperRotamer sr1 = new SuperRotamer(supRot1[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]);
//		SuperRotamer sr2 = new SuperRotamer(supRot2[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]);
		
		SuperRotamer sr1 = new SuperRotamer(globalRots[i[0]][i[1]][i[2]]);
		SuperRotamer sr2 = new SuperRotamer(globalRots[i[3]][i[4]][i[5]]);
		
		if(!doDih)
		return new RotamerPairEntry(i[0], sr1, i[3], sr2, E[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]], pruned[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]);
		
		return new RotamerPairEntry(i[0], sr1, i[3], sr2, E[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]], maxE[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]],
				pruned[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]],rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]],rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]);
		
	}


	public RotamerPairEntry getTerm(int p1, int a1, int r1, int p2,	int a2, int r2, int[][][][] globalRots) {
		
//		SuperRotamer sr1 = new SuperRotamer(supRot1[p1][a1][r1][p2][a2][r2]);
//		SuperRotamer sr2 = new SuperRotamer(supRot2[p1][a1][r1][p2][a2][r2]);
		
		SuperRotamer sr1 = new SuperRotamer(globalRots[p1][a1][r1]);
		SuperRotamer sr2 = new SuperRotamer(globalRots[p2][a2][r2]);
		
		if(!doDih)
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], pruned[p1][a1][r1][p2][a2][r2]);
		
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], maxE[p1][a1][r1][p2][a2][r2],
				pruned[p1][a1][r1][p2][a2][r2],rotDih1[p1][a1][r1][p2][a2][r2],rotDih2[p1][a1][r1][p2][a2][r2]);
		
	}
	
	public RotamerPairEntry getTerm(int[] i, int[][][][] globalRots) {
		
		int p1 = i[0];
		int a1 = i[1];
		int r1 = i[2];
		int p2 = i[3];
		int a2 = i[4];
		int r2 = i[5];
		
//		SuperRotamer sr1 = new SuperRotamer(supRot1[p1][a1][r1][p2][a2][r2]);
//		SuperRotamer sr2 = new SuperRotamer(supRot2[p1][a1][r1][p2][a2][r2]);
		
		SuperRotamer sr1 = new SuperRotamer(globalRots[p1][a1][r1]);
		SuperRotamer sr2 = new SuperRotamer(globalRots[p2][a2][r2]);
		
		if(!doDih)
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], pruned[p1][a1][r1][p2][a2][r2]);
		
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], maxE[p1][a1][r1][p2][a2][r2],
				pruned[p1][a1][r1][p2][a2][r2],rotDih1[p1][a1][r1][p2][a2][r2],rotDih2[p1][a1][r1][p2][a2][r2]);
		
	}
	
	public RotamerPairEntry getTerm(int[] i1, int[] i2,int[][][][] globalRots) {
		
		int p1 = i1[0];
		int a1 = i1[1];
		int r1 = i1[2];
		int p2 = i2[0];
		int a2 = i2[1];
		int r2 = i2[2];
		
//		SuperRotamer sr1 = new SuperRotamer(supRot1[p1][a1][r1][p2][a2][r2]);
//		SuperRotamer sr2 = new SuperRotamer(supRot2[p1][a1][r1][p2][a2][r2]);
		
		SuperRotamer sr1 = new SuperRotamer(globalRots[p1][a1][r1]);
		SuperRotamer sr2 = new SuperRotamer(globalRots[p2][a2][r2]);
		
		if(!doDih)
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], pruned[p1][a1][r1][p2][a2][r2]);
		
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], maxE[p1][a1][r1][p2][a2][r2],
				pruned[p1][a1][r1][p2][a2][r2],rotDih1[p1][a1][r1][p2][a2][r2],rotDih2[p1][a1][r1][p2][a2][r2]);
		
	}
	
	public RotamerPairEntry getTerm(Index3 i1, Index3 i2,int[][][][] globalRots) {
		
		int p1 = i1.pos;
		int a1 = i1.aa;
		int r1 = i1.rot;
		int p2 = i2.pos;
		int a2 = i2.aa;
		int r2 = i2.rot;
		
//		SuperRotamer sr1 = new SuperRotamer(supRot1[p1][a1][r1][p2][a2][r2]);
//		SuperRotamer sr2 = new SuperRotamer(supRot2[p1][a1][r1][p2][a2][r2]);
		
		SuperRotamer sr1 = new SuperRotamer(globalRots[p1][a1][r1]);
		SuperRotamer sr2 = new SuperRotamer(globalRots[p2][a2][r2]);
		
		if(!doDih)
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], pruned[p1][a1][r1][p2][a2][r2]);
		
		return new RotamerPairEntry(p1, sr1, p2, sr2, E[p1][a1][r1][p2][a2][r2], maxE[p1][a1][r1][p2][a2][r2], 
				pruned[p1][a1][r1][p2][a2][r2],rotDih1[p1][a1][r1][p2][a2][r2],rotDih2[p1][a1][r1][p2][a2][r2]);
		
	}

	/*
	 * Shrink the dim=6 of the pairs matrix
	 */
	public void shrinkDim6(int[] ind, boolean[] tmpPruned) {
		
		int numLeft = 0;
		for(int i=0; i<tmpPruned.length;i++)
			if(!tmpPruned[i])
				numLeft++;
		
		double[] newE = new double[numLeft];
		double[] newMaxE = new double[numLeft];
		boolean[] newPruned = new boolean[numLeft];
		//First 6 dim are indices, last dim are rotamers
		double[][] newRotDih1 = null;
		double[][] newRotDih2 = null;
		if(doDih){
			newRotDih1 = new double[numLeft][];
			newRotDih2 = new double[numLeft][];
		}
		
		
		int ctr=0;
		for(int i=0; i<tmpPruned.length;i++){
			if(!tmpPruned[i]){
				newE[ctr] = E[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
				
				newPruned[ctr] = pruned[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
				if(doDih){
					newMaxE[ctr] = maxE[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
					newRotDih1[ctr] = rotDih1[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
					newRotDih2[ctr] = rotDih2[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
				}
//				newSupRot2[ctr] = supRot2[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][i];
				ctr++;
			}
				
		}
				
		E[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newE;
		
		pruned[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newPruned;
		if(doDih){
			maxE[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newMaxE;
			rotDih1[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newRotDih1;
			rotDih2[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newRotDih2;
		}
//		supRot2[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = newSupRot2;
		
		
	}


	public void shrinkDim3(int[] ind, boolean[] tmpPruned) {
		int numLeft = 0;
		for(int i=0; i<tmpPruned.length;i++)
			if(!tmpPruned[i])
				numLeft++;
		
		double[][][][] newE = new double[numLeft][][][];
		double[][][][] newMaxE = new double[numLeft][][][];
		boolean[][][][] newPruned = new boolean[numLeft][][][];
		//First 6 dim are indices, last dim are rotamers
		double[][][][][] newRotDih1 = null;
		double[][][][][] newRotDih2 = null;
		if(doDih){
		newRotDih1 = new double[numLeft][][][][];
		newRotDih2 = new double[numLeft][][][][];
		}
		//int[][][][][] newSupRot2 = new int[numLeft][][][][];
		
		
		int ctr=0;
		for(int i=0; i<tmpPruned.length;i++){
			if(!tmpPruned[i]){
				newE[ctr] = E[ind[0]][ind[1]][i];
				
				newPruned[ctr] = pruned[ind[0]][ind[1]][i];
				if(doDih){
					newMaxE[ctr] = maxE[ind[0]][ind[1]][i];
				newRotDih1[ctr] = rotDih1[ind[0]][ind[1]][i];
				newRotDih2[ctr] = rotDih2[ind[0]][ind[1]][i];
				}
//				newSupRot2[ctr] = supRot2[ind[0]][ind[1]][i];
				ctr++;	
			}
			
		}
				
		E[ind[0]][ind[1]] = newE;
		
		pruned[ind[0]][ind[1]] = newPruned;
		if(doDih){
			maxE[ind[0]][ind[1]] = newMaxE;
		rotDih1[ind[0]][ind[1]] = newRotDih1;
		rotDih2[ind[0]][ind[1]] = newRotDih2;
		}
//		supRot2[ind[0]][ind[1]] = newSupRot2;
		
	}


	public void write(String fileName) {
		KSParser.outputObject(E,fileName+".pairsE");
		
		KSParser.outputObject(pruned,fileName+".pairsPruned");
		//First 6 dim are indices, last dim are rotamers
//		KSParser.outputObject(supRot1,fileName+".pairsSupRot1");
//		KSParser.outputObject(supRot2,fileName+".pairsSupRot2");
		
		if(doDih){
			KSParser.outputObject(maxE,fileName+".pairsMaxE");
			KSParser.outputObject(rotDih1,fileName+".pairsRotDih1");
			KSParser.outputObject(rotDih2,fileName+".pairsRotDih2");
		}
	}

	public static PairMats read(String fileName, boolean doDih){
		
		PairMats pm = new PairMats(doDih);
		
		try{
//			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName+".pairsE"));
			pm.E = (double[][][][][][])KSParser.loadObject(fileName+".pairsE");//in.readObject();
//			in.close();
//			in = new ObjectInputStream(new FileInputStream(fileName+".pairsPruned"));
			pm.pruned = (boolean[][][][][][])KSParser.loadObject(fileName+".pairsPruned");//in.readObject();
//			in.close();
//			in = new ObjectInputStream(new FileInputStream(fileName+".pairsSupRot1"));
//			pm.supRot1 = (int[][][][][][][])in.readObject();
//			in.close();
//			in = new ObjectInputStream(new FileInputStream(fileName+".pairsSupRot2"));
//			pm.supRot2 = (int[][][][][][][])in.readObject();
//			in.close();
			if(doDih){
//				in = new ObjectInputStream(new FileInputStream(fileName+".pairsMaxE"));
				pm.maxE = (double[][][][][][])KSParser.loadObject(fileName+".pairsMaxE");//in.readObject();
//				in.close();
//			in = new ObjectInputStream(new FileInputStream(fileName+".pairsRotDih1"));
			pm.rotDih1 = (double[][][][][][][])KSParser.loadObject(fileName+".pairsRotDih1");//in.readObject();
//			in.close();
//			in = new ObjectInputStream(new FileInputStream(fileName+".pairsRotDih2"));
			pm.rotDih2 = (double[][][][][][][])KSParser.loadObject(fileName+".pairsRotDih2");//in.readObject();
//			in.close();
		}
		}
		catch (Exception e){
			//e.printStackTrace();
			System.out.println("Couldn't load/read file: "+fileName);
			return null;
		}		
		
		return pm;
		
		
		
	}


	 public void extendDim(int[] i, int toExtendBy){
	        int dim = i.length;
	        switch(dim){
	            /*case 0:
	                E = new double[length][][][][][];
	                pruned = new boolean[length][][][][][];
	                supRot1 = new int[length][][][][][][];
	                supRot2 = new int[length][][][][][][];
	                break;
	            case 1:
	                E[i[0]] = new double[length][][][][];
	                pruned[i[0]] = new boolean[length][][][][];
	                supRot1[i[0]] = new int[length][][][][][];
	                supRot2[i[0]] = new int[length][][][][][];
	                break;*/
	            case 2:
	                double[][][][] tmpE = E[i[0]][i[1]];
	                boolean[][][][] tmpPruned = pruned[i[0]][i[1]];
	                //int[][][][][] tmpSupRot1 = supRot1[i[0]][i[1]];
	                //int[][][][][] tmpSupRot2 = supRot2[i[0]][i[1]];
	                E[i[0]][i[1]] = new double[E[i[0]][i[1]].length+toExtendBy][][][];
	                pruned[i[0]][i[1]] = new boolean[E[i[0]][i[1]].length+toExtendBy][][][];
	                if(doDih){
	                	double[][][][] tmpMaxE = maxE[i[0]][i[1]];
		                maxE[i[0]][i[1]] = new double[maxE[i[0]][i[1]].length+toExtendBy][][][];
		                System.arraycopy(tmpMaxE, 0, maxE[i[0]][i[1]], 0, tmpMaxE.length);
		                
	                	double[][][][][] tmpRotDih1 = rotDih1[i[0]][i[1]];
		                double[][][][][] tmpRotDih2 = rotDih2[i[0]][i[1]];
		                rotDih1[i[0]][i[1]] = new double[E[i[0]][i[1]].length+toExtendBy][][][][];
		                rotDih2[i[0]][i[1]] = new double[E[i[0]][i[1]].length+toExtendBy][][][][];
		                System.arraycopy(tmpRotDih1, 0, rotDih1[i[0]][i[1]], 0, tmpRotDih1.length);
		                System.arraycopy(tmpRotDih2, 0, rotDih2[i[0]][i[1]], 0, tmpRotDih2.length);
	                }
	                //supRot2[i[0]][i[1]] = new int[E[i[0]][i[1]].length+toExtendBy][][][][];

	                System.arraycopy(tmpE, 0, E[i[0]][i[1]], 0, tmpE.length);
	                System.arraycopy(tmpPruned, 0, pruned[i[0]][i[1]], 0, tmpPruned.length);
	                //System.arraycopy(tmpSupRot1, 0, supRot1[i[0]][i[1]], 0, tmpSupRot1.length);
	                //System.arraycopy(tmpSupRot2, 0, supRot2[i[0]][i[1]], 0, tmpSupRot2.length);

	                break;
	            /*case 3:
	                E[i[0]][i[1]][i[2]] = new double[length][][];
	                pruned[i[0]][i[1]][i[2]] = new boolean[length][][];
	                supRot1[i[0]][i[1]][i[2]] = new int[length][][][];
	                supRot2[i[0]][i[1]][i[2]] = new int[length][][][];
	                break;
	            case 4:
	                E[i[0]][i[1]][i[2]][i[3]] = new double[length][];
	                pruned[i[0]][i[1]][i[2]][i[3]] = new boolean[length][];
	                supRot1[i[0]][i[1]][i[2]][i[3]] = new int[length][][];
	                supRot2[i[0]][i[1]][i[2]][i[3]] = new int[length][][];
	                break;*/
	            case 5:
	                double[] tmpE_5 = E[i[0]][i[1]][i[2]][i[3]][i[4]];
	                boolean[] tmpPruned_5 = pruned[i[0]][i[1]][i[2]][i[3]][i[4]];
	                //int[][] tmpSupRot1_5 = supRot1[i[0]][i[1]][i[2]][i[3]][i[4]];
	                //int[][] tmpSupRot2_5 = supRot2[i[0]][i[1]][i[2]][i[3]][i[4]];

	                E[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy];
	                pruned[i[0]][i[1]][i[2]][i[3]][i[4]] = new boolean[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy];
	                //supRot1[i[0]][i[1]][i[2]][i[3]][i[4]] = new int[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy][];
	                //supRot2[i[0]][i[1]][i[2]][i[3]][i[4]] = new int[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy][];

	                System.arraycopy(tmpE_5, 0, E[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpE_5.length);
	                System.arraycopy(tmpPruned_5, 0, pruned[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpPruned_5.length);
	                //System.arraycopy(tmpSupRot1_5, 0, supRot1[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpSupRot1_5.length);
	                //System.arraycopy(tmpSupRot2_5, 0, supRot2[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpSupRot2_5.length);
	                if(doDih){
	                	double[] tmpMaxE_5 = maxE[i[0]][i[1]][i[2]][i[3]][i[4]];
		                maxE[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[maxE[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy];
		                System.arraycopy(tmpE_5, 0, maxE[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpMaxE_5.length);
	                	
	                	double[][] tmpRotDih1_5 = rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]];
		                double[][] tmpRotDih2_5 = rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]];
	                	rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy][];
		                rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]] = new double[E[i[0]][i[1]][i[2]][i[3]][i[4]].length+toExtendBy][];
	                	System.arraycopy(tmpRotDih1_5, 0, rotDih1[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpRotDih1_5.length);
		                System.arraycopy(tmpRotDih2_5, 0, rotDih2[i[0]][i[1]][i[2]][i[3]][i[4]], 0, tmpRotDih2_5.length);
	                }


	                break;
	            default:
	                System.out.println("Index too high for pair mats.");
	                break;
	        }
	    }


	public void setSupRot(int p1, int a1, int r1, int p2, int a2, int r2,
			int[] supRot1, int[] supRot2) {
		//this.supRot1[p1][a1][r1][p2][a2][r2] = new int[supRot1.length];
		//System.arraycopy(supRot1, 0, this.supRot1[p1][a1][r1][p2][a2][r2], 0, supRot1.length);
		
		//this.supRot2[p1][a1][r1][p2][a2][r2] = new int[supRot2.length];
		//System.arraycopy(supRot2, 0, this.supRot2[p1][a1][r1][p2][a2][r2], 0, supRot2.length);
		
	}



	

	


	
	
	
}
