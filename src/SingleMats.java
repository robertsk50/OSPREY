import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.Serializable;


public class SingleMats implements Serializable {
	//EMatrixEntry[][][] intraE;
	double[][][] E;
	double[][][] maxE;
	boolean[][][] pruned;
	//First 6 dim are indices, last dim are rotamers
	int[][][][] supRot;
	double[][][][] rotDih;
	boolean doDih;
	
	SingleMats(boolean doDih){
		this.doDih = doDih;
	}
	
	
	public void setDihedrals(int[] i, double[] diheds){
		this.rotDih[i[0]][i[1]][i[2]] = diheds;
	}
	
	
	public void setE(int[] i, double E){
		this.E[i[0]][i[1]][i[2]] = E;
	}
	
	public void setMaxE(int[] i, double E){
		this.maxE[i[0]][i[1]][i[2]] = E;
	}
	
	public void addE(int[] i, double E){
		this.E[i[0]][i[1]][i[2]] += E;
	}
	
	public void setE(EMatrixEntrySlim eme){
		this.E[eme.index[0]][eme.index[1]][eme.index[2]] = eme.minE;	
	}
	
	public void setMaxE(EMatrixEntrySlim eme){
		this.maxE[eme.index[0]][eme.index[1]][eme.index[2]] = eme.maxE;	
	}
	
	public void setDihed(EMatrixEntrySlim eme){
		this.rotDih[eme.index[0]][eme.index[1]][eme.index[2]] = eme.rotDih1;	
	}
	
	public void setAll(int pos, int aa, int rot, double E, boolean pruned, int[] supRot){
		this.E[pos][aa][rot] = E;
		this.pruned[pos][aa][rot] = pruned;
		this.supRot[pos][aa][rot] = supRot;
	}
	
	public void setAll(int pos, int aa, int rot, RotamerEntry re){
		this.E[pos][aa][rot] = re.minE();
		this.pruned[pos][aa][rot] = re.isPruned();
		this.supRot[pos][aa][rot] = re.r.rotamers;
	}
	
	public void setAll(int[] i, RotamerEntry re){
		int pos=i[0];
		int aa =i[1];
		int rot=i[2];
		
		this.E[pos][aa][rot] = re.minE();
		this.pruned[pos][aa][rot] = re.isPruned();
		this.supRot[pos][aa][rot] = re.r.rotamers;
	}
	
	
	public void addDim(int[] i, int length){
		int dim = i.length;
		switch(dim){
			case 0:
				E = new double[length][][];
				pruned = new boolean[length][][];
				supRot = new int[length][][][];
				if(doDih){
					rotDih = new double[length][][][];
					maxE = new double[length][][];
				}
				break;
			case 1:
				E[i[0]] = new double[length][];
				pruned[i[0]] = new boolean[length][];
				supRot[i[0]] = new int[length][][];
				if(doDih){
					rotDih[i[0]] = new double[length][][];
					maxE[i[0]] = new double[length][];
				}
				break;
			case 2:
				E[i[0]][i[1]] = new double[length];
				pruned[i[0]][i[1]] = new boolean[length];
				supRot[i[0]][i[1]] = new int[length][];
				if(doDih){
					rotDih[i[0]][i[1]] = new double[length][];
					maxE[i[0]][i[1]] = new double[length];
				}
				break;
			default:
				System.out.println("Index too high for single mats.");
				break;
		}
	}
	
	public void copy(int[] newI, SingleMats singles, int[] oldI){
		this.E[newI[0]][newI[1]][newI[2]] = singles.E[oldI[0]][oldI[1]][oldI[2]];
		this.pruned[newI[0]][newI[1]][newI[2]] = singles.pruned[oldI[0]][oldI[1]][oldI[2]];
		this.supRot[newI[0]][newI[1]][newI[2]] = new int[singles.supRot[oldI[0]][oldI[1]][oldI[2]].length];
		System.arraycopy(singles.supRot[oldI[0]][oldI[1]][oldI[2]], 0, this.supRot[newI[0]][newI[1]][newI[2]], 0, singles.supRot[oldI[0]][oldI[1]][oldI[2]].length);
		
		if(doDih){
			this.maxE[newI[0]][newI[1]][newI[2]] = singles.maxE[oldI[0]][oldI[1]][oldI[2]];
			this.rotDih[newI[0]][newI[1]][newI[2]] = new double[singles.rotDih[oldI[0]][oldI[1]][oldI[2]].length];
			System.arraycopy(singles.rotDih[oldI[0]][oldI[1]][oldI[2]], 0, this.rotDih[newI[0]][newI[1]][newI[2]], 0, singles.rotDih[oldI[0]][oldI[1]][oldI[2]].length);
	}
	}
	
	public void setSupRot(int p1, int a1, int r1, int[] supRot) {
		this.supRot[p1][a1][r1] = new int[supRot.length];
		System.arraycopy(supRot, 0, this.supRot[p1][a1][r1], 0, supRot.length);
	}
	
	public int[] getRot(int[] i){
		return supRot[i[0]][i[1]][i[2]];
	}
	
	public int[] getRot(int pos, int aa, int rot){
		return supRot[pos][aa][rot];
	}
	
	public int[] getRot(Index3 i){
		return supRot[i.pos][i.aa][i.rot];
	}
	
	
	public RotamerEntry getTerm(int p1, int a1, int r1) {
		
		SuperRotamer sr1 = new SuperRotamer(supRot[p1][a1][r1]);
		
		if(rotDih == null)
			return new RotamerEntry(p1, sr1, E[p1][a1][r1], Double.MAX_VALUE,pruned[p1][a1][r1]);
		
		return new RotamerEntry(p1, sr1, E[p1][a1][r1], maxE[p1][a1][r1],pruned[p1][a1][r1],rotDih[p1][a1][r1]);
		
	}
	
	public RotamerEntry getTerm(int[] i) {
		
		SuperRotamer sr1 = new SuperRotamer(supRot[i[0]][i[1]][i[2]]);
		
		if(doDih)
			return new RotamerEntry(i[0], sr1, E[i[0]][i[1]][i[2]], maxE[i[0]][i[1]][i[2]],pruned[i[0]][i[1]][i[2]],rotDih[i[0]][i[1]][i[2]]);
		else
			return new RotamerEntry(i[0], sr1, E[i[0]][i[1]][i[2]], Double.MAX_VALUE, pruned[i[0]][i[1]][i[2]]);
		
	}
	
	public RotamerEntry getTerm(Index3 i) {
		
		SuperRotamer sr1 = new SuperRotamer(supRot[i.pos][i.aa][i.rot]);
		
		if(doDih)
			return new RotamerEntry(i.pos, sr1, E[i.pos][i.aa][i.rot], maxE[i.pos][i.aa][i.rot],pruned[i.pos][i.aa][i.rot],rotDih[i.pos][i.aa][i.rot]);
		else
		return new RotamerEntry(i.pos, sr1, E[i.pos][i.aa][i.rot], pruned[i.pos][i.aa][i.rot]);
		
	}
	
	
		public void write(String fileName) {
			KSParser.outputObject(E,fileName+".singlesE");
			
			KSParser.outputObject(pruned,fileName+".singlesPruned");
			//First 6 dim are indices, last dim are rotamers
			KSParser.outputObject(supRot,fileName+".singlesSupRot");
			if(doDih){
				KSParser.outputObject(rotDih,fileName+".singlesRotDih");
				KSParser.outputObject(maxE,fileName+".singlesMaxE");
			}
			
			
}
		public static SingleMats read(String fileName, boolean doDih){
			
			SingleMats sm = new SingleMats(doDih);
			
			try{
//				ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName+".singlesE"));
				sm.E = (double[][][])KSParser.loadObject(fileName+".singlesE");//in.readObject();
//				in.close();
				
//				in = new ObjectInputStream(new FileInputStream(fileName+".singlesPruned"));
				sm.pruned = (boolean[][][])KSParser.loadObject(fileName+".singlesPruned");//in.readObject();
//				in.close();
//				in = new ObjectInputStream(new FileInputStream(fileName+".singlesSupRot"));
				sm.supRot = (int[][][][])KSParser.loadObject(fileName+".singlesSupRot");//in.readObject();
//				in.close();
				if(doDih){
//					in = new ObjectInputStream(new FileInputStream(fileName+".singlesMaxE"));
					sm.maxE = (double[][][])KSParser.loadObject(fileName+".singlesMaxE");//in.readObject();
//					in.close();
				
//					in = new ObjectInputStream(new FileInputStream(fileName+".singlesRotDih"));
					sm.rotDih = (double[][][][])KSParser.loadObject(fileName+".singlesRotDih");//in.readObject();
//					in.close();
				}
				
			}
			catch (Exception e){
				//e.printStackTrace();
				System.out.println("Could not find/read file: "+fileName);
				return null;
			}		
			
			return sm;
			
			
			
		}

		public void extendDim(int[] i, int toExtendBy){
	        int dim = i.length;
	        switch(dim){
	            /*case 0:
	                E = new double[length][][];
	                pruned = new boolean[length][][];
	                supRot = new int[length][][][];
	                break;
	            case 1:
	                E[i[0]] = new double[length][];
	                pruned[i[0]] = new boolean[length][];
	                supRot[i[0]] = new int[length][][];
	                break;*/
	            case 2:
	                double[] tmpE = E[i[0]][i[1]];  
	                boolean[] tmpPruned = pruned[i[0]][i[1]];
	                int[][] tmpSupRot = supRot[i[0]][i[1]];
	                E[i[0]][i[1]] = new double[E[i[0]][i[1]].length+toExtendBy];
	               
	                pruned[i[0]][i[1]] = new boolean[pruned[i[0]][i[1]].length+toExtendBy];
	                supRot[i[0]][i[1]] = new int[supRot[i[0]][i[1]].length+toExtendBy][];

	                System.arraycopy(tmpE, 0, E[i[0]][i[1]], 0, tmpE.length);
	               
	                System.arraycopy(tmpPruned, 0, pruned[i[0]][i[1]], 0, tmpPruned.length);
	                System.arraycopy(tmpSupRot, 0, supRot[i[0]][i[1]], 0, tmpSupRot.length);
	                
	                if(doDih){
	                	double[] tmpMaxE = maxE[i[0]][i[1]];
	                	maxE[i[0]][i[1]] = new double[maxE[i[0]][i[1]].length+toExtendBy];
		                System.arraycopy(tmpMaxE, 0, maxE[i[0]][i[1]], 0, tmpMaxE.length);
		                double[][] tmpRotDih = rotDih[i[0]][i[1]];
		                rotDih[i[0]][i[1]] = new double[rotDih[i[0]][i[1]].length+toExtendBy][];
		                System.arraycopy(tmpRotDih, 0, rotDih[i[0]][i[1]], 0, tmpRotDih.length);
	                }

	                break;
	            default:
	                System.out.println("Index too high for single mats.");
	                break;
	        }
	    }


		
	}

