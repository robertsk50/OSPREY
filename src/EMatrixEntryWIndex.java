import java.io.Serializable;


public class EMatrixEntryWIndex implements Serializable {
	
	EMatrixEntry eme;
	int[] index;
	
	public EMatrixEntryWIndex(EMatrixEntry ematent,int[] ind) {
		eme = ematent;
		index = ind;
	}
	
	public EMatrixEntryWIndex(EMatrixEntry ematent,Index3 ind) {
		eme = ematent;
		index = new int[3];
		index[0] = ind.pos; index[1]=ind.aa;index[2]=ind.rot;
	}
	
	/**
	 * This constructor makes a new PairMatrix entry from two intra matrix entries
	 * @param emat
	 * @param eme1
	 * @param eme2
	 */
	public EMatrixEntryWIndex(Emat emat, EMatrixEntryWIndex eme1,EMatrixEntryWIndex eme2){
		eme = emat.pairs.getTerm(eme1.index, eme2.index, emat.singles.supRot);
		index = new int[eme1.index.length+eme2.index.length];
		for(int i=0; i<index.length;i++){
			if(i>=eme1.index.length){
				index[i] = eme2.index[i%eme1.index.length];
			}else{
				index[i] = eme1.index[i];
			}
		}
	}

	public int pos1(){
		return index[0];
	}
	
	public int aa1(){
		return index[1];
	}
	
	public int rot1(){
		return index[2];
	}
	
	public int pos2(){
		return index[3];
	}
	
	public int aa2(){
		return index[4];
	}
	
	public int rot2(){
		return index[5];
	}
	
	public int[] rot1index(){
		int[] rot1 = {index[0],index[1],index[2]};
		return rot1;
	}
	
	public Index3 rot1index3(){
		Index3 rot1 = new Index3(index[0],index[1],index[2]);
		return rot1;
	}
	
	public int[] rot2index(){
		int[] rot1 = {index[3],index[4],index[5]};
		return rot1;
	}
	
	public Index3 rot2index3(){
		Index3 rot2 = new Index3(index[3],index[4],index[5]);
		return rot2;
	}
	
	public String toString(){
		String s = "";
		
		s += index[0]+"_"+index[1]+"_"+index[2];
		if(index.length > 3)
			s += " "+index[3]+"_"+index[4]+"_"+index[5];
		
		return s;
	}
	
	
}
