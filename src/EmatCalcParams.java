import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeSet;


public class EmatCalcParams implements Serializable{
	int pos1;
	int pos2;
	
	TreeSet<Integer> AAs1 = null;
	TreeSet<Integer> AAs2 = null;
	
	ArrayList<Index3> rotamers = null;
	
	EmatCalcParams(int p1,int p2, TreeSet<Integer> AAs1, TreeSet<Integer> AAs2){
		pos1 = p1;
		pos2 = p2;
		this.AAs1 = AAs1;
		this.AAs2 = AAs2;
	}
	
	EmatCalcParams(int p1,int p2, ArrayList<Index3> rotamers){
		pos1 = p1;
		pos2 = p2;
		this.rotamers = rotamers;
	}
	
	EmatCalcParams(int p1,int p2){
		pos1 = p1;
		pos2 = p2;
	
	}
	
}
