import java.io.Serializable;
import java.util.ArrayList;



class Index3wVal implements Comparable<Index3wVal>,Serializable{

	ArrayList<Index3> i3s;
	Double val;
	int pos;

	Index3wVal(Index3[] i3s, double val ){
		this.i3s = new ArrayList<Index3>();
		for(int i=0; i<i3s.length;i++)
			this.i3s.add(i3s[i]);
		
		//this.i3s = i3s;
		this.val = val;
	}
	
	Index3wVal(Index3[] i3s, double val,int pos ){
		this.i3s = new ArrayList<Index3>();
		for(int i=0; i<i3s.length;i++)
			this.i3s.add(i3s[i]);
		
		//this.i3s = i3s;
		this.val = val;
		this.pos = pos;
	}
	
	Index3wVal(ArrayList<Index3> i3s, double val ){
		
		this.i3s = i3s;
		this.val = val;
	}
	Index3wVal(ArrayList<Index3> i3s, double val, int pos ){
		
		this.i3s = i3s;
		this.val = val;
		this.pos = pos;
	}

	@Override
	public int compareTo(Index3wVal arg) {
		//Descending order;
		return (-1 * val.compareTo(arg.val));			
	}
	
	public boolean containsPos(int pos){
		for(int i=0; i<i3s.size();i++ )
			if(i3s.get(i).pos == pos)
				return true;
		
		return false;
	}

}