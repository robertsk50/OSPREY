import java.io.Serializable;
import java.util.StringTokenizer;

@SuppressWarnings("serial")
class Index3 implements Comparable<Index3>,Serializable {
	int pos;
	int aa;
	int rot;
	
	
	
	public Index3(int curPos, int curAA, int curRot) {
		this.pos = curPos;
		this.aa = curAA;
		this.rot = curRot;
	}
	
	public Index3(int[] i) {
		this.pos = i[0];
		this.aa = i[1];
		this.rot = i[2];
	}
	
	public Index3(String s){
		//Trim parens
		String sSub = s.substring(1, s.length()-1);
		StringTokenizer st = new StringTokenizer(sSub, ",");
		pos = new Integer(st.nextToken()); 
		aa = new Integer(st.nextToken());
		rot = new Integer(st.nextToken());
	}
	
	@Override
	public boolean equals (Object o){
		if(o instanceof Index3 ){
			Index3 b = (Index3) o;
			if(pos == b.pos && aa == b.aa && rot == b.rot)
				return true;
		}
		return false;
	}
	
	
	
	public String toString(){
		return "("+pos+","+aa+","+rot+")";
	}

	@Override
	public int compareTo(Index3 o) {
		if(pos < o.pos){
			return -1;
		}
		else if(pos == o.pos){//pos either = or great
			if(aa < o.aa){
				return -1;
			}
			else if(aa == o.aa){//pos either = or great
				if(rot < o.rot){
					return -1;
				}
				else if(rot == o.rot){//pos either = or great
					return 0;
				}
				else
					return 1;
			}
			else
				return 1;
		}
		else
			return 1;
		
		
	}
	
	@Override
	public int hashCode() { 
	    int hashCode = 23;
	    
	    hashCode += hashCode*31 + new Integer(pos).hashCode();
	    hashCode += hashCode*31 + new Integer(aa).hashCode();
	    hashCode += hashCode*31 + new Integer(rot).hashCode();
	    //System.out.println("in hashCode() [" + Integer.toString(hashCode) + "]");
	    return(hashCode);
	  }
	
	

	
}
