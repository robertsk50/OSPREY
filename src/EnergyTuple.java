import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Vector;




class EnergyTuple implements Comparable<EnergyTuple>,Serializable{
//		double[][][] E;
		double intraE;
		Index3[] rots;
		//Integer[] indexes;
	
		ArrayList<EnergyTuple> children;
				
		EnergyTuple(ArrayList<Index3> rots){
			
			this.rots = rots.toArray(new Index3[0]);
			children = new ArrayList<EnergyTuple>();
		}
		
		@Override
		public boolean equals (Object o){
			if(o instanceof EnergyTuple ){
				EnergyTuple p1 = (EnergyTuple)o;
				if(rots.length != p1.rots.length)
					return false;
				for(int i=0; i<rots.length;i++){
					if(rots[i] != p1.rots[i])
						return false;
				}
				
				return true;
			}
			return false;
		}
		
		@Override
		public int compareTo(EnergyTuple pair2) {
			for(int i=0; i< rots.length;i++){
				if(rots[i].compareTo(pair2.rots[i]) < 0)
					return -1;
				else if(rots[i].compareTo(pair2.rots[i]) > 0)
					return 1;
			}
			//All equal
			return 0;
			
		}
		
		@Override
		public int hashCode() { 
		    int hashCode = 23;
		    for(int i=0; i<rots.length;i++)
		    	hashCode += hashCode*31 + rots[i].hashCode();
		    //System.out.println("in hashCode() [" + Integer.toString(hashCode) + "]");
		    return(hashCode);
		  }

		public boolean isParent(EnergyTuple child) {
			if(child.rots.length != rots.length+1)
				return false;
			
			int rotsNonMatching = 0;
			for(Index3 r1: child.rots){
				boolean matchedRot = false;
				for(Index3 r2:rots){
					if(r1.equals(r2)){
						matchedRot = true;
					}
				}
				if(!matchedRot)
					rotsNonMatching++;
			}
			if(rotsNonMatching == 1)
				return true;
			
			return false;
			
		}

		public boolean noSharedPositions(EnergyTuple tuple) {
		for(Index3 r1: rots){
			for(Index3 r2: tuple.rots){
				if(r1.pos == r2.pos)
					return false;
			}
		}
		return true;
	}

		public boolean containsRot(Index3 r2) {
			for(Index3 r1:rots ){
				if(r1.equals(r2))
					return true;
			}
			
			return false;
		}

		public ArrayList<Index3> sharedRots(EnergyTuple tuple) {
			ArrayList<Index3> rots = new ArrayList<Index3>();
			for(Index3 r1:rots ){
				for(Index3 r2: tuple.rots)
				if(r1.equals(r2))
					rots.add(r2);
			}
			return rots;
		}

		//KER: Returns rots in tuple that are not in "this"
		public ArrayList<Index3> nonSharedRots(EnergyTuple tuple) {
			ArrayList<Index3> rots = new ArrayList<Index3>();
			
			for(Index3 r1:tuple.rots ){
				boolean nonSharedRot = true;
				for(Index3 r2: this.rots){
					if(r1.equals(r2))
						nonSharedRot = false;
				}
				if(nonSharedRot)
					rots.add(r1);
			}
			return rots;
		}

		public boolean shareRots(EnergyTuple tuple) {
			for(Index3 r1:rots ){
				for(Index3 r2: tuple.rots)
				if(r1.equals(r2))
					return true;
			}
			return false;
		}
		
		public boolean shareRots(ArrayList<Index3> rots) {
			for(Index3 r1:this.rots ){
				for(Index3 r2: rots)
				if(r1.equals(r2))
					return true;
			}
			return false;
		}

		public boolean allSharedPosMatch(EnergyTuple tuple) {
			boolean sharePos = false;
			for(Index3 r1:rots ){
				for(Index3 r2: tuple.rots)
				if(r1.pos == r2.pos){
					sharePos = true;
					if(!r1.equals(r2))
						return false;
				}
			}
			if(sharePos)
				return true;
			else
				return false;
		}

		
		public boolean containsAllRot(EnergyTuple tuple) {
			if(rots.length < tuple.rots.length)
				return false;
			
			for(Index3 r1:tuple.rots ){
				boolean foundPos = false;
				for(Index3 r2: rots){
					if(r1.pos == r2.pos){
						if(!r1.equals(r2))
							return false;
						else
							foundPos = true;
					}
				}
				if(!foundPos)
					return false;
			}
			
			return true;
			
		}
		
		
	}
