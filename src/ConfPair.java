
public class ConfPair implements Comparable{
	EMatrixEntryWIndex[] conf;
	double[] energy;
	public ConfPair(EMatrixEntryWIndex[] conformation, double[] e){
		conf = new EMatrixEntryWIndex[conformation.length];
		for(int i=0; i<conformation.length;i++)
			conf[i] = conformation[i];
		energy = new double[e.length];
		for(int i=0; i<e.length;i++)
			energy[i] = e[i];
		
	}
	
	@Override
	public int compareTo(Object o) throws ClassCastException {
		// TODO Auto-generated method stub
		if(!(o instanceof ConfPair))
			throw new ClassCastException("Another confPair was expected.");
		double otherE = ((ConfPair) o).energy[0];
		if(otherE >= energy[0])
			return 1;
		else
			return -1;
		
	}

}
