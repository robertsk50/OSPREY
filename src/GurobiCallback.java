import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;


public class GurobiCallback extends GRBCallback {

	double objVal = Double.NEGATIVE_INFINITY;
	double iterCount = 0;
	long maxIter;
	boolean didAbort = false;
	
	public GurobiCallback(long maxIter){
		this.maxIter = maxIter;
	}
	
	@Override
	protected void callback() {
		if(where == GRB.Callback.SIMPLEX){
			try {
				
				iterCount = getDoubleInfo(GRB.Callback.SPX_ITRCNT);
				if(iterCount > maxIter){
					objVal = getDoubleInfo(GRB.Callback.SPX_OBJVAL);
					didAbort = true;
					abort();
				}
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
			

 	}
	
	public boolean didAbort(){
		return didAbort;
	}
	
	public double getObjVal(){
		return objVal;
	}

}
