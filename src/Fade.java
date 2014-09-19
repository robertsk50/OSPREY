
public class Fade {

	String fadeintervalname;
	double min0;
	double max0;
	double fmin;
	double fmax;
	double dfademin;
	double dfademax;
	final boolean smooth;
	
	public Fade(String fadeintervalname, double min0, double fmin,
			double fmax, double max0, boolean smoothed) {
		this.fadeintervalname = fadeintervalname;
		this.min0 = min0;
		this.fmin = fmin;
		this.fmax = fmax;
		this.max0 = max0;
		this.smooth = smoothed;
		this.dfademin = 1.0/(fmin-min0);
		this.dfademax = 1.0/(max0-fmax);
	}


	public double getSmoothVal(double x){
		double val;
	
		
		if (x <= fmax) {
			if (x <= min0) { return 0.0; }       // in (-\infty, min0]
			if (x >= fmin) { return 1.0; }  // in [fmin, fmax]
			if (smooth){
				final double z = ((x - min0)* dfademin);
				val = z*z*(3-2*z);
				//deriv = -6*z*(z-1)*dfademin;
			} else {
				//deriv = dfademin; 
				val = (x - min0) * dfademin;  // in (min0,fmin)
			}
		} else {
			if (x >= max0) { return 0.0; }       // in [max0, \infty)
			if (smooth){
				final double z = ((x - fmax) * dfademax);
				val = z*z*(2*z-3) + 1;
				//deriv = 6*z*(z-1)*dfademax;
			} else {
				//deriv = -dfademax; 
				val = (max0 - x) * dfademax; // in (fmax,max0)
			}
		}
		
		return val;
		
//		//Smooth Value Using 10% and 90% for cutoffs
//		if(variable < min0 || variable > max0)
//			smoothVal = 0.0;
//		else if(variable < fmin)
//			smoothVal = (variable - xmin)/ba;
//		else if(variable > fmax)
//			smoothVal = (xmax - variable)/dc;
//		
//		return smoothVal;
	}
	
}
