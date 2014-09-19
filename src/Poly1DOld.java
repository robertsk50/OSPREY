import java.util.StringTokenizer;


public class Poly1DOld {

	static final int VALUE = 0;
	static final int DERIV = 1;
	
	//Coefficients for energy function started with largest degree
	double[] coefficients;
	int degree; //degree of the polynomial
	double xmin; //Min and max are the boundaries where the polynomial is valid
	double xmax;
	
	double fmin; //Smoothing Interval Terms
	double fmax; 
	double b_a;
	double d_c;
	
	double min_val = Double.POSITIVE_INFINITY; //Values of energy function if outside polynomial range
	double max_val = Double.POSITIVE_INFINITY;
	
	Poly1DOld(String s){
		
		StringTokenizer st = new StringTokenizer(s);
		
		xmin = new Double(st.nextToken());
		xmax = new Double(st.nextToken());
		
		if(xmin > xmax){
			double tmp = xmin;
			xmin = xmax;
			xmax = tmp;
		}
		
		degree = new Integer(st.nextToken()) - 1;
		
		coefficients=new double[degree+1];
		for(int i=0; i<coefficients.length;i++){
			coefficients[i] = new Double(st.nextToken());
		}
		
		b_a = (xmax-xmin)*0.10;
		d_c = b_a;
		
		fmin = b_a + xmin;
		fmax = xmax - d_c;
		
		
		
	}
	
	void printFunction(){
		double incr = (xmax-xmin)/180;
		
		for(double i=xmin; i <=xmax; i+= incr){
			System.out.println(i+" "+getVal(i));
		}
		
	}
	
	
	////////////////////////////////////////////////////////////////////////////////
	/// //KER: Taken from Rosetta 3.4 polynomial.cc
	///	@begin operator()
	///
	/// @brief evaluate the polynomial and its derivative.
	///
	/// @detailed
	///
	/// @param  variable - [in] - evaluate polynomial(value)
	/// @param  value - [out] - returned output
	/// @param  deriv - [out] - returned output
	///
	/// @global_read
	///
	/// @global_write
	///
	/// @remarks
	///  Note the coefficients must be in reverse order: low to high
	///
	///  Polynomial value and derivative using Horner's rule
	///  value = Sum_(i = 1,...,N) [ coeff_i * variable^(i-1) ]
	///  deriv = Sum_(i = 2,...,N) [ ( i - 1 ) * coeff_i * variable^(i-2) ]
	///  JSS: Horner's rule for evaluating polynomials is based on rewriting the polynomial as:
	///  JSS: p(x)  = a0 + x*(a1 + x*(a2 + x*(...  x*(aN)...)))
	///  JSS: or value_k = a_k + x*value_k+1 for k = N-1 to 0
	///  JSS: and the derivative is
	///  JSS: deriv_k = value_k+1 + deriv_k+1 for k = N-1 to 1
	///
	/// @references
	///
	/// @authors Jack Snoeyink
	/// @authors Matthew O'Meara
	///
	/// @last_modified Matthew O'Meara
	/////////////////////////////////////////////////////////////////////////////////
	void getValueAndDeriv(double variable, double[] retArr){
		
		if(variable <= xmin){
			retArr[VALUE] = min_val;
			retArr[DERIV] = 0.0;
			return;
		}
		if(variable >= xmax){
			retArr[VALUE] = max_val;
			retArr[DERIV] = 0.0;
			return;
		}
		
		
		double value = coefficients[0];
		double deriv = 0.0;
		for(int i=1; i <= degree; i++){
			deriv *= variable;
			deriv += value;
			value *= variable;
			value += coefficients[i];
		}
		retArr[VALUE] = value;
		retArr[DERIV] = deriv;
		
	}
	
	double getVal(double variable){
		
		if(variable <= xmin){
			return min_val;
		}
		if(variable >= xmax){
			return max_val;
			
		}
		
		
		double value = coefficients[0];
		for(int i=1; i <= degree; i++){
			value *= variable;
			value += coefficients[i];
		}
		
		
		return value;
		
		
	}
	
	public double getSmoothVal(double variable){
		double smoothVal = 1.0;
	
		//Smooth Value Using 10% and 90% for cutoffs
		if(variable < xmin || variable > xmax)
			smoothVal = 0.0;
		else if(variable < fmin)
			smoothVal = (variable - xmin)/b_a;
		else if(variable > fmax)
			smoothVal = (xmax - variable)/d_c;
		
		return smoothVal;
	}
		
	
}
