
public class MathUtils {

	public static double[] cross(double[] v1, double[] v2){
		double[] ret = new double[3];
		ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
		ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
		return ret;
	}
	
}
