
public class Util {
	// Returns the min of an array
	public static double min(double array []){
		int minIndex = 0;
		for (int i = 1; i < array.length; i++){
			if(array[minIndex] > array[i]){
				minIndex = i;
			}			
		}
		return array[minIndex];
	}
	
	//KER: helper function that returns the min of the array in pos 0
	//KER: and the index of the Min at pos 1
	public static double[] indexOfMin(double[] arr){
		double min = Double.POSITIVE_INFINITY;
		double[] index = new double[2];
		for(int i=0; i<arr.length;i++){
			if(arr[i] < min){
				min = arr[i];
				index[0] = min;
				index[1] = i;
			}
		}
		return index;
	}

	//KER: helper function that returns the min of the array in pos 0
	//KER: and the index of the Min at pos 1
	public static int[] indexOfMin(int[] arr){
		int min = Integer.MAX_VALUE;
		int[] index = new int[2];
		for(int i=0; i<arr.length;i++){
			if(arr[i] < min){
				min = arr[i];
				index[0] = min;
				index[1] = i;
			}
		}
		return index;
	}
	
}
