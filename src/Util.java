
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
		// Returns the index of the queue node that has the minimum fscore in an array.	
	public static int indexOfMinQueueNodeInArray(PGQueueNode nodeArray[]){
		double minimumFscore = Double.POSITIVE_INFINITY;
		int indexOfMin = -1;
		for(int i = 0; i < nodeArray.length; i++){
			if(nodeArray[i].fScore< minimumFscore){
				minimumFscore = nodeArray[i].fScore;
				indexOfMin = i;
			}
			
		}
		return indexOfMin;
	}
	// Compute an RMS between two arrays of dihedrals.  Corrects the angle if on the edges (-180 and +180) 
	public static double RootMeanSquareDihedrals(double array1[], double array2[]){
		double squared_sum = 0;
		if (array1.length != array2.length){
			System.out.println("Computing root mean square of two arrays with different sizes.");
			System.exit(1);
		}
		for (int i = 0; i < array1.length; i++){
			double minElem, maxElem;
			if(array1[i] < array2[i]){
				minElem = array1[i];
				maxElem = array2[i];
			}
			else{
				minElem = array2[i];
				maxElem = array1[i];
			}
			double distance = maxElem - minElem;
			// Handle the "border" cases where dihedral values are close to 180. 
			if((minElem+360 - maxElem) < distance){
				distance = minElem+360 - maxElem;
			}			
			squared_sum += distance*distance;
		}
		return Math.sqrt(squared_sum/((double)array1.length));
	}
}
