import java.util.Comparator;


public class ArrayIndexComparator implements Comparator<Integer>{

	private final double[] array;
	
	public ArrayIndexComparator(double[] array){
		this.array = array;
	}
	
	public ArrayIndexComparator(int[] array){
		this.array = new double[array.length];
		for(int i=0; i<this.array.length;i++){
			this.array[i] = array[i];
		}
	}
	
	public Integer[] createIndexArray(){
		Integer[] indexes = new Integer[array.length];
		for(int i=0; i<array.length;i++){
			indexes[i] = i;
		}
		return indexes;
	}
	
	@Override
	public int compare(Integer index1, Integer index2){
		if(array[index1] > array[index2])
			return -1;
		else if(array[index1] < array[index2])
			return 1;
		else return 0;
	}
	
}
