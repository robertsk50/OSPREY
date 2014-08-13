import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;

/*
 * Created on 06-25-2014
 *
 */

/**
 * @author Kyle Roberts
 *
 * 
 */
public class MatrixRead {

	public static void main(String[] args) {


		
		
		PairwiseEnergyMatrix a1 = new PairwiseEnergyMatrix();
		PairwiseEnergyMatrix a2 = new PairwiseEnergyMatrix();
		

        String f1 = "/usr/xtmp/kroberts/tutorials/OSPREY/DHFR_noEpic_working/PEM/PEMmin_COM.dat";
		String f2 = "/usr/xtmp/kroberts/tutorials/OSPREY/DHFR_noEpic_master/PEM/PEMmin_COM.dat";
		


		//String fOut = "1amuAspSCPEMmax_full.dat";

                try{
			ObjectInputStream in1 = new ObjectInputStream(new FileInputStream(f1));
            ObjectInputStream in2 = new ObjectInputStream(new FileInputStream(f2));
        
			a1.eMatrix = (double [][][][][][])in1.readObject();
			a2.eMatrix = (double [][][][][][])in2.readObject();
			
			in1.close();
            in2.close();
  		}
		catch (Exception e){
                    e.printStackTrace();
                }

                int same = 0;
                int different = 0;
                for(int i=0; i<a1.eMatrix.length;i++){
                    if(a1.eMatrix[i] != null)
                    for(int j=0; j<a1.eMatrix[i].length;j++){
                        if(a1.eMatrix[i][j] != null)
                        for(int k=0; k<a1.eMatrix[i][j].length;k++){
                            if(a1.eMatrix[i][j][k] != null)
                            for(int l=0; l<a1.eMatrix[i][j][k].length;l++){
                                if(a1.eMatrix[i][j][k][l] != null)
                                for(int m=0; m<a1.eMatrix[i][j][k][l].length;m++){
                                    if(a1.eMatrix[i][j][k][l][m] != null)
                                    for(int n=0; n<a1.eMatrix[i][j][k][l][m].length;n++){
                                    	if(i== 21 && j==6 && k==7 && l == 22 && m==6 && n==1){
                                    		System.out.println(a1.eMatrix[i][j][k][l][m][n] +" "+a2.eMatrix[i][j][k][l][m][n] );
                                    	}
                                        if(Math.abs(a1.eMatrix[i][j][k][l][m][n] - a2.eMatrix[i][j][k][l][m][n]) < 0.1){
                                            same++;
                                            //System.out.println("Same: "+i+" "+j+" "+k+" "+l+" "+m+" "+n+" : "+a1.eMatrix[i][j][k][l][m][n] +" "+a2.eMatrix[i][j][k][l][m][n] );
                                        }
                                        else{
                                            different++;
                                           System.out.println("Diff: "+i+" "+j+" "+k+" "+l+" "+m+" "+n+" : "+a1.eMatrix[i][j][k][l][m][n] +" "+a2.eMatrix[i][j][k][l][m][n] );
                                        }
                                    }
                                }
                            }

                        }
                    }

                }
                
                System.out.println("SAME: "+same);
                System.out.println("DIFFERENT: "+different);
                
			
		/*for (int i=0; i<matrixSize; i++){
			for (int j=0; j<matrixSize; j++){
				if ((i<1369)&&(j<1369))
					aOut[i][j] = a1[i][j];
				else
					aOut[i][j] = a2[i][j];
			}
		}
		
		writeO(fOut,aOut);*/
		
		/*final int matrixSize = 1374;
		float a1[][] = new float [matrixSize][matrixSize];
		float a2[][] = new float [matrixSize][matrixSize];
		float a3[][] = new float [matrixSize][matrixSize];
		float aOut[][] = new float [matrixSize][matrixSize];
		
		String f1 = "1amuLeuBBPEMmax_pt1.dat";
		String f2 = "1amuLeuBBPEMmax_pt2.dat";
		String f3 = "1amuLeuBBPEMmax_pt3.dat";
		String fOut = "1amuLeuBBPEMmax.dat";
		
		a1 = readO(f1,a1);
		a2 = readO(f2,a2);
		a3 = readO(f3,a3);
		
		for (int i=0; i<matrixSize; i++){
			for (int j=0; j<matrixSize; j++){
				if (a1[i][j]!=0.0f)
					aOut[i][j] = a1[i][j];
				else if (a2[i][j]!=0.0f)
					aOut[i][j] = a2[i][j];
				else if (a3[i][j]!=0.0f)
					aOut[i][j] = a3[i][j];
			}
		}
		
		int num = 0;
		for (int i=0; i<matrixSize; i++){
			for (int j=0; j<matrixSize; j++){
				if (aOut[i][j]!=0.0f)
					num++;
			}
		}
		
		writeO(fOut,aOut);
		
		System.out.println(num);*/
	}
	
	public static float [][] readO(String fname, float o[][]){
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fname));
			o = (float [][])in.readObject();
			in.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while reading energy matrices");
			System.exit(0);
		}
		return o;
	}
	
	public static void writeO(String fname, float o[][]){
		try{
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(fname));
			out.writeObject(o);
			out.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing file");
			System.exit(0);
		}
	}
}
