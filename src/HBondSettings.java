
public class HBondSettings {
	boolean doHbondE = false;
	double hbondScale = 0.0;
	String dsspFile = "";
	
	HBondSettings(double hbondScale, String dsspFile){
		this.hbondScale = hbondScale;
		if(this.hbondScale > 0){
			doHbondE = true;
			this.dsspFile = dsspFile;
		}
		
	}
	
}
