import java.io.Serializable;


public class DEEsettings implements Serializable {

	boolean gold;
	boolean split1;
	boolean split2;
	boolean mb;
	
	
	DEEsettings(boolean gold,boolean split1,boolean split2,boolean mb){
		this.gold = gold;
		this.split1 = split1;
		this.split2 = split2;
		this.mb = mb;
	}
	
}
