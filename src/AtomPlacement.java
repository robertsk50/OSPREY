
public class AtomPlacement {
	String a1Name;
	String a2Name;
	String a3Name;
	
	double len;
	double ang;
	double dihe;
	
	public AtomPlacement(String a1Name,String a2Name, String a3Name,double len, double ang, double dihe) {
		this.a1Name = a1Name;
		this.a2Name = a2Name;
		this.a3Name = a3Name;
		
		this.len = len;
		this.ang = ang;
		this.dihe = dihe;
	}
	
	public AtomPlacement(Atom a1, Atom a2, Atom a3, Atom a4) {
		this.a1Name = a1.name;
		this.a2Name = a2.name;
		this.a3Name = a3.name;
		
		this.len = a3.distance(a4);
		this.ang = a4.angle(a2, a3);
		this.dihe = a4.torsion(a1, a2, a3);
	}
}
