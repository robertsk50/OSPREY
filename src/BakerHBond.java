import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

public class BakerHBond {

	public static final boolean debug = false;

	public static final byte SP2 = 1;
	public static final byte SP3 = 2;
	public static final byte RING = 2;

	//KER: HBond Poly1DOld Functions
	Poly1DOld delta_bb_helix;
	Poly1DOld delta_bb_sheet;
	Poly1DOld delta_bb_other;
	Poly1DOld delta_sc_sp2;
	Poly1DOld delta_sc_sp3;

	Poly1DOld chi_bb_helix;
	Poly1DOld chi_bb_sheet;
	Poly1DOld chi_bb_other;
	Poly1DOld phi_bb_helix;
	Poly1DOld phi_bb_sheet;
	Poly1DOld phi_bb_other;
	Poly1DOld theta_bb_helix;
	Poly1DOld theta_bb_sheet;
	Poly1DOld theta_bb_other;

	Poly1DOld chi_sc_sp2_short;
	Poly1DOld chi_sc_sp2_long;
	Poly1DOld phi_sc_sp2_short;
	Poly1DOld phi_sc_sp2_long;
	Poly1DOld theta_sc_sp2_short;
	Poly1DOld theta_sc_sp2_long;

	Poly1DOld phi_sc_sp3_short;
	Poly1DOld theta_sc_sp3_short;
	Poly1DOld phi_sc_sp3_long;
	Poly1DOld theta_sc_sp3_long;

	ArrayList<HbondPair> hbondTerms; //Hbond Pairs
	ArrayList<HbondPair>[] partHBlist;

	double hbondScale;

	BakerHBond(double hbondScaleFactor, String dsspFile, Molecule m){
		this.hbondScale = hbondScaleFactor;
		try {
			readDSSP(dsspFile,m);
		}
		catch(Exception e ){
			System.out.println("ERROR: An error occurred while reading DSSP file: "+dsspFile);
			System.exit(0);
		}


		readPoly1DOld();

	}


	private void readDSSP(String dsspFile, Molecule m) throws Exception {
		FileInputStream is = new FileInputStream( dsspFile );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;

		// Skip over the first line of header info
		curLine = bufread.readLine();
		while(!(getToken(curLine,2).equals("RESIDUE"))){
			curLine = bufread.readLine();
		}

		ArrayList<SSPair> residues = new ArrayList<SSPair>();

		// 1. Read Residue secondary structure
		curLine = bufread.readLine();		
		// Until we're at a blank line (or until we've read numAtomTypes)
		while (curLine != null) {
			String res = curLine.substring(5,10);
			res = res.trim();
			String SStype = curLine.substring(16, 17);
			residues.add(new SSPair(res,SStype));
			curLine = bufread.readLine();
		}
		bufread.close();

		Iterator<SSPair> iter = residues.iterator();

		int ctr=0;
		String nextPDBresNum = m.residue[ctr].getResNumberString();
		while(iter.hasNext()){
			SSPair ssp = iter.next();
			if(ssp.residue.equals(nextPDBresNum)){
				if(ssp.SStype.equals("H"))
					m.residue[ctr].SStype = Molecule.HELIX; 
				else if(ssp.SStype.equals("E"))
					m.residue[ctr].SStype = Molecule.STRAND;
				else
					m.residue[ctr].SStype = Molecule.OTHER;

				ctr++;
				if(ctr < m.residue.length){
					nextPDBresNum = m.residue[ctr].getResNumberString();
				}
			}
		}


	}

	//Each string represents a polynomial function
	//mindomain, maxdomain, degree, coefficients
	private void readPoly1DOld(){

		delta_bb_helix = new Poly1DOld("1.70 2.60 7 76.92576075872044327753 -920.31423728022900832002 4582.98159388960266369395 -12230.47471070210121979471 18585.45580048343617818318 -15372.59121156383480411023 5445.07452913820725370897 ");
		delta_bb_sheet = new Poly1DOld("1.70 2.60 7 144.34922152025424679778 -1884.15932310602352117712 10261.76991717469900322612 -29877.53999254936570650898 49086.73306215617776615545 -43168.39618908205738989636 15875.40397483624838059768 ");
		delta_bb_other = new Poly1DOld("1.70 2.60 7 196.42741354417250931874 -2484.78990744927750711213 13070.16612191789863572922 -36624.76264512752823065966 57723.94278221482090884820 -48573.53827961719071026891 17067.48005813081908854656 ");
		delta_sc_sp2 = new Poly1DOld("1.45 3.00 7 -15.65765464358041292314 215.33240577975831797630 -1213.59836768308628052182 3577.06985723827438050648 -5794.21906707604011899093 4869.95323566881415899843 -1651.50066037367196258856 ");
		delta_sc_sp3 = new Poly1DOld("1.60 3.00 7 -42.31970036824986181045 590.48118211852386139071 -3392.85116566297801909968 10263.53084035706160648260 -17216.76429296930655254982 15164.63863552352813712787 -5472.86050937056279508397 ");
		chi_bb_helix = new Poly1DOld("-180.00 180.00 5 0.00000000961038168880 -0.00000113871638696979 -0.00037868589809405862 0.03800930299233887255 4.66385932829456351101 ");
		phi_bb_helix = new Poly1DOld("0.26 -1.00 7 33.21247794981988477048 -37.68933459537822727725 -62.63560337026657975912 60.45920219352220925657 62.92834384918143086907 8.74080551597517541040 -0.28531776573979594769 ");
		theta_bb_helix = new Poly1DOld("0.09 -1.00 7 342.42851535798280337985 560.09726474317949396209 -5.92491183262715370716 -385.60022370399616420400 -169.36291868437217544852 1.55364191702288212404 7.51796669968754560642 ");
		chi_bb_sheet = new Poly1DOld("-180.00 180.00 9 -0.00000000000000002410 -0.00000000000000094448 0.00000000000169221646 0.00000000003755034428 -0.00000003943090905246 -0.00000029074458397694 0.00035081894740413851 0.00026085932547542504 -0.69743806191339363654 ");
		phi_bb_sheet = new Poly1DOld("0.26 -1.00 7 -30.60491834465477012373 -131.40114346116789079133 -147.82544216565378292216 -38.11140874654714849612 12.78070615011911570491 8.48280508773128438804 2.77747808337573687254 ");
		theta_bb_sheet = new Poly1DOld("0.09 -1.00 4 1.57503870188756001092 -1.06484333943499343000 7.16935628756164167186 6.74927432214103362185 ");
		chi_bb_other = new Poly1DOld("-180.00 180.00 8 0.00000000000000187922 -0.00000000000012371591 -0.00000000012264568621 0.00000000996420554609 0.00000210756502485180 -0.00020282841545786161 -0.00215309358577880439 1.11783093883408080060 ");
		//chi_bb_other.printFunction();
		phi_bb_other = new Poly1DOld("0.26 -1.00 7 -7.74234267695908595641 -87.14824460418199691958 -135.42189508877558523636 -50.14595702491587303484 16.12583756394147371793 11.80465464213827608830 0.51275398347199618865 ");
		theta_bb_other = new Poly1DOld("0.09 -1.00 5 -14.83314401229022116979 -24.89844884059162666290 -8.43802158761153009436 8.76727282946360553240 4.47130231637154551549 ");
		chi_sc_sp2_short = new Poly1DOld("0.00 180.00 6 0.00000000004304423759 0.00000001119712218952 -0.00000708716113943940 0.00068371691994235397 -0.00647416044380971518 -0.05735940642119304045 ");
		//chi_sc_sp2_short.printFunction();
		phi_sc_sp2_short = new Poly1DOld("0.34 -1.00 5 -7.87625495668602582100 -12.25817564525566183420 4.11346097819066613965 8.62576263592866787633 0.06509380929670126326 ");
		theta_sc_sp2_short = new Poly1DOld("-0.17 -1.00 5 -11.77820902027137428547 -36.43126824079791248323 -39.53971358030011629126 -8.54491815494485962290 3.45575596246188432303 ");
		chi_sc_sp2_long = new Poly1DOld("0.00 180.00 8 0.00000000000005893090 -0.00000000003909603506 0.00000001017599946156 -0.00000130600611716664 0.00008493675404931176 -0.00253493923961143871 0.01723449624692728610 0.60254004576696240658 ");
		phi_sc_sp2_long = new Poly1DOld("0.50 -1.00 4 0.85822499622818793075 3.02208096366379708186 1.85514165923868112884 -1.20125735479177331300 ");
		theta_sc_sp2_long = new Poly1DOld("0.17 -0.98 4 -1.77513856757575716472 -1.34801805450061595870 2.56092069703684765969 0.27895441899056011570 ");
		phi_sc_sp3_short = new Poly1DOld("0.17 -1.00 4 -0.92188267159473979984 9.47562101080037422207 9.03284152160650499752 -0.24658458159896870510 ");
		theta_sc_sp3_short = new Poly1DOld("-0.17 -1.00 5 6.30244240590163684601 11.10939560974528106385 5.97747799525861189807 8.43782213557545901494 4.59307407681592039239 ");
		phi_sc_sp3_long = new Poly1DOld("0.50 -1.00 4 -0.02329171617200254119 3.52753356229713066483 3.40660535785749329918 -0.88323596956894478982 ");
		theta_sc_sp3_long = new Poly1DOld("0.17 -1.00 3 1.79892178902436628007 3.18903183717485250170 -0.21928164239458927676 ");



		/*delta_bb_helix = new Poly1DOld("1.70 2.60 7 76.92576076 -920.31423728 4582.98159389 -12230.47471070 18585.45580048 -15372.59121156 5445.07452914 ");
			delta_bb_sheet = new Poly1DOld("1.70 2.60 7 144.34922152 -1884.15932311 10261.76991717 -29877.53999255 49086.73306216 -43168.39618908 15875.40397484 ");
			delta_bb_other = new Poly1DOld("1.70 2.60 7 196.42741354 -2484.78990745 13070.16612192 -36624.76264513 57723.94278221 -48573.53827962 17067.48005813 ");
			delta_sc_sp2 = new Poly1DOld("1.45 3.00 7 -15.65765464 215.33240578 -1213.59836768 3577.06985724 -5794.21906708 4869.95323567 -1651.50066037 ");
			delta_sc_sp3 = new Poly1DOld("1.60 3.00 7 -42.31970037 590.48118212 -3392.85116566 10263.53084036 -17216.76429297 15164.63863552 -5472.86050937 ");
			chi_bb_helix = new Poly1DOld("-180.00 180.00 5 0.00000001 -0.00000114 -0.00037869 0.03800930 4.66385933 ");
			phi_bb_helix = new Poly1DOld("0.26 -1.00 7 33.21247795 -37.68933460 -62.63560337 60.45920219 62.92834385 8.74080552 -0.28531777 ");
			theta_bb_helix = new Poly1DOld("0.09 -1.00 7 342.42851536 560.09726474 -5.92491183 -385.60022370 -169.36291868 1.55364192 7.51796670 ");
			chi_bb_sheet = new Poly1DOld("-180.00 180.00 9 -0.00000000 -0.00000000 0.00000000 0.00000000 -0.00000004 -0.00000029 0.00035082 0.00026086 -0.69743806 ");
			phi_bb_sheet = new Poly1DOld("0.26 -1.00 7 -30.60491834 -131.40114346 -147.82544217 -38.11140875 12.78070615 8.48280509 2.77747808 ");
			theta_bb_sheet = new Poly1DOld("0.09 -1.00 4 1.57503870 -1.06484334 7.16935629 6.74927432 ");
			chi_bb_other = new Poly1DOld("-180.00 180.00 13 0.00000000 -0.00000000 -0.00000000 0.00000000 -0.00000000 -0.00000000 0.00000000 0.00000000 -0.00000008 0.00000108 0.00023680 -0.00043897 0.78625239 ");
			phi_bb_other = new Poly1DOld("0.26 -1.00 7 -7.74234268 -87.14824460 -135.42189509 -50.14595702 16.12583756 11.80465464 0.51275398 ");
			theta_bb_other = new Poly1DOld("0.09 -1.00 5 -14.83314401 -24.89844884 -8.43802159 8.76727283 4.47130232 ");
			chi_sc_sp2_short = new Poly1DOld("0.00 180.00 6 0.00000000 0.00000001 -0.00000709 0.00068372 -0.00647416 -0.05735941 ");
			phi_sc_sp2_short = new Poly1DOld("0.34 -1.00 5 -7.87625496 -12.25817565 4.11346098 8.62576264 0.06509381 ");
			theta_sc_sp2_short = new Poly1DOld("-0.17 -1.00 5 -11.77820902 -36.43126824 -39.53971358 -8.54491815 3.45575596 ");
			chi_sc_sp2_long = new Poly1DOld("0.00 180.00 8 0.00000000 -0.00000000 0.00000001 -0.00000131 0.00008494 -0.00253494 0.01723450 0.60254005 ");
			phi_sc_sp2_long = new Poly1DOld("0.50 -1.00 4 0.85822500 3.02208096 1.85514166 -1.20125735 ");
			theta_sc_sp2_long = new Poly1DOld("0.17 -0.98 4 -1.77513857 -1.34801805 2.56092070 0.27895442 ");
			phi_sc_sp3_short = new Poly1DOld("0.17 -1.00 4 -0.92188267 9.47562101 9.03284152 -0.24658458 ");
			theta_sc_sp3_short = new Poly1DOld("-0.17 -1.00 5 6.30244241 11.10939561 5.97747800 8.43782214 4.59307408 ");
			phi_sc_sp3_long = new Poly1DOld("0.50 -1.00 4 -0.02329172 3.52753356 3.40660536 -0.88323597 ");
			theta_sc_sp3_long = new Poly1DOld("0.17 -1.00 3 1.79892179 3.18903184 -0.21928164 ");*/



	}


	public void calculateHBondEnergy(double[] coordinates, int curIndex,
			double[] energyTerms, Molecule m,boolean updateTerms,Amber96ext amb96ff) {

		int atomix3, atomjx3, atomi, atomj, atomk, atomkx3;//, atoml, atomlx3, atomm, atommx3;
		int atomb=0, atombx3=0, atomb2=0, atomb2x3=0; //atom bases
		double rij, rij2;
		double cx, cy, cz, cmag;
		double dx, dy, dz, ex, ey, ez, dmag, emag;
		double fx, fy, fz, fmag;
		double ux, uy, uz, vx, vy, vz, umag, vmag;
		double costheta;
		double cosPhi, Ehb;
		byte acceptType;
		byte donorSS, acceptSS;

		double chi=0;



		ArrayList<HbondPair> hbonds;

		if(curIndex == -1){
			hbonds = hbondTerms;
		}
		else{
			hbonds = partHBlist[curIndex];
		}

		Iterator<HbondPair> iter = hbonds.iterator();

		while(iter.hasNext()){
			HbondPair hbp = iter.next();

			/*if(hbp.donor.isBBatom || hbp.accept.isBBatom)
					continue;

				if(hbp.donor.strandNumber == hbp.accept.strandNumber)
					continue;*/

			//determine donor acceptor type
			//donorType = getHybridization(hbp.donor.forceFieldType);
			acceptType = getHybridization(hbp.accept.forceFieldType);

			donorSS = m.residue[hbp.donor.moleculeResidueNumber].SStype;
			acceptSS = m.residue[hbp.accept.moleculeResidueNumber].SStype;




			atomi = hbp.donor.moleculeAtomNumber;
			atomj = hbp.hydro.moleculeAtomNumber;
			atomk = hbp.accept.moleculeAtomNumber;


			atomix3 = atomi * 3;
			atomjx3 = atomj * 3;
			atomkx3 = atomk * 3;


			//H-A
			dx = coordinates[atomjx3] - coordinates[atomkx3];
			dy = coordinates[atomjx3 + 1] - coordinates[atomkx3 + 1];
			dz = coordinates[atomjx3 + 2] - coordinates[atomkx3 + 2];
			rij2 = dx*dx + dy*dy + dz*dz;
			dmag = Math.sqrt(rij2);

			rij = dmag;

			if(dmag < 1.4 || dmag > 3.0){
				continue; //Only calculate if 1.4<= R <= 3.0 angstroms
			}

			//OH shouldn't be an acceptor unless the lone pair is pointing toward the donor H.
			boolean badPair = false;
			if(hbp.accept.elementType.equals("O")){ 
				for(int i1 = 0; i1<hbp.accept.bond.length;i1++){
					if((m.atom[hbp.accept.bond[i1]].elementType).equals("H")){
						//Bad if donor is closer to H than to acceptor
						if(hbp.donor.distance(m.atom[hbp.accept.bond[i1]]) < hbp.donor.distance(hbp.accept)){
							badPair = true;
						}
					}
				}
			}
			if(badPair)
				continue;

			double Edelta = 0.0;
			double sDelta = 0.0; //Smoothing term

			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					Edelta = delta_bb_helix.getVal(rij);
					sDelta = delta_bb_helix.getSmoothVal(rij);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					Edelta = delta_bb_sheet.getVal(rij);
					sDelta = delta_bb_sheet.getSmoothVal(rij);
				}
				else{
					Edelta = delta_bb_other.getVal(rij);
					sDelta = delta_bb_other.getSmoothVal(rij);
				}
			}
			else if(acceptType == SP2 || acceptType == RING){
				Edelta = delta_sc_sp2.getVal(rij);
				sDelta = delta_sc_sp2.getSmoothVal(rij);
			}
			else if(acceptType == SP3){
				Edelta = delta_sc_sp3.getVal(rij);
				sDelta = delta_sc_sp3.getSmoothVal(rij);
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(Edelta == Double.POSITIVE_INFINITY)
				continue;


			//rij6 = rij2 * rij2 * rij2;
			//rij10 = rij2*rij2*rij6;
			//rij12 = rij6 * rij6;



			//D-H
			cx = coordinates[atomix3] - coordinates[atomjx3];
			cy = coordinates[atomix3 + 1] - coordinates[atomjx3 + 1];
			cz = coordinates[atomix3 + 2] - coordinates[atomjx3 + 2];
			cmag = Math.sqrt(cx*cx + cy*cy + cz*cz);


			//angle = (float)Math.toDegrees(Math.acos((dx*ex + dy*ey + dz*ez) / (dmag * -emag)));
			costheta = (cx*dx + cy*dy + cz*dz) / (cmag * -dmag);
			//theta = Math.toDegrees(Math.acos(costheta));

			double Etheta = 0.0;
			double sTheta = 1.0; //Smoothing term

			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					Etheta = theta_bb_helix.getVal(costheta);
					sTheta = theta_bb_helix.getSmoothVal(costheta);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					Etheta = theta_bb_sheet.getVal(costheta);
					sTheta = theta_bb_sheet.getSmoothVal(costheta);
				}
				else{
					Etheta = theta_bb_other.getVal(costheta);
					sTheta = theta_bb_other.getSmoothVal(costheta);
				}
			}
			else if(rij < 2.1){
				if(acceptType == SP2 || acceptType == RING){
					Etheta = theta_sc_sp2_short.getVal(costheta);
					sTheta = theta_sc_sp2_short.getSmoothVal(costheta);
				}
				else if(acceptType == SP3){
					Etheta = theta_sc_sp3_short.getVal(costheta);
					sTheta = theta_sc_sp3_short.getSmoothVal(costheta);
				}
			}
			else if(rij >= 2.1){
				if(acceptType == SP2 || acceptType == RING){
					Etheta = theta_sc_sp2_long.getVal(costheta);
					sTheta = theta_sc_sp2_long.getSmoothVal(costheta);
				}
				else if(acceptType == SP3){
					Etheta = theta_sc_sp3_long.getVal(costheta);
					sTheta = theta_sc_sp3_long.getSmoothVal(costheta);
				}
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(Etheta == Double.POSITIVE_INFINITY)
				continue;

			//Determine Phi
			//Hydro-Accept-Base
			Atom base = null;
			Atom base2 = null;

			for(int i=0; i<hbp.accept.bond.length;i++){
				if(!m.atom[hbp.accept.bond[i]].elementType.equals("H")){ //Want a heavy atom
					if(base == null){
						base = m.atom[hbp.accept.bond[i]];
						atomb = base.moleculeAtomNumber;
						atombx3 = atomb*3;
					}
					else{
						base2 = m.atom[hbp.accept.bond[i]];
						atomb2 = base2.moleculeAtomNumber;
						atomb2x3 = atomb2*3;
					}
				}else if(base2 == null){
					base2 = m.atom[hbp.accept.bond[i]];
					atomb2 = base2.moleculeAtomNumber;
					atomb2x3 = atomb2*3;
				}
			}


			//Accept-Base
			if(hbp.accept.hybridization == RING){ //If acceptor is a ring, then we average the two base atoms
				ex = coordinates[atomkx3] - 0.5*(coordinates[atombx3]+coordinates[atomb2x3]);
				ey = coordinates[atomkx3 + 1] - 0.5*(coordinates[atombx3 + 1]+coordinates[atomb2x3 + 1]);
				ez = coordinates[atomkx3 + 2] - 0.5*(coordinates[atombx3 + 2]+coordinates[atomb2x3 + 2]);
			} else{
				ex = coordinates[atomkx3] - coordinates[atombx3];
				ey = coordinates[atomkx3 + 1] - coordinates[atombx3 + 1];
				ez = coordinates[atomkx3 + 2] - coordinates[atombx3 + 2];
			}
			emag = Math.sqrt(ex*ex + ey*ey + ez*ez);

			cosPhi = (dx*ex + dy*ey + dz*ez) / (dmag * -emag);
			//phi = Math.toDegrees(Math.acos( cosPhi ));

			double Ephi = 0.0;
			double sPhi = 1.0; //Smoothing term

			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					Ephi = phi_bb_helix.getVal(cosPhi);
					sPhi = phi_bb_helix.getSmoothVal(cosPhi);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					Ephi = phi_bb_sheet.getVal(cosPhi);
					sPhi = phi_bb_sheet.getSmoothVal(cosPhi);
				}
				else{
					Ephi = phi_bb_other.getVal(cosPhi);
					sPhi = phi_bb_other.getSmoothVal(cosPhi);
				}
			}
			else if(rij < 2.1){
				if(acceptType == SP2 || acceptType == RING){
					Ephi = phi_sc_sp2_short.getVal(cosPhi);
					sPhi = phi_sc_sp2_short.getSmoothVal(cosPhi);
				}
				else if(acceptType == SP3){
					Ephi = phi_sc_sp3_short.getVal(cosPhi);
					sPhi = phi_sc_sp3_short.getSmoothVal(cosPhi);
				}
			}
			else if(rij >= 2.1){
				if(acceptType == SP2 || acceptType == RING){
					Ephi = phi_sc_sp2_long.getVal(cosPhi);
					sPhi = phi_sc_sp2_long.getSmoothVal(cosPhi);
				}
				else if(acceptType == SP3){
					Ephi = phi_sc_sp3_long.getVal(cosPhi);
					sPhi = phi_sc_sp3_long.getSmoothVal(cosPhi);
				}
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(Ephi == Double.POSITIVE_INFINITY)
				continue;	

			double Echi = 0.0;

			//Determine chi
			//Not calculating chi right now
			if(acceptType == SP2  || acceptType == RING){
				if(base2 == null){
					for(int i=0; i<base.bond.length;i++){
						if(m.atom[base.bond[i]].elementType.equals("C") && m.atom[base.bond[i]].moleculeResidueNumber == base.moleculeResidueNumber 
								&& m.atom[base.bond[i]].moleculeAtomNumber != atomk){ //Doesn't equal the acceptor
							base2 = m.atom[base.bond[i]];
						}
					}	
					if(base2 == null){
						for(int i=0; i<base.bond.length;i++){
							if(m.atom[base.bond[i]].elementType.equals("N") && m.atom[base.bond[i]].moleculeResidueNumber == base.moleculeResidueNumber
									&& m.atom[base.bond[i]].moleculeAtomNumber != atomk){
								base2 = m.atom[base.bond[i]];
							}
						}
					}
					if(base2 == null){ //if it's still equal to null we find the connected atom with the most bonds
						int maxBonds = -1; //if tied we should recurse, but I don't do that yet
						int maxAtom = -1;
						for(int i=0; i<base.bond.length;i++){
							if(m.atom[base.bond[i]].bond.length > maxBonds){
								maxBonds = m.atom[base.bond[i]].bond.length;
								maxAtom = base.bond[i];
							}
						}
						base2 = m.atom[maxAtom];

					}
					atomb2 = base2.moleculeAtomNumber;
					atomb2x3 = atomb2*3;
					}

				//Accept-Base
				fx = coordinates[atombx3] - coordinates[atomb2x3];
				fy = coordinates[atombx3 + 1] - coordinates[atomb2x3 + 1];
				fz = coordinates[atombx3 + 2] - coordinates[atomb2x3 + 2];
				fmag = Math.sqrt(fx*fx + fy*fy + fz*fz);


				// Cross product: a x b = (aybz-azby, -axbz+azbx, axby-aybx)
				// 'u' and 'v' are normals to planes
				// u = f x e, v = e x d
				ux = fy*ez - fz*ey;
				uy = fz*ex - fx*ez;
				uz = fx*ey - fy*ex;
				umag = Math.sqrt(ux*ux + uy*uy + uz*uz);
				vx = ey*dz - ez*dy;
				vy = ez*dx - ex*dz;
				vz = ex*dy - ey*dx;
				vmag = Math.sqrt(vx*vx + vy*vy + vz*vz);
				// Dot product again
				chi = Math.toDegrees(Math.acos((ux*vx + uy*vy + uz*vz) / (umag * vmag)));

				// BUT, that doesn't solve the handedness (sign) problem for the dihedral!
				// To do that, we look at the angle between 'f' and 'v'
				// Dot product again
				if( Math.toDegrees(Math.acos((fx*vx + fy*vy + fz*vz) / (fmag * vmag))) > 90.0 )
				{ chi = -chi; }

				if(hbp.donor.isBBatom && hbp.accept.isBBatom){
					if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
						Echi = chi_bb_helix.getVal(chi);
					}
					else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
						Echi = chi_bb_sheet.getVal(chi);
					}
					else{
						Echi = chi_bb_other.getVal(chi);
					}
				}
				else if(rij < 2.1){
					if(chi < 0){chi = -chi;}
					if(acceptType == SP2 || acceptType == RING){
						Echi = chi_sc_sp2_short.getVal(chi);
					}	
				}
				else if(rij >= 2.1){
					if(chi < 0){chi = -chi;}
					if(acceptType == SP2 || acceptType == RING){
						Echi = chi_sc_sp2_long.getVal(chi);
					}
				}
				else{
					System.out.println("Unrecognized Histogram");
					System.exit(0);
				}

				if(Echi == Double.POSITIVE_INFINITY)
					continue;


			}


			Ehb = hbp.multiplier * sDelta * sTheta * sPhi * (Edelta+Etheta+Ephi+Echi);

			//Ehb = Echi;

			//				if(hbp.donor.strandNumber != hbp.accept.strandNumber){
			if(debug){
				System.out.print("Ehb: "+Ehb+" "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+hbp.donor.name+" "+m.residue[hbp.accept.moleculeResidueNumber].fullName+" "+hbp.accept.name+" ");
				System.out.println(Edelta+" ( "+rij+" ) "+Etheta+" ( "+costheta+" ) "+Ephi+" ( "+cosPhi+" ) "+Echi+" ( "+chi+" ) ");
			}
			//				}
			if(updateTerms)
				amb96ff.updatePairTerms(hbondScale*Ehb,hbp.donor.moleculeResidueNumber,hbp.accept.moleculeResidueNumber);
			energyTerms[4] += Ehb; 
		}

		energyTerms[4] *= hbondScale;


	}

	public void calculateHBondGradient(int curIndex, Molecule m) {



		int atomDx3, atomHx3, atomD, atomH, atomA, atomAx3, atomC, atomCx3;
		int atomB=0, atomBx3=0, atomB2=0, atomB2x3=0;
		double rij, rij2;
		double DHx, DHy, DHz, DHmag;
		double BCx, BCy, BCz, BCmag;
		double HAx, HAy, HAz, ABx, ABy, ABz, HAmag, ABmag;
		double cosTheta;
		double[] A = new double[3];
		double[] B = new double[3];
		double Amag,Bmag;
		double cosPhi;
		byte acceptType;
		byte donorSS, acceptSS;
		double[] DH = new double[3];
		double[] HA = new double[3];
		double[] AB = new double[3];

		double chi=0;

		double[] HAvalDeriv = new double[2];
		double[] cosThetaValDeriv = new double[2];
		double[] cosPhiValDeriv = new double[2];
		double[] chiValDeriv = new double[2];
		double[] dCosThetaDR = new double[6];
		double[] dCosPhiDR = new double[6];
		double[] dChiDR = new double[12];


		ArrayList<HbondPair> hbonds;

		if(curIndex == -1){
			hbonds = hbondTerms;
		}
		else{
			hbonds = partHBlist[curIndex];
		}

		Iterator<HbondPair> iter = hbonds.iterator(); 
		while(iter.hasNext()){
			HbondPair hbp = iter.next();

			/*if(hbp.donor.isBBatom || hbp.accept.isBBatom)
				continue;

			if(hbp.donor.strandNumber == hbp.accept.strandNumber)
				continue;*/

			//determine donor acceptor type
			//donorType = getHybridization(hbp.donor.forceFieldType);
			acceptType = getHybridization(hbp.accept.forceFieldType);

			donorSS = m.residue[hbp.donor.moleculeResidueNumber].SStype;
			acceptSS = m.residue[hbp.accept.moleculeResidueNumber].SStype;


			atomD = hbp.donor.moleculeAtomNumber;
			atomH = hbp.hydro.moleculeAtomNumber;
			atomA = hbp.accept.moleculeAtomNumber;


			atomDx3 = atomD * 3;
			atomHx3 = atomH * 3;
			atomAx3 = atomA * 3;


			//H-A
			HAx = m.actualCoordinates[atomHx3] - m.actualCoordinates[atomAx3];
			HAy = m.actualCoordinates[atomHx3 + 1] - m.actualCoordinates[atomAx3 + 1];
			HAz = m.actualCoordinates[atomHx3 + 2] - m.actualCoordinates[atomAx3 + 2];
			rij2 = HAx*HAx + HAy*HAy + HAz*HAz;
			HAmag = Math.sqrt(rij2);

			rij = HAmag;

			if(HAmag < 1.4 || HAmag > 3.0){
				continue; //Only calculate if 1.4<= R <= 3.0 angstroms
			}

			//OH shouldn't be an acceptor unless the lone pair is pointing toward the donor H.
			boolean badPair = false;
			if(hbp.accept.elementType.equals("O")){ 
				for(int i1 = 0; i1<hbp.accept.bond.length;i1++){
					if((m.atom[hbp.accept.bond[i1]].elementType).equals("H")){
						//Bad if donor is closer to H than to acceptor
						if(hbp.donor.distance(m.atom[hbp.accept.bond[i1]]) < hbp.donor.distance(hbp.accept)){
							badPair = true;
						}
					}
				}
			}
			if(badPair)
				continue;

			//double Edelta = 0.0;
			double sDelta = 0.0; //Smoothing term

			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					delta_bb_helix.getValueAndDeriv(rij,HAvalDeriv);
					sDelta = delta_bb_helix.getSmoothVal(rij);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					delta_bb_sheet.getValueAndDeriv(rij,HAvalDeriv);
					sDelta = delta_bb_sheet.getSmoothVal(rij);
				}
				else{
					delta_bb_other.getValueAndDeriv(rij,HAvalDeriv);
					sDelta = delta_bb_other.getSmoothVal(rij);
				}
			}
			else if(acceptType == SP2 || acceptType == RING){
				delta_sc_sp2.getValueAndDeriv(rij,HAvalDeriv);
				sDelta = delta_sc_sp2.getSmoothVal(rij);
			}
			else if(acceptType == SP3){
				delta_sc_sp3.getValueAndDeriv(rij,HAvalDeriv);
				sDelta = delta_sc_sp3.getSmoothVal(rij);
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(HAvalDeriv[Poly1DOld.DERIV] == 0)
				continue;







			//D-H
			DHx = m.actualCoordinates[atomDx3] - m.actualCoordinates[atomHx3];
			DHy = m.actualCoordinates[atomDx3 + 1] - m.actualCoordinates[atomHx3 + 1];
			DHz = m.actualCoordinates[atomDx3 + 2] - m.actualCoordinates[atomHx3 + 2];
			DHmag = Math.sqrt(DHx*DHx + DHy*DHy + DHz*DHz);


			//angle = (float)Math.toDegrees(Math.acos((dx*ex + dy*ey + dz*ez) / (dmag * -emag)));
			cosTheta = dot(DHx,DHy,DHz,HAx,HAy,HAz) / (DHmag * -HAmag);

			double sTheta = 1.0; //Smoothing term


			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					theta_bb_helix.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_bb_helix.getSmoothVal(cosTheta);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					theta_bb_sheet.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_bb_sheet.getSmoothVal(cosTheta);
				}
				else{
					theta_bb_other.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_bb_other.getSmoothVal(cosTheta);
				}
			}
			else if(rij < 2.1){
				if(acceptType == SP2 || acceptType == RING){
					theta_sc_sp2_short.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_sc_sp2_short.getSmoothVal(cosTheta);
				}
				else if(acceptType == SP3){
					theta_sc_sp3_short.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_sc_sp3_short.getSmoothVal(cosTheta);
				}
			}
			else if(rij >= 2.1){
				if(acceptType == SP2 || acceptType == RING){
					theta_sc_sp2_long.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_sc_sp2_long.getSmoothVal(cosTheta);
				}
				else if(acceptType == SP3){
					theta_sc_sp3_long.getValueAndDeriv(cosTheta,cosThetaValDeriv);
					sTheta = theta_sc_sp3_long.getSmoothVal(cosTheta);
				}
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(cosThetaValDeriv[Poly1DOld.DERIV] == 0.0)
				continue;


			DH[0] = DHx;DH[1] = DHy;DH[2] = DHz;
			HA[0] = HAx;HA[1] = HAy;HA[2] = HAz;


			dCosThetaDR = cosThetaDR1(DH,DHmag,HA,HAmag);




			//Determine Phi
			//Hydro-Accept-Base
			Atom base = null;
			Atom base2 = null;

			for(int i=0; i<hbp.accept.bond.length;i++){
				if(!m.atom[hbp.accept.bond[i]].elementType.equals("H")){ //Want a heavy atom
					if(base == null){
						base = m.atom[hbp.accept.bond[i]];
						atomB = base.moleculeAtomNumber;
						atomBx3 = atomB*3;
					}
					else{
						base2 = m.atom[hbp.accept.bond[i]];
						atomB2 = base2.moleculeAtomNumber;
						atomB2x3 = atomB2*3;
					}
				}else if(base2 == null){
					base2 = m.atom[hbp.accept.bond[i]];
					atomB2 = base2.moleculeAtomNumber;
					atomB2x3 = atomB2*3;
				}
			}

			atomB = base.moleculeAtomNumber;
			atomBx3 = atomB*3;

			//Accept-Base
			if(hbp.accept.hybridization == RING){
				ABx = m.actualCoordinates[atomAx3] - 0.5*(m.actualCoordinates[atomBx3]+m.actualCoordinates[atomB2x3]);
				ABy = m.actualCoordinates[atomAx3 + 1] - 0.5*(m.actualCoordinates[atomBx3 + 1]+m.actualCoordinates[atomB2x3 + 1]);
				ABz = m.actualCoordinates[atomAx3 + 2] - 0.5*(m.actualCoordinates[atomBx3 + 2]+m.actualCoordinates[atomB2x3 + 2]);
			} else{
				ABx = m.actualCoordinates[atomAx3] - m.actualCoordinates[atomBx3];
				ABy = m.actualCoordinates[atomAx3 + 1] - m.actualCoordinates[atomBx3 + 1];
				ABz = m.actualCoordinates[atomAx3 + 2] - m.actualCoordinates[atomBx3 + 2];
			}
			
			//Accept-Base
			ABmag = Math.sqrt(ABx*ABx + ABy*ABy + ABz*ABz);


			cosPhi = (HAx*ABx + HAy*ABy + HAz*ABz) / (HAmag * -ABmag) ;

			double sPhi = 1.0; //Smoothing term

			if(hbp.donor.isBBatom && hbp.accept.isBBatom){
				if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
					phi_bb_helix.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_bb_helix.getSmoothVal(cosPhi);
				}
				else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
					phi_bb_sheet.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_bb_sheet.getSmoothVal(cosPhi);
				}
				else{
					phi_bb_other.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_bb_other.getSmoothVal(cosPhi);
				}
			}
			else if(rij < 2.1){
				if(acceptType == SP2 || acceptType == RING){
					phi_sc_sp2_short.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_sc_sp2_short.getSmoothVal(cosPhi);
				}
				else if(acceptType == SP3){
					phi_sc_sp3_short.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_sc_sp3_short.getSmoothVal(cosPhi);
				}
			}
			else if(rij >= 2.1){
				if(acceptType == SP2 || acceptType == RING){
					phi_sc_sp2_long.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_sc_sp2_long.getSmoothVal(cosPhi);
				}
				else if(acceptType == SP3){
					phi_sc_sp3_long.getValueAndDeriv(cosPhi,cosPhiValDeriv);
					sPhi = phi_sc_sp3_long.getSmoothVal(cosPhi);
				}
			}
			else{
				System.out.println("Unrecognized Histogram");
				System.exit(0);
			}

			if(cosPhiValDeriv[Poly1DOld.DERIV] == 0.0)
				continue;	


			AB[0] = ABx;AB[1] = ABy;AB[2] = ABz;


			dCosPhiDR = cosThetaDR1(HA,HAmag,AB,ABmag);



			//Precompute hbondScale*multiplier
			double newMult = hbondScale*hbp.multiplier*sDelta*sTheta*sPhi;


			//Determine chi
			//Using derivative formulation from 
			if(acceptType == SP2 || acceptType == RING ){
				Atom R1 = null;
				for(int i=0; i<base.bond.length;i++){
					if(m.atom[base.bond[i]].elementType.equals("C") && m.atom[base.bond[i]].moleculeResidueNumber == base.moleculeResidueNumber){
						R1 = m.atom[base.bond[i]];
					}
				}	
				if(R1 == null){
					for(int i=0; i<base.bond.length;i++){
						if(m.atom[base.bond[i]].elementType.equals("N") && m.atom[base.bond[i]].moleculeResidueNumber == base.moleculeResidueNumber){
							R1 = m.atom[base.bond[i]];
						}
					}
				}
				if(R1 == null){ //if it's still equal to null we find the connected atom with the most bonds
					int maxBonds = -1; //if tied we should recurse, but I don't do that yet
					int maxAtom = -1;
					for(int i=0; i<base.bond.length;i++){
						if(m.atom[base.bond[i]].bond.length > maxBonds){
							maxBonds = m.atom[base.bond[i]].bond.length;
							maxAtom = base.bond[i];
						}
					}
					R1 = m.atom[maxAtom];

				}

				atomC = R1.moleculeAtomNumber;
				atomCx3 = atomC*3;

				//Base(B) - Base2(C)
				BCx = m.actualCoordinates[atomBx3] - m.actualCoordinates[atomCx3];
				BCy = m.actualCoordinates[atomBx3 + 1] - m.actualCoordinates[atomCx3 + 1];
				BCz = m.actualCoordinates[atomBx3 + 2] - m.actualCoordinates[atomCx3 + 2];
				BCmag = Math.sqrt(BCx*BCx + BCy*BCy + BCz*BCz);


				// Cross product: a x b = (aybz-azby, -axbz+azbx, axby-aybx)
				// 'u' and 'v' are normals to planes
				// u = f x e, v = e x d

				double[] BC = {BCx,BCy,BCz};
				double[] CB = {-BCx,-BCy,-BCz};
				double[] BA = {-ABx,-ABy,-ABz};

				A = cross(CB,BA);
				B = cross(HA,BA);

				Amag = Math.sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);

				Bmag = Math.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
				// Dot product again
				chi = Math.toDegrees(Math.acos(dot(A,B) / (Amag * Bmag)));

				// BUT, that doesn't solve the handedness (sign) problem for the dihedral!
				// To do that, we look at the angle between 'f' and 'v'
				// Dot product again
				if( Math.toDegrees(Math.acos(dot(B,BC) / (BCmag * Bmag))) > 90.0 )
				{ chi = -chi; }

				if(hbp.donor.isBBatom && hbp.accept.isBBatom){
					if(donorSS == Molecule.HELIX && acceptSS == Molecule.HELIX){
						chi_bb_helix.getValueAndDeriv(chi,chiValDeriv);
					}
					else if(donorSS == Molecule.STRAND && acceptSS == Molecule.STRAND){
						chi_bb_sheet.getValueAndDeriv(chi,chiValDeriv);
					}
					else{
						chi_bb_other.getValueAndDeriv(chi,chiValDeriv);
					}
				}
				else if(rij < 2.1){
					if(chi < 0){chi = -chi;}
					if(acceptType == SP2 || acceptType == RING){
						chi_sc_sp2_short.getValueAndDeriv(chi,chiValDeriv);
					}	
				}
				else if(rij >= 2.1){
					if(chi < 0){chi = -chi;}
					if(acceptType == SP2 || acceptType == RING){
						chi_sc_sp2_long.getValueAndDeriv(chi,chiValDeriv);
					}
				}
				else{
					System.out.println("Unrecognized Histogram");
					System.exit(0);
				}

				if(chiValDeriv[Poly1DOld.DERIV] == 0.0)
					continue;


				//Compute dChi/dR for all four atoms
				double Asquare = dot(A,A);
				double Bsquare = dot(B,B);
				double term1 = ABmag/Asquare;  //-|G|/A^2
				dChiDR[0] = -term1*A[0];
				dChiDR[1] = -term1*A[1];
				dChiDR[2] = -term1*A[2];

				double term2 = dot(CB,BA)/(Asquare*ABmag);
				double term3 = dot(HA,BA)/(Bsquare*ABmag); 

				dChiDR[3] = dChiDR[0]+term2 *A[0]-term3*B[0];
				dChiDR[4] = dChiDR[1]+term2 *A[1]-term3*B[1];
				dChiDR[5] = dChiDR[2]+term2 *A[2]-term3*B[2];

				double term4 = ABmag/Bsquare;

				dChiDR[6] = term3*B[0]-term2*A[0]-term4*B[0];
				dChiDR[7] = term3*B[1]-term2*A[1]-term4*B[1];
				dChiDR[8] = term3*B[2]-term2*A[2]-term4*B[2];

				dChiDR[9] = term4*B[0];
				dChiDR[10] = term4*B[1];
				dChiDR[11] = term4*B[2];


				//System.out.println("Hbond: "+m.residue[m.atom[atomH].moleculeResidueNumber].fullName+" "+ m.atom[atomH].name+" "+m.residue[m.atom[atomA].moleculeResidueNumber].fullName+" "+ m.atom[atomA].name);


				//Add chi gradient
				m.gradient[atomCx3] += dChiDR[0]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomCx3+1] += dChiDR[1]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomCx3+2] += dChiDR[2]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomBx3] += dChiDR[3]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomBx3+1] += dChiDR[4]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomBx3+2] += dChiDR[5]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomAx3] += dChiDR[6]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomAx3+1] += dChiDR[7]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomAx3+2] += dChiDR[8]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomHx3] += dChiDR[9]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomHx3+1] += dChiDR[10]*chiValDeriv[Poly1DOld.DERIV]*newMult;
				m.gradient[atomHx3+2] += dChiDR[1]*chiValDeriv[Poly1DOld.DERIV]*newMult;


			
			}

			//Add AHdist gradient
			m.gradient[atomHx3] += HAx*HAvalDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomHx3+1] += HAy*HAvalDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomHx3+2] += HAz*HAvalDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3] += -HAx*HAvalDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3+1] += -HAy*HAvalDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3+2] += -HAz*HAvalDeriv[Poly1DOld.DERIV]*newMult;

			//Add cosTheta gradient
			m.gradient[atomDx3] += dCosThetaDR[0]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomDx3+1] += dCosThetaDR[1]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomDx3+2] += dCosThetaDR[2]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3] += dCosThetaDR[3]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3+1] += dCosThetaDR[4]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomAx3+2] += dCosThetaDR[5]*cosThetaValDeriv[Poly1DOld.DERIV]*newMult;


			//Add cosPhi gradient
			m.gradient[atomHx3] += dCosPhiDR[0]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomHx3+1] += dCosPhiDR[1]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomHx3+2] += dCosPhiDR[2]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomBx3] += dCosPhiDR[3]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomBx3+1] += dCosPhiDR[4]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;
			m.gradient[atomBx3+2] += dCosPhiDR[5]*cosPhiValDeriv[Poly1DOld.DERIV]*newMult;


			//if(Ehb<-1 && hbp.donor.strandNumber != hbp.accept.strandNumber)
			//	System.out.println("Ehb: "+Ehb+" "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+hbp.donor.name+" "+m.residue[hbp.accept.moleculeResidueNumber].fullName+" "+hbp.accept.name); 
		}


	}

	//ffType: forcefield type
	private byte getHybridization(String ffType) {
		//Check SP2
		if(ffType.equals("N") || ffType.equals("O") || ffType.equals("O2") || ffType.equals("NA") ||
				ffType.equals("NC") || ffType.equals("N*") || ffType.equals("N2") || ffType.equals("NT")){
			return SP2;
		}
		//Check SP3
		else if(ffType.equals("N3") || ffType.equals("OH") || ffType.equals("OW") || ffType.equals("OS") ){
			return SP3;
		}
		else if(ffType.equals("NB")){
			return RING;
		}
		else{
			System.out.println("Do not recognize hybridization type of atom: "+ffType);
			System.exit(0);
			return -1;
		}
	}

	private class HbondPair{
		Atom donor;
		Atom hydro;
		Atom accept;
		double multiplier;

		HbondPair(Atom d, Atom h, Atom a,double multiplier){
			donor = d;
			hydro = h;
			accept = a;
			this.multiplier = multiplier;
		}
	}

	private class SSPair{
		String residue;
		String SStype;


		SSPair(String res, String type){
			residue = res;
			SStype = type;
		}
	}

	/******************************/
	// This function returns the xth token in string s
	private String getToken(String s, int x) {

		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
				st.nextToken();
			else {
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	public static double dot(double v1x, double v1y, double v1z, 
			double v2x, double v2y, double v2z){

		return v1x*v2x+v1y*v2y+v1z*v2z;

	}

	public static double dot(double[] a, double[] b) {
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}

	public static double[] cross(double[] v1, double[] v2){
		double[] ret = new double[3];
		ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
		ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
		return ret;
	}

	public static float[] cross(float[] v1, float[] v2){
		float[] ret = new float[3];
		ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
		ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
		return ret;
	}

	//Calculate the gradient of cosTheta for atom 1 and atom 2
	//See Tuzun Noid and Sumpter 2000
	private double[] cosThetaDR1(double[] r12,double r12mag, double[] r23, double r23mag) {
		double[] deriv = new double[6];

		double[] r21hat = new double[3];
		double[] r32hat = new double[3];

		for(int i=0; i<3;i++){
			r21hat[i] = -r12[i]/r12mag;
			r32hat[i] = -r23[i]/r23mag;
		}

		double[] num1 = cross(r32hat, cross(r21hat,r32hat));
		double[] num2 = cross(r21hat, cross(r32hat,r21hat));


		for(int i=0; i<3;i++){
			deriv[i] = num1[i]/r12mag;
			deriv[i+3] = num2[i]/r23mag;
		}


		return deriv;
	}


	public void initializeHbondTerms() {
		hbondTerms = new ArrayList<HbondPair>();
	}


	public void addIfValid(int atom1, int atom2,int res1, int res2,Molecule m, double multiplier) {
		boolean badPair = false;
		Atom hydro;
		Atom accept;
		Atom donor;
		//KER: I'm going to require that the atoms be from different residues. Not sure if that's correct
		//KER: But it seems unlikely that a residue can h-bond to itself
		if(m.atom[atom1].moleculeResidueNumber == m.atom[atom2].moleculeResidueNumber)
			return;

		if(m.atom[atom1].elementType.equals("H")){
			hydro = m.atom[atom1];
			accept = m.atom[atom2];
		}
		else if(m.atom[atom2].elementType.equals("H")){
			hydro = m.atom[atom2];
			accept = m.atom[atom1];
		}
		else{return;}

		//check donor
		if(m.atom[hydro.bond[0]].elementType.equals("O") || m.atom[hydro.bond[0]].elementType.equals("N")){
			donor = m.atom[hydro.bond[0]];
		} else{return;}

		//check accept
		if(accept.elementType.equals("O") || accept.elementType.equals("N")){
			//N acceptor won't have H bound, O can always be acceptor
			if(accept.elementType.equals("N")){
				for(int i1 = 0; i1<accept.bond.length;i1++){
					if((m.atom[accept.bond[i1]].elementType).equals("H")){
						badPair = true;
					}
				}
			}
		} else{return;}

		//If we've reached here then this is a valid Hbond pair
		if(!badPair )
			hbondTerms.add(new HbondPair(donor,hydro,accept,1.0)); //1.0 interface multiplier since not in this version

	}


	public void setupPartialHBondArrays(int numRows, int maxNumColumns,
			int[][] atomList, int[] numColumns, Molecule m) {

		partHBlist = new ArrayList[numRows];

		for(int q=0; q<numRows;q++){
			partHBlist[q] = new ArrayList<HbondPair>();
			int[] tempAtomList = new int[m.numberOfAtoms];
			for(int i=0;i<numColumns[q];i++)
				tempAtomList[atomList[q][i]] = 1;

			Iterator<HbondPair> iter = hbondTerms.iterator();
			while(iter.hasNext()){
				HbondPair hbp = iter.next();
				if(tempAtomList[hbp.hydro.moleculeAtomNumber] + tempAtomList[hbp.accept.moleculeAtomNumber] > 0){
					partHBlist[q].add(hbp);
				}
			}

		}

	}

}
