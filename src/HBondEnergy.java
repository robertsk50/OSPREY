import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;

public class HBondEnergy {

	public static final byte SP2_HYBRID = 1;
	public static final byte SP3_HYBRID = 2;
	public static final byte RING_HYBRID = 3;
	public static final byte UNKNOWN_HYBRID = 4;

	final static int HB_EVAL_TYPE_COUNT = (HBDonChemType.hbdon_MAX.ordinal()-1)*(HBAccChemType.hbacc_MAX.ordinal()-1)*(HBSeqSep.seq_sep_MAX.ordinal()-1);
	final static double MAX_HB_ENERGY = 0;
	
	static HBEvalType[][][] hbe;

	public int[] numNeighborsHB;
	ArrayList<HbondPair> hbondTerms; //Hbond Pairs
	ArrayList<HbondPair>[] partHBlist;
	double hbondScale;
	
	enum HBondWeightType {
		hbw_NONE,
		hbw_SR_BB,
		hbw_LR_BB,
		hbw_SR_BB_SC,
		hbw_LR_BB_SC,
		hbw_SC,
		hbw_MAX
	}

	enum HBDonChemType {
		hbdon_NONE,
		hbdon_PBA, // hbdon_PROTEIN_BB_AMIDE
		hbdon_CXA, // hbdon_CARBOXAMIDE
		hbdon_IMD, // hbdon_IMIDAZOL_DELTA
		hbdon_IME, // hbdon_IMIDAZOL_EPSILON
		hbdon_IND, // hbdon_INDOL
		hbdon_AMO, // hbdon_AMINO
		hbdon_GDE, // hbdon_GUANIDINIUM_EPSILON
		hbdon_GDH, // hbdon_DIHYDRO_GUANIDINIUM
		hbdon_AHX, // hbdon_AROMATIC_HYDROXYL
		hbdon_HXL, // hbdon_HYDROXYL
		hbdon_H2O, // hbdon_WATER
		hbdon_GENERIC_BB,
		hbdon_GENERIC_SC,
		hbdon_MAX
	}

	enum HBAccChemType {
		hbacc_NONE,
		hbacc_PBA, // hbacc_PROTEIN_BB_AMIDE
		hbacc_CXA, // hbacc_CARBOXAMIDE
		hbacc_CXL, // hbacc_CARBOXYL
		hbacc_IMD, // hbacc_IMIDAZOL_DELTA
		hbacc_IME, // hbacc_IMIDAZOL_EPSILON
		hbacc_AHX, // hbacc_AROMATIC_HYDROXYL
		hbacc_HXL, // hbacc_HYDROXY
		hbacc_PCA_DNA, // hbacc_PHOSPHATE_CARBONYL_DNA
		hbacc_PES_DNA, // hbacc_PHOSPHATE_ESTER_DNA
		hbacc_RRI_DNA, // hbacc_RIBOSE_RING_DNA
		hbacc_PCA_RNA, // hbacc_PHOSPHATE_CARBONYL_RNA
		hbacc_PES_RNA, // hbacc_PHOSPHATE_ESTER_RNA
		hbacc_RRI_RNA, // hbacc_RIBOSE_RING_RNA
		hbacc_H2O, // hbacc_WATER
		hbacc_GENERIC_SP2BB,
		hbacc_GENERIC_SP2SC,
		hbacc_GENERIC_SP3BB,
		hbacc_GENERIC_SP3SC,
		hbacc_GENERIC_RINGBB,
		hbacc_GENERIC_RINGSC,
		hbacc_MAX
	}

	enum HBSeqSep{
		seq_sep_other, // // all other sequence separation not specified
		seq_sep_M4, // // acc_rsd.seqpos() - don_rsd.seqpos() = -4
		seq_sep_M3, // // acc_rsd.seqpos() - don_rsd.seqpos() = -3
		seq_sep_M2, // // acc_rsd.seqpos() - don_rsd.seqpos() = -2
		seq_sep_PM1, // // abs(acc_rsd.seqpos() - don_rsd.seqpos()) = 1
		seq_sep_P2, // // acc_rsd.seqpos() - don_rsd.seqpos() = 2
		seq_sep_P3, // // acc_rsd.seqpos() - don_rsd.seqpos() = -3
		seq_sep_P4, // // acc_rsd.seqpos() - don_rsd.seqpos() = 4
		seq_sep_MAX
	}

	enum HBGeoDimType {
		hbgd_NONE,

		// distance from the acceptor atom to the hydrogen
		hbgd_AHdist,

		// cosine of the base-acceptor-hydrogen angle
		// the base is acceptor hybridization dependent:
		//    sp2 hybrid -> base = res.atom_base(atm_num)
		//    sp3 hybrid -> base = res.abase2(atm_num)
		//    ring hybrid-> base = (res.atom_base(atm_num) + res.abase2(atm_num))/2
		// let BAunit = unit vector from base to the acceptor
		// let AHunit = unit vector from acceptor to hydrogen
		// cosBAH = BAunit <dot> AHunit
		hbgd_cosBAH,

		// cosine of the acceptor-hydrogen-donor angle
		// let AHunit = unit vector from acceptor to the hydrogen
		// let HDunit = unit vector from hydrogen to donor
		// cosAHD = AHDunit <dot> HDunit
		hbgd_cosAHD,


		// the angle formed by the acceptor-hydrogen-donor
		// this the interior angle measured in radians.

		// In score12, the hydrogen bond score function evaluated the cosine
		// of exterior BAH and AHD angles rather than the angles
		// themselves. This was done for two reasons:
		//
		//    1) The cosine of the exterior angle is easy to evaluate, just
		//    take the dot product of the noralized bond vectors.
		//
		//    2) When projecting uniform density density over cartesian
		//    space onto the theta angle in spherical coordinates, the
		//    resulting density is not uniform. This happens because the
		//    change in the volume of the conic section per unit angle
		//    depends on the angle itself. It turns out that the
		//    distribution of the cosine of the angle is uniform. Therefore
		//    when estimating distributions, as is done for knowledge based
		//    potentials, one should normalize the distribution by computing
		//    the distrbution in "cosine" space.
		//
		// Because density estimation should be done in cosine space and the
		// cosine of the angles is easy to evaluate, the polynomials in the
		// hydrogen bond score function were defined as functions of the
		// cosine of the angles.
		//
		// A limitation of this parametrization is that the dynamic range
		// from optimal AHD angle (180 degrees) to decent AHD angle (~160
		// degrees) is compressed. This can be seen by noticing that the
		// acos(x) around zero is steep, so a small change in x results in
		// large change in acos(x). To create polynomials that have such a
		// tight distribution requires them to be relatively high degree.
		//
		// As an alternative, the AHD angle can be used directly in the
		// parametrization.
		hbgd_AHD,

		// Torsional angle about the base-acceptor bond vector
		// Not yet implemented (11/09) but coming soon...
		hbgd_chi,
		hbgd_MAX

	}

	enum HBEvalType{
		hbe_NONE,
		hbe_dPBAaPBAsepM4helix,
		hbe_dPBAaPBAsepM3turn,
		hbe_dPBAaPBAsepM2turn,
		hbe_dPBAaPBAsepPM1,
		hbe_dPBAaPBAsepP2turn,
		hbe_dPBAaPBAsepP3turn,
		hbe_dPBAaPBAsepP4helix,
		hbe_dPBAaPBAsepother,
		hbe_dCXAaPBAsepPM1,
		hbe_dIMDaPBAsepPM1,
		hbe_dIMEaPBAsepPM1,
		hbe_dINDaPBAsepPM1,
		hbe_dAMOaPBAsepPM1,
		hbe_dGDEaPBAsepPM1,
		hbe_dGDHaPBAsepPM1,
		hbe_dAHXaPBAsepPM1,
		hbe_dHXLaPBAsepPM1,
		hbe_dCXAaPBAsepother,
		hbe_dIMDaPBAsepother,
		hbe_dIMEaPBAsepother,
		hbe_dINDaPBAsepother,
		hbe_dAMOaPBAsepother,
		hbe_dGDEaPBAsepother,
		hbe_dGDHaPBAsepother,
		hbe_dAHXaPBAsepother,
		hbe_dHXLaPBAsepother,
		hbe_dH2OaPBA,
		hbe_dPBAaCXAsepPM1,
		hbe_dPBAaCXAsepother,
		hbe_dCXAaCXA,
		hbe_dIMDaCXA,
		hbe_dIMEaCXA,
		hbe_dINDaCXA,
		hbe_dAMOaCXA,
		hbe_dGDEaCXA,
		hbe_dGDHaCXA,
		hbe_dAHXaCXA,
		hbe_dHXLaCXA,
		hbe_dH2OaCXA,
		hbe_dPBAaCXLsepPM1,
		hbe_dPBAaCXLsepother,
		hbe_dCXAaCXL,
		hbe_dIMDaCXL,
		hbe_dIMEaCXL,
		hbe_dINDaCXL,
		hbe_dAMOaCXL,
		hbe_dGDEaCXL,
		hbe_dGDHaCXL,
		hbe_dAHXaCXL,
		hbe_dHXLaCXL,
		hbe_dH2OaCXL,
		hbe_dPBAaIMDsepPM1,
		hbe_dPBAaIMDsepother,
		hbe_dCXAaIMD,
		hbe_dIMDaIMD,
		hbe_dIMEaIMD,
		hbe_dINDaIMD,
		hbe_dAMOaIMD,
		hbe_dGDEaIMD,
		hbe_dGDHaIMD,
		hbe_dAHXaIMD,
		hbe_dHXLaIMD,
		hbe_dH2OaIMD,
		hbe_dPBAaIMEsepPM1,
		hbe_dPBAaIMEsepother,
		hbe_dCXAaIME,
		hbe_dIMDaIME,
		hbe_dIMEaIME,
		hbe_dINDaIME,
		hbe_dAMOaIME,
		hbe_dGDEaIME,
		hbe_dGDHaIME,
		hbe_dAHXaIME,
		hbe_dHXLaIME,
		hbe_dH2OaIME,
		hbe_dPBAaAHXsepPM1,
		hbe_dPBAaAHXsepother,
		hbe_dCXAaAHX,
		hbe_dIMDaAHX,
		hbe_dIMEaAHX,
		hbe_dINDaAHX,
		hbe_dAMOaAHX,
		hbe_dGDEaAHX,
		hbe_dGDHaAHX,
		hbe_dAHXaAHX,
		hbe_dHXLaAHX,
		hbe_dH2OaAHX,
		hbe_dPBAaHXLsepPM1,
		hbe_dPBAaHXLsepother,
		hbe_dCXAaHXL,
		hbe_dIMDaHXL,
		hbe_dIMEaHXL,
		hbe_dINDaHXL,
		hbe_dAMOaHXL,
		hbe_dGDEaHXL,
		hbe_dGDHaHXL,
		hbe_dAHXaHXL,
		hbe_dHXLaHXL,
		hbe_dH2OaHXL,
		hbe_dPBAaPCA_DNAsepPM1,
		hbe_dPBAaPCA_DNAsepother,
		hbe_dCXAaPCA_DNA,
		hbe_dIMDaPCA_DNA,
		hbe_dIMEaPCA_DNA,
		hbe_dINDaPCA_DNA,
		hbe_dAMOaPCA_DNA,
		hbe_dGDEaPCA_DNA,
		hbe_dGDHaPCA_DNA,
		hbe_dAHXaPCA_DNA,
		hbe_dHXLaPCA_DNA,
		hbe_dH2OaPCA_DNA,
		hbe_dPBAaPCA_RNAsepPM1,
		hbe_dPBAaPCA_RNAsepother,
		hbe_dCXAaPCA_RNAsepPM1,
		hbe_dCXAaPCA_RNAsepother,
		hbe_dIMDaPCA_RNAsepPM1,
		hbe_dIMDaPCA_RNAsepother,
		hbe_dIMEaPCA_RNAsepPM1,
		hbe_dIMEaPCA_RNAsepother,
		hbe_dINDaPCA_RNAsepPM1,
		hbe_dINDaPCA_RNAsepother,
		hbe_dAMOaPCA_RNAsepPM1,
		hbe_dAMOaPCA_RNAsepother,
		hbe_dGDEaPCA_RNAsepPM1,
		hbe_dGDEaPCA_RNAsepother,
		hbe_dGDHaPCA_RNAsepPM1,
		hbe_dGDHaPCA_RNAsepother,
		hbe_dAHXaPCA_RNAsepPM1,
		hbe_dAHXaPCA_RNAsepother,
		hbe_dHXLaPCA_RNAsepPM1,
		hbe_dHXLaPCA_RNAsepother,
		hbe_dH2OaPCA_RNA,
		hbe_dPBAaPES_DNAsepPM1,
		hbe_dPBAaPES_DNAsepother,
		hbe_dCXAaPES_DNA,
		hbe_dIMDaPES_DNA,
		hbe_dIMEaPES_DNA,
		hbe_dINDaPES_DNA,
		hbe_dAMOaPES_DNA,
		hbe_dGDEaPES_DNA,
		hbe_dGDHaPES_DNA,
		hbe_dAHXaPES_DNA,
		hbe_dHXLaPES_DNA,
		hbe_dH2OaPES_DNA,
		hbe_dPBAaPES_RNAsepPM1,
		hbe_dPBAaPES_RNAsepother,
		hbe_dCXAaPES_RNAsepPM1,
		hbe_dCXAaPES_RNAsepother,
		hbe_dIMDaPES_RNAsepPM1,
		hbe_dIMDaPES_RNAsepother,
		hbe_dIMEaPES_RNAsepPM1,
		hbe_dIMEaPES_RNAsepother,
		hbe_dINDaPES_RNAsepPM1,
		hbe_dINDaPES_RNAsepother,
		hbe_dAMOaPES_RNAsepPM1,
		hbe_dAMOaPES_RNAsepother,
		hbe_dGDEaPES_RNAsepPM1,
		hbe_dGDEaPES_RNAsepother,
		hbe_dGDHaPES_RNAsepPM1,
		hbe_dGDHaPES_RNAsepother,
		hbe_dAHXaPES_RNAsepPM1,
		hbe_dAHXaPES_RNAsepother,
		hbe_dHXLaPES_RNAsepPM1,
		hbe_dHXLaPES_RNAsepother,
		hbe_dH2OaPES_RNA,
		hbe_dPBAaRRI_DNAsepPM1,
		hbe_dPBAaRRI_DNAsepother,
		hbe_dCXAaRRI_DNA,
		hbe_dIMDaRRI_DNA,
		hbe_dIMEaRRI_DNA,
		hbe_dINDaRRI_DNA,
		hbe_dAMOaRRI_DNA,
		hbe_dGDEaRRI_DNA,
		hbe_dGDHaRRI_DNA,
		hbe_dAHXaRRI_DNA,
		hbe_dHXLaRRI_DNA,
		hbe_dH2OaRRI_DNA,
		hbe_dPBAaRRI_RNAsepPM1,
		hbe_dPBAaRRI_RNAsepother,
		hbe_dCXAaRRI_RNAsepPM1,
		hbe_dCXAaRRI_RNAsepother,
		hbe_dIMDaRRI_RNAsepPM1,
		hbe_dIMDaRRI_RNAsepother,
		hbe_dIMEaRRI_RNAsepPM1,
		hbe_dIMEaRRI_RNAsepother,
		hbe_dINDaRRI_RNAsepPM1,
		hbe_dINDaRRI_RNAsepother,
		hbe_dAMOaRRI_RNAsepPM1,
		hbe_dAMOaRRI_RNAsepother,
		hbe_dGDEaRRI_RNAsepPM1,
		hbe_dGDEaRRI_RNAsepother,
		hbe_dGDHaRRI_RNAsepPM1,
		hbe_dGDHaRRI_RNAsepother,
		hbe_dAHXaRRI_RNAsepPM1,
		hbe_dAHXaRRI_RNAsepother,
		hbe_dHXLaRRI_RNAsepPM1,
		hbe_dHXLaRRI_RNAsepother,
		hbe_dH2OaRRI_RNA,
		hbe_dPBAaH2O,
		hbe_dCXAaH2O,
		hbe_dIMDaH2O,
		hbe_dIMEaH2O,
		hbe_dINDaH2O,
		hbe_dAMOaH2O,
		hbe_dGDEaH2O,
		hbe_dGDHaH2O,
		hbe_dAHXaH2O,
		hbe_dHXLaH2O,
		hbe_dH2OaH2O,
		hbe_GENERIC_SP2BB_SR,
		hbe_GENERIC_SP2BB_LR,
		hbe_GENERIC_SP3BB_SR,
		hbe_GENERIC_SP3BB_LR,
		hbe_GENERIC_RINGBB_SR,
		hbe_GENERIC_RINGBB_LR,
		hbe_GENERIC_SP2BSC_SR,
		hbe_GENERIC_SP2BSC_LR,
		hbe_GENERIC_SP3BSC_SR,
		hbe_GENERIC_SP3BSC_LR,
		hbe_GENERIC_RINGBSC_SR,
		hbe_GENERIC_RINGBSC_LR,
		hbe_GENERIC_SP2SCSC_SR,
		hbe_GENERIC_SP2SCSC_LR,
		hbe_GENERIC_SP3SCSC_SR,
		hbe_GENERIC_SP3SCSC_LR,
		hbe_GENERIC_RINGSCSC_SR,
		hbe_GENERIC_RINGSCSC_LR,
		hbe_MAX

	}


	public static void initHBElookup(){


		hbe = new HBEvalType[HBDonChemType.hbdon_MAX.ordinal()][HBAccChemType.hbacc_MAX.ordinal()][HBSeqSep.seq_sep_MAX.ordinal()];

		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_NONE.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_NONE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_NONE;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_M4.ordinal()] = HBEvalType.hbe_dPBAaPBAsepM4helix;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_M3.ordinal()] = HBEvalType.hbe_dPBAaPBAsepM3turn;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_M2.ordinal()] = HBEvalType.hbe_dPBAaPBAsepM2turn;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaPBAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_P2.ordinal()] = HBEvalType.hbe_dPBAaPBAsepP2turn;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_P3.ordinal()] = HBEvalType.hbe_dPBAaPBAsepP3turn;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_P4.ordinal()] = HBEvalType.hbe_dPBAaPBAsepP4helix;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaPBAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dCXAaPBAsepPM1;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMDaPBAsepPM1;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMEaPBAsepPM1;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dINDaPBAsepPM1;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAMOaPBAsepPM1;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDEaPBAsepPM1;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDHaPBAsepPM1;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAHXaPBAsepPM1;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dHXLaPBAsepPM1;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaPBAsepother;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaPBAsepother;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaPBAsepother;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaPBAsepother;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaPBAsepother;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaPBAsepother;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaPBAsepother;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaPBAsepother;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaPBAsepother;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaPBA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaCXAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaCXAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaCXA;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaCXA;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaCXA;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaCXA;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaCXA;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaCXA;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaCXA;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaCXA;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaCXA;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaCXA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaCXLsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaCXLsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaCXL;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaCXL;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaCXL;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaCXL;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaCXL;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaCXL;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaCXL;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaCXL;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaCXL;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaCXL;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaIMDsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaIMDsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaIMD;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaIMD;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaIMD;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaIMD;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaIMD;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaIMD;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaIMD;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaIMD;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaIMD;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaIMD;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaIMEsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaIMEsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaIME;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaIME;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaIME;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaIME;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaIME;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaIME;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaIME;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaIME;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaIME;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaIME;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaAHXsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaAHXsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaAHX;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaAHX;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaAHX;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaAHX;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaAHX;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaAHX;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaAHX;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaAHX;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaAHX;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaAHX;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaHXLsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaHXLsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaHXL;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaHXL;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaHXL;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaHXL;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaHXL;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaHXL;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaHXL;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaHXL;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaHXL;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaHXL;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaPCA_DNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaPCA_DNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaPCA_DNA;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaPCA_DNA;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaPCA_DNA;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaPCA_DNA;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaPCA_DNA;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaPCA_DNA;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaPCA_DNA;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaPCA_DNA;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaPCA_DNA;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaPCA_DNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dCXAaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMDaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMEaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dINDaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAMOaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDEaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDHaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAHXaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dHXLaPCA_RNAsepPM1;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaPCA_RNAsepother;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaPCA_RNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaPES_DNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaPES_DNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaPES_DNA;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaPES_DNA;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaPES_DNA;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaPES_DNA;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaPES_DNA;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaPES_DNA;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaPES_DNA;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaPES_DNA;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaPES_DNA;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaPES_DNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dCXAaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMDaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMEaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dINDaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAMOaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDEaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDHaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAHXaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dHXLaPES_RNAsepPM1;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaPES_RNAsepother;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaPES_RNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaRRI_DNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaRRI_DNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaRRI_DNA;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaRRI_DNA;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaRRI_DNA;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaRRI_DNA;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaRRI_DNA;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaRRI_DNA;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaRRI_DNA;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaRRI_DNA;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaRRI_DNA;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaRRI_DNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dPBAaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dCXAaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMDaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dIMEaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dINDaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAMOaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDEaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dGDHaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dAHXaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_dHXLaRRI_RNAsepPM1;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaRRI_RNAsepother;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaRRI_RNA;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dPBAaH2O;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dCXAaH2O;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMDaH2O;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dIMEaH2O;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dINDaH2O;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAMOaH2O;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDEaH2O;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dGDHaH2O;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dAHXaH2O;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dHXLaH2O;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_dH2OaH2O;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBB_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBB_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BB_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_BB.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PBA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PES_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PES_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP2BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_RRI_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_RRI_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP3BB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_RINGBB.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_CXA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_CXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PCA_DNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_PCA_RNA.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2BSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP2SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP2SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_AHX.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_HXL.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_H2O.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3BSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_SP3SC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_SP3SCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_IMD.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_IME.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_GENERIC_SC.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_SR;
		hbe[HBDonChemType.hbdon_PBA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGBSC_LR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_CXA.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_IMD.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_IME.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_IND.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_AMO.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_GDE.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_GDH.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_AHX.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_HXL.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_PM1.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_SR;
		hbe[HBDonChemType.hbdon_H2O.ordinal()][HBAccChemType.hbacc_GENERIC_RINGSC.ordinal()][HBSeqSep.seq_sep_other.ordinal()] = HBEvalType.hbe_GENERIC_RINGSCSC_LR;
	}

	static HashMap<String, Poly1D> polyFunctions;
	static HashMap<String, Fade>	fadeFunctions;

	static HBondWeightType[] weight_type_lookup;

	//Polynomials
	static Poly1D[] AHdist_poly_lookup; 
	static Poly1D[] cosBAH_short_poly_lookup;
	static Poly1D[] cosBAH_long_poly_lookup;
	static Poly1D[] cosAHD_short_poly_lookup;
	static Poly1D[] cosAHD_long_poly_lookup;

	//Fade intervals
	static Fade[] AHdist_short_fade_lookup;
	static Fade[] AHdist_long_fade_lookup;
	static Fade[] cosBAH_fade_lookup; 
	static Fade[] cosAHD_fade_lookup;


	public HBondEnergy(Molecule m, double scale) {
		init(m);
		hbondScale = scale;
	}
	
	public void init(Molecule m){
		//Init Neighbor list
		initHBNeighbors(m);

		if(hbe != null)
			return;

		initHBElookup();
		loadPoly();
		loadFade();
		loadHBEval();

	}


	static void loadPoly(){

		polyFunctions = new HashMap<String, Poly1D>();

		FileInputStream is;
		try {
			is = new FileInputStream( EnvironmentVars.getDataDir().concat("hbonds/HBPoly1D.csv") );

			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;
			String[] tokens;
			String polynomial_name;
			String geo_dim_name;
			HBGeoDimType geometric_dimension;
			double xmin,xmax,min_val,max_val;
			double root1=0.0,root2=0.0;
			int id,degree;
			double[] coefficients;

			curLine = bufread.readLine();
			while( curLine != null ){
				tokens = curLine.split(",");
				int i=0;
				id = Integer.parseInt(tokens[i++]);
				polynomial_name = tokens[i++];
				i++; // classic name field
				geo_dim_name = tokens[i++];
				geometric_dimension = HBGeoDimType.valueOf(geo_dim_name);
				xmin = Double.parseDouble(tokens[i++]);
				xmax = Double.parseDouble(tokens[i++]);
				min_val = Double.parseDouble(tokens[i++]);
				max_val = Double.parseDouble(tokens[i++]);
				i++;i++;
				//root1 = Double.parseDouble(tokens[i++]);
				//root2 = Double.parseDouble(tokens[i++]);
				degree = Integer.parseInt(tokens[i++]);
				coefficients = new double[degree];
				for(int j=0;j<tokens.length-i;j++){
					coefficients[j] = Double.parseDouble(tokens[j+i]);
				}

				Poly1D p = new Poly1D(polynomial_name, geometric_dimension, 
						xmin, xmax, min_val, max_val, 
						root1, root2, degree, coefficients);

				polyFunctions.put(polynomial_name, p);
				curLine = bufread.readLine();

			}
		} catch (Exception e) {
			System.out.println("Problem reading HBPoly1D.csv");
			e.printStackTrace();
		}

	}

	static void loadFade(){

		fadeFunctions = new HashMap<String, Fade>();

		FileInputStream is;
		try {
			is = new FileInputStream( EnvironmentVars.getDataDir().concat("hbonds/HBFadeIntervals.csv") );

			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;
			String[] tokens;
			String fade_interval_name;
			boolean smoothed;
			double min0,fmin,fmax,max0;
			int id,degree;
			double[] coefficients;
			String junction_type;
			curLine = bufread.readLine();
			while( curLine != null ){
				tokens = curLine.split(",");
				int i=0;
				id = Integer.parseInt(tokens[i++]);
				fade_interval_name = tokens[i++];
				junction_type = tokens[i++];
				if(junction_type.equals("smoothed"))
					smoothed = true;
				else if(junction_type.equals("piecewise_linear"))
					smoothed = false;
				else{
					System.out.println("Junction needs to be smooth or linear");
					smoothed = false;
				}

				min0 = Double.parseDouble(tokens[i++]);
				fmin = Double.parseDouble(tokens[i++]);
				fmax = Double.parseDouble(tokens[i++]);
				max0 = Double.parseDouble(tokens[i++]);


				Fade fade_interval = new Fade(fade_interval_name,min0,fmin,fmax,max0,smoothed);

				fadeFunctions.put(fade_interval_name, fade_interval);
				curLine = bufread.readLine();
			}
		} catch (Exception e) {
			System.out.println("Problem reading HBPoly1D.csv");
			e.printStackTrace();
		}

	}

	static void loadHBEval(){

		//Polynomials
		AHdist_poly_lookup = new Poly1D[HB_EVAL_TYPE_COUNT]; 
		cosBAH_short_poly_lookup = new Poly1D[HB_EVAL_TYPE_COUNT];
		cosBAH_long_poly_lookup = new Poly1D[HB_EVAL_TYPE_COUNT];
		cosAHD_short_poly_lookup = new Poly1D[HB_EVAL_TYPE_COUNT];
		cosAHD_long_poly_lookup = new Poly1D[HB_EVAL_TYPE_COUNT];

		//Fade intervals
		AHdist_short_fade_lookup = new Fade[HB_EVAL_TYPE_COUNT];
		AHdist_long_fade_lookup = new Fade[HB_EVAL_TYPE_COUNT];
		cosBAH_fade_lookup = new Fade[HB_EVAL_TYPE_COUNT]; 
		cosAHD_fade_lookup = new Fade[HB_EVAL_TYPE_COUNT];

		weight_type_lookup = new HBondWeightType[HB_EVAL_TYPE_COUNT];

		FileInputStream is;
		try {
			is = new FileInputStream( EnvironmentVars.getDataDir().concat("hbonds/HBEval.csv") );

			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
			String curLine = null;
			String[] tokens;
			String AHdist_poly_name, cosBAH_short_poly_name, cosBAH_long_poly_name,
			cosAHD_short_poly_name, cosAHD_long_poly_name, chi_poly_name,
			don_chem_type_name, acc_chem_type_name, seq_sep_type_name,
			AHdist_short_fade_name, AHdist_long_fade_name,
			cosBAH_fade_name, cosAHD_fade_name;
			HBDonChemType don_chem_type;
			HBAccChemType acc_chem_type;
			HBSeqSep seq_sep_type;
			Poly1D AHdist_poly, cosBAH_short_poly, cosBAH_long_poly,
			cosAHD_short_poly, cosAHD_long_poly, chi_poly;
			Fade AHdist_short_fade, AHdist_long_fade, cosBAH_fade, cosAHD_fade;
			String weight_type_name;
			HBondWeightType weight_type;

			curLine = bufread.readLine();
			while( curLine != null ){
				tokens = curLine.split(",");
				int i=0;
				don_chem_type_name = tokens[i++];
				don_chem_type = HBDonChemType.valueOf(don_chem_type_name);

				acc_chem_type_name = tokens[i++];
				acc_chem_type = HBAccChemType.valueOf(acc_chem_type_name);

				seq_sep_type_name = tokens[i++];
				seq_sep_type = HBSeqSep.valueOf(seq_sep_type_name);

				HBEvalType hbe_type = hbe[don_chem_type.ordinal()][acc_chem_type.ordinal()][seq_sep_type.ordinal()];

				AHdist_short_fade_name = tokens[i++];
				AHdist_short_fade = fadeFunctions.get(AHdist_short_fade_name);
				AHdist_short_fade_lookup[hbe_type.ordinal()] = AHdist_short_fade;

				AHdist_long_fade_name = tokens[i++];
				AHdist_long_fade = fadeFunctions.get(AHdist_long_fade_name);
				AHdist_long_fade_lookup[hbe_type.ordinal()] = AHdist_long_fade;

				cosBAH_fade_name = tokens[i++];
				cosBAH_fade = fadeFunctions.get(cosBAH_fade_name);
				cosBAH_fade_lookup[hbe_type.ordinal()] = cosBAH_fade;

				cosAHD_fade_name = tokens[i++];
				cosAHD_fade = fadeFunctions.get(cosAHD_fade_name);
				cosAHD_fade_lookup[hbe_type.ordinal()] = cosAHD_fade;

				i++; //fade for chi dimension

				AHdist_poly_name = tokens[i++];
				AHdist_poly = polyFunctions.get(AHdist_poly_name);
				AHdist_poly_lookup[hbe_type.ordinal()] = AHdist_poly;

				cosBAH_short_poly_name = tokens[i++];
				cosBAH_short_poly = polyFunctions.get(cosBAH_short_poly_name);
				cosBAH_short_poly_lookup[hbe_type.ordinal()] = cosBAH_short_poly;

				cosBAH_long_poly_name = tokens[i++];
				cosBAH_long_poly = polyFunctions.get(cosBAH_long_poly_name);
				cosBAH_long_poly_lookup[hbe_type.ordinal()] = cosBAH_long_poly;


				cosAHD_short_poly_name = tokens[i++];
				cosAHD_short_poly = polyFunctions.get(cosAHD_short_poly_name);
				cosAHD_short_poly_lookup[hbe_type.ordinal()] = cosAHD_short_poly;

				cosAHD_long_poly_name = tokens[i++];
				cosAHD_long_poly = polyFunctions.get(cosAHD_long_poly_name);
				cosAHD_long_poly_lookup[hbe_type.ordinal()] = cosAHD_long_poly;

				weight_type_name = tokens[i++];
				weight_type = HBondWeightType.valueOf(weight_type_name);
				weight_type_lookup[hbe_type.ordinal()] = weight_type;

				curLine = bufread.readLine();
			}
		} catch (Exception e) {
			System.out.println("Problem reading HBPoly1D.csv");
			e.printStackTrace();
		}

	}


	public static HBEvalType getEvalType(Atom donor, Residue donorRes, Atom accept, Residue acceptRes) {
		HBDonChemType don_type = donor.hbdon;
		HBAccChemType acc_type = accept.hbacc;
		HBSeqSep      seq_sep = getSeqSep(don_type,acc_type,donorRes,acceptRes);
		HBEvalType    eval_type = hbe[don_type.ordinal()][acc_type.ordinal()][seq_sep.ordinal()];
		return eval_type;

	}


	private static HBSeqSep getSeqSep(HBDonChemType don_type,
			HBAccChemType acc_type, Residue donorRes, Residue acceptRes) {

		int sep = Integer.MAX_VALUE;
		if(donorRes.strandNumber == acceptRes.strandNumber)
			sep = acceptRes.moleculeResidueNumber - donorRes.moleculeResidueNumber;

		switch(don_type){
		case hbdon_NONE:
		case hbdon_H2O:
			return HBSeqSep.seq_sep_other;
		case hbdon_PBA:
			switch(acc_type){
			case hbacc_NONE:
			case hbacc_H2O:
				return HBSeqSep.seq_sep_other;
			case hbacc_PBA:
				switch(sep){
				case -4: return HBSeqSep.seq_sep_M4; 
				case -3: return HBSeqSep.seq_sep_M3; 
				case -2: return HBSeqSep.seq_sep_M2; 
				case -1:
				case 1: return HBSeqSep.seq_sep_PM1; 
				case 2: return HBSeqSep.seq_sep_P2; 
				case 3: return HBSeqSep.seq_sep_P3; 
				case 4: return HBSeqSep.seq_sep_P4; 
				default: return HBSeqSep.seq_sep_other; 
				}
			default:
				if (sep == 1 || sep == -1) { return HBSeqSep.seq_sep_PM1;}
				else { return HBSeqSep.seq_sep_other; }

			} 
		case hbdon_CXA:
		case hbdon_IMD:
		case hbdon_IME:
		case hbdon_IND:
		case hbdon_AMO:
		case hbdon_GDE:
		case hbdon_GDH:
		case hbdon_AHX:
		case hbdon_HXL:
			switch(acc_type){
			case hbacc_NONE:
			case hbacc_CXA:
			case hbacc_CXL:
			case hbacc_IMD:
			case hbacc_IME:
			case hbacc_AHX:
			case hbacc_HXL:
			case hbacc_PCA_DNA:
			case hbacc_PES_DNA:
			case hbacc_RRI_DNA:
			case hbacc_H2O:
				return HBSeqSep.seq_sep_other; 
			default:
				if (sep == 1 || sep == -1) { return HBSeqSep.seq_sep_PM1;}
				else { return HBSeqSep.seq_sep_other; } 
			}
		default:
			if (sep == 1 || sep == -1) { return HBSeqSep.seq_sep_PM1;}
			else { return HBSeqSep.seq_sep_other; } 
		}

	}


	///////////////////////////////////////////////////////////////////////////////
	static HBDonChemType getDonType(Residue res, Atom a, boolean isProtein)
	{

		if (a.isBBatom){
			if (isProtein) {
				return HBDonChemType.hbdon_PBA;
			} else {
				return HBDonChemType.hbdon_GENERIC_BB;
			}
		} else {
			switch(res.name){
			case "ASN": case "GLN": return HBDonChemType.hbdon_CXA; 
			case "HIS": case "HID": case "HIE": case "HIP":
				if (a.name.equals("ND1")){
					return HBDonChemType.hbdon_IMD;
				} else {
					//assert( a.name.equals("NE2"));
					return HBDonChemType.hbdon_IME;
				} 
			case "TRP":
				return HBDonChemType.hbdon_IND; 
			case "LYS":
				return HBDonChemType.hbdon_AMO; 
			case "ARG": 
				if (a.name.equals("NE")){
					return HBDonChemType.hbdon_GDE;
				} else {
					assert(a.name.equals("NH1") || a.name.equals("NH2"));
					return HBDonChemType.hbdon_GDH;
				} 
			case "TYR":
				return HBDonChemType.hbdon_AHX; 
			case "SER":
			case "THR":
				return HBDonChemType.hbdon_HXL; 
			case "ALA": case "CYS": case "ASP": case "GLU": case "PHE":
			case "GLY": case "ILE": case "LEU": case "MET": case "PRO": case "VAL":
				return HBDonChemType.hbdon_NONE; 
			case "HOH":
				return HBDonChemType.hbdon_H2O; 
			default:
				//tr << "WARNING: Unknown Hydrogen Bond donor type for: " + don_rsd.name1() + I(3, don_rsd.seqpos()) + " " + don_rsd.atom_name( datm) + ".  Using hbdon_GENERIC_SC.";
				return HBDonChemType.hbdon_GENERIC_SC; 
			}
		}
		//return HBDonChemType.hbdon_NONE;
	}

	///////////////////////////////////////////////////////////////////////////////
	static HBAccChemType getAccType(Residue res, Atom a, boolean isProtein)
	{
		if( a.isBBatom){
			if( isProtein ) {
				return HBAccChemType.hbacc_PBA;
			} else {
				// generic types; for backwards compatibility; prefer functional group based chem type
				switch (a.hybridization){
				case SP2_HYBRID:
					//				tr << "WARNING: Unknown hydrogen bond acceptor type for: " + acc_rsd.name1() + I(3, acc_rsd.seqpos()) + " " + acc_rsd.atom_name(aatm) + ".  Using hbdon_GENERIC_SP2BB.";
					return HBAccChemType.hbacc_GENERIC_SP2BB; 
				case SP3_HYBRID:
					//tr << "WARNING: Unknown hydrogen bond acceptor type for: " + acc_rsd.name1() + I(3, acc_rsd.seqpos()) + " " + acc_rsd.atom_name(aatm) + ".  Using hbdon_GENERIC_SP3BB.";
					return HBAccChemType.hbacc_GENERIC_SP3BB; 
				case RING_HYBRID:
					//tr << "WARNING: Unknown hydrogen bond acceptor type for: " + acc_rsd.name1() + I(3, acc_rsd.seqpos()) + " " + acc_rsd.atom_name(aatm) + ".  Using hbdon_GENERIC_RINGBB.";
					return HBAccChemType.hbacc_GENERIC_RINGBB;
				case UNKNOWN_HYBRID:
					//tr << "WARNING: Unknown hydrogen bond acceptor type for: " + acc_rsd.name1() + I(3, acc_rsd.seqpos()) + " " + acc_rsd.atom_name(aatm) + ".  Using hbdon_GENERIC_RINGBB.";
					return HBAccChemType.hbacc_NONE; 
				}
			}
		} else {
			switch(res.name){
			case "ASN": case "GLN": return HBAccChemType.hbacc_CXA; 
			case "ASP": case "GLU": return HBAccChemType.hbacc_CXL; 
			case "HIS": case "HID": case "HIE": case "HIP":
				if (a.name.equals("ND1")){
					return HBAccChemType.hbacc_IMD;
				} else {
					//assert( a.name.equals("NE2"));
					return HBAccChemType.hbacc_IME;
				} 
			case "ALA": case "CYS": case "PHE": case "GLY": case "ILE": case "LEU":
			case "MET": case "PRO": case "VAL": case "TYR":
				return HBAccChemType.hbacc_AHX; 
			case "SER": case "THR": return HBAccChemType.hbacc_HXL; 
			case "LYS": case "ARG": case "TRP": 
				return HBAccChemType.hbacc_NONE;
			case "HOH":
				return HBAccChemType.hbacc_H2O; 
			default:
				// generic types; for backwards compatibility; prefer functional group based chem type
				switch(a.hybridization){
				case SP2_HYBRID:
					return HBAccChemType.hbacc_GENERIC_SP2SC; 
				case SP3_HYBRID:
					return HBAccChemType.hbacc_GENERIC_SP3SC; 
				case RING_HYBRID:
					return HBAccChemType.hbacc_GENERIC_RINGSC;
				case UNKNOWN_HYBRID:
					return HBAccChemType.hbacc_NONE; 
				}
			}

		}
		return HBAccChemType.hbacc_NONE;
	}
	
	public void initializeHbondTerms() {
		hbondTerms = new ArrayList<HbondPair>();
	}
	
	void addIfValid(int atom1, int atom2, int res1, int res2, Molecule m, double multiplier){
		Atom hydro=null;
		Atom accept;
		Atom donor;
		Residue donorRes;
		Residue acceptRes;
		//KER: I'm going to require that the atoms be from different residues. Not sure if that's correct
		//KER: But it seems unlikely that a residue can h-bond to itself
		if(m.atom[atom1].moleculeResidueNumber == m.atom[atom2].moleculeResidueNumber)
			return;
		
		if(m.atom[atom1].hbacc != null && m.atom[atom2].hbdon != null){
			accept = m.atom[atom1];
			donor = m.atom[atom2];
			acceptRes = m.residue[res1];
			donorRes = m.residue[res2];
		}else if(m.atom[atom2].hbacc != null && m.atom[atom1].hbdon != null){
			accept = m.atom[atom2];
			donor = m.atom[atom1];
			acceptRes = m.residue[res2];
			donorRes = m.residue[res1];
		}else{
			return;
		}
		
		
		//Find H
		for(int q=0;q<donor.bond.length;q++)
			if(m.atom[donor.bond[q]].elementType.equals("H")){
				hydro = m.atom[donor.bond[q]];
				//identify hbe type
				HBondEnergy.HBEvalType hbe = HBondEnergy.getEvalType(donor,donorRes,accept,acceptRes);
				hbondTerms.add(new HbondPair(donor,hydro,accept,multiplier,hbe));
			}
	}

	void calculateHBondEnergy(double[] coordinates, int curIndex,
			double[] energyTerms, Molecule m) {



		int atomix3, atomjx3, atomi, atomj, atomk, atomkx3, atomb=0,atombx3 = 0,atomb2=0,atomb2x3=0,atoml, atomlx3, atomm, atommx3;
		int ix4;
		double rij, rij2, rij6, rij10, rij12;
		double cx, cy, cz, cmag;
		double dx, dy, dz, ex, ey, ez, dmag, emag;
		double theta=0.0, thetadeg;
		double fx, fy, fz, fmag, gx,gy,gz,gmag;
		double ux, uy, uz, vx, vy, vz, umag, vmag;
		double dihed;
		double negcostheta;
		double cos2theta,cos4theta, negcosPhi, cos2phi, Ehb;
		byte donorType, acceptType;
		byte donorSS, acceptSS;

		double phi=0;
		double chi=0;



		Iterator<HbondPair> iter = hbondTerms.iterator();

		while(iter.hasNext()){
			HbondPair hbp = iter.next();

			HBEvalType hbe = hbp.hbe;
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

			double Edelta = 0.0;

			double sDeltaShort = 0.0; //Smoothing term
			double sDeltaLong = 0.0; //Smoothing term

			Edelta = AHdist_poly_lookup[hbe.ordinal()].getVal(rij); 
			sDeltaShort = AHdist_short_fade_lookup[hbe.ordinal()].getSmoothVal(rij); 
			sDeltaLong = AHdist_long_fade_lookup[hbe.ordinal()].getSmoothVal(rij);

			//D-H
			cx = coordinates[atomix3] - coordinates[atomjx3];
			cy = coordinates[atomix3 + 1] - coordinates[atomjx3 + 1];
			cz = coordinates[atomix3 + 2] - coordinates[atomjx3 + 2];
			cmag = Math.sqrt(cx*cx + cy*cy + cz*cz);

			//angle = (float)Math.toDegrees(Math.acos((dx*ex + dy*ey + dz*ez) / (dmag * -emag)));
			//costheta = (cx*dx + cy*dy + cz*dz) / (cmag * -dmag);
			negcostheta = (cx*dx + cy*dy + cz*dz) / (cmag * dmag);
			//theta = Math.toDegrees(Math.acos(costheta));

			boolean use_cosAHD = (cosBAH_short_poly_lookup[hbe.ordinal()].geoType == HBGeoDimType.hbgd_AHD);

			if(! use_cosAHD)
				theta = Math.PI - Math.acos(negcostheta);

			double EthetaShort = 0.0;
			double EthetaLong = 0.0;
			double sTheta = 1.0; //Smoothing term

			if(use_cosAHD){
				EthetaShort = cosAHD_short_poly_lookup[hbe.ordinal()].getVal(negcostheta);
				EthetaLong = cosAHD_long_poly_lookup[hbe.ordinal()].getVal(negcostheta);
			}else{
				EthetaShort = cosAHD_short_poly_lookup[hbe.ordinal()].getVal(theta);
				EthetaLong = cosAHD_long_poly_lookup[hbe.ordinal()].getVal(theta);
			}
			sTheta = cosAHD_fade_lookup[hbe.ordinal()].getSmoothVal(negcostheta);



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
			if(hbp.accept.hybridization == RING_HYBRID){
				ex = coordinates[atomkx3] - 0.5*(coordinates[atombx3]+coordinates[atomb2x3]);
				ey = coordinates[atomkx3 + 1] - 0.5*(coordinates[atombx3 + 1]+coordinates[atomb2x3 + 1]);
				ez = coordinates[atomkx3 + 2] - 0.5*(coordinates[atombx3 + 2]+coordinates[atomb2x3 + 2]);
			} else{
				ex = coordinates[atomkx3] - coordinates[atombx3];
				ey = coordinates[atomkx3 + 1] - coordinates[atombx3 + 1];
				ez = coordinates[atomkx3 + 2] - coordinates[atombx3 + 2];
			}
			emag = Math.sqrt(ex*ex + ey*ey + ez*ez);

			negcosPhi = (dx*ex + dy*ey + dz*ez) / (dmag * emag);
			//phi = Math.toDegrees(Math.acos( cosPhi ));

			double EphiShort = 0.0;
			double EphiLong = 0.0;
			double sPhi = 1.0; //Smoothing term

			EphiShort = cosBAH_short_poly_lookup[hbe.ordinal()].getVal(negcosPhi);
			EphiLong = cosBAH_long_poly_lookup[hbe.ordinal()].getVal(negcosPhi);
			sPhi = cosBAH_fade_lookup[hbe.ordinal()].getSmoothVal(negcosPhi);






			//energy = Pr*FxD*FxH + FSr*(PSxD*FxH + FxD*PSxH) + FLr*(PLxD*FxH + FxD*PLxH);


			//energy = Pr*FxD*FxH + FSr*(PSxD*FxH + FxD*PSxH) + FLr*(PLxD*FxH + FxD*PLxH);

			Ehb = hbp.multiplier * (Edelta*sTheta*sPhi + sDeltaShort*(EthetaShort*sPhi + sTheta*EphiShort) + sDeltaLong*(EthetaLong*sPhi + sTheta*EphiLong));


 			double Echi = 0.0;

			//Determine chi
			if(hbp.accept.hybridization == SP2_HYBRID ||
					(hbp.accept.hbacc == HBAccChemType.hbacc_AHX || 
					hbp.accept.hbacc == HBAccChemType.hbacc_HXL )){
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

//				atomm = R1.moleculeAtomNumber;
//				atommx3 = atomm*3;

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
				chi = Math.acos((ux*vx + uy*vy + uz*vz) / (umag * vmag));

				// BUT, that doesn't solve the handedness (sign) problem for the dihedral!
				// To do that, we look at the angle between 'f' and 'v'
				// Dot product again
				if( Math.acos((fx*vx + fy*vy + fz*vz) / (fmag * vmag)) > 1.57079633 )
				{ chi = -chi; }

				double chiPenalty = 0.0;
				if(hbp.accept.hybridization == SP2_HYBRID)
					chiPenalty = bah_chi_compute_energy_sp2(negcosPhi,chi);
				else if((hbp.accept.hbacc == HBAccChemType.hbacc_AHX || 
					hbp.accept.hbacc == HBAccChemType.hbacc_HXL ))
					chiPenalty = bah_chi_compute_energy_sp3(chi);
				
				Ehb += chiPenalty;


			}

			Ehb = fade_energy(Ehb);

			if(Ehb >= MAX_HB_ENERGY)
				continue;


			double environmentScale = getNeighborScale(hbe,m,hbp.donor,hbp.accept);

			if(Amber96ext.debug)
				System.out.println("Acc: " + m.residue[hbp.accept.moleculeResidueNumber].fullName+" Don: "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+Ehb+" "+environmentScale+" "+
					numNeighborsHB[hbp.accept.moleculeResidueNumber]+" "+numNeighborsHB[hbp.donor.moleculeResidueNumber]);
			
			Ehb *= environmentScale;
			
			
			//Ehb = hbp.multiplier * sDelta * sTheta * sPhi * (Edelta+Etheta+Ephi+Echi);

			//Ehb = Echi;

			//			if(hbp.donor.strandNumber != hbp.accept.strandNumber){
			//			if(debug){
			//				System.out.print("Ehb: "+Ehb+" "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+hbp.donor.name+" "+m.residue[hbp.accept.moleculeResidueNumber].fullName+" "+hbp.accept.name+" ");
			//				System.out.println(Edelta+" ( "+rij+" ) "+Etheta+" ( "+costheta+" ) "+Ephi+" ( "+cosPhi+" ) "+Echi+" ( "+chi+" ) ");
			//			}
			//			}
			energyTerms[4] += Ehb; 
		}

		energyTerms[4] *= hbondScale;


	}

	
	//Copy of previous function, except we add to the AmberPairs terms
	void calculateHBondEnergyUpd(ArrayList<HbondPair> hbonds,float[] coordinates, int curIndex,
			double[] energyTerms, Molecule m,double hbondScale, Amber96ext amb96ff ) {



		int atomix3, atomjx3, atomi, atomj, atomk, atomkx3, atomb=0,atombx3 = 0,atomb2=0,atomb2x3=0,atoml, atomlx3, atomm, atommx3;
		int ix4;
		double rij, rij2, rij6, rij10, rij12;
		double cx, cy, cz, cmag;
		double dx, dy, dz, ex, ey, ez, dmag, emag;
		double theta=0.0, thetadeg;
		double fx, fy, fz, fmag, gx,gy,gz,gmag;
		double ux, uy, uz, vx, vy, vz, umag, vmag;
		double dihed;
		double negcostheta;
		double cos2theta,cos4theta, negcosPhi, cos2phi, Ehb;
		byte donorType, acceptType;
		byte donorSS, acceptSS;

		double phi=0;
		double chi=0;



		Iterator<HbondPair> iter = hbonds.iterator();

		while(iter.hasNext()){
			HbondPair hbp = iter.next();

			/*if(hbp.donor.isBBatom || hbp.accept.isBBatom)
				continue;

			if(hbp.donor.strandNumber == hbp.accept.strandNumber)
				continue;*/

			//determine donor acceptor type
			//donorType = getHybridization(hbp.donor.forceFieldType);
			//acceptType = getHybridization(hbp.accept.forceFieldType);

			//donorSS = m.residue[hbp.donor.moleculeResidueNumber].SStype;
			//acceptSS = m.residue[hbp.accept.moleculeResidueNumber].SStype;

			HBEvalType hbe = hbp.hbe;
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
//			boolean badPair = false;
//			if(hbp.accept.elementType.equals("O")){ 
//				for(int i1 = 0; i1<hbp.accept.bond.length;i1++){
//					if((m.atom[hbp.accept.bond[i1]].elementType).equals("H")){
//						//Bad if donor is closer to H than to acceptor
//						if(hbp.donor.distance(m.atom[hbp.accept.bond[i1]]) < hbp.donor.distance(hbp.accept)){
//							badPair = true;
//						}
//					}
//				}
//			}
//			if(badPair)
//				continue;

			double Edelta = 0.0;

			double sDeltaShort = 0.0; //Smoothing term
			double sDeltaLong = 0.0; //Smoothing term

			Edelta = AHdist_poly_lookup[hbe.ordinal()].getVal(rij); 
			sDeltaShort = AHdist_short_fade_lookup[hbe.ordinal()].getSmoothVal(rij); 
			sDeltaLong = AHdist_long_fade_lookup[hbe.ordinal()].getSmoothVal(rij);

			//D-H
			cx = coordinates[atomix3] - coordinates[atomjx3];
			cy = coordinates[atomix3 + 1] - coordinates[atomjx3 + 1];
			cz = coordinates[atomix3 + 2] - coordinates[atomjx3 + 2];
			cmag = Math.sqrt(cx*cx + cy*cy + cz*cz);

			//angle = (float)Math.toDegrees(Math.acos((dx*ex + dy*ey + dz*ez) / (dmag * -emag)));
			//costheta = (cx*dx + cy*dy + cz*dz) / (cmag * -dmag);
			negcostheta = (cx*dx + cy*dy + cz*dz) / (cmag * dmag);
			//theta = Math.toDegrees(Math.acos(costheta));

			boolean use_cosAHD = (cosBAH_short_poly_lookup[hbe.ordinal()].geoType == HBGeoDimType.hbgd_AHD);

			if(! use_cosAHD)
				theta = Math.PI - Math.acos(negcostheta);

			double EthetaShort = 0.0;
			double EthetaLong = 0.0;
			double sTheta = 1.0; //Smoothing term

			if(use_cosAHD){
				EthetaShort = cosAHD_short_poly_lookup[hbe.ordinal()].getVal(negcostheta);
				EthetaLong = cosAHD_long_poly_lookup[hbe.ordinal()].getVal(negcostheta);
			}else{
				EthetaShort = cosAHD_short_poly_lookup[hbe.ordinal()].getVal(theta);
				EthetaLong = cosAHD_long_poly_lookup[hbe.ordinal()].getVal(theta);
			}
			sTheta = cosAHD_fade_lookup[hbe.ordinal()].getSmoothVal(negcostheta);



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
			if(hbp.accept.hybridization == RING_HYBRID){
				ex = coordinates[atomkx3] - 0.5*(coordinates[atombx3]+coordinates[atomb2x3]);
				ey = coordinates[atomkx3 + 1] - 0.5*(coordinates[atombx3 + 1]+coordinates[atomb2x3 + 1]);
				ez = coordinates[atomkx3 + 2] - 0.5*(coordinates[atombx3 + 2]+coordinates[atomb2x3 + 2]);
			} else{
				ex = coordinates[atomkx3] - coordinates[atombx3];
				ey = coordinates[atomkx3 + 1] - coordinates[atombx3 + 1];
				ez = coordinates[atomkx3 + 2] - coordinates[atombx3 + 2];
			}
			emag = Math.sqrt(ex*ex + ey*ey + ez*ez);

			negcosPhi = (dx*ex + dy*ey + dz*ez) / (dmag * emag);
			//phi = Math.toDegrees(Math.acos( cosPhi ));

			double EphiShort = 0.0;
			double EphiLong = 0.0;
			double sPhi = 1.0; //Smoothing term

			EphiShort = cosBAH_short_poly_lookup[hbe.ordinal()].getVal(negcosPhi);
			EphiLong = cosBAH_long_poly_lookup[hbe.ordinal()].getVal(negcosPhi);
			sPhi = cosBAH_fade_lookup[hbe.ordinal()].getSmoothVal(negcosPhi);






			//energy = Pr*FxD*FxH + FSr*(PSxD*FxH + FxD*PSxH) + FLr*(PLxD*FxH + FxD*PLxH);


			//energy = Pr*FxD*FxH + FSr*(PSxD*FxH + FxD*PSxH) + FLr*(PLxD*FxH + FxD*PLxH);

			Ehb = hbp.multiplier * (Edelta*sTheta*sPhi + sDeltaShort*(EthetaShort*sPhi + sTheta*EphiShort) + sDeltaLong*(EthetaLong*sPhi + sTheta*EphiLong));


 			double Echi = 0.0;

			//Determine chi
			if(hbp.accept.hybridization == SP2_HYBRID ||
					(hbp.accept.hbacc == HBAccChemType.hbacc_AHX || 
					hbp.accept.hbacc == HBAccChemType.hbacc_HXL )){
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

//				atomm = R1.moleculeAtomNumber;
//				atommx3 = atomm*3;

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
				chi = Math.acos((ux*vx + uy*vy + uz*vz) / (umag * vmag));

				// BUT, that doesn't solve the handedness (sign) problem for the dihedral!
				// To do that, we look at the angle between 'f' and 'v'
				// Dot product again
				if( Math.acos((fx*vx + fy*vy + fz*vz) / (fmag * vmag)) > 1.57079633 )
				{ chi = -chi; }

				double chiPenalty = 0.0;
				if(hbp.accept.hybridization == SP2_HYBRID)
					chiPenalty = bah_chi_compute_energy_sp2(negcosPhi,chi);
				else if((hbp.accept.hbacc == HBAccChemType.hbacc_AHX || 
					hbp.accept.hbacc == HBAccChemType.hbacc_HXL ))
					chiPenalty = bah_chi_compute_energy_sp3(chi);
				
				Ehb += chiPenalty;


			}

			Ehb = fade_energy(Ehb);

			if(Ehb >= MAX_HB_ENERGY)
				continue;


			double environmentScale = getNeighborScale(hbe,m,hbp.donor,hbp.accept);

			if(Amber96ext.debug)
				System.out.println("Acc: " + m.residue[hbp.accept.moleculeResidueNumber].fullName+" Don: "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+Ehb+" "+environmentScale+" "+
					numNeighborsHB[hbp.accept.moleculeResidueNumber]+" "+numNeighborsHB[hbp.donor.moleculeResidueNumber]);
			
			Ehb *= environmentScale;
			
			
			//Ehb = hbp.multiplier * sDelta * sTheta * sPhi * (Edelta+Etheta+Ephi+Echi);

			//Ehb = Echi;

			//			if(hbp.donor.strandNumber != hbp.accept.strandNumber){
			//			if(debug){
			//				System.out.print("Ehb: "+Ehb+" "+m.residue[hbp.donor.moleculeResidueNumber].fullName+" "+hbp.donor.name+" "+m.residue[hbp.accept.moleculeResidueNumber].fullName+" "+hbp.accept.name+" ");
			//				System.out.println(Edelta+" ( "+rij+" ) "+Etheta+" ( "+costheta+" ) "+Ephi+" ( "+cosPhi+" ) "+Echi+" ( "+chi+" ) ");
			//			}
			//			}
			amb96ff.updatePairTerms(hbondScale*Ehb,hbp.donor.moleculeResidueNumber,hbp.accept.moleculeResidueNumber);
			energyTerms[4] += Ehb; 
		}

		energyTerms[4] *= hbondScale;


	}

	void setupPartialHBondArrays(int numRows, int maxNumColumns,
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
	
	///@brief Fade the energy smoothly to zero over the energy range [-0.1, 0.1]
///@detail Because of the additive functional form, in order to make
///derivative continuous at the boundary of definition, we fade the
///energy function smoothly to zero.
///
/// Check that f(x) = -0.025 + 0.5x - 2.5x^2 satisfies
///     f(-.1) = -0.025 +   0.5*(-.1) - 2.5*(-.1)^2 = -.1
///     f( .1) = -0.025 +    0.5*(.1) -  2.5*(.1)^2 = 0
///    f'(-.1) =  0.5   - 2.5*2*(-.1)               = 1
///     f'(.1) =  0.5   -  2.5*2*(.1)               = 0
public static double fade_energy(double energy) {
	
	if(energy > 0.1){
		return 0.0;
	} else if (energy > -0.1){
		return (-0.025 + 0.5*energy - 2.5*energy*energy);
	}
	return energy;
}


	private static double bah_chi_compute_energy_sp3(double chi) {
		// just add in a penalty directly to the energy sum; the chi-penalty
		// is only multiplied in for the sp2 term.
		double max_penalty = 0.25;
		double cos2ChiShifted = max_penalty * ( 1 + Math.cos(chi)) / 2;
		return cos2ChiShifted;
	}

	static final double d = 0.75;
	static final double m = 1.6;
	static final double l = .357;
	static double bah_chi_compute_energy_sp2(double xH, double chi) {


		double PI_minus_BAH = Math.acos(xH);
		double BAH = Math.PI - ( PI_minus_BAH );

		double H = ((Math.cos(2*chi) + 1) * 0.5);
		double F=0,G=0;

		if ( BAH >= Math.PI * 2/3 ) {
			F = d/2 * Math.cos(3 * PI_minus_BAH) + d/2 - 0.5;
			G = d - 0.5;
		} else if ( BAH >= Math.PI * (2/3 - l)) {
			double outer_rise = (Math.cos(Math.PI - (Math.PI*2/3 -  BAH)/l));
			F = m/2 * outer_rise + m/2 - 0.5;
			G = (m - d)/2 * outer_rise + (m - d)/2 + d - 0.5;
		} else {
			F = m-0.5;
			G = m-0.5;
		}

		return (H*F + (1-H)*G);
	}

	private double getNeighborScale(HBEvalType hbe, Molecule m, Atom donor, Atom accept) {



		double weight =  1.0;
		if ( hbe_is_SC_type(hbe) ) {
			int donorNB = numNeighborsHB[donor.moleculeResidueNumber];
			int acceptNB = numNeighborsHB[accept.moleculeResidueNumber];
			weight = hb_env_dep_burial_lin( acceptNB, donorNB );
		}
		return weight;

	}

	private static double hb_env_dep_burial_lin(int acceptNB, int donorNB) {
		return (burial_weight(acceptNB) + burial_weight(donorNB));
	}

	static double burial_weight(int nb)
	{
		if ( nb < 7 ) return 0.1;
		if ( nb > 24 ) return 0.5;
		return (nb-2.75)*(0.5/21.25);
	}

	static boolean	hbe_is_BB_type( HBEvalType hbe )
	{
		switch(hbe){
		case hbe_NONE: return false; 
		case hbe_dPBAaPBAsepM4helix: return true; 
		case hbe_dPBAaPBAsepM3turn: return true; 
		case hbe_dPBAaPBAsepM2turn: return true; 
		case hbe_dPBAaPBAsepPM1: return true; 
		case hbe_dPBAaPBAsepP2turn: return true; 
		case hbe_dPBAaPBAsepP3turn: return true; 
		case hbe_dPBAaPBAsepP4helix: return true; 
		case hbe_dPBAaPBAsepother: return true; 
		case hbe_dCXAaPBAsepPM1: return false; 
		case hbe_dIMDaPBAsepPM1: return false; 
		case hbe_dIMEaPBAsepPM1: return false; 
		case hbe_dINDaPBAsepPM1: return false; 
		case hbe_dAMOaPBAsepPM1: return false; 
		case hbe_dGDEaPBAsepPM1: return false; 
		case hbe_dGDHaPBAsepPM1: return false; 
		case hbe_dAHXaPBAsepPM1: return false; 
		case hbe_dHXLaPBAsepPM1: return false; 
		case hbe_dCXAaPBAsepother: return false; 
		case hbe_dIMDaPBAsepother: return false; 
		case hbe_dIMEaPBAsepother: return false; 
		case hbe_dINDaPBAsepother: return false; 
		case hbe_dAMOaPBAsepother: return false; 
		case hbe_dGDEaPBAsepother: return false; 
		case hbe_dGDHaPBAsepother: return false; 
		case hbe_dAHXaPBAsepother: return false; 
		case hbe_dHXLaPBAsepother: return false; 
		case hbe_dH2OaPBA: return false; 
		case hbe_dPBAaCXAsepPM1: return false; 
		case hbe_dPBAaCXAsepother: return false; 
		case hbe_dCXAaCXA: return false; 
		case hbe_dIMDaCXA: return false; 
		case hbe_dIMEaCXA: return false; 
		case hbe_dINDaCXA: return false; 
		case hbe_dAMOaCXA: return false; 
		case hbe_dGDEaCXA: return false; 
		case hbe_dGDHaCXA: return false; 
		case hbe_dAHXaCXA: return false; 
		case hbe_dHXLaCXA: return false; 
		case hbe_dH2OaCXA: return false; 
		case hbe_dPBAaCXLsepPM1: return false; 
		case hbe_dPBAaCXLsepother: return false; 
		case hbe_dCXAaCXL: return false; 
		case hbe_dIMDaCXL: return false; 
		case hbe_dIMEaCXL: return false; 
		case hbe_dINDaCXL: return false; 
		case hbe_dAMOaCXL: return false; 
		case hbe_dGDEaCXL: return false; 
		case hbe_dGDHaCXL: return false; 
		case hbe_dAHXaCXL: return false; 
		case hbe_dHXLaCXL: return false; 
		case hbe_dH2OaCXL: return false; 
		case hbe_dPBAaIMDsepPM1: return false; 
		case hbe_dPBAaIMDsepother: return false; 
		case hbe_dCXAaIMD: return false; 
		case hbe_dIMDaIMD: return false; 
		case hbe_dIMEaIMD: return false; 
		case hbe_dINDaIMD: return false; 
		case hbe_dAMOaIMD: return false; 
		case hbe_dGDEaIMD: return false; 
		case hbe_dGDHaIMD: return false; 
		case hbe_dAHXaIMD: return false; 
		case hbe_dHXLaIMD: return false; 
		case hbe_dH2OaIMD: return false; 
		case hbe_dPBAaIMEsepPM1: return false; 
		case hbe_dPBAaIMEsepother: return false; 
		case hbe_dCXAaIME: return false; 
		case hbe_dIMDaIME: return false; 
		case hbe_dIMEaIME: return false; 
		case hbe_dINDaIME: return false; 
		case hbe_dAMOaIME: return false; 
		case hbe_dGDEaIME: return false; 
		case hbe_dGDHaIME: return false; 
		case hbe_dAHXaIME: return false; 
		case hbe_dHXLaIME: return false; 
		case hbe_dH2OaIME: return false; 
		case hbe_dPBAaAHXsepPM1: return false; 
		case hbe_dPBAaAHXsepother: return false; 
		case hbe_dCXAaAHX: return false; 
		case hbe_dIMDaAHX: return false; 
		case hbe_dIMEaAHX: return false; 
		case hbe_dINDaAHX: return false; 
		case hbe_dAMOaAHX: return false; 
		case hbe_dGDEaAHX: return false; 
		case hbe_dGDHaAHX: return false; 
		case hbe_dAHXaAHX: return false; 
		case hbe_dHXLaAHX: return false; 
		case hbe_dH2OaAHX: return false; 
		case hbe_dPBAaHXLsepPM1: return false; 
		case hbe_dPBAaHXLsepother: return false; 
		case hbe_dCXAaHXL: return false; 
		case hbe_dIMDaHXL: return false; 
		case hbe_dIMEaHXL: return false; 
		case hbe_dINDaHXL: return false; 
		case hbe_dAMOaHXL: return false; 
		case hbe_dGDEaHXL: return false; 
		case hbe_dGDHaHXL: return false; 
		case hbe_dAHXaHXL: return false; 
		case hbe_dHXLaHXL: return false; 
		case hbe_dH2OaHXL: return false; 
		case hbe_dPBAaPCA_DNAsepPM1: return false; 
		case hbe_dPBAaPCA_DNAsepother: return false; 
		case hbe_dCXAaPCA_DNA: return false; 
		case hbe_dIMDaPCA_DNA: return false; 
		case hbe_dIMEaPCA_DNA: return false; 
		case hbe_dINDaPCA_DNA: return false; 
		case hbe_dAMOaPCA_DNA: return false; 
		case hbe_dGDEaPCA_DNA: return false; 
		case hbe_dGDHaPCA_DNA: return false; 
		case hbe_dAHXaPCA_DNA: return false; 
		case hbe_dHXLaPCA_DNA: return false; 
		case hbe_dH2OaPCA_DNA: return false; 
		case hbe_dPBAaPCA_RNAsepPM1: return true; 
		case hbe_dPBAaPCA_RNAsepother: return true; 
		case hbe_dCXAaPCA_RNAsepPM1: return false; 
		case hbe_dCXAaPCA_RNAsepother: return false; 
		case hbe_dIMDaPCA_RNAsepPM1: return false; 
		case hbe_dIMDaPCA_RNAsepother: return false; 
		case hbe_dIMEaPCA_RNAsepPM1: return false; 
		case hbe_dIMEaPCA_RNAsepother: return false; 
		case hbe_dINDaPCA_RNAsepPM1: return false; 
		case hbe_dINDaPCA_RNAsepother: return false; 
		case hbe_dAMOaPCA_RNAsepPM1: return false; 
		case hbe_dAMOaPCA_RNAsepother: return false; 
		case hbe_dGDEaPCA_RNAsepPM1: return false; 
		case hbe_dGDEaPCA_RNAsepother: return false; 
		case hbe_dGDHaPCA_RNAsepPM1: return false; 
		case hbe_dGDHaPCA_RNAsepother: return false; 
		case hbe_dAHXaPCA_RNAsepPM1: return false; 
		case hbe_dAHXaPCA_RNAsepother: return false; 
		case hbe_dHXLaPCA_RNAsepPM1: return false; 
		case hbe_dHXLaPCA_RNAsepother: return false; 
		case hbe_dH2OaPCA_RNA: return false; 
		case hbe_dPBAaPES_DNAsepPM1: return false; 
		case hbe_dPBAaPES_DNAsepother: return false; 
		case hbe_dCXAaPES_DNA: return false; 
		case hbe_dIMDaPES_DNA: return false; 
		case hbe_dIMEaPES_DNA: return false; 
		case hbe_dINDaPES_DNA: return false; 
		case hbe_dAMOaPES_DNA: return false; 
		case hbe_dGDEaPES_DNA: return false; 
		case hbe_dGDHaPES_DNA: return false; 
		case hbe_dAHXaPES_DNA: return false; 
		case hbe_dHXLaPES_DNA: return false; 
		case hbe_dH2OaPES_DNA: return false; 
		case hbe_dPBAaPES_RNAsepPM1: return true; 
		case hbe_dPBAaPES_RNAsepother: return true; 
		case hbe_dCXAaPES_RNAsepPM1: return false; 
		case hbe_dCXAaPES_RNAsepother: return false; 
		case hbe_dIMDaPES_RNAsepPM1: return false; 
		case hbe_dIMDaPES_RNAsepother: return false; 
		case hbe_dIMEaPES_RNAsepPM1: return false; 
		case hbe_dIMEaPES_RNAsepother: return false; 
		case hbe_dINDaPES_RNAsepPM1: return false; 
		case hbe_dINDaPES_RNAsepother: return false; 
		case hbe_dAMOaPES_RNAsepPM1: return false; 
		case hbe_dAMOaPES_RNAsepother: return false; 
		case hbe_dGDEaPES_RNAsepPM1: return false; 
		case hbe_dGDEaPES_RNAsepother: return false; 
		case hbe_dGDHaPES_RNAsepPM1: return false; 
		case hbe_dGDHaPES_RNAsepother: return false; 
		case hbe_dAHXaPES_RNAsepPM1: return false; 
		case hbe_dAHXaPES_RNAsepother: return false; 
		case hbe_dHXLaPES_RNAsepPM1: return false; 
		case hbe_dHXLaPES_RNAsepother: return false; 
		case hbe_dH2OaPES_RNA: return false; 
		case hbe_dPBAaRRI_DNAsepPM1: return false; 
		case hbe_dPBAaRRI_DNAsepother: return false; 
		case hbe_dCXAaRRI_DNA: return false; 
		case hbe_dIMDaRRI_DNA: return false; 
		case hbe_dIMEaRRI_DNA: return false; 
		case hbe_dINDaRRI_DNA: return false; 
		case hbe_dAMOaRRI_DNA: return false; 
		case hbe_dGDEaRRI_DNA: return false; 
		case hbe_dGDHaRRI_DNA: return false; 
		case hbe_dAHXaRRI_DNA: return false; 
		case hbe_dHXLaRRI_DNA: return false; 
		case hbe_dH2OaRRI_DNA: return false; 
		case hbe_dPBAaRRI_RNAsepPM1: return false; 
		case hbe_dPBAaRRI_RNAsepother: return false; 
		case hbe_dCXAaRRI_RNAsepPM1: return false; 
		case hbe_dCXAaRRI_RNAsepother: return false; 
		case hbe_dIMDaRRI_RNAsepPM1: return false; 
		case hbe_dIMDaRRI_RNAsepother: return false; 
		case hbe_dIMEaRRI_RNAsepPM1: return false; 
		case hbe_dIMEaRRI_RNAsepother: return false; 
		case hbe_dINDaRRI_RNAsepPM1: return false; 
		case hbe_dINDaRRI_RNAsepother: return false; 
		case hbe_dAMOaRRI_RNAsepPM1: return false; 
		case hbe_dAMOaRRI_RNAsepother: return false; 
		case hbe_dGDEaRRI_RNAsepPM1: return false; 
		case hbe_dGDEaRRI_RNAsepother: return false; 
		case hbe_dGDHaRRI_RNAsepPM1: return false; 
		case hbe_dGDHaRRI_RNAsepother: return false; 
		case hbe_dAHXaRRI_RNAsepPM1: return false; 
		case hbe_dAHXaRRI_RNAsepother: return false; 
		case hbe_dHXLaRRI_RNAsepPM1: return false; 
		case hbe_dHXLaRRI_RNAsepother: return false; 
		case hbe_dH2OaRRI_RNA: return false; 
		case hbe_dPBAaH2O: return false; 
		case hbe_dCXAaH2O: return false; 
		case hbe_dIMDaH2O: return false; 
		case hbe_dIMEaH2O: return false; 
		case hbe_dINDaH2O: return false; 
		case hbe_dAMOaH2O: return false; 
		case hbe_dGDEaH2O: return false; 
		case hbe_dGDHaH2O: return false; 
		case hbe_dAHXaH2O: return false; 
		case hbe_dHXLaH2O: return false; 
		case hbe_dH2OaH2O: return false; 
		case hbe_GENERIC_SP2BB_SR: return true; 
		case hbe_GENERIC_SP2BB_LR: return true; 
		case hbe_GENERIC_SP3BB_SR: return true; 
		case hbe_GENERIC_SP3BB_LR: return true; 
		case hbe_GENERIC_RINGBB_SR: return true; 
		case hbe_GENERIC_RINGBB_LR: return true; 
		case hbe_GENERIC_SP2BSC_SR: return false; 
		case hbe_GENERIC_SP2BSC_LR: return false; 
		case hbe_GENERIC_SP3BSC_SR: return false; 
		case hbe_GENERIC_SP3BSC_LR: return false; 
		case hbe_GENERIC_RINGBSC_SR: return false; 
		case hbe_GENERIC_RINGBSC_LR: return false; 
		case hbe_GENERIC_SP2SCSC_SR: return false; 
		case hbe_GENERIC_SP2SCSC_LR: return false; 
		case hbe_GENERIC_SP3SCSC_SR: return false; 
		case hbe_GENERIC_SP3SCSC_LR: return false; 
		case hbe_GENERIC_RINGSCSC_SR: return false; 
		case hbe_GENERIC_RINGSCSC_LR: return false; 
		}
		return false; // <- to assure gcc a bool will be returned
	}
	static boolean hbe_is_SC_type( HBEvalType hbe ){return !hbe_is_BB_type(hbe);}

	static final double MAX_ATOMIC_INTERACTION = 4.35;
	public void initHBNeighbors(Molecule m) {
		if(numNeighborsHB != null)
			return;
		
		//get CB atoms
		Atom[] CBatoms = new Atom[m.residue.length];
		double[] nbRadii = new double[m.residue.length];
		for(int i=0; i<m.residue.length;i++){
			String resName = m.residue[i].name;
			String nbrAtomName = "CB";
			if(resName.equals("GLY"))
				nbrAtomName = "CA";
			for(Atom a : m.residue[i].atom){
				if(a.name.equals(nbrAtomName)){
					CBatoms[i] = a;
					break;
				}
			}
			switch(m.residue[i].name){
			case "ALA": nbRadii[i] = 3.4473;break;
			case "ARG": nbRadii[i] = 6.1209;break;
			case "ASN": nbRadii[i] = 3.4473;break;
			case "ASP": nbRadii[i] = 3.4473;break;
			case "CYD": nbRadii[i] = 3.4473;break;
			case "CYS": nbRadii[i] = 3.4473;break;
			case "CYV": nbRadii[i] = 3.4473;break;
			case "CYZ": nbRadii[i] = 3.4473;break;
			case "GLN": nbRadii[i] = 3.7348;break;
			case "GLU": nbRadii[i] = 3.6154;break;
			case "GLY": nbRadii[i] = 3.4473;break;
			case "HID": nbRadii[i] = 5.437512;break;
			case "HIE": nbRadii[i] = 3.6718;break;
			case "ILE": nbRadii[i] = 3.4473;break;
			case "LEU": nbRadii[i] = 3.4473;break;
			case "LYS": nbRadii[i] = 5.0084;break;
			case "MET": nbRadii[i] = 4.1766;break;
			case "PHE": nbRadii[i] = 4.2832;break;
			case "PRO": nbRadii[i] = 3.4473;break;
			case "SER": nbRadii[i] = 3.4473;break;
			case "THR": nbRadii[i] = 3.4473;break;
			case "TRP": nbRadii[i] = 5.3639;break;
			case "TYR": nbRadii[i] = 5.6694;break;
			case "VAL": nbRadii[i] = 4.137511;break;
			}
		}
		
		numNeighborsHB = new int[m.residue.length];
		//get number of neighbors
		for(int i=0; i<m.residue.length;i++){
			int numNBR = 0;
			for(int j=0; j<m.residue.length;j++){
				double dist = CBatoms[i].distance(CBatoms[j]);
				if(dist*dist < 100){//(nbRadii[i]+MAX_ATOMIC_INTERACTION+nbRadii[j])*(nbRadii[i]+MAX_ATOMIC_INTERACTION+nbRadii[j])){
					numNBR++;
				}
			}
			numNeighborsHB[i] = numNBR;
		}
		
		
	}
	
	class HbondPair{
		Atom donor;
		Atom hydro;
		Atom accept;
		double multiplier;
		HBondEnergy.HBEvalType hbe;

		HbondPair(Atom d, Atom h, Atom a,double multiplier,HBondEnergy.HBEvalType hbe){
			donor = d;
			hydro = h;
			accept = a;
			this.multiplier = multiplier;
			this.hbe = hbe;
		}
	}

	public void calculateHBondGradient(int curIndex, Molecule m2) {
		System.out.println("Gradient not implemented yet for this Hbond Function.");
		System.out.println("Exiting...");
		System.exit(0);
		
	}

	
}
