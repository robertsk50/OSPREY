import java.io.File;
import java.io.Serializable;
import java.math.BigInteger;

public class Settings {

	public enum DEEMETHOD {
		GOLDSTEIN, SPLITFLAGS1, SPLITFLAGS2, MBPAIRS, FULLPAIRS
	}
	
	public enum ASTARMETHOD{
		ORIG,PGREORDER,ASGUROBI,MIN,LPGUROBI,ASMPLP,BYSEQ,WCSP,BYSUBROT,ASWCSP,ASGUROBIREORDER,ASWCSPREORDER, BYSEQREORDER
	}
	
	public enum CONTRACTMETHOD {
		LEASTPAIRS, PERCENTPRUNED, CLOSESTRES, LARGESTDIFF, PERCENTLEAST
	}
	
	public static String getRunName(ParamSet sParams){
		return ((String)sParams.getValue("RUNNAME"));
	}
	
	
	public class Minimization{
		boolean doMinimize;
		boolean minimizeBB;
		boolean doBackrubs;
		String backrubFile;
		// DEEPer parameters
		boolean doPerturbations;
		boolean pertScreen;
		String pertFile;
		boolean minimizePerts;
		String screenOutFile;
		boolean selectPerturbations;
		
		boolean useCCD;
		
		RotamerSearch.MINIMIZATIONSCHEME minScheme;
		
		Minimization(ParamSet sParams){
			doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
			minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB", "false"))).booleanValue();
			doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
			backrubFile = "";
			if(doBackrubs){
				backrubFile = ((String)sParams.getValue("BACKRUBFILE"));
			}
			// DEEPer parameters
			doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
			pertScreen = (new Boolean((String)sParams.getValue("PERTURBATIONSCREEN","false"))).booleanValue();//Triggers perturbation screen: pruning-only run with rigid perturbations
			pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
			minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
			screenOutFile = ((String)sParams.getValue("SCREENOUTFILE","screenOutFileDefaultName.pert"));//Name of file for outputting results of screen (same format as PERTURBATIONFILE)
			selectPerturbations = (new Boolean((String)sParams.getValue("SELECTPERTURBATIONS","false"))).booleanValue();//Should perturbations be automatically selected?
			Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();
	
			useCCD = (new Boolean((String)sParams.getValue("USECCD","true"))).booleanValue();//CCD minimization
			CCDMinimizer.EConvTol = (new Double((String)sParams.getValue("ECONVTOL","0.01"))).doubleValue();
			
			if (!doMinimize) //no minimization
				minimizeBB = false;
			if (!minimizeBB) //not backbone minimization
				doBackrubs = false;
			
			minScheme = RotamerSearch.MINIMIZATIONSCHEME.valueOf(sParams.getValue("MINIMIZATIONSCHEME","PAIRWISE").toUpperCase());
		}
		
	}
	
	public class DEE implements Serializable {

		int algOption;
		boolean doDACS;
		boolean useFlags;
		boolean distrDACS;
		boolean distrDEE;
		boolean preprocPairs;
		boolean scaleInt;
		double initEw;
		double pruningE;
		double stericE;
		double pairSt;
		double maxIntScale = 0;
		double minRatioDiff=0;
		int initDepth=0;
		int subDepth=0;
		int diffFact=0;
		boolean typeDep;
		int maxFullPairs;
		int maxDEELoopNum;
		boolean doGold;
		boolean doSplit1;
		boolean doSplit2;
		boolean doMB;
		DEEsettings deeSettings;
		boolean useTriples;
		boolean magicBulletTriples;
		int magicBulletNumTriples;
		// 2010: Use energy window MinDEE method.  If this is set to true,
		//   MinDEE will use traditional DEE with an energy window (initEw) 
		//   for pruning.  Max terms will be ignored and only the min terms for pruning and 
		boolean useMinDEEPruningEw;
		double Ival;
		double interval;
		
		
		boolean gold;
		boolean split1;
		boolean split2;
		boolean mb;
		
		
//		DEE(boolean gold,boolean split1,boolean split2,boolean mb){
//			this.gold = gold;
//			this.split1 = split1;
//			this.split2 = split2;
//			this.mb = mb;
//		}
		
		DEE(ParamSet sParams){
			algOption = (new Integer((String)sParams.getValue("ALGOPTION", "3"))).intValue();
			doDACS = (new Boolean((String)sParams.getValue("DODACS", "false"))).booleanValue();
			useFlags = (new Boolean((String)sParams.getValue("SPLITFLAGS", "false"))).booleanValue();
			distrDACS = (new Boolean((String)sParams.getValue("DISTRDACS", "false"))).booleanValue();
			distrDEE = (new Boolean((String)sParams.getValue("DISTRDEE","false"))).booleanValue();
			preprocPairs = (new Boolean((String)sParams.getValue("PREPROCPAIRS", "true"))).booleanValue();
			scaleInt = (new Boolean((String)sParams.getValue("SCALEINT", "false"))).booleanValue();		
			initEw = (new Double((String)sParams.getValue("INITEW","0"))).doubleValue();
			pruningE = (new Double((String)sParams.getValue("PRUNINGE", "100.0"))).doubleValue();
			stericE = (new Double((String)sParams.getValue("STERICE","30.0"))).doubleValue();
			pairSt = (new Double((String)sParams.getValue("PAIRST", "100.0"))).doubleValue();
			maxIntScale =0;
			if(scaleInt){
				maxIntScale = (new Double((String)sParams.getValue("MAXINTSCALE","0"))).doubleValue();
			}
			
			if(doDACS){
				minRatioDiff = (new Double((String)sParams.getValue("MINRATIODIFF"))).doubleValue();
				initDepth = (new Integer((String)sParams.getValue("INITDEPTH"))).intValue();
				subDepth = (new Integer((String)sParams.getValue("SUBDEPTH"))).intValue();
				diffFact = (new Integer((String)sParams.getValue("DIFFFACT"))).intValue();	
			}  
			typeDep = (new Boolean((String)sParams.getValue("TYPEDEP","false"))).booleanValue();
			//DEE Settings
			maxFullPairs = (new Integer((String)sParams.getValue("MAXFULLPAIRS", "-1"))).intValue();
			maxDEELoopNum = (new Integer((String)sParams.getValue("maxDEELoopNum", "1000000"))).intValue();
			doGold = (new Boolean((String)sParams.getValue("DOGOLD", "true"))).booleanValue();
			doSplit1 = (new Boolean((String)sParams.getValue("DOSPLIT1", "false"))).booleanValue();
			doSplit2 = (new Boolean((String)sParams.getValue("DOSPLIT2", "false"))).booleanValue();
			doMB = (new Boolean((String)sParams.getValue("DOMB", "false"))).booleanValue();
			deeSettings = new DEEsettings(doGold,doSplit1,doSplit2,doMB);
			useTriples = (new Boolean((String)sParams.getValue("USETRIPLES","false"))).booleanValue();
			magicBulletTriples = (new Boolean((String)sParams.getValue("MAGICBULLETTRIPLES","true"))).booleanValue();//Use only "magic bullet" competitor triples
			magicBulletNumTriples = (new Integer((String)sParams.getValue("MAGICBULLETNUMTRIPLES","5"))).intValue();//Number of magic bullet triples to use
			// 2010: Use energy window MinDEE method.  If this is set to true,
			//   MinDEE will use traditional DEE with an energy window (initEw) 
			//   for pruning.  Max terms will be ignored and only the min terms for pruning and 
			useMinDEEPruningEw = (new Boolean((String)sParams.getValue("imindee", "true"))).booleanValue();
			Ival = 0.0f;
			interval = 0;
			if(useMinDEEPruningEw)
				Ival = (new Double((String)sParams.getValue("IVAL", "5.0"))).doubleValue();
		}
		
		int bak_maxFullPairs;
		void setMaxFullPairs(int maxFullPairs){
			bak_maxFullPairs = this.maxFullPairs;
			this.maxFullPairs = maxFullPairs;
		}
		
		void resetMaxFullPairs(){
			this.maxFullPairs = bak_maxFullPairs;
		}
		
		int bak_maxDEELoopNum;
		void setMaxDEELoopNum(int maxDEELoopNum){
			bak_maxDEELoopNum = this.maxDEELoopNum;
			this.maxDEELoopNum = maxDEELoopNum;
		}
		
		void resetMaxDEELoopNum(){
			this.maxDEELoopNum = bak_maxDEELoopNum;
		}
		
	}
	
	public class Enum{
		
		boolean approxMinGMEC;
		double lambda;
		boolean useFlagsAStar;
		ASTARMETHOD asMethod;
		int numMaxMut;
		double bestE;
		
		Enum(ParamSet sParams){
			EnvironmentVars.useMPLP = (new Boolean((String)sParams.getValue("USEMPLP", "false"))).booleanValue(); 		//from Pablo
			approxMinGMEC = (new Boolean((String)sParams.getValue("APPROXMINGMEC", "false"))).booleanValue();
			lambda = (new Double((String)sParams.getValue("LAMBDA", "0"))).doubleValue();
			useFlagsAStar = (new Boolean((String)sParams.getValue("USEFLAGSASTAR","false"))).booleanValue();
			asMethod = ASTARMETHOD.valueOf(sParams.getValue("ASTARMETHOD","ASWCSPREORDER").toUpperCase());
			numMaxMut = (new Integer((String)sParams.getValue("NUMMAXMUT", "1000"))).intValue();
			bestE = (new Double((String)sParams.getValue("BESTE","10000000"))).doubleValue();
		}
		
	}
	
	public class Emat{
		
		String runNameEMatrixMin;
		boolean useEref;
		boolean addWTRotsSomehow;
		boolean addOrigRots;
		boolean addWTRot;
		
		Emat(ParamSet sParams, String runName, boolean doPerturbations){
			runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
			File tmpFile = new File(runNameEMatrixMin);
			File ematDir = tmpFile.getParentFile();
			if(ematDir != null && !ematDir.exists()){
				ematDir.mkdirs();
			}
			useEref = (new Boolean((String)sParams.getValue("USEEREF","true"))).booleanValue();
			addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
			addOrigRots = false;
			addWTRot = false;
			if( addWTRotsSomehow ){
				if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
					addWTRot = true;
				else//Otherwise just add them to the rotamer library
					addOrigRots = true;
			}
		}
		
	}
	
	public class InteractionGraph{
		
		boolean genInteractionGraph;
		double distCutoff=0;
		double eInteractionCutoff=0;
		boolean neighborList;
		
		
		InteractionGraph(ParamSet sParams){
			genInteractionGraph = (new Boolean((String)sParams.getValue("GENINTERACTIONGRAPH","false"))).booleanValue();
			distCutoff=0;
			eInteractionCutoff=0;
			if(genInteractionGraph){
				distCutoff = (new Double((String)sParams.getValue("DISTCUTOFF"))).doubleValue();
				eInteractionCutoff = (new Double((String)sParams.getValue("EINTERACTIONCUTOFF"))).doubleValue();
			}
			neighborList = new Boolean(sParams.getValue("NEIGHBORLIST","false"));
			if(neighborList){
				distCutoff = new Float(sParams.getValue("DISTCUTOFF"));
			}
			
		}
		
		
	}
	
	public class Output{
		String outputConfInfo;
		String outputPruneInfo;
		String pdbOutDir;
		
		Output(ParamSet sParams, String runName){
			outputConfInfo = (String)(sParams.getValue("OUTPUTCONFINFO","c_"+runName));
			outputPruneInfo = (String)(sParams.getValue("OUTPUTPRUNEINFO","p_"+runName));
			pdbOutDir = sParams.getValue("pdbOutDir","pdbs");
		}
	}
	
	public class KStar{
		
		int numMutations;
		String mutFileName;
		boolean repeatSearch;
		double targetVol;
		double volWindow;
		boolean resumeSearch;
		String resumeFilename;
		double gamma;
		double epsilon;
		boolean saveTopConfs;
		boolean printTopConfs;
		int numTopConfs;
		boolean useMaxKSconfs;
		BigInteger maxKSconfs;
		
		KStar(ParamSet sParams, String runName){
			numMutations = (new Integer((String)sParams.getValue("NUMMUTATIONS", "1"))).intValue();
			mutFileName = (String)sParams.getValue("MUTFILENAME",runName+".mut");
			repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH","true"))).booleanValue();

			targetVol = (new Double((String)sParams.getValue("TARGETVOLUME","0.0"))).doubleValue();
			volWindow = (new Double((String)sParams.getValue("VOLUMEWINDOW","50000000"))).doubleValue();
			resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
			resumeFilename = "";
			if(resumeSearch){
				resumeFilename = (String)sParams.getValue("RESUMEFILENAME");
			}
			gamma = (new Double((String)sParams.getValue("GAMMA", "0"))).doubleValue();
			epsilon = (new Double((String)sParams.getValue("EPSILON","0.3"))).doubleValue();
			
			saveTopConfs = (new Boolean((String)sParams.getValue("SAVETOPCONFSASPDB","false"))).booleanValue();
			printTopConfs = (new Boolean((String)sParams.getValue("SAVETOPCONFSROTS","false"))).booleanValue();
			numTopConfs = (new Integer((String)sParams.getValue("NUMTOPCONFSTOSAVE","0"))).intValue();

			if(printTopConfs || saveTopConfs){
				//KER: make directory for the confs to be printed to.
				File ksConfDir = new File(EnvironmentVars.ksConfDir);
				if(!ksConfDir.exists())
					ksConfDir.mkdir();
			}
			
			useMaxKSconfs = new Boolean((String)sParams.getValue("useMaxKSconfs","false")).booleanValue();
			maxKSconfs = BigInteger.ZERO;
			if(useMaxKSconfs)
				maxKSconfs = new BigInteger(sParams.getValue("maxKSconfs"));
			
		}
		
		
	}
	
	public class ImprovedBounds{
	
		boolean superRotamers;
		boolean subRotamers;
		boolean doTuples;
		boolean readTuples;
		String tupleFile;
		boolean keepAStree;
		boolean addTuplesByDistance;
		
		CONTRACTMETHOD contractMethod;
//		EnvironmentVars.setContractMethod(contractMethod);
		
		public ImprovedBounds(ParamSet sParams) {
		
			superRotamers = (new Boolean((String)sParams.getValue("SUPERROTAMERS", "false"))).booleanValue();
			subRotamers = (new Boolean((String)sParams.getValue("SUBROTAMERS", "false"))).booleanValue();
			doTuples = (new Boolean((String)sParams.getValue("tuples", "false"))).booleanValue();
			readTuples = (new Boolean((String)sParams.getValue("readtuples", "false"))).booleanValue();
			tupleFile = "";
			if(readTuples)
				tupleFile = sParams.getValue("tupleFile");
			keepAStree = (new Boolean((String)sParams.getValue("keepAStree", "false"))).booleanValue();
			addTuplesByDistance = new Boolean(sParams.getValue("ADDTUPLESBYDISTANCE","false"));
			String contractMethodStr = sParams.getValue("CONTRACTMETHOD","LEASTPAIRS");
			contractMethod = CONTRACTMETHOD.valueOf(contractMethodStr);
		}
	}
}
