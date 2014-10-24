import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.PriorityBlockingQueue;




public class WCSPOptimization {

	static final long MAXUPPER = 1537228672809129301L;

	static final String[] wcsp2preprocLoc = {EnvironmentVars.dataDir+"/exe/toulbar2"};
	static final String[] wcsp2optLoc     = {EnvironmentVars.dataDir+"/exe/toulbar2"};


	//	static final String[] wcsp2preprocLoc = {"/usr/project/dlab/Users/kroberts/toulbar/toulbar2.0.9.6.0-kro/build_release/bin/Linux/toulbar2"};
	//"/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build_release/bin/Linux/toulbar2"};
	//	static final String[] wcsp2optLoc = {"/usr/project/dlab/Users/kroberts/toulbar/toulbar2.0.9.6.0-kro/build_release/bin/Linux/toulbar2"};
	//"/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build_release/bin/Linux/toulbar2"};

	//	static final String wcsp2preprocLocBack = "/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build_release/bin/Linux/toulbar2";
	//	static final String wcsp2optLocBack = "/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build_release/bin/Linux/toulbar2";
	//static final String wcsp2optLoc = "/usr/project/dlab/Users/kroberts/Troubleshooting/toulbar2/SpeedUp/example.1MJC/bin/toulbar2";
	//static final String wcsp2preprocLoc = "/usr/project/dlab/Users/kroberts/Troubleshooting/toulbar2/SpeedUp/example.1MJC/bin/toulbar2";
	String outDir = EnvironmentVars.localDir;
	String filename = "";
	long upperBound = MAXUPPER;
	//map of pos and index to arrayIndex
	Index3 rots[][];
	int osprey2wcsp[][][];
	Emat emat;
	int numTotalNodes;
	double lowestE;
	int confLength;
	int numConstraints;

	RotConf bestConf; 

	PriorityBlockingQueue<RotConf> allConfs;

	//For by sequence nodes
	public WCSPOptimization(PGQueueNode node,
			Emat emat,int[][] seqIndicesPerLevel,
			int numRotRemainingBySeq[][],
			int numNodesForLevel[],double upperE) {

		this.emat = emat;
		allConfs = new PriorityBlockingQueue<RotConf>();
		this.numTotalNodes = numTotalNodes;
		this.outDir += System.getProperty("user.dir")+File.separator+"dat";
		int random = new Random().nextInt(1000000);
		this.filename = "toulbar_"+random+"_"+System.currentTimeMillis()+".wcsp";


		//Find Min E and get Problem Constraints
		int numVars = node.confSoFar.length;
		confLength = numVars;
		int maxDom = -1;

		lowestE = Double.POSITIVE_INFINITY;
		rots = new Index3[node.confSoFar.length][];
		osprey2wcsp = new int[emat.singles.E.length][][];
		for(int p1=0; p1<node.confSoFar.length;p1++){

			if(node.confSoFar[p1] >= 0){ //Is not empty
				//Set rots to have the same number of rotamer entries as the assigned sequence does
				rots[p1] = new Index3[numRotRemainingBySeq[p1][node.confSoFar[p1]]];
				osprey2wcsp[p1] = new int[emat.singles.E[p1].length][];
				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;
				int ctr = 0;
				for(int r1=0; r1<emat.singles.E[p1][seqIndicesPerLevel[p1][node.confSoFar[p1]]].length;r1++){
					Index3 rot = new Index3(p1,seqIndicesPerLevel[p1][node.confSoFar[p1]],r1);// nodeIndexOffset[p1]+seqIndexOffset[p1][node.confSoFar[p1]] + r1;
					if(!emat.getSinglePruned(rot)){
						rots[p1][r1] = rot;
						if(osprey2wcsp[rot.pos][rot.aa]== null )
							osprey2wcsp[rot.pos][rot.aa]= new int[emat.singles.E[rot.pos][rot.aa].length]; 
						osprey2wcsp[rot.pos][rot.aa][rot.rot] = ctr;
						double minIndVoxE = emat.getSingleMinE(rot);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						if(minIndVoxE < lowestE)
							lowestE = minIndVoxE;
						ctr++;
					}	

				}
			}
			else{
				rots[p1] = new Index3[numNodesForLevel[p1]];
				osprey2wcsp[p1] = new int[emat.singles.E[p1].length][];
				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;

				SinglesIterator siter = emat.singlesIterator(p1);
				int r1 = 0;
				while(siter.hasNext()){
					EMatrixEntryWIndex emeWI = siter.next();
					if(!emeWI.eme.isPruned()){
						Index3 rot = new Index3(p1,emeWI.aa1(),emeWI.rot1());
						rots[p1][r1] = rot;
						if(osprey2wcsp[rot.pos][rot.aa]== null )
							osprey2wcsp[rot.pos][rot.aa]= new int[emat.singles.E[rot.pos][rot.aa].length];
						osprey2wcsp[rot.pos][rot.aa][rot.rot] = r1;
						double minIndVoxE = emat.getSingleMinE(rot);
						if(minIndVoxE < lowestE)
							lowestE = minIndVoxE;

						r1++;
					}
				}
			}
		}

		//Go through pairs
		int numPairs = 0;
		for(int p1=0; p1<node.confSoFar.length;p1++){
			for(int p2=p1+1;p2<node.confSoFar.length;p2++){
				if(p1!=p2 && emat.areNeighbors(p1, p2)){
					numPairs++;
					for(int r1=0;r1<rots[p1].length;r1++){
						Index3 rot1 = rots[p1][r1];
						for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
							Index3 rot2 = rots[p2][r2];

							if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
								continue;

							double E = emat.getPairMinE(rot1, rot2);
							if(E < lowestE)
								lowestE = E;

						}
					}
				}
			}
		}
		numConstraints = numVars+numPairs;

		if(upperE == Double.POSITIVE_INFINITY || upperE >= MAXUPPER)
			upperBound = MAXUPPER;
		else
			upperBound = getCostForFullConf(upperE, MAXUPPER);

		//Make dat dir if it doesn't exist
		File datDir = new File(outDir);
		if(!datDir.exists()){
			datDir.mkdirs();
		}

		try{  
			FileWriter fstream = new FileWriter(outDir+File.separator+filename, false);
			BufferedWriter out = new BufferedWriter(fstream);

			KSParser.outputObject(lowestE,outDir+File.separator+filename+".minE");

			//Write header
			out.write("toulbar_osprey "+numVars+" "+maxDom+" "+numConstraints+" "+upperBound+"\n");
			for(Index3 i[]:rots )
				out.write(i.length+" ");
			out.write("\n");


			for(int p1=0; p1<node.confSoFar.length;p1++){
				int ctr=0;
				out.write("1 "+p1+" "+upperBound+" "+rots[p1].length+"\n");

				for(int r1=0; r1<rots[p1].length;r1++){
					Index3 rot = rots[p1][r1];
					double E = emat.getSingleMinE(rot);
					String cost = getCost(E,upperBound);
					out.write(r1+" "+cost+"\n");
				}

			}

			//Go through pairs

			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int p2=p1+1;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						String[] s = new String[rots[p1].length*rots[p2].length];
						numPairs = 0;
						for(int r1=0;r1<rots[p1].length;r1++){
							Index3 rot1 = rots[p1][r1];
							for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 rot2 = rots[p2][r2];

								if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
									continue;

								double E = emat.getPairMinE(rot1, rot2);
								String cost = getCost(E,upperBound);

								s[numPairs] = (r1+" "+r2+" "+cost+"\n"); 
								numPairs++;

							}
						}
						out.write("2 "+p1+" "+p2+" "+upperBound+" "+(numPairs)+"\n");
						for(int i=0; i<numPairs;i++)
							out.write(s[i]);
					}
				}
			}

			out.close(); 
		}catch(Exception E){
			E.printStackTrace();
		}

	}

	//For subrotamer nodes
	public WCSPOptimization(PGQueueNode node,
			Emat emat,int[] numParentRotPerLvl,
			int[][] parentRotIndexPerLvl, int[][] numSubRotPerParentRot,
			ArrayList<HashMap<Integer,ArrayList<Index3>>> subRotsPerLvlPerParent,
			int numNodesForLevel[],double upperE) {

		this.emat = emat;
		allConfs = new PriorityBlockingQueue<RotConf>();
		this.numTotalNodes = numTotalNodes;
		this.outDir += System.getProperty("user.dir")+File.separator+"dat";
		int random = new Random().nextInt(1000000);
		this.filename = "toulbar_"+random+"_"+System.currentTimeMillis()+".wcsp";


		//Find Min E and get Problem Constraints
		int numVars = node.confSoFar.length;
		confLength = numVars;
		int maxDom = -1;

		lowestE = Double.POSITIVE_INFINITY;
		rots = new Index3[node.confSoFar.length][];
		osprey2wcsp = new int[emat.singles.E.length][][];
		for(int p1=0; p1<node.confSoFar.length;p1++){

			if(node.confSoFar[p1] >= 0){ //Is not empty
				//Set rots to have the same number of rotamer entries as the assigned sequence does
				rots[p1] = new Index3[numSubRotPerParentRot[p1][node.confSoFar[p1]]];
				osprey2wcsp[p1] = new int[emat.singles.E[p1].length][];
				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;
				ArrayList<Index3> subRotamers = subRotsPerLvlPerParent.get(p1).get(parentRotIndexPerLvl[p1][node.confSoFar[p1]]);
				for(int r1=0; r1<subRotamers.size();r1++){
					Index3 rot = subRotamers.get(r1);
					rots[p1][r1] = rot;
					if(osprey2wcsp[rot.pos][rot.aa]== null )
						osprey2wcsp[rot.pos][rot.aa]= new int[emat.singles.E[rot.pos][rot.aa].length]; 
					osprey2wcsp[rot.pos][rot.aa][rot.rot] = r1;
					double minIndVoxE = emat.getSingleMinE(rot);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
					if(minIndVoxE < lowestE)
						lowestE = minIndVoxE;

				}
			}
			else{
				rots[p1] = new Index3[numNodesForLevel[p1]];
				osprey2wcsp[p1] = new int[emat.singles.E[p1].length][];
				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;

				SinglesIterator siter = emat.singlesIterator(p1);
				int r1 = 0;
				while(siter.hasNext()){
					EMatrixEntryWIndex emeWI = siter.next();
					if(!emeWI.eme.isPruned()){
						Index3 rot = new Index3(p1,emeWI.aa1(),emeWI.rot1());
						rots[p1][r1] = rot;
						if(osprey2wcsp[rot.pos][rot.aa]== null )
							osprey2wcsp[rot.pos][rot.aa]= new int[emat.singles.E[rot.pos][rot.aa].length];
						osprey2wcsp[rot.pos][rot.aa][rot.rot] = r1;
						double minIndVoxE = emat.getSingleMinE(rot);
						if(minIndVoxE < lowestE)
							lowestE = minIndVoxE;

						r1++;
					}
				}
			}
		}

		//Go through pairs
		int numPairs = 0;
		for(int p1=0; p1<node.confSoFar.length;p1++){
			for(int p2=p1+1;p2<node.confSoFar.length;p2++){
				if(p1!=p2 && emat.areNeighbors(p1, p2)){
					numPairs++;
					for(int r1=0;r1<rots[p1].length;r1++){
						Index3 rot1 = rots[p1][r1];
						for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
							Index3 rot2 = rots[p2][r2];

							if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
								continue;

							double E = emat.getPairMinE(rot1, rot2);
							if(E < lowestE)
								lowestE = E;

						}
					}
				}
			}
		}
		numConstraints = numVars+numPairs;

		if(upperE == Double.POSITIVE_INFINITY)
			upperBound = MAXUPPER;
		else
			upperBound = getCostForFullConf(upperE, MAXUPPER);

		//Make dat dir if it doesn't exist
		File datDir = new File(outDir);
		if(!datDir.exists()){
			datDir.mkdirs();
		}

		try{  
			FileWriter fstream = new FileWriter(outDir+File.separator+filename, false);
			BufferedWriter out = new BufferedWriter(fstream);

			KSParser.outputObject(lowestE,outDir+File.separator+filename+".minE");

			//Write header
			out.write("toulbar_osprey "+numVars+" "+maxDom+" "+numConstraints+" "+upperBound+"\n");
			for(Index3 i[]:rots )
				out.write(i.length+" ");
			out.write("\n");


			for(int p1=0; p1<node.confSoFar.length;p1++){
				int ctr=0;
				out.write("1 "+p1+" "+upperBound+" "+rots[p1].length+"\n");

				for(int r1=0; r1<rots[p1].length;r1++){
					Index3 rot = rots[p1][r1];
					double E = emat.getSingleMinE(rot);
					String cost = getCost(E,upperBound);
					out.write(r1+" "+cost+"\n");
				}

			}

			//Go through pairs

			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int p2=p1+1;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						String[] s = new String[rots[p1].length*rots[p2].length];
						numPairs = 0;
						for(int r1=0;r1<rots[p1].length;r1++){
							Index3 rot1 = rots[p1][r1];
							for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 rot2 = rots[p2][r2];

								if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
									continue;

								double E = emat.getPairMinE(rot1, rot2);
								String cost = getCost(E,upperBound);

								s[numPairs] = (r1+" "+r2+" "+cost+"\n"); 
								numPairs++;

							}
						}
						out.write("2 "+p1+" "+p2+" "+upperBound+" "+(numPairs)+"\n");
						for(int i=0; i<numPairs;i++)
							out.write(s[i]);
					}
				}
			}

			out.close(); 
		}catch(Exception E){
			E.printStackTrace();
		}

	}

	//For Normal A* Bounds
	public WCSPOptimization(PGQueueNode node,
			Emat emat, int numNodesForLevel[],double upperE,Index3[][] twoDTo3D) {

		this.emat = emat;
		allConfs = new PriorityBlockingQueue<RotConf>();
		this.numTotalNodes = numTotalNodes;
		this.outDir += System.getProperty("user.dir")+File.separator+"dat";
		int random = new Random().nextInt(1000000);
		this.filename = "toulbar_"+random+"_"+System.currentTimeMillis()+".wcsp";


		//Find Min E and get Problem Constraints
		int numVars = node.confSoFar.length;
		confLength = numVars;
		int maxDom = -1;

		lowestE = Double.POSITIVE_INFINITY;
		rots = new Index3[node.confSoFar.length][];

		for(int p1=0; p1<node.confSoFar.length;p1++){

			if(node.confSoFar[p1] >= 0){ //Is not empty
				//Set rots to have the single rotamer entry
				rots[p1] = new Index3[1];

				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;
				Index3 rot = twoDTo3D[p1][node.confSoFar[p1]];
				rots[p1][0] = rot;
				double minIndVoxE = emat.getSingleMinE(rot);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
				if(minIndVoxE < lowestE)
					lowestE = minIndVoxE;


			}
			else{
				rots[p1] = new Index3[numNodesForLevel[p1]];

				if(rots[p1].length > maxDom)
					maxDom = rots[p1].length;

				SinglesIterator siter = emat.singlesIterator(p1);
				int r1 = 0;
				while(siter.hasNext()){
					EMatrixEntryWIndex emeWI = siter.next();
					if(!emeWI.eme.isPruned()){
						Index3 rot = new Index3(p1,emeWI.aa1(),emeWI.rot1());
						rots[p1][r1] = rot;

						double minIndVoxE = emat.getSingleMinE(rot);
						if(minIndVoxE < lowestE)
							lowestE = minIndVoxE;

						r1++;
					}
				}
			}
		}

		//Go through pairs
		int numPairs = 0;
		for(int p1=0; p1<node.confSoFar.length;p1++){
			for(int p2=p1+1;p2<node.confSoFar.length;p2++){
				if(p1!=p2 && emat.areNeighbors(p1, p2)){
					numPairs++;
					for(int r1=0;r1<rots[p1].length;r1++){
						Index3 rot1 = rots[p1][r1];
						for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
							Index3 rot2 = rots[p2][r2];

							if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
								continue;

							double E = emat.getPairMinE(rot1, rot2);
							if(E < lowestE)
								lowestE = E;

						}
					}
				}
			}
		}
		numConstraints = numVars+numPairs;

		if(upperE == Double.POSITIVE_INFINITY)
			upperBound = MAXUPPER;
		else
			upperBound = getCostForFullConf(upperE, MAXUPPER);

		//Make dat dir if it doesn't exist
		File datDir = new File(outDir);
		if(!datDir.exists()){
			datDir.mkdirs();
		}

		try{  
			FileWriter fstream = new FileWriter(outDir+File.separator+filename, false);
			BufferedWriter out = new BufferedWriter(fstream);

			KSParser.outputObject(lowestE,outDir+File.separator+filename+".minE");

			//Write header
			out.write("toulbar_osprey "+numVars+" "+maxDom+" "+numConstraints+" "+upperBound+"\n");
			for(Index3 i[]:rots )
				out.write(i.length+" ");
			out.write("\n");


			for(int p1=0; p1<node.confSoFar.length;p1++){
				int ctr=0;
				out.write("1 "+p1+" "+upperBound+" "+rots[p1].length+"\n");

				for(int r1=0; r1<rots[p1].length;r1++){
					Index3 rot = rots[p1][r1];
					double E = emat.getSingleMinE(rot);
					String cost = getCost(E,upperBound);
					out.write(r1+" "+cost+"\n");
				}

			}

			//Go through pairs

			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int p2=p1+1;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						String[] s = new String[rots[p1].length*rots[p2].length];
						numPairs = 0;
						for(int r1=0;r1<rots[p1].length;r1++){
							Index3 rot1 = rots[p1][r1];
							for(int r2=0;r2<rots[p2].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 rot2 = rots[p2][r2];

								if(emat.getPairPruned(rot1, rot2)) //Will be null if not neighbors
									continue;

								double E = emat.getPairMinE(rot1, rot2);
								String cost = getCost(E,upperBound);

								s[numPairs] = (r1+" "+r2+" "+cost+"\n"); 
								numPairs++;

							}
						}
						out.write("2 "+p1+" "+p2+" "+upperBound+" "+(numPairs)+"\n");
						for(int i=0; i<numPairs;i++)
							out.write(s[i]);
					}
				}
			}

			out.close(); 
		}catch(Exception E){
			E.printStackTrace();
		}

	}


	public double optimize(String[] additionalCommands){
		//String[] commands = {"/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build/bin/Linux/toulbar2","-s",outDir+File.separator+filename};

		int commandLength = 0;
		if(additionalCommands != null)
			commandLength = additionalCommands.length;

		String[] commands = new String[5+commandLength];

		for(int i=0; i<commandLength;i++){
			commands[i+1] = additionalCommands[i];
		}
		commands[commandLength+1] = "-s";
		commands[commandLength+2] = "-e:";
		commands[commandLength+3] = "-f:";
		commands[commandLength+4] = outDir+File.separator+filename;

		String sol = "";
		boolean solFound = false;
		int ctr=0;
		while(!solFound && ctr<wcsp2optLoc.length){	
			commands[0] = wcsp2optLoc[ctr++];
			try {
				Process p = Runtime.getRuntime().exec(commands);
				// any error message?
				//			StreamGobbler errorGobbler = new 
				//                StreamGobbler(p.getErrorStream(), "ERROR");            

				// any output?
				//            StreamGobbler outputGobbler = new 
				//                StreamGobbler(p.getInputStream(), "OUTPUT");

				// kick them off
				//            errorGobbler.start();
				//            outputGobbler.start();

				ArrayList<String> output = new ArrayList<String>();
				InputStreamReader isr = new InputStreamReader(p.getInputStream());
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				boolean newSolution = false;
				while( (line = br.readLine()) != null){
					output.add(line);

					if(newSolution){
						sol = line;
						solFound = true;
					}

					if(line.contains("New solution")){
						newSolution = true;
					}else{
						newSolution = false;
					}

					//				if(line.startsWith("Optimum: ")){
					//					String[] aLine = line.split(" ");
					//					optimum = new Long(aLine[1]);
					//				}


				}


				// any error???
				int exitVal = p.waitFor();
				//            System.out.println("ExitValue: " + exitVal); 

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		//Find Solution
		if(sol.equals(""))
			return Double.POSITIVE_INFINITY;
		else{
			return wcspSol2ospreyConf(sol);
		}
	}



	public double getBound(String[] additionalCommands){
		//String[] commands = {"/home/home1/kroberts/Downloads/toulbar2.0.9.5.0-Release-sources/build/bin/Linux/toulbar2","-s",outDir+File.separator+filename};

		int commandLength = 0;
		if(additionalCommands != null)
			commandLength = additionalCommands.length;

		String[] commands = new String[6+commandLength];

		for(int i=0; i<commandLength;i++){
			commands[i+1] = additionalCommands[i];
		}
		commands[commandLength+1] = "-s";
		commands[commandLength+2] = "-e:";
		commands[commandLength+3] = "-f:";
		commands[commandLength+4] = "-preproconly";
		commands[commandLength+5] = outDir+File.separator+filename;

		String sol = "";
		long lowerBound = Long.MAX_VALUE;
		boolean solFound = false;
		int ctr = 0;
		while(!solFound && ctr < wcsp2preprocLoc.length){
			commands[0] = wcsp2preprocLoc[ctr++];
			try {
				Process p = Runtime.getRuntime().exec(commands);
				// any error message?
				//			StreamGobbler errorGobbler = new 
				//                StreamGobbler(p.getErrorStream(), "ERROR");            

				// any output?
				//            StreamGobbler outputGobbler = new 
				//                StreamGobbler(p.getInputStream(), "OUTPUT");

				// kick them off
				//            errorGobbler.start();
				//            outputGobbler.start();

				ArrayList<String> output = new ArrayList<String>();
				InputStreamReader isr = new InputStreamReader(p.getInputStream());
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				boolean newSolution = false;
				while( (line = br.readLine()) != null){
					output.add(line);

					//Initial lower and upper bounds: [10,1537228672809129301[
					if(line.contains("Initial lower and upper bounds")){
						String[] aLine = line.split("[ \\[,]");
						lowerBound = new Long( aLine[6] );
						solFound = true;
					}




				}


				// any error???
				int exitVal = p.waitFor();
				//            System.out.println("ExitValue: " + exitVal); 

			} catch (IOException e) {
				if(e.getCause().getMessage().equalsIgnoreCase("error=2, No such file or directory")){
					System.out.println("Could not find the wcsp executable: "+wcsp2preprocLoc[0] );
					System.out.println("Please make sure it exists.");
					System.out.println("Exiting...");
					System.exit(0);
				}
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		//Find Solution
		return getEfromCost(lowerBound, numConstraints);
	}


	public void getAllConfs(double Ew, double bestE){

		long costL = getCostForFullConf(Ew+bestE, upperBound);
		String cost = (new Long(costL)).toString();

		String[] commands = {wcsp2optLoc[0],"-s","-a","-e:","-f:","-ub="+cost,outDir+File.separator+filename};

		String sol = "";
		boolean solFound = false;
		int ctr=0;
		while(!solFound && ctr < wcsp2optLoc.length){
			commands[0] = wcsp2optLoc[ctr++];
			try {
				Process p = Runtime.getRuntime().exec(commands);
				// any error message?
				//			StreamGobbler errorGobbler = new 
				//                StreamGobbler(p.getErrorStream(), "ERROR");            

				// any output?
				//            StreamGobbler outputGobbler = new 
				//                StreamGobbler(p.getInputStream(), "OUTPUT");

				// kick them off
				//            errorGobbler.start();
				//            outputGobbler.start();

				ArrayList<String> output = new ArrayList<String>();
				InputStreamReader isr = new InputStreamReader(p.getInputStream());
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				while( (line = br.readLine()) != null){
					output.add(line);

					if(line.contains("solution:")){
						//Found Solution
						sol = line.substring(line.indexOf(":")+2, line.length());
						wcspSol2ospreyConf(sol);
						solFound = true;
					}


				}


				// any error???
				int exitVal = p.waitFor();
				//            System.out.println("ExitValue: " + exitVal); 

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}


	//
	//	private void outputToulbarWCSP(ParamSet sParams,Molecule m, Emat emat) {
	//	
	//		String runName = ((String)sParams.getValue("RUNNAME"));
	//		String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
	//		boolean useEref = (new Boolean((String)sParams.getValue("USEEREF","true"))).booleanValue();
	//
	//		
	//		
	//		
	//		
	//		if(useEref){
	//			if(emat.eRef == null){
	//				System.out.println("Reference energies turned on but not calculated");
	//				System.out.println("Exiting....");
	//				System.exit(0);
	//			}
	//			if( !emat.hasEref)
	//				RotamerSearch.addEref(emat,m,emat.eRef, -1, null);
	//		}
	//		if(EnvironmentVars.useEntropy && !emat.hasEntropy)
	//			RotamerSearch.addEntropyTerm(emat, m, -1);
	//		
	//		
	////		//KER: Always do as much goldstein as possible
	////		boolean goldDone = false;
	////		while(!goldDone){
	////
	////			int unprunedBeforeGold = countUnprunedRot(emat);
	////
	////			System.out.println("NonPrunedRot: "+ unprunedBeforeGold);	
	////
	////			double initEw = 0.0;
	////			boolean doMinimize = false;
	////			boolean useFlags = false;
	////			boolean minimizeBB = false;
	////			boolean localUseMinDEEPruningEw = true;
	////			double Ival = 0.0;
	////			boolean removeRot = true;
	////			
	////			//Depending on the chosen algorithm option, apply the corresponding pruning criteria;
	////			if(numProc > 2 && emat.numMutPos() > 2 )
	////				doDistrDEEMaster(emat, sParams, initEw, doMinimize, 
	////						minimizeBB, false, 0.0f,Ival,false,false,DEEMETHOD.GOLDSTEIN);
	////			else
	////				RotamerSearch.DoDEEGoldstein(emat,initEw, doMinimize, useFlags, minimizeBB,false,
	////						localUseMinDEEPruningEw,Ival,false,null,null,removeRot);
	////
	////			if(removeRot)
	////				emat.removePrunedRotReducedMem(false);
	////
	////			//check how many rotamers/pairs are pruned this run
	////			int unprunedAfterGold = countUnprunedRot(emat);
	////
	////			System.out.println("NumUnprunedRot: "+unprunedAfterGold );
	////
	////			if(unprunedBeforeGold == unprunedAfterGold)
	////				goldDone = true;
	////
	////		}
	////		emat.write("toulbarMat.dat",aaRotLib);
	//		
	//		int[][][] rotIndices = new int[emat.numMutPos()][][];
	//		
	//		
	//		//Get problem constraints
	//		int numVars = emat.numMutPos();
	//		int[] numRotPerPos = emat.numRotPerPos();
	//		int maxDom  = -1; 
	//		for(int i: numRotPerPos)
	//			if(i>maxDom)
	//				maxDom = i;
	//		int numPairs = 0;
	//		for(int i=0; i<numVars;i++)
	//			for(int j=i+1; j<numVars;j++)
	//				if(emat.areNeighbors(i, j))
	//					numPairs++;
	//		int numConstraints = numVars+numPairs;
	//		long upperBound = 1537228672809129301L;
	//		
	//		//Find Min E
	//		
	//		double minE = Double.POSITIVE_INFINITY;
	//		SinglesIterator tmpiter = emat.singlesIterator();
	//		while(tmpiter.hasNext()){
	//			EMatrixEntryWIndex emeWI = tmpiter.next();
	//			if(emeWI.eme.minE() < minE)
	//				minE = emeWI.eme.minE();
	//		}
	//		PairsIterator piter = emat.pairsIterator();
	//		while(piter.hasNext()){
	//			EMatrixEntryWIndex emeWI = piter.next();
	//			if(emeWI.eme.minE() < minE)
	//				minE = emeWI.eme.minE();
	//		}
	//		
	//		
	//		//Make dat dir if it doesn't exist
	//		File datDir = new File("dat");
	//		if(!datDir.exists()){
	//			datDir.mkdirs();
	//		}
	//		
	//		
	//		try{  
	//		FileWriter fstream = new FileWriter("dat/"+runName+".wcsp", false);
	//        BufferedWriter out = new BufferedWriter(fstream);
	//        
	//        FileWriter fstream2 = new FileWriter("dat/"+runName+".rots", false);
	//        BufferedWriter rotOut = new BufferedWriter(fstream2);
	//        
	//        KSParser.outputObject(minE,"dat/"+runName+".minE");
	//        
	//        
	//        //Write header
	//        out.write(runName+" "+numVars+" "+maxDom+" "+numConstraints+" "+upperBound+"\n");
	//        for(int i: numRotPerPos)
	//        	out.write(i+" ");
	//        out.write("\n");
	//	    
	//        for(int p=0; p<emat.numMutPos();p++){
	//        	int ctr=0;
	//        	
	//        	rotIndices[p] = new int[emat.singles.E[p].length][];
	//        	
	//        	int numRot = emat.numRot(p);
	//		    SinglesIterator siter = emat.singlesIterator(p);
	//		    out.write("1 "+p+" "+upperBound+" "+numRot+"\n");
	//			while(siter.hasNext()){
	//				
	//				EMatrixEntryWIndex emeWI = siter.next();
	//				Residue r = m.residue[emat.resByPos.get(emeWI.pos1()).get(0)];
	//				Rotamer rot = r.rl.getRot(emat.singles.getRot(emeWI.index)[0]);
	//				int rotNum = emeWI.rot1()+1;
	//				
	//				String cost = getCost(minE,emeWI.eme.minE(),upperBound);
	//				  
	////	          out.write(ctr+" ""ROT::BKB\t"+r.getResNumberString()+"-"+rot.aaType.name+"-"+rotNum+" "+emeWI.eme.minE()+" "+emeWI.eme.minE()+" "+rot.aaIndex+" "+rot.rlIndex+"\n");
	//	          out.write(ctr+" "+cost+"\n");
	//				
	//	          if(rotIndices[p][emeWI.aa1()] == null){
	//	        	  rotIndices[p][emeWI.aa1()] = new int[emat.singles.E[p][emeWI.aa1()].length];
	//	          }
	//	          rotIndices[p][emeWI.aa1()][emeWI.rot1()] = ctr;
	//	          
	//	          rotOut.write(p+" "+ctr+" "+rot.aaType.name+" "+rot.aaIndex+" "+rot.rlIndex+"\n");
	//	          
	//              ctr++;
	////              out.write("ROT::SELF\t"+r.getResNumberString()+"-"+rot.aaType.name+"-"+rotNum+" "+minE+" "+minE+" "+rot.aaIndex+" "+rot.rlIndex+"\n");
	//			    
	//			}
	//        }
	//		
	//        for(int p1=0;p1<emat.numMutPos();p1++){
	//        	for(int p2=p1+1; p2<emat.numMutPos();p2++){
	//        		if(emat.areNeighbors(p1, p2)){
	//		piter = emat.pairsIterator(p1,p2);
	//		out.write("2 "+p1+" "+p2+" "+upperBound+" "+(emat.numRot(p1)*emat.numRot(p2))+"\n");
	//		while(piter.hasNext()){
	//			EMatrixEntryWIndex emeWI = piter.next();
	//			
	//			Residue r1 = m.residue[emat.resByPos.get(emeWI.pos1()).get(0)];
	//			Rotamer rot1 = r1.rl.getRot(emat.singles.getRot(emeWI.rot1index())[0]);
	//			
	//			Residue r2 = m.residue[emat.resByPos.get(emeWI.pos2()).get(0)];
	//			Rotamer rot2 = r2.rl.getRot(emat.singles.getRot(emeWI.rot2index())[0]);
	//			
	//			String my_resNum1 = r1.getResNumberString();
	//			String my_resNum2 = r2.getResNumberString();
	//			String my_resname1 = rot1.aaType.name;
	//			String my_resname2 = rot2.aaType.name;
	//			int my_rotNum1 = emeWI.rot1()+1;
	//			int my_rotNum2 = emeWI.rot2()+1;
	//			
	//			double curEnergy = emeWI.eme.minE(); 
	//			String cost = getCost(minE, curEnergy, upperBound);
	//			//### Print out rot-rot interaction energies
	//			int rotIndex1 = rotIndices[emeWI.pos1()][emeWI.aa1()][emeWI.rot1()];
	//			int rotIndex2 = rotIndices[emeWI.pos2()][emeWI.aa2()][emeWI.rot2()];
	//            //out.write("ROT::ROT\t"+my_resNum1+"-"+my_resname1+"-"+my_rotNum1+"::"+my_resNum2+"-"+my_resname2+"-"+my_rotNum2+" "+curEnergy+" "+curEnergy+"\n");
	//			out.write(rotIndex1+" "+rotIndex2+" "+cost+"\n"); 
	//		
	//			
	//			}}}}
	//		
	//		
	//		     
	//		     out.close();
	//		     rotOut.close();
	//		}catch (Exception e){//Catch exception if any
	//		       System.err.println("Error: " + e.getMessage());
	//		       e.printStackTrace();
	//		     }
	//		
	//		
	//	}
	//	


	/**
	 * 
	 * Given the wcsp solution create the OSPREY Conf file
	 * 
	 * Input: System.cfg MutSearch.cfg Toulbar.sol 
	 * 
	 */
	private double wcspSol2ospreyConf(String sol) {


		//			BufferedReader bufread = null;
		//			String curLine = null;
		//			boolean done = false;

		// First attempt to open and read the config file
		//			try{
		//				FileInputStream is = new FileInputStream("sol");
		//				bufread = new BufferedReader(new InputStreamReader(is));
		//
		//				//Only need to read the first line since that's where the solution is
		//				curLine = bufread.readLine();
		//
		//				bufread.close();
		//			}catch (Exception e) {
		//				e.printStackTrace();
		//			}

		String[] solution = sol.split(" ");

		EMatrixEntryWIndex[] curConf = new EMatrixEntryWIndex[confLength]; 
		Index3 gmecRots[] = new Index3[confLength]; 

		//Find the GMEC
		//			System.out.print("1 ");
		for(int p=0;p<solution.length-1;p++){
			int index = new Integer(solution[p+1]);
			gmecRots[p] = rots[p][index];
			int[] rot = {gmecRots[p].pos,gmecRots[p].aa,gmecRots[p].rot};
			curConf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(gmecRots[p]),rot);
		}


		//Calculate minE
		double minE = 0;//pairwiseMinEnergyMatrix.templ_E;
		for(int p1=0;p1<gmecRots.length;p1++){
			minE += emat.getSingleMinE(gmecRots[p1]);
			for(int p2=p1+1;p2<gmecRots.length;p2++){
				if(emat.areNeighbors(p1, p2))
					minE += emat.getPairMinE(gmecRots[p1], gmecRots[p2]);
			}
		}

		RotConf curRotConf = new RotConf(curConf,minE);

		if(bestConf == null || minE < bestConf.E){
			bestConf = curRotConf;
		}

		allConfs.add(curRotConf);
		//			System.out.println("unMinE: "+minE+" minE: "+minE+" bestE: "+minE);

		return minE;

	}

	public void cleanUp(){
		File f = new File(outDir+File.separator+filename);
		f.delete();
		f = new File(outDir+File.separator+filename+".minE");
		f.delete();
	}


	private String getCost(double E, long upperBound) {
		long Pshift = 100000000L;

		double cost;

		if((E-lowestE) >= upperBound/Pshift){
			cost = upperBound;
		}else{
			cost = E-lowestE;
			cost *= Pshift;

		}
		long costL = (long) cost;

		if(costL < 0)
			System.out.println("Cost should never be 0");

		String s = new Long(costL).toString();
		return s;
	}

	private long getCostForFullConf(double E, long upperBound) {
		long Pshift = 100000000L;

		double cost;

		//		if((E-lowestE) >= upperBound/Pshift){
		//			cost = upperBound;
		//		}else{
		//			cost = E-lowestE;
		//			cost *= Pshift;
		//			
		//		}

		cost = E - (lowestE*numConstraints);
		cost *= Pshift;
		long costL = (long) cost;

		if(costL < 0)
			System.out.println("Cost should never be 0");

		return costL;
	}

	private double getEfromCost(long cost, int numConstrs) {
		double Pshift = 100000000;

		double costD = cost/Pshift;

		double E = costD + (lowestE*numConstrs);


		return E;
	}



	//Given the postprocessed arity 0-1 wcsp, prune the energy matrix based
	//on the initEw value and lowerbound
	public void pruneEmat(String wcspFile, double pruneE){


		long pruneCost = getCostForFullConf(pruneE, upperBound);

		long zeroCost = 0;
		long[][] costs = new long[rots.length][];

		BufferedReader bufread = null;
		String curLine = null;

		// First attempt to open and read the wcsp file
		try{
			FileInputStream is = new FileInputStream(wcspFile);
			bufread = new BufferedReader(new InputStreamReader(is));

			curLine = bufread.readLine(); //Header
			curLine = bufread.readLine(); //Domains
			curLine = bufread.readLine(); //First actual line

			while(curLine != null){
				int arity = new Integer(KSParser.getToken(curLine,1));

				if(arity == 0){ //Add to the template E
					zeroCost = new Long(KSParser.getToken(curLine,2));
				}else if(arity == 1){ //Singles terms
					String[] aLine = curLine.split(" ");
					int pos = new Integer(aLine[1]);
					costs[pos] = new long[rots[pos].length];
					long defaultVal = new Long(aLine[2]);
					for(int i=0; i<costs[pos].length;i++){costs[pos][i] = defaultVal;} //Default is MAX
					int numRots = new Integer(aLine[3]);
					for(int i=0; i<numRots;i++){
						curLine = bufread.readLine();
						int index = new Integer(KSParser.getToken(curLine,1));
						//						if(index!=i){
						//							System.out.println("Something wrong with wcsp loading");
						//						}
						long cost = new Long(KSParser.getToken(curLine,2));
						costs[pos][index] = cost;
					}
				}
				//				else if(arity == 2){ //Pair terms
				//					int p1 = new Integer(getToken(curLine,2));
				//					int p2 = new Integer(getToken(curLine,3));
				//					int numPairs = new Integer(getToken(curLine,5));
				//					for(int i=0; i<numPairs;i++){
				//						curLine = bufread.readLine();
				//						int i1 = new Integer(getToken(curLine,1));
				//						int i2 = new Integer(getToken(curLine,2));
				//						String cost = getToken(curLine,3);
				//						double E = getEfromCost(minE,cost,upperBound);
				//						Index3 r1 = emat.getGlobalRot(emat.resByPos.get(p1).get(0), rots[p1][i1]);
				//						Index3 r2 = emat.getGlobalRot(emat.resByPos.get(p2).get(0), rots[p2][i2]);
				//						
				//						emat.pairs.E[r1.pos][r1.aa][r1.rot][r2.pos][r2.aa][r2.rot] = E;
				//						emat.pairs.E[r2.pos][r2.aa][r2.rot][r1.pos][r1.aa][r1.rot] = E;
				//					}
				//				}


				curLine = bufread.readLine();	
			}

			bufread.close();
		}catch (Exception e) {
			e.printStackTrace();
		}

		pruneCost -= zeroCost;
		int numPruned = 0;
		//Look through the rotamers and prune
		for(int p=0; p<rots.length;p++){
			if(costs[p] == null)
				continue;
			for(int i=0; i<rots[p].length;i++){
				if(costs[p][i] > pruneCost){
					numPruned++;
					Index3 rot = rots[p][i];
					emat.setSinglePruned(rots[p][i],true);
				}
			}
		}

		System.out.println("Pruned "+numPruned+" rotamers with local consistency.");
		emat.removePrunedRotReducedMem(false);

	}

	//Given the postprocessed arity 0-1-2 wcsp, update the energy matrix to match the 
	//costs from the wcsp local consistency enforcement
	public void updateEmat(String wcspFile){

		long zeroCost = 0;
		long[][] costs = new long[rots.length][];

		BufferedReader bufread = null;
		String curLine = null;

		//First set all energies to max and then only reset the ones we read in
		SinglesIterator siter = emat.singlesIterator();
		while(siter.hasNext()){
			EMatrixEntryWIndex emeWI = siter.next();
			emat.setE(emeWI.index, Double.POSITIVE_INFINITY);
		}
		PairsIterator piter = emat.pairsIterator();
		while(piter.hasNext()){
			EMatrixEntryWIndex emeWI = piter.next();
			emat.setE(emeWI.index, Double.POSITIVE_INFINITY);
		}
		
		// First attempt to open and read the wcsp file
		try{
			FileInputStream is = new FileInputStream(wcspFile);
			bufread = new BufferedReader(new InputStreamReader(is));

			curLine = bufread.readLine(); //Header
			String[] aLine = curLine.split(" ");
			int numConstrs = new Integer(aLine[3]);
			curLine = bufread.readLine(); //Domains
			curLine = bufread.readLine(); //First actual line
			double Pshift = 1000000000.0;
			while(curLine != null){
				int arity = new Integer(KSParser.getToken(curLine,1));

				if(arity == 0){ //Add to the template E
					zeroCost = new Long(KSParser.getToken(curLine,2));
					emat.templ_E += getEfromCost(zeroCost, numConstraints);
//					emat.templ_E += zeroCost/100000000.0;
				}else if(arity == 1){ //Singles terms
					aLine = curLine.split(" ");
					int pos = new Integer(aLine[1]);
					costs[pos] = new long[rots[pos].length];
					long defaultVal = new Long(aLine[2]);
					for(int i=0; i<costs[pos].length;i++){costs[pos][i] = defaultVal;} //Default is MAX
					int numRots = new Integer(aLine[3]);
					Index3 r1 = null;
					for(int i=0; i<numRots;i++){
						curLine = bufread.readLine();
						int index = new Integer(KSParser.getToken(curLine,1));
						long cost = new Long(KSParser.getToken(curLine,2));
						double E = getEfromCost(cost, 1);
						r1 = rots[pos][index];
						emat.singles.E[r1.pos][r1.aa][r1.rot] = cost/Pshift;
					}
					if(numRots == 1){ //If there is only one rotamer the pair will not be listed, so we need to make it's pairs viable
						double E = getEfromCost(0, 1);
						for(int p2=0;p2<emat.singles.E.length;p2++){
							if(pos == p2 || !emat.areNeighbors(pos, p2))
								continue;
							for(int a2=0; a2<emat.singles.E[p2].length;a2++){
								for(int r2=0; r2<emat.singles.E[p2][a2].length;r2++){
									emat.pairs.E[r1.pos][r1.aa][r1.rot][p2][a2][r2] = 0;
									emat.pairs.E[p2][a2][r2][r1.pos][r1.aa][r1.rot] = 0;
								}
							}
						}
					}
				}
				else if(arity == 2){ //Pair terms
					aLine = curLine.split(" ");
					int p1 = new Integer(aLine[1]);
					int p2 = new Integer(aLine[2]);
					int numPairs = new Integer(aLine[4]);
					for(int i=0; i<numPairs;i++){
						//do nothing for now
						curLine = bufread.readLine();
						int i1 = new Integer(KSParser.getToken(curLine,1));
						int i2 = new Integer(KSParser.getToken(curLine,2));
						String cost = KSParser.getToken(curLine,3);
						double E = getEfromCost(new Long(cost), 1);
						Index3 r1 = rots[p1][i1];
						Index3 r2 = rots[p2][i2];
						
						emat.pairs.E[r1.pos][r1.aa][r1.rot][r2.pos][r2.aa][r2.rot] = new Long(cost)/Pshift;//E;
						emat.pairs.E[r2.pos][r2.aa][r2.rot][r1.pos][r1.aa][r1.rot] = new Long(cost)/Pshift;//E;
					}
				}


				curLine = bufread.readLine();	
			}

			bufread.close();
		}catch (Exception e) {
			e.printStackTrace();
		}

	}




}
