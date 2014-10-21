import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.StringTokenizer;

import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBExpr;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;


public class GurobiOptimization {

	GRBVar[][][] singleVars;
	GRBVar[][][][][][] pairVars;

	GRBConstr[] singleConstraints;
	GRBConstr[][][][] pairConstraints;

	GRBLinExpr[] singleExprs;
	GRBLinExpr[][][][] pairExprs;

	HashMap<ArrayList<Index3>,TupleConstraint> tupleConstraints;

	final static boolean debug = false;

	GRBEnv    env;
	GRBModel  model;

	GRBConstr badConstr = null;

	int excludedConformations = 0;
	int excludedSeqs = 0;

	public GurobiOptimization(Emat emat,boolean doILP) {
		singleVars = new GRBVar[emat.singles.E.length][][];
		pairVars   = new GRBVar[emat.pairs.E.length][][][][][];

		char varType = GRB.BINARY;
		if (!doILP)
			varType = GRB.CONTINUOUS;

		double[][][][][][] mat = emat.pairs.E;

		for(int p1=0;p1<mat.length;p1++){
			singleVars[p1] = new GRBVar[mat[p1].length][];
			pairVars[p1] = new GRBVar[mat[p1].length][][][][];
			for(int a1=0;a1<mat[p1].length;a1++){
				if(mat[p1][a1]!=null){
					singleVars[p1][a1] = new GRBVar[mat[p1][a1].length];
					pairVars[p1][a1] = new GRBVar[mat[p1][a1].length][][][];
					for(int r1=0;r1<mat[p1][a1].length;r1++){
						pairVars[p1][a1][r1] = new GRBVar[mat[p1][a1][r1].length][][];
						for(int p2=0;p2<mat[p1][a1][r1].length;p2++){
							if(p1!=p2 && mat[p1][a1][r1][p2] != null){
								pairVars[p1][a1][r1][p2] = new GRBVar[mat[p1][a1][r1][p2].length][];
								for(int a2=0;a2<mat[p1][a1][r1][p2].length;a2++){
									if(mat[p1][a1][r1][p2][a2] != null){
										pairVars[p1][a1][r1][p2][a2] = new GRBVar[mat[p1][a1][r1][p2][a2].length];

									}}}}}}}}


		try{

			env   = new GRBEnv("pemMin.log");
			model = new GRBModel(env);



			SinglesIterator s_iter = emat.singlesIterator();
			while(s_iter.hasNext()){
				EMatrixEntryWIndex rot = s_iter.next();
				if(rot.eme.isPruned())
					continue;
				double E = rot.eme.minE();
				if(Math.abs(E) < 0.00)
					E = 0.0;
				else if(E > 1000)
					E = 1000;
				try{
					singleVars[rot.pos1()][rot.aa1()][rot.rot1()] = model.addVar(0.0, 1.0, E, varType, "s_"+rot.pos1()+"_"+rot.aa1()+"_"+rot.rot1());
				}catch(Exception error){
					error.printStackTrace();
				}
			}

			PairsIterator p_iter = emat.pairsIterator(); 
			while(p_iter.hasNext()){
				EMatrixEntryWIndex rot = p_iter.next();
				if(rot.pos1() >= rot.pos2() || rot.eme.isPruned() || emat.getSinglePruned(rot.rot1index()) || emat.getSinglePruned(rot.rot2index()))
					continue;
				double E = rot.eme.minE();
				if(Math.abs(E) < 0.00)
					E = 0.0;
				else if(E > 1000)
					E = 1000;
				pairVars[rot.pos1()][rot.aa1()][rot.rot1()][rot.pos2()][rot.aa2()][rot.rot2()] = model.addVar(0.0, 1.0, E, varType, "p_"+rot.pos1()+"_"+rot.aa1()+"_"+rot.rot1()+"_"+rot.pos2()+"_"+rot.aa2()+"_"+rot.rot2());
				pairVars[rot.pos2()][rot.aa2()][rot.rot2()][rot.pos1()][rot.aa1()][rot.rot1()] =pairVars[rot.pos1()][rot.aa1()][rot.rot1()][rot.pos2()][rot.aa2()][rot.rot2()]; 
			}

			// Integrate new variables
			model.update();

			int constrCtr = 0;
			singleConstraints = new GRBConstr[emat.numMutPos()];
			singleExprs = new GRBLinExpr[emat.numMutPos()];
			for(int p1 =0; p1< emat.numMutPos();p1++){
				SinglesIterator s_pos_iter = emat.singlesIterator(p1);
				GRBLinExpr expr = new GRBLinExpr();
				while(s_pos_iter.hasNext()){
					EMatrixEntryWIndex rot = s_pos_iter.next(); 
					if(rot.eme.isPruned())
						continue;
					expr.addTerm(1.0, singleVars[rot.pos1()][rot.aa1()][rot.rot1()]);
				}
				singleExprs[p1] = expr;
				singleConstraints[p1] = model.addConstr(expr, GRB.EQUAL, 1.0, "s_constr_"+p1);
				constrCtr++;
			}

			int pairconstrCtr = 0;
			pairConstraints = new GRBConstr[mat.length][][][];
			pairExprs = new GRBLinExpr[mat.length][][][];
			for(int p1=0;p1<mat.length;p1++){
				pairConstraints[p1] = new GRBConstr[mat[p1].length][][];
				pairExprs[p1] = new GRBLinExpr[mat[p1].length][][];
				for(int a1=0;a1<mat[p1].length;a1++)
					if(mat[p1][a1]!=null){
						pairConstraints[p1][a1] = new GRBConstr[mat[p1][a1].length][];
						pairExprs[p1][a1] = new GRBLinExpr[mat[p1][a1].length][];
						for(int r1=0;r1<mat[p1][a1].length;r1++){
							pairConstraints[p1][a1][r1] = new GRBConstr[mat[p1][a1][r1].length];
							pairExprs[p1][a1][r1] = new GRBLinExpr[mat[p1][a1][r1].length];
							if(emat.getSinglePruned(p1, a1, r1))
								continue;
							for(int p2=0;p2<mat[p1][a1][r1].length;p2++)
								if(p1!=p2 && mat[p1][a1][r1][p2] != null){
									GRBLinExpr expr = new GRBLinExpr();
									for(int a2=0;a2<mat[p1][a1][r1][p2].length;a2++){
										if(mat[p1][a1][r1][p2][a2] != null)
											for(int r2=0;r2<mat[p1][a1][r1][p2][a2].length;r2++){
												if(emat.getPairPruned(p1,a1,r1,p2,a2,r2) || emat.getSinglePruned(p1, a1, r1) || emat.getSinglePruned(p2, a2, r2))
													continue;
												expr.addTerm(1.0, pairVars[p1][a1][r1][p2][a2][r2]);
											}
									}

									expr.addTerm(-1.0, singleVars[p1][a1][r1]);
									pairExprs[p1][a1][r1][p2] = expr;
									pairConstraints[p1][a1][r1][p2] = model.addConstr(expr, GRB.EQUAL, 0.0, "p_constr"+p1+"_"+a1+"_"+r1+"_"+p2);
									pairconstrCtr++;
								}}}}

			model.update();
			//			model.write("gurobiModel.lp");



		}catch(Exception E){
			System.out.println("Something went wrong");
			E.printStackTrace();
		}

	}


	public GurobiOptimization(PGQueueNode node, Emat emat, int[] numNodesForLevel, Index3[][] twoDTo3D, int numTotalNodes, int numThreads){
		singleVars = new GRBVar[numNodesForLevel.length][1][];
		pairVars   = new GRBVar[numNodesForLevel.length][1][][][][];



		for(int p1=0; p1<node.confSoFar.length;p1++){
			if(node.confSoFar[p1] == -1){//Empty
				singleVars[p1][0] = new GRBVar[numNodesForLevel[p1]];
				pairVars[p1][0] = new GRBVar[numNodesForLevel[p1]][][][];
			}else{	
				singleVars[p1][0] = new GRBVar[1];
				pairVars[p1][0] = new GRBVar[1][][][];
			}
			for(int r1=0;r1<singleVars[p1][0].length;r1++){
				pairVars[p1][0][r1] = new GRBVar[numNodesForLevel.length][1][];
				for(int p2=0;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						if(node.confSoFar[p2] == -1){
							pairVars[p1][0][r1][p2][0] = new GRBVar[numNodesForLevel[p2]];
						}else{
							pairVars[p1][0][r1][p2][0] = new GRBVar[1];
						}
					}}}}


		try{

			env   = new GRBEnv();
			env.set( GRB.IntParam.OutputFlag, 0 );
			env.set( GRB.IntParam.Threads, numThreads);
			model = new GRBModel(env);


			for(int p1=0; p1<node.confSoFar.length;p1++){
				if(node.confSoFar[p1] >= 0){ //Is not empty
					Index3 index1 = twoDTo3D[p1][node.confSoFar[p1]];
					double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
					singleVars[p1][0][0] = model.addVar(0.0, 1.0, minIndVoxE, GRB.CONTINUOUS, "s_"+p1+"_0_0");

				}
				else{
					for(int r1=0;r1<singleVars[p1][0].length;r1++){
						Index3 index1 = twoDTo3D[p1][r1];
						double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						singleVars[p1][0][r1] = model.addVar(0.0, 1.0, minIndVoxE, GRB.CONTINUOUS, "s_"+p1+"_0_"+r1);
					}
				}
			}


			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					for(int p2=p1+1;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){
							for(int r2=0;r2<singleVars[p2][0].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][node.confSoFar[p1]];
								}else{
									index1 = twoDTo3D[p1][r1];
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][node.confSoFar[p2]];
								}else{
									index2 = twoDTo3D[p2][r2];
								}
								if(!emat.areNeighbors(index1.pos, index2.pos) || emat.getPairPruned(index1, index2)) //Will be null if not neighbors
									continue;

								double E = emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
								pairVars[p1][0][r1][p2][0][r2] = model.addVar(0.0, 1.0, E, GRB.CONTINUOUS, "p_"+p1+"_0_"+r1+"_"+p2+"_0_"+r2);
								pairVars[p2][0][r2][p1][0][r1] =pairVars[p1][0][r1][p2][0][r2];	
							}
						}
					}
				}
			}

			// Integrate new variables
			model.update();

			int constrCtr = 0;
			singleConstraints = new GRBConstr[node.confSoFar.length];
			singleExprs = new GRBLinExpr[node.confSoFar.length];
			for(int p1 =0; p1< node.confSoFar.length;p1++){

				GRBLinExpr expr = new GRBLinExpr();
				for(int r1 = 0; r1<singleVars[p1][0].length;r1++){
					expr.addTerm(1.0, singleVars[p1][0][r1]);
				}
				singleExprs[p1] = expr;
				singleConstraints[p1] = model.addConstr(expr, GRB.EQUAL, 1.0, "s_constr_"+p1);
				constrCtr++;
			}

			int pairconstrCtr = 0;
			pairConstraints = new GRBConstr[node.confSoFar.length][1][][];
			pairExprs = new GRBLinExpr[node.confSoFar.length][1][][];
			for(int p1=0;p1<node.confSoFar.length;p1++){
				pairConstraints[p1][0] = new GRBConstr[singleVars[p1][0].length][];
				pairExprs[p1][0] = new GRBLinExpr[singleVars[p1][0].length][];
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					pairConstraints[p1][0][r1] = new GRBConstr[node.confSoFar.length];
					pairExprs[p1][0][r1] = new GRBLinExpr[node.confSoFar.length];
					for(int p2=0;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){ //need to check they are neighbors
							GRBLinExpr expr = new GRBLinExpr();
							for(int r2=0;r2<singleVars[p2][0].length;r2++){
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][node.confSoFar[p1]];//nodeIndexOffset[p1]+node.confSoFar[p1];
								}else{
									index1 = twoDTo3D[p1][r1];//nodeIndexOffset[p1]+r1;
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][node.confSoFar[p2]];//nodeIndexOffset[p2]+node.confSoFar[p2];
								}else{
									index2 = twoDTo3D[p2][r2];//nodeIndexOffset[p2]+r2;
								}
								if(emat.getPairPruned(index1, index2))//[index1][index2].eme.isPruned()) //If they aren't neighbors don't include
									continue;

								expr.addTerm(1.0, pairVars[p1][0][r1][p2][0][r2]);
							}

							expr.addTerm(-1.0, singleVars[p1][0][r1]);
							pairExprs[p1][0][r1][p2] = expr;
							pairConstraints[p1][0][r1][p2] = model.addConstr(expr, GRB.EQUAL, 0.0, "p_constr"+p1+"_0_"+r1+"_"+p2);
							pairconstrCtr++;
						}
					}}}

			model.update();
			//model.write("gurobiModel.lp");



		}catch(Exception E){
			System.out.println("Something went wrong");
			E.printStackTrace();
		}


	}

	public GurobiOptimization(PGQueueNode node, Emat emat, int[] numNodesForLevel, 
			Index3[][] twoDTo3D, int[][] numRotRemainingBySeq, int[][] seqIndexOffset, int numTotalNodes,boolean doInteger){
		singleVars = new GRBVar[numNodesForLevel.length][1][];
		pairVars   = new GRBVar[numNodesForLevel.length][1][][][][];


		char varMode;

		if(doInteger)
			varMode = GRB.BINARY;
		else
			varMode = GRB.CONTINUOUS;

		for(int p1=0; p1<node.confSoFar.length;p1++){
			if(node.confSoFar[p1] == -1){//Empty
				singleVars[p1][0] = new GRBVar[numNodesForLevel[p1]];
				pairVars[p1][0] = new GRBVar[numNodesForLevel[p1]][][][];
			}else{	
				singleVars[p1][0] = new GRBVar[numRotRemainingBySeq[p1][node.confSoFar[p1]]];
				pairVars[p1][0] = new GRBVar[numRotRemainingBySeq[p1][node.confSoFar[p1]]][][][];
			}
			for(int r1=0;r1<singleVars[p1][0].length;r1++){
				pairVars[p1][0][r1] = new GRBVar[numNodesForLevel.length][1][];
				for(int p2=0;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						if(node.confSoFar[p2] == -1){
							pairVars[p1][0][r1][p2][0] = new GRBVar[numNodesForLevel[p2]];
						}else{
							pairVars[p1][0][r1][p2][0] = new GRBVar[numRotRemainingBySeq[p2][node.confSoFar[p2]]];
						}
					}}}}


		try{

			env   = new GRBEnv();
			env.set( GRB.IntParam.OutputFlag, 0 );
			model = new GRBModel(env);


			for(int p1=0; p1<node.confSoFar.length;p1++){
				if(node.confSoFar[p1] >= 0){ //Is not empty
					for(int r1=0; r1<singleVars[p1][0].length;r1++){
						Index3 index1 = twoDTo3D[p1][seqIndexOffset[p1][node.confSoFar[p1]] + r1];
						double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						singleVars[p1][0][r1] = model.addVar(0.0, 1.0, minIndVoxE, varMode, "s_"+p1+"_0_"+r1);
					}

				}
				else{
					for(int r1=0;r1<singleVars[p1][0].length;r1++){
						Index3 index1 = twoDTo3D[p1][r1];
						double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						singleVars[p1][0][r1] = model.addVar(0.0, 1.0, minIndVoxE, varMode, "s_"+p1+"_0_"+r1);
					}
				}
			}


			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					for(int p2=p1+1;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){
							for(int r2=0;r2<singleVars[p2][0].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][seqIndexOffset[p1][node.confSoFar[p1]]+r1];
								}else{
									index1 = twoDTo3D[p1][r1];
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][seqIndexOffset[p2][node.confSoFar[p2]]+r2];
								}else{
									index2 = twoDTo3D[p2][r2];
								}
								if(!emat.areNeighbors(p1, p2) || emat.getPairPruned(index1, index2)) //Will be null if not neighbors
									continue;

								double E = emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
								pairVars[p1][0][r1][p2][0][r2] = model.addVar(0.0, 1.0, E, varMode, "p_"+p1+"_0_"+r1+"_"+p2+"_0_"+r2);
								pairVars[p2][0][r2][p1][0][r1] =pairVars[p1][0][r1][p2][0][r2];	
							}
						}
					}
				}
			}

			// Integrate new variables
			model.update();

			int constrCtr = 0;
			singleConstraints = new GRBConstr[node.confSoFar.length];
			singleExprs = new GRBLinExpr[node.confSoFar.length];
			for(int p1 =0; p1< node.confSoFar.length;p1++){

				GRBLinExpr expr = new GRBLinExpr();
				for(int r1 = 0; r1<singleVars[p1][0].length;r1++){
					expr.addTerm(1.0, singleVars[p1][0][r1]);
				}
				singleExprs[p1] = expr;
				singleConstraints[p1] = model.addConstr(expr, GRB.EQUAL, 1.0, "s_constr_"+p1);
				constrCtr++;
			}

			int pairconstrCtr = 0;
			pairConstraints = new GRBConstr[node.confSoFar.length][1][][];
			pairExprs = new GRBLinExpr[node.confSoFar.length][1][][];
			for(int p1=0;p1<node.confSoFar.length;p1++){
				pairConstraints[p1][0] = new GRBConstr[singleVars[p1][0].length][];
				pairExprs[p1][0] = new GRBLinExpr[singleVars[p1][0].length][];
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					pairConstraints[p1][0][r1] = new GRBConstr[node.confSoFar.length];
					pairExprs[p1][0][r1] = new GRBLinExpr[node.confSoFar.length];
					for(int p2=0;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){ //need to check they are neighbors
							GRBLinExpr expr = new GRBLinExpr();
							for(int r2=0;r2<singleVars[p2][0].length;r2++){
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][seqIndexOffset[p1][node.confSoFar[p1]]+r1];
								}else{
									index1 = twoDTo3D[p1][r1];
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][seqIndexOffset[p2][node.confSoFar[p2]]+r2];
								}else{
									index2 = twoDTo3D[p2][r2];
								}
								
								if(!emat.areNeighbors(p1, p2) || emat.getPairPruned(index1, index2)) //If they aren't neighbors don't include
									continue;
								

								expr.addTerm(1.0, pairVars[p1][0][r1][p2][0][r2]);
							}

							expr.addTerm(-1.0, singleVars[p1][0][r1]);
							pairExprs[p1][0][r1][p2] = expr;
							pairConstraints[p1][0][r1][p2] = model.addConstr(expr, GRB.EQUAL, 0.0, "p_constr"+p1+"_0_"+r1+"_"+p2);
							pairconstrCtr++;
						}
					}}}

			model.update();
			//model.write("gurobiModel.lp");



		}catch(Exception E){
			System.out.println("Something went wrong");
			E.printStackTrace();
		}


	}


	public GurobiOptimization(PGQueueNode node, Emat emat, int[] numNodesForLevel, Index3[][] twoDTo3D, int numTotalNodes,boolean min){
		singleVars = new GRBVar[numNodesForLevel.length][1][];
		pairVars   = new GRBVar[numNodesForLevel.length][1][][][][];



		for(int p1=0; p1<node.confSoFar.length;p1++){
			if(node.confSoFar[p1] == -1){//Empty
				singleVars[p1][0] = new GRBVar[numNodesForLevel[p1]];
				pairVars[p1][0] = new GRBVar[numNodesForLevel[p1]][][][];
			}else{	
				singleVars[p1][0] = new GRBVar[1];
				pairVars[p1][0] = new GRBVar[1][][][];
			}
			for(int r1=0;r1<singleVars[p1][0].length;r1++){
				pairVars[p1][0][r1] = new GRBVar[numNodesForLevel.length][1][];
				for(int p2=0;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						if(node.confSoFar[p2] == -1){
							pairVars[p1][0][r1][p2][0] = new GRBVar[numNodesForLevel[p2]];
						}else{
							pairVars[p1][0][r1][p2][0] = new GRBVar[1];
						}
					}}}}


		try{

			env   = new GRBEnv();
			env.set( GRB.IntParam.OutputFlag, 0 );
			model = new GRBModel(env);


			for(int p1=0; p1<node.confSoFar.length;p1++){
				if(node.confSoFar[p1] >= 0){ //Is not empty
					//int index1 = nodeIndexOffset[p1]+node.confSoFar[p1];
					//double minIndVoxE = pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
					//KER: put 0 here since we are now minimizing the full terms
					singleVars[p1][0][0] = model.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS, "s_"+p1+"_0_0");

				}
				else{
					for(int r1=0;r1<singleVars[p1][0].length;r1++){
						Index3 index1 = twoDTo3D[p1][r1];//nodeIndexOffset[p1]+r1;
						double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						singleVars[p1][0][r1] = model.addVar(0.0, 1.0, minIndVoxE, GRB.CONTINUOUS, "s_"+p1+"_0_"+r1);
					}
				}
			}


			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					for(int p2=p1+1;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){
							for(int r2=0;r2<singleVars[p2][0].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 index1,index2;
								boolean valid1 = false;
								boolean valid2 = false;
								if(node.confSoFar[p1] >= 0){
									valid1 = true;
									index1 = twoDTo3D[p1][node.confSoFar[p1]];//nodeIndexOffset[p1]+node.confSoFar[p1];
								}else{
									index1 = twoDTo3D[p1][r1];//nodeIndexOffset[p1]+r1;
								}
								if(node.confSoFar[p2] >= 0){
									valid2 = true;
									index2 = twoDTo3D[p2][node.confSoFar[p2]];//nodeIndexOffset[p2]+node.confSoFar[p2];
								}else{
									index2 = twoDTo3D[p2][r2];//nodeIndexOffset[p2]+r2;
								}
								if(!emat.areNeighbors(index1.pos, index2.pos) || emat.getPairPruned(index1, index2)) //Will be null if not neighbors
									continue;

								double E = 0.0;
								if(!valid1 || !valid2) //If both valid we leave energy as 0.0
									E = emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
								pairVars[p1][0][r1][p2][0][r2] = model.addVar(0.0, 1.0, E, GRB.CONTINUOUS, "p_"+p1+"_0_"+r1+"_"+p2+"_0_"+r2);
								pairVars[p2][0][r2][p1][0][r1] =pairVars[p1][0][r1][p2][0][r2];	
							}
						}
					}
				}
			}

			// Integrate new variables
			model.update();

			int constrCtr = 0;
			singleConstraints = new GRBConstr[node.confSoFar.length];
			singleExprs = new GRBLinExpr[node.confSoFar.length];
			for(int p1 =0; p1< node.confSoFar.length;p1++){

				GRBLinExpr expr = new GRBLinExpr();
				for(int r1 = 0; r1<singleVars[p1][0].length;r1++){
					expr.addTerm(1.0, singleVars[p1][0][r1]);
				}
				singleExprs[p1] = expr;
				singleConstraints[p1] = model.addConstr(expr, GRB.EQUAL, 1.0, "s_constr_"+p1);
				constrCtr++;
			}

			int pairconstrCtr = 0;
			pairConstraints = new GRBConstr[node.confSoFar.length][1][][];
			pairExprs = new GRBLinExpr[node.confSoFar.length][1][][];
			for(int p1=0;p1<node.confSoFar.length;p1++){
				pairConstraints[p1][0] = new GRBConstr[singleVars[p1][0].length][];
				pairExprs[p1][0] = new GRBLinExpr[singleVars[p1][0].length][];
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					pairConstraints[p1][0][r1] = new GRBConstr[node.confSoFar.length];
					pairExprs[p1][0][r1] = new GRBLinExpr[node.confSoFar.length];
					for(int p2=0;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){ //need to check they are neighbors
							GRBLinExpr expr = new GRBLinExpr();
							for(int r2=0;r2<singleVars[p2][0].length;r2++){
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][node.confSoFar[p1]];//nodeIndexOffset[p1]+node.confSoFar[p1];
								}else{
									index1 = twoDTo3D[p1][r1];//nodeIndexOffset[p1]+r1;
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][node.confSoFar[p2]];//nodeIndexOffset[p2]+node.confSoFar[p2];
								}else{
									index2 = twoDTo3D[p2][r2];//nodeIndexOffset[p2]+r2;
								}
								if(emat.getPairPruned(index1, index2))//[index1][index2].eme.isPruned()) //If they aren't neighbors don't include
									continue;

								expr.addTerm(1.0, pairVars[p1][0][r1][p2][0][r2]);
							}

							expr.addTerm(-1.0, singleVars[p1][0][r1]);
							pairExprs[p1][0][r1][p2] = expr;
							pairConstraints[p1][0][r1][p2] = model.addConstr(expr, GRB.EQUAL, 0.0, "p_constr"+p1+"_0_"+r1+"_"+p2);
							pairconstrCtr++;
						}
					}}}

			model.update();
			//model.write("gurobiModel.lp");



		}catch(Exception E){
			System.out.println("Something went wrong");
			E.printStackTrace();
		}


	}

	//For tuples calculation, removes g score from the LP step
	public GurobiOptimization(boolean tuples, PGQueueNode node, Emat emat, int[] numNodesForLevel, Index3[][] twoDTo3D, int numTotalNodes,int verbosity){
		singleVars = new GRBVar[numNodesForLevel.length][1][];
		pairVars   = new GRBVar[numNodesForLevel.length][1][][][][];



		for(int p1=0; p1<node.confSoFar.length;p1++){
			if(node.confSoFar[p1] == -1){//Empty
				singleVars[p1][0] = new GRBVar[numNodesForLevel[p1]];
				pairVars[p1][0] = new GRBVar[numNodesForLevel[p1]][][][];
			}else{	
				singleVars[p1][0] = new GRBVar[1];
				pairVars[p1][0] = new GRBVar[1][][][];
			}
			for(int r1=0;r1<singleVars[p1][0].length;r1++){
				pairVars[p1][0][r1] = new GRBVar[numNodesForLevel.length][1][];
				for(int p2=0;p2<node.confSoFar.length;p2++){
					if(p1!=p2 && emat.areNeighbors(p1, p2)){
						if(node.confSoFar[p2] == -1){
							pairVars[p1][0][r1][p2][0] = new GRBVar[numNodesForLevel[p2]];
						}else{
							pairVars[p1][0][r1][p2][0] = new GRBVar[1];
						}
					}}}}


		try{

			env   = new GRBEnv();
			env.set( GRB.IntParam.OutputFlag, verbosity );
			model = new GRBModel(env);


			for(int p1=0; p1<node.confSoFar.length;p1++){
				if(node.confSoFar[p1] >= 0){ //Is not empty
					Index3 index1 = twoDTo3D[p1][node.confSoFar[p1]];
					double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
					singleVars[p1][0][0] = model.addVar(0.0, 1.0, 0, GRB.CONTINUOUS, "s_"+p1+"_0_0");

				}
				else{
					for(int r1=0;r1<singleVars[p1][0].length;r1++){
						Index3 index1 = twoDTo3D[p1][r1];
						double minIndVoxE = emat.getSingleMinE(index1);//pairwiseMinEnergyMatrix[index1][numTotalNodes].eme.minE();
						singleVars[p1][0][r1] = model.addVar(0.0, 1.0, minIndVoxE, GRB.CONTINUOUS, "s_"+p1+"_0_"+r1);
					}
				}
			}


			for(int p1=0; p1<node.confSoFar.length;p1++){
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					for(int p2=p1+1;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){
							for(int r2=0;r2<singleVars[p2][0].length;r2++){ //The options here is that either p1 or p2 is already defined and I need to use confSoFar
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][node.confSoFar[p1]];
								}else{
									index1 = twoDTo3D[p1][r1];
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][node.confSoFar[p2]];
								}else{
									index2 = twoDTo3D[p2][r2];
								}
								if(!emat.areNeighbors(index1.pos, index2.pos) || emat.getPairPruned(index1, index2)) //Will be null if not neighbors
									continue;

								double E;
								if(node.confSoFar[p1] >= 0 && node.confSoFar[p2] >= 0)
									E = 0;
								else
									E = emat.getPairMinE(index1, index2);//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
								pairVars[p1][0][r1][p2][0][r2] = model.addVar(0.0, 1.0, E, GRB.CONTINUOUS, "p_"+p1+"_0_"+r1+"_"+p2+"_0_"+r2);
								pairVars[p2][0][r2][p1][0][r1] =pairVars[p1][0][r1][p2][0][r2];	
							}
						}
					}
				}
			}

			// Integrate new variables
			model.update();

			int constrCtr = 0;
			singleConstraints = new GRBConstr[node.confSoFar.length];
			singleExprs = new GRBLinExpr[node.confSoFar.length];
			for(int p1 =0; p1< node.confSoFar.length;p1++){

				GRBLinExpr expr = new GRBLinExpr();
				for(int r1 = 0; r1<singleVars[p1][0].length;r1++){
					expr.addTerm(1.0, singleVars[p1][0][r1]);
				}
				singleExprs[p1] = expr;
				singleConstraints[p1] = model.addConstr(expr, GRB.EQUAL, 1.0, "s_constr_"+p1);
				constrCtr++;
			}

			int pairconstrCtr = 0;
			pairConstraints = new GRBConstr[node.confSoFar.length][1][][];
			pairExprs = new GRBLinExpr[node.confSoFar.length][1][][];
			for(int p1=0;p1<node.confSoFar.length;p1++){
				pairConstraints[p1][0] = new GRBConstr[singleVars[p1][0].length][];
				pairExprs[p1][0] = new GRBLinExpr[singleVars[p1][0].length][];
				for(int r1=0;r1<singleVars[p1][0].length;r1++){
					pairConstraints[p1][0][r1] = new GRBConstr[node.confSoFar.length];
					pairExprs[p1][0][r1] = new GRBLinExpr[node.confSoFar.length];
					for(int p2=0;p2<node.confSoFar.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1, p2)){ //need to check they are neighbors
							GRBLinExpr expr = new GRBLinExpr();
							for(int r2=0;r2<singleVars[p2][0].length;r2++){
								Index3 index1,index2;
								if(node.confSoFar[p1] >= 0){
									index1 = twoDTo3D[p1][node.confSoFar[p1]];//nodeIndexOffset[p1]+node.confSoFar[p1];
								}else{
									index1 = twoDTo3D[p1][r1];//nodeIndexOffset[p1]+r1;
								}
								if(node.confSoFar[p2] >= 0){
									index2 = twoDTo3D[p2][node.confSoFar[p2]];//nodeIndexOffset[p2]+node.confSoFar[p2];
								}else{
									index2 = twoDTo3D[p2][r2];//nodeIndexOffset[p2]+r2;
								}
								if(emat.getPairPruned(index1, index2))//[index1][index2].eme.isPruned()) //If they aren't neighbors don't include
									continue;

								expr.addTerm(1.0, pairVars[p1][0][r1][p2][0][r2]);
							}

							expr.addTerm(-1.0, singleVars[p1][0][r1]);
							pairExprs[p1][0][r1][p2] = expr;
							pairConstraints[p1][0][r1][p2] = model.addConstr(expr, GRB.EQUAL, 0.0, "p_constr"+p1+"_0_"+r1+"_"+p2);
							pairconstrCtr++;
						}
					}}}

			model.update();
			//model.write("gurobiModel.lp");



		}catch(Exception E){
			System.out.println("Something went wrong");
			E.printStackTrace();
		}


	}


	public double optimize(){

		try{

			model.optimize();

//						GRBVar[] vars = model.getVars();
//						for(int i=0; i<vars.length;i++){
//							if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){
//			
////								if(vars[i].get(GRB.StringAttr.VarName).startsWith("tup"))
//									System.out.println(vars[i].get(GRB.StringAttr.VarName)
//											+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));
//			
			//					String varName = vars[i].get(GRB.StringAttr.VarName);
			//					if(varName.startsWith("s")){
			//						StringTokenizer st = new StringTokenizer(varName,"_");
			//						st.nextElement(); //remove the "s" at the beginning
			//						int p = new Integer(st.nextToken()).intValue();
			//						int a = new Integer(st.nextToken()).intValue();
			//						int r = new Integer(st.nextToken()).intValue();
			//						int[] index = {p,a,r};
			//						conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
			//					}
			//					else if(varName.startsWith("tup_s")){
			//						StringTokenizer st = new StringTokenizer(varName,"_");
			//						st.nextElement(); //remove the "tup" at the beginning
			//						st.nextElement(); //remove the "s"
			//						ArrayList<Index3> rots = new ArrayList<Index3>();
			//						while(st.hasMoreTokens()){
			//							int p = new Integer(st.nextToken()).intValue();
			//							int a = new Integer(st.nextToken()).intValue();
			//							int r = new Integer(st.nextToken()).intValue();
			//							int[] index = {p,a,r};
			//							rots.add(new Index3(p,a,r));
			//							conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
			//						}
			//						confTuples.add(energyTuples.get(rots));
			//					}
			
//							}
//						}


			return model.get(GRB.DoubleAttr.ObjVal);
		}catch(Exception E){
			//			try {
			//							model.computeIIS();
			//							model.write("IIS.ilp");
			//						} catch (GRBException e) {
			//							e.printStackTrace();
			//						}
			//			E.printStackTrace();
			//KER: Presumably the model was infeasible because all possible pairs for two positions were pruned
			return Double.POSITIVE_INFINITY;
		}


	}

	public double getObjVal(){
		try{
			return model.get(GRB.DoubleAttr.ObjVal);
		}catch(Exception E){
			return Double.POSITIVE_INFINITY;
		}
	}

	public void printParams(){
		GRBVar[] vars = model.getVars();
		for(int i=0; i<vars.length;i++){
			try {
				if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){

					System.out.println(vars[i].get(GRB.StringAttr.VarName)
							+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));

					//				String varName = vars[i].get(GRB.StringAttr.VarName);
					//				if(varName.startsWith("s")){
					//					StringTokenizer st = new StringTokenizer(varName,"_");
					//					st.nextElement(); //remove the "s" at the beginning
					//					int p = new Integer(st.nextToken()).intValue();
					//					int a = new Integer(st.nextToken()).intValue();
					//					int r = new Integer(st.nextToken()).intValue();
					//					int[] index = {p,a,r};
					//					conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
					//				}
					//				else if(varName.startsWith("tup_s")){
					//					StringTokenizer st = new StringTokenizer(varName,"_");
					//					st.nextElement(); //remove the "tup" at the beginning
					//					st.nextElement(); //remove the "s"
					//					ArrayList<Index3> rots = new ArrayList<Index3>();
					//					while(st.hasMoreTokens()){
					//						int p = new Integer(st.nextToken()).intValue();
					//						int a = new Integer(st.nextToken()).intValue();
					//						int r = new Integer(st.nextToken()).intValue();
					//						int[] index = {p,a,r};
					//						rots.add(new Index3(p,a,r));
					//						conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
					//					}
					//					confTuples.add(energyTuples.get(rots));
					//				}

				}
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}


	public GurobiConf getMostProbable(Emat emat){
		EMatrixEntryWIndex[] conf = new EMatrixEntryWIndex[singleVars.length];
		double[] confCoeff = new double[singleVars.length];
		LinkedList<EnergyTuple> confTuples = new LinkedList<EnergyTuple>();

		try{

			GRBVar[] vars = model.getVars();
			for(int i=0; i<vars.length;i++){
				if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){

					//System.out.println(vars[i].get(GRB.StringAttr.VarName)
					//		+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));

					String varName = vars[i].get(GRB.StringAttr.VarName);
					if(varName.startsWith("s")){
						StringTokenizer st = new StringTokenizer(varName,"_");
						st.nextElement(); //remove the "s" at the beginning
						int p = new Integer(st.nextToken()).intValue();
						int a = new Integer(st.nextToken()).intValue();
						int r = new Integer(st.nextToken()).intValue();
						int[] index = {p,a,r};
						if(confCoeff[p] < vars[i].get(GRB.DoubleAttr.X)){
							conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
							confCoeff[p] = vars[i].get(GRB.DoubleAttr.X);
						}
					}
				}
			}
		} catch(Exception E){
			System.out.println("Something went wrong, presumably we're out of conformations");
			return new GurobiConf(null,null);
		}

		return new GurobiConf(conf, confTuples);
	}


	public GurobiConf optimize(Emat emat, HashMap<ArrayList<Index3>,EnergyTuple> energyTuples) {

		EMatrixEntryWIndex[] conf = new EMatrixEntryWIndex[singleVars.length];
		LinkedList<EnergyTuple> confTuples = new LinkedList<EnergyTuple>();

		//Emat emat = new Emat("PEMmin_COM.dat");
		try{
			//model.write("gurobiModel.lp");
			model.optimize();

			//double myTotal = emat.templ_E;
			//			double myTotal = 0;
			//			for(int p1=0;p1<singleVars.length;p1++)
			//				for(int a1=0;a1<singleVars[p1].length;a1++)
			//					if(singleVars[p1][a1]!=null)
			//						for(int r1=0;r1<singleVars[p1][a1].length;r1++){
			//							if(singleVars[p1][a1][r1].get(GRB.DoubleAttr.X) != 0.0){
			//								//System.out.print(singleVars[p1][a1][r1].get(GRB.StringAttr.VarName)
			//								//		+ " " +singleVars[p1][a1][r1].get(GRB.DoubleAttr.X));
			//								//System.out.println(" "+emat.singles.E[p1][a1][r1]);
			//								myTotal += emat.singles.E[p1][a1][r1];
			//								int[] index = {p1,a1,r1};
			//								conf[p1] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
			//								//if(singleVars[p1][a1][r1].get(GRB.DoubleAttr.X) < 1){
			//								//	System.out.println("Non-Integer");
			//								//}
			//							}
			//							/*for(int p2=p1+1;p2<pairVars[p1][a1][r1].length;p2++)
			//								if(p1!=p2 && pairVars[p1][a1][r1][p2] != null){
			//									for(int a2=0;a2<pairVars[p1][a1][r1][p2].length;a2++){
			//										if(pairVars[p1][a1][r1][p2][a2] != null)
			//											for(int r2=0;r2<pairVars[p1][a1][r1][p2][a2].length;r2++){
			//												try{
			//													if(pairVars[p1][a1][r1][p2][a2][r2] != null && pairVars[p1][a1][r1][p2][a2][r2].get(GRB.DoubleAttr.X) != 0.0){
			//													System.out.print(pairVars[p1][a1][r1][p2][a2][r2].get(GRB.StringAttr.VarName)
			//															+ " " +pairVars[p1][a1][r1][p2][a2][r2].get(GRB.DoubleAttr.X));
			//													System.out.println(" "+emat.pairs.E[p1][a1][r1][p2][a2][r2]);
			//													myTotal += emat.pairs.E[p1][a1][r1][p2][a2][r2];
			//													}
			//													}catch(Exception E){
			//														System.out.println("WHY!?!");
			//													}
			//											}
			//									}
			//								}*/
			//						}
			//			
			//			
			//
			//			//System.out.println("GurobiVars"); 

			GRBVar[] vars = model.getVars();
			for(int i=0; i<vars.length;i++){
				if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){

					if(debug)
					System.out.println(vars[i].get(GRB.StringAttr.VarName)
							+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));

					String varName = vars[i].get(GRB.StringAttr.VarName);
					if(varName.startsWith("s")){
						StringTokenizer st = new StringTokenizer(varName,"_");
						st.nextElement(); //remove the "s" at the beginning
						int p = new Integer(st.nextToken()).intValue();
						int a = new Integer(st.nextToken()).intValue();
						int r = new Integer(st.nextToken()).intValue();
						int[] index = {p,a,r};
						conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
					}
					else if(varName.startsWith("tup_s")){
						StringTokenizer st = new StringTokenizer(varName,"_");
						st.nextElement(); //remove the "tup" at the beginning
						st.nextElement(); //remove the "s"
						ArrayList<Index3> rots = new ArrayList<Index3>();
						while(st.hasMoreTokens()){
							int p = new Integer(st.nextToken()).intValue();
							int a = new Integer(st.nextToken()).intValue();
							int r = new Integer(st.nextToken()).intValue();
							int[] index = {p,a,r};
							rots.add(new Index3(p,a,r));
							conf[p] = new EMatrixEntryWIndex(emat.singles.getTerm(index),index);
						}
						confTuples.add(energyTuples.get(rots));
					}

				}
			}

			//System.out.println("TotalE: "+myTotal);
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

		}catch(Exception E){
			//			try {
			//				model.computeIIS();
			//				model.write("IIS.ilp");
			//			} catch (GRBException e) {
			//				e.printStackTrace();
			//			}
			//			System.out.println("Something went wrong");
			//			E.printStackTrace();
			//			return null;
			System.out.println("Something went wrong, presumably we're out of conformations");
			return new GurobiConf(null,null);
		}

		return new GurobiConf(conf, confTuples);

	}

	public void removeConf(EMatrixEntryWIndex[] conf){
		GRBLinExpr expr = new GRBLinExpr();


		for(int i=0; i<conf.length;i++){
			expr.addTerm(1.0, singleVars[conf[i].pos1()][conf[i].aa1()][conf[i].rot1()]);
		}
		try {
			model.addConstr(expr, GRB.LESS_EQUAL, conf.length-1, "exclude_constr"+excludedConformations);
			model.update();
		} catch (GRBException e) {

			e.printStackTrace();
		}

		excludedConformations++;
	}
//
//	public void addTuple(EnergyTuple child, ArrayList<EnergyTuple> parents) {
//		if(tupleConstraints == null){
//			tupleConstraints = new ArrayList<TupleConstraint>();
//		}
//
//		//KER: We need to add variables for the tuple and all of the tuple pairs
//
//		try{
//			//KER: First add the tuple variable
//			String tupPrefSingle = "tup_s";
//			String tupPrefPair = "tup_p";
//			String tupPrefPairConstr = "tup_p_constr";
//			for(Index3 rot: child.rots){
//				tupPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//				tupPrefPair += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//				tupPrefPairConstr += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//			}
//			double E = child.intraE;
//			if(Math.abs(E) < 0.00)
//				E = 0.0;
//			else if(E > 1000)
//				E = 1000;
//
//			GRBVar tupVar = model.addVar(0.0, 1.0, E, GRB.BINARY, tupPrefSingle);
//
//			//KER: now add all of the tupPair variables
//			GRBVar[][][] tupPairVars = new GRBVar[singleVars.length][][];
//			int posCtr=-1;
//			//int posOffset = 0;
//			for(int p=0; p<child.E.length;p++){
//				posCtr++;
//				if(child.E[p] != null){
//					while(posInTuple(posCtr,child.rots)){
//						posCtr++;
//					}
//					tupPairVars[posCtr] = new GRBVar[child.E[p].length][];
//					for(int a=0; a<child.E[p].length;a++){
//						if(child.E[p][a] != null){
//							tupPairVars[posCtr][a] = new GRBVar[child.E[p][a].length];
//							for(int r=0; r<child.E[p][a].length;r++){
//								E = child.E[p][a][r];
//								if(Math.abs(E) < 0.00)
//									E = 0.0;
//								else if(E > 1000)
//									E = 1000;
//								tupPairVars[posCtr][a][r] = model.addVar(0.0, 1.0, E, GRB.BINARY, tupPrefPair+"_"+posCtr+"_"+a+"_"+r);
//							}
//						}}}
//			}
//
//			model.update();
//
//			//KER: We need to add new tuple variable to single contraints
//			for(Index3 rot: child.rots){
//				model.remove(singleConstraints[rot.pos]);
//				singleExprs[rot.pos].addTerm(1.0, tupVar);
//				singleConstraints[rot.pos] = model.addConstr(singleExprs[rot.pos], GRB.EQUAL, 1.0, "s_constr_"+rot.pos);
//
//				//We also add the variable to all of the pairs that contain this rotamer
//				//							for(int i=0; i<pairExprs[rot.pos][rot.aa][rot.rot].length;i++){
//				//								model.remove(pairConstraints[rot.pos][rot.aa][rot.rot][i]);
//				//								pairExprs[rot.pos][rot.aa][rot.rot][i].addTerm(-1, tupVar);
//				//								if(posInTuple(i,child.rots))//For pair constraints that have both positions defined by the rotamer add the tuple variable to the LHS
//				//									pairExprs[rot.pos][rot.aa][rot.rot][i].addTerm(1,tupVar);
//				//								pairConstraints[rot.pos][rot.aa][rot.rot][i] = model.addConstr(pairExprs[rot.pos][rot.aa][rot.rot][i], GRB.EQUAL, 0.0, "p_constr"+rot.pos+"_"+rot.aa+"_"+rot.rot+"_"+i);
//				//							}
//				//							model.update();
//				//							
//				//							//For pair constraints that have both positions defined by the rotamer add the tuple variable to the LHS
//				//							for(Index3 rot2:child.rots){
//				//								if(rot.pos == rot2.pos)
//				//									continue;
//				//								model.remove(pairConstraints[rot.pos][rot.aa][rot.rot][i]);
//				//								pairExprs[rot.pos][rot.aa][rot.rot][i].addTerm(-1, tupVar);
//				//								pairConstraints[rot.pos][rot.aa][rot.rot][i] = model.addConstr(pairExprs[rot.pos][rot.aa][rot.rot][i], GRB.EQUAL, 0.0, "p_constr"+rot.pos+"_"+rot.aa+"_"+rot.rot+"_"+i);
//				//							}
//				model.update();
//				posCtr=-1;
//				for(int p=0; p<child.E.length;p++){
//					posCtr++;
//					if(child.E[p] != null){
//						while(posInTuple(posCtr,child.rots)){
//							/*model.remove(pairConstraints[rot.pos][rot.aa][rot.rot][posCtr]);
//							pairExprs[rot.pos][rot.aa][rot.rot][posCtr].addTerm(1.0, tupVar);*/
//							posCtr++;
//						}
//						try{	
//							model.remove(pairConstraints[rot.pos][rot.aa][rot.rot][posCtr]);
//						}catch(Exception E2){
//							System.out.println("DELETE ME!!!");
//						}
//						for(int a=0; a<child.E[p].length;a++){
//							if(child.E[p][a] != null){
//
//								for(int r=0; r<child.E[p][a].length;r++){
//									try{
//										pairExprs[rot.pos][rot.aa][rot.rot][posCtr].addTerm(1.0, tupPairVars[posCtr][a][r]);
//									}catch(Exception E5){
//										E5.printStackTrace();
//									}
//
//								}
//							}}
//						pairExprs[rot.pos][rot.aa][rot.rot][posCtr].addTerm(-1.0, tupVar);
//						pairConstraints[rot.pos][rot.aa][rot.rot][posCtr] = model.addConstr(pairExprs[rot.pos][rot.aa][rot.rot][posCtr],GRB.EQUAL,0,"p_constr"+rot.pos+"_"+rot.aa+"_"+rot.rot+"_"+posCtr);
//					}
//				}
//
//			}
//			model.update();
//
//			posCtr=-1;
//
//			for(int p=0; p<child.E.length;p++){
//				posCtr++;
//				if(child.E[p] != null){
//					while(posInTuple(posCtr,child.rots)){
//						/*model.remove(pairConstraints[rot.pos][rot.aa][rot.rot][posCtr]);
//						pairExprs[rot.pos][rot.aa][rot.rot][posCtr].addTerm(1.0, tupVar);*/
//						posCtr++;
//					}
//					for(int a=0; a<tupPairVars[posCtr].length;a++){
//						for(int r=0; r<child.E[p][a].length;r++){
//
//							for(Index3 rot:child.rots){
//								try{
//									if(pairConstraints[posCtr][a][r][rot.pos] == null)
//										continue;
//								}catch(Exception E4){
//									System.out.println("DELETE ME!!!");
//								}
//
//								model.remove(pairConstraints[posCtr][a][r][rot.pos]);
//								pairExprs[posCtr][a][r][rot.pos].addTerm(1.0, tupPairVars[posCtr][a][r]);
//								pairConstraints[posCtr][a][r][rot.pos] = model.addConstr(pairExprs[posCtr][a][r][rot.pos], GRB.EQUAL, 0, "p_constr"+posCtr+"_"+a+"_"+r+"_"+rot.pos);
//							}
//
//						}
//					}
//
//				}
//
//			}
//
//
//			posCtr = -1;
//			//KER: Now add the new pair term constraints
//
//			for(int p=0; p<child.E.length;p++){
//				posCtr++;
//				if(child.E[p] != null){
//					GRBLinExpr expr = new GRBLinExpr();
//					while(posInTuple(posCtr,child.rots)){
//						posCtr++;
//					}
//					for(int a=0; a<tupPairVars[posCtr].length;a++){
//						for(int r=0; r<child.E[p][a].length;r++){
//							expr.addTerm(1.0, tupPairVars[posCtr][a][r]);
//						}
//					}
//					expr.addTerm(-1.0, tupVar);
//
//					//KER: don't add the constraint yet because it will have to be updated below
//					model.addConstr(expr, GRB.EQUAL, 0, tupPrefPairConstr+"_"+posCtr);
//					//tupleExprs[posCtr] = expr;
//
//				}
//
//			}
//
//			//KER: Add constraint that disallows the independent rotamers (this should be done for all previous tuples)
//			GRBLinExpr expr = new GRBLinExpr();
//			for(Index3 rot: child.rots){
//				expr.addTerm(1,singleVars[rot.pos][rot.aa][rot.rot]);
//			}
//
//			model.addConstr(expr, GRB.LESS_EQUAL, child.rots.length-1, "constr_"+tupPrefSingle);
//
//			for(EnergyTuple parent: parents){
//
//				String parentTupPrefSingle = "tup_s";
//				for(Index3 rot: parent.rots){
//					parentTupPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//				}
//
//				expr = new GRBLinExpr();
//				//Find the rotamer that's in child but not in the parent
//				Index3 addedRot = null;
//				for(Index3 rot: child.rots){
//					if(!posInTuple(rot.pos, parent.rots))
//						addedRot = rot;
//				}
//				try{
//					expr.addTerm(1, singleVars[addedRot.pos][addedRot.aa][addedRot.rot]);
//					expr.addTerm(1, model.getVarByName(parentTupPrefSingle));
//					model.addConstr(expr, GRB.LESS_EQUAL, 1, "constr_par_"+tupPrefSingle+"_"+addedRot.pos+"_"+addedRot.aa+"_"+addedRot.rot);
//				}catch(Exception E3){
//					System.out.println("DELETE ME!!!");
//				}
//
//			}
//
//
//			//KER: We need to make sure all of the tuple-tuple pairs are set up correctly
//			//KER: For every tuple already added if the two tuples don't have any rotamers in common
//			//KER: then add the pairwise variables to the pairwise constraint for the other tuple
//			GRBLinExpr[][] tupleExprs = new GRBLinExpr[3][child.rots.length];
//			GRBConstr[][] tupleConstrs = new GRBConstr[3][child.rots.length];
//			int ctr=0;
//			double bigNum = 5000;
//			for(Index3 rot:child.rots){
//				GRBLinExpr expr1 = new GRBLinExpr();
//				GRBLinExpr expr2 = new GRBLinExpr();
//				GRBLinExpr expr3 = new GRBLinExpr();
//				GRBVar binary = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "binary"+rot.pos+"_"+rot.aa+"_"+rot.rot+"_"+tupPrefPair);
//				model.update();
//				expr1.addTerm(bigNum, binary); //A >= 0 - Ly
//				expr2.addTerm(-bigNum, binary); // B >= 0 - M*(1-y)
//				expr3.addTerm(bigNum, binary); // B <= 0 + M*(1-y)
//				expr1.addTerm(-1,tupVar);
//				expr1.addTerm(-1, singleVars[rot.pos][rot.aa][rot.rot]);
//				for(TupleConstraint tc: tupleConstraints){
//					if(child.noSharedPositions(tc.tuple)){
//						expr2.addTerm(-1, tc.tupleVar);
//						try{
//							expr2.addTerm(1, tc.vars[rot.pos][rot.aa][rot.rot]);		
//						}catch(Exception E7){
//							E7.printStackTrace();
//						}
//						expr3.addTerm(-1, tc.tupleVar);
//						expr3.addTerm(1, tc.vars[rot.pos][rot.aa][rot.rot]);
//					}
//					else if(tc.tuple.containsRot(rot)){
//						expr1.addTerm(-1, tc.tupleVar);
//					}
//					//Add all ways (singles) that rot could be turned on
//				}
//				tupleConstrs[0][ctr] = model.addConstr(expr1, GRB.GREATER_EQUAL, 0, "if1_"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tupPrefPair);
//				tupleConstrs[1][ctr] = model.addConstr(expr2, GRB.GREATER_EQUAL, -bigNum, "if2a"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tupPrefPair);
//				tupleConstrs[2][ctr] = model.addConstr(expr3, GRB.LESS_EQUAL, bigNum, "if2b"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tupPrefPair);
//				tupleExprs[0][ctr] = expr1;
//				tupleExprs[1][ctr] = expr2;
//				tupleExprs[2][ctr] = expr2;
//				ctr++;
//			}
//
//			//KER: we also need to make sure we add the current tuple to the tuples that have already
//			//KER: been added
//			for(TupleConstraint tc: tupleConstraints){
//
//				String tmpPrefPair = "tup_p";
//				for(Index3 rot: tc.tuple.rots){
//					tmpPrefPair += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//				}
//
//				ctr=0;
//				for(Index3 rot: tc.tuple.rots){
//					if(tc.tuple.noSharedPositions(child)){
//						model.remove(tc.constrs[1][ctr]);
//						model.remove(tc.constrs[2][ctr]);
//						tc.exprs[1][ctr].addTerm(-1, tupVar);
//						try{
//							tc.exprs[1][ctr].addTerm(1, tupPairVars[rot.pos][rot.aa][rot.rot]);
//						}catch(Exception E6){
//							E6.printStackTrace();
//						}
//						tc.exprs[2][ctr].addTerm(-1, tupVar);
//						tc.exprs[2][ctr].addTerm(1, tupPairVars[rot.pos][rot.aa][rot.rot]);
//						tc.constrs[1][ctr] = model.addConstr(tc.exprs[1][ctr], GRB.GREATER_EQUAL, -bigNum, "if2a"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tmpPrefPair);
//						tc.constrs[2][ctr] = model.addConstr(tc.exprs[2][ctr], GRB.LESS_EQUAL, bigNum, "if2b"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tmpPrefPair);
//					}
//					else if(child.containsRot(rot)){
//						model.remove(tc.constrs[0][ctr]);
//						tc.exprs[0][ctr].addTerm(-1, tupVar);
//						tc.constrs[0][ctr] = model.addConstr(tc.exprs[0][ctr], GRB.GREATER_EQUAL, 0, "if1"+"_"+rot.pos+"_"+rot.aa+"_"+rot.rot+tmpPrefPair);
//					}
//					ctr++;
//				}
//			}
//
//			//KER: Force the largest tuple for the conformation to be used
//			for(TupleConstraint tc: tupleConstraints){
//
//				//If at least one rotamer is in common, then require only the longer one to be populated
//				if(tc.tuple.rots.length < child.rots.length && tc.tuple.allSharedPosMatch(child)){
//					String tmpPrefSingle = "tup_s";
//					for(Index3 rot: tc.tuple.rots){
//						tmpPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//					}
//
//					ArrayList<Index3> nonSharedRots = tc.tuple.nonSharedRots(child);
//					expr = new GRBLinExpr();
//					ctr=0; 
//					expr.addTerm(1.0, tc.tupleVar);
//					for(Index3 rot: nonSharedRots){
//						expr.addTerm(1.0, singleVars[rot.pos][rot.aa][rot.rot]);
//						ctr++;
//					}
//
//					model.addConstr(expr, GRB.LESS_EQUAL, ctr, "lrgTup_"+tupPrefSingle+"_"+tmpPrefSingle);
//				}
//			}
//
//
//			//			System.out.println("DELETE ME");
//			//			expr = new GRBLinExpr();
//			//			expr.addTerm(1.0, tupVar);
//			//			model.addConstr(expr, GRB.EQUAL, 1, "stupidConstr");
//			//			///END DELETE ME
//
//			//DELETE ME!!!!
//			/*if(child.rots.length == 2){
//			expr = new GRBLinExpr();
//			expr.addTerm(1, singleVars[3][0][8]);
//			badConstr = model.addConstr(expr, GRB.EQUAL, 1, "fixPos3");
//			model.update();
//			}*/
//
//			model.update();
//			//model.write("gurobiModel.lp");
//
//			tupleConstraints.add(new TupleConstraint(child,tupVar,tupPairVars,tupleExprs,tupleConstrs));
//
//		}
//		catch(Exception E){
//			System.out.println("SOMETHING IS WRONG!!!");
//			E.printStackTrace();
//		}
//	}


	public void addTuple(EnergyTuple child, Emat emat) {
		if(tupleConstraints == null){
			tupleConstraints = new HashMap<ArrayList<Index3>,TupleConstraint>();
		}

		//KER: We need to add variables for the tuple and all of the tuple pairs

		ArrayList<Index3> tmpRots = new ArrayList<Index3>(child.rots.length);
		for(Index3 i3 : child.rots )
			tmpRots.add(i3);
		
		if(!tupleConstraints.containsKey(tmpRots)){
		
		try{
			//KER: First add the tuple variable
			String tupPrefSingle = "tup";
//			String tupPrefPair = "tup_p";
			String tupPrefSingleConstr = "tup_constr";
			for(Index3 rot: child.rots){
				tupPrefSingle += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
//				tupPrefPair += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
				tupPrefSingleConstr += "_"+rot.pos+"_"+rot.aa+"_"+rot.rot;
			}
			double cost;// = child.intraE;
			
			
			//Get all subTuples
			LinkedList<TupleConstraint> allTups = getSubTuples(child);
			
			
			double curCost = child.intraE - getTupleE(child,null,emat);
			for(TupleConstraint tc: allTups)
				curCost -= tc.cost; 
			
			if(curCost < 0)
				curCost = 0;
			
			GRBVar tupVar = model.addVar(0.0, 1.0, curCost, GRB.BINARY, tupPrefSingle);

			model.update();

			TupleConstraint newTupleConstr = new TupleConstraint(child, tupVar, curCost);
			

			//Add tuple constraints
			GRBLinExpr expr = new GRBLinExpr();
			for(Index3 rot : child.rots){
				expr.addTerm(1.0, singleVars[rot.pos][rot.aa][rot.rot]);
			}
			expr.addTerm(-1.0, tupVar);
			model.addConstr(expr, GRB.LESS_EQUAL, child.rots.length-1 , tupPrefSingleConstr+"_all");
			
			//Add tuples <= x(rot) for each rot
			for(Index3 rot : child.rots){
				expr = new GRBLinExpr();
				expr.addTerm(1.0, tupVar);
				expr.addTerm(-1.0, singleVars[rot.pos][rot.aa][rot.rot]);
				
				model.addConstr(expr, GRB.LESS_EQUAL, 0.0 , tupPrefSingleConstr+"_p_"+rot.pos);
			}
			
			tupleConstraints.put(tmpRots, newTupleConstr);
			
			
			//Need to check if we need to update any of the tupleVar costs
			for(TupleConstraint tc: tupleConstraints.values()){
				if(tc.tuple.rots.length > child.rots.length && tc.tuple.containsAllRot(child)){
					tc.cost -= curCost;
					
					if(tc.cost < 0)
						tc.cost = 0;
					
					tc.tupleVar.set(GRB.DoubleAttr.Obj, tc.cost);
					
				}
					
			}

			model.update();
			//model.write("gurobiModel.lp");

		}
		catch(Exception E){
			System.out.println("SOMETHING IS WRONG!!!");
			E.printStackTrace();
		}
		}
	}


	private boolean posInTuple(int pos, Index3[] rots){
		for(Index3 r: rots){
			if(r.pos == pos)
				return true;
		}

		return false;
	}

	// This function returns the xth token in string s
	public static String getToken(String s, int x) {

		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s,"_");

		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
				st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	public class GurobiConf{
		EMatrixEntryWIndex[] conf;
		LinkedList<EnergyTuple> tuples = null;

		GurobiConf(EMatrixEntryWIndex[] c, LinkedList<EnergyTuple> curTuples){
			conf = c;
			tuples = curTuples;
		}

	}

	public class TupleConstraint{
		EnergyTuple tuple;
		GRBVar tupleVar;
		double cost;
		
		public TupleConstraint(EnergyTuple t, GRBVar tv, double cost) {
			tuple = t;
			tupleVar = tv;
			this.cost = cost;
		}
	}

	/*
	 * Dissallow the sequence of the current conf from the search
	 */
	public void removeSeq(EMatrixEntryWIndex[] conf) {
		GRBLinExpr expr = new GRBLinExpr();


		for(int i=0; i<conf.length;i++){
			for(int j=0; j<singleVars[conf[i].pos1()][conf[i].aa1()].length;j++)
				expr.addTerm(1.0, singleVars[conf[i].pos1()][conf[i].aa1()][j]);
		}
		try {
			model.addConstr(expr, GRB.LESS_EQUAL, conf.length-1, "exclude_seq"+excludedSeqs);
			model.update();
		} catch (GRBException e) {

			e.printStackTrace();
		}

		excludedSeqs++;
	}


	public int[] getConfFromSeq(int[] confSoFar, int[][] seqIndexOffset) {
		int conf[] = new int[confSoFar.length];

		GRBVar[] vars = model.getVars();
		for(int i=0; i<vars.length;i++){
			try {
				if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){

//					System.out.println(vars[i].get(GRB.StringAttr.VarName)
//							+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));

					String varName = vars[i].get(GRB.StringAttr.VarName);
					if(varName.startsWith("s")){
						StringTokenizer st = new StringTokenizer(varName,"_");
						st.nextElement(); //remove the "s" at the beginning
						int p = new Integer(st.nextToken()).intValue();
						int a = new Integer(st.nextToken()).intValue();
						int r = new Integer(st.nextToken()).intValue();
						conf[p] = seqIndexOffset[p][confSoFar[p]]+r;
					}

				}
			} catch (GRBException e) {
				e.printStackTrace();
				return null;
			}
		}	

		return conf;
	}

	//This is the slow way, but is useful just to check that everything works.
	private LinkedList<TupleConstraint> getSubTuples(EnergyTuple tuple){
		LinkedList<TupleConstraint> retTups = new LinkedList<TupleConstraint>();
		
		for(TupleConstraint tc : tupleConstraints.values()){
			if(tuple.containsAllRot(tc.tuple)){
				retTups.add(tc);
			}
		}
		
		return retTups;
		
	}
	
	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	public double getTupleE(EnergyTuple tuple, LinkedList<EnergyTuple> curTuples, Emat emat) {

		double bestE = 0.0;

		boolean[] excludeLevel = new boolean[singleVars.length];
		for(int i=0; i<excludeLevel.length;i++)
			excludeLevel[i] = true;
		
		for(Index3 rot : tuple.rots){
			excludeLevel[rot.pos] = false;
		}
		
		if(curTuples != null && curTuples.size() > 0){
			for(EnergyTuple curTuple:curTuples){

				for(int i=0; i<curTuple.rots.length;i++){
					excludeLevel[curTuple.rots[i].pos] = true;
				}

				//Intra E for tuple
				bestE += curTuple.intraE;

				//Pair E for tuple
				for (Index3 i3: tuple.rots){
					if(!excludeLevel[i3.pos]){
						
						for(Index3 i3_2: curTuple.rots){
							double pairEtoAdd = emat.getPairMinE(i3, i3_2);
							bestE += pairEtoAdd;//pairwiseMinEnergyMatrix[index1][index2].eme.minE();
						}
					}
				}
			}
		}

		for (Index3 i1: tuple.rots){

			if(!excludeLevel[i1.pos]){
				double tmpSingleE = emat.getSingleMinE(i1); //Add the intra-rotamer energy
				
				bestE += tmpSingleE;
				
				for(Index3 i2: tuple.rots){
					if(i1.pos < i2.pos && !excludeLevel[i2.pos] && emat.areNeighbors(i1.pos, i2.pos)){
						double tmpPairE = emat.getPairMinE(i1, i2);
						bestE += tmpPairE;
					}
				}
			}

		}

		return bestE;
	}



	public void calcBoundPerPos(Index3wVal[] boundPerPos) {
		GRBVar[] vars = model.getVars();
		for(int i=0; i<vars.length;i++){
			try {
				if(vars[i].get(GRB.DoubleAttr.X) > 0.0001){// && vars[i].get(GRB.StringAttr.VarName).startsWith("s")){

//				System.out.println(vars[i].get(GRB.StringAttr.VarName)
//						+ " " +vars[i].get(GRB.DoubleAttr.X)+" "+vars[i].get(GRB.DoubleAttr.Obj));

					String varName = vars[i].get(GRB.StringAttr.VarName);
					double E = vars[i].get(GRB.DoubleAttr.Obj);
					
					ArrayList<Index3> rots = new ArrayList<Index3>();
					StringTokenizer st = new StringTokenizer(varName,"_");
					st.nextElement(); //remove the prefix at the beginning
					while(st.hasMoreTokens()){
						int p = new Integer(st.nextToken()).intValue();
						int a = new Integer(st.nextToken()).intValue();
						int r = new Integer(st.nextToken()).intValue();
						rots.add(new Index3(p,a,r));
					}
					for(Index3 rot: rots){
						boundPerPos[rot.pos].val += E/(double)rots.size(); 
					}
					

				}
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (GRBException e) {
				e.printStackTrace();
			}
		}
	}


}
