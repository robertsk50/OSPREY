/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	DEEIndirect.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Performs indirect DEE rotamer and pairs pruning
 * considering the entire protein as the system being pruned.
 * Uses rigid DEE or iMinDEE
 */


public class DEEIndirect extends DEE {

	//useFlags and useMinDEEPruningEw are not used
	//because we always use flags and iMinDEE is assumed if doing minimization

	private double pairK[][][][][];//Precomputed K(i_r,i_t) matrix
	boolean inZ[];//Indicates which of the flexible residue positions are in the pruning zone

	int prunedSinglesCurRun;
	int prunedPairsCurRun;

	//constructor
	DEEIndirect (Emat arpMatrix, MutableResParams strMut, double initEw,
			StrandRotamers strandLRot[], boolean residueMut[],
			boolean doMin, boolean mb, boolean dDEE, boolean minBB,
			boolean typeDep, boolean aIMinDEE,
			double aIval, boolean tripFlags[][][][][][][][][], boolean doPerts, boolean inZ[]) {

		init(arpMatrix, null, arpMatrix.numMutPos(),
				strMut, initEw, strandLRot, null, doMin, null, null,
				arpMatrix.pairs.pruned, true, minBB, null, null, typeDep, aIMinDEE, aIval,
				mb, dDEE, residueMut, tripFlags, doPerts);

		this.inZ = inZ;


		if(doMinimize && (!doIMinDEE) ){
			System.err.println("ERROR: iMinDEE required for indirect pruning and thus for algOption >= 4");
			System.exit(1);
		}

		//Now allocate space for pairK
		pairK = new double[numMutable][][][][];





		for(int pos=0; pos<numMutable; pos++){

			if(inZ[pos]){
				pairK[pos] = new double[arpMatrix.pairs.E[pos].length][][][];

				for(int AA=0; AA<pairK[pos].length; AA++){
					if(arpMatrix.pairs.E[pos][AA] == null)
						continue;

					pairK[pos][AA] = new double[arpMatrix.pairs.E[pos][AA].length][][];

					for(int rot=0; rot<pairK[pos][AA].length; rot++){
						pairK[pos][AA][rot] = new double[arpMatrix.pairs.E[pos][AA][rot].length][];

						for(int AA2=0; AA2<arpMatrix.singles.E[pos].length; AA2++){
							if(arpMatrix.singles.E[pos][AA2] == null)
								continue;

							pairK[pos][AA][rot][AA2] = new double[arpMatrix.singles.E[pos][AA2].length];
						}
					}
				}
			}
		}
	}


	//Compute the conformations that can be eliminated
	//Return a boolean array in which an element is true if
	//the corresponding r at i can be eliminated, and false otherwise
	//Also saves "true" to pairwiseMinEnergyMatrix.pairs.pruned for pairs that can be eliminated
	public void ComputeEliminatedRotConf(){

		boolean done = false;
		numRuns = 1;

		while (!done){

			prunedSinglesCurRun = 0;
			prunedPairsCurRun = 0;

			System.out.println("Current run: "+numRuns);

			//First precompute K(i_r,i_t) for each rotamer pair in the pruning zone
			for (int curPos=0; curPos<numMutable; curPos++){

				if( inZ[curPos] ){

					System.out.print("Starting precomputation: residue " + curPos);

					//There are two types of indices for amino acids:
					//those among the allowed AAs at some residue
					//(variable names here for this are "AA" or "AA[extra letter or number]")
					//and those among all the 20 possible AAs
					//(used in energy matrix, pairK, etc; names are like curAA or AANumAtPos)

					for (int AA=0; AA<pairwiseMinEnergyMatrix.singles.E[curPos].length; AA++){

						System.out.print(".");
						for(int curRot=0; curRot<pairwiseMinEnergyMatrix.singles.E[curPos][AA].length; curRot++){//This is i_r

							if ( !pairwiseMinEnergyMatrix.getSinglePruned(curPos, AA, curRot) ){//not already pruned

								//Find i_t now
								for(int AAt=0; AAt<pairwiseMinEnergyMatrix.singles.E[curPos].length; AAt++ ){

									for(int tRot=0; tRot<pairwiseMinEnergyMatrix.singles.E[curPos][AAt].length; tRot++){
										if(!pairwiseMinEnergyMatrix.getSinglePruned(curPos, AAt, tRot))
											pairK[curPos][AA][curRot][AAt][tRot] = precomputePairK(curPos, AA, curRot, AAt, tRot);
									}
								}
							}
						}
					}
					System.out.println("done");
				}
			}


			//Now try to prune singles.
			for (int curPos=0; curPos<numMutable; curPos++){

				if( inZ[curPos] ){

					System.out.print("Starting residue "+curPos);

					for (int AA=0; AA<pairwiseMinEnergyMatrix.singles.E[curPos].length; AA++){

						System.out.print(".");

						for(int curRot=0; curRot<pairwiseMinEnergyMatrix.singles.E[curPos][AA].length; curRot++){

							if (!pairwiseMinEnergyMatrix.getSinglePruned(curPos, AA, curRot)){//not already pruned

								if (CanEliminate(curPos, AA, curRot)){
									pairwiseMinEnergyMatrix.setSinglePruned(curPos, AA, curRot, true);
									prunedSinglesCurRun++;
								}
								else
									pairwiseMinEnergyMatrix.setSinglePruned(curPos, AA, curRot, false);
							}
						}
					}
					System.out.println("done");

				}
			}


			//Now try to prune pairs

			for (int curPos1=0; curPos1<numMutable; curPos1++){

				if ( inZ[curPos1] && ( (magicBullet)||(!distrDEE)||(resInPair[curPos1]) ) ){ //mb-pairs or not distrDEE or cur res is in distr pair; also, this res must be in the pruning zone


					System.out.print("Starting residue "+curPos1);
					System.out.print("..");

					for (int curPos2=curPos1+1; curPos2<numMutable; curPos2++){

						if ( inZ[curPos2] && ( (magicBullet)||(!distrDEE)||(resInPair[curPos2]) ) ){ //mb-pairs or not distrDEE or cur res is in distr pair


							for (int AA1=0; AA1<pairwiseMinEnergyMatrix.singles.E[curPos1].length; AA1++){

								for(int curRot1=0; curRot1<pairwiseMinEnergyMatrix.singles.E[curPos1][AA1].length; curRot1++){

									if (!pairwiseMinEnergyMatrix.getSinglePruned(curPos1, AA1, curRot1)){//not already pruned

										for (int AA2=0; AA2<pairwiseMinEnergyMatrix.singles.E[curPos2].length; AA2++){
											for(int curRot2=0; curRot2<pairwiseMinEnergyMatrix.singles.E[curPos2][AA2].length; curRot2++){

												if (!pairwiseMinEnergyMatrix.getSinglePruned(curPos2, AA2, curRot2)){//not already pruned

													if (!pairwiseMinEnergyMatrix.pairs.pruned[curPos1][AA1][curRot1][curPos2][AA2][curRot2]){ //rotamer pair not already pruned

														if (CanEliminate(curPos1, AA1, curRot1, curPos2, AA2, curRot2)){
															pairwiseMinEnergyMatrix.pairs.pruned[curPos1][AA1][curRot1][curPos2][AA2][curRot2] = true;
															pairwiseMinEnergyMatrix.pairs.pruned[curPos2][AA2][curRot2][curPos1][AA1][curRot1] = true;
															prunedPairsCurRun++;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
					System.out.println("done");
				}
			}

			//Higher-order indirect pruning would go here

			System.out.println("Number of rotamers pruned this run: "+prunedSinglesCurRun);
			System.out.println("Number of pairs pruned this run: "+prunedPairsCurRun);
			System.out.println("DEE: The minimum difference is "+minDiff);
			System.out.println();

			if ( prunedSinglesCurRun==0 && prunedPairsCurRun==0 ) //no rotamers or pairs pruned this run, so done
				done = true;
			else
				numRuns++;
		}

	}

	//Called only by ComputeEliminatedRotConf(.)
	private boolean CanEliminate (int posNum, int AANumAtPos, int rotNumAtPos){

		double checkSum;

		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)

		if ((!pairwiseMinEnergyMatrix.getSinglePruned(posNum, AANumAtPos, rotNumAtPos))){ //not pruned yet

			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning

			for (int altAA=0; altAA<pairwiseMinEnergyMatrix.singles.E[posNum].length; altAA++){

				if( (!typeDependent) || (altAA==AANumAtPos) ){

					for (int altRot=0; altRot<pairwiseMinEnergyMatrix.singles.E[posNum][altAA].length; altRot++){

						//if t and r are not actually the same rotamer of the same AA
						if (!((altAA==AANumAtPos)&&(altRot==rotNumAtPos))){

							//at this point, we know what r at i and t at i are

							//if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
							if (!pairwiseMinEnergyMatrix.getSinglePruned(posNum, altAA, altRot)){ //not pruned

								checkSum = 0;

								int s_AAtypes[] = new int[numMutable];//AA types and rotamers used to prune
								int s_rots[] = new int[numMutable];

								for(int pos2=0; pos2<numMutable; pos2++){

									if(pos2 == posNum){
										checkSum += pairK[posNum][AANumAtPos][rotNumAtPos][altAA][altRot];
										continue;
									}
									else if( inZ[pos2] ){

										double Kmaxmin = Double.NEGATIVE_INFINITY;

										for(int sAA=0; sAA<pairwiseMinEnergyMatrix.singles.E[pos2].length; sAA++){
											
											for(int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[pos2][sAA].length; sRot++){

												if( pairwiseMinEnergyMatrix.getSinglePruned(pos2, sAA, sRot)
														|| pairwiseMinEnergyMatrix.pairs.pruned[pos2][sAA][sRot][posNum][altAA][altRot] )
													continue;//Don't use this j_s

													double sKmin = Double.POSITIVE_INFINITY;

													for(int uAA=0; uAA<pairwiseMinEnergyMatrix.singles.E[pos2].length; uAA++){

														for(int uRot=0; uRot<pairwiseMinEnergyMatrix.singles.E[pos2][uAA].length; uRot++){

															if( pairwiseMinEnergyMatrix.getSinglePruned(pos2, uAA, uRot)
																	|| pairwiseMinEnergyMatrix.pairs.pruned[pos2][uAA][uRot][posNum][AANumAtPos][rotNumAtPos] )
																continue;

															if( pairK[pos2][uAA][uRot][sAA][sRot] < sKmin )
																sKmin = pairK[pos2][uAA][uRot][sAA][sRot];
														}
													}

													if(sKmin == Double.POSITIVE_INFINITY)
														return true;//Rotamer r can be pruned because it is incompatible with all unpruned rotamers at pos2

													if(sKmin > Kmaxmin){
														boolean isCompatible = true;//is index_s compatible with the rest of the s-indices used?
														for(int a=0;a<pos2;a++){
															if( inZ[a] && (a != posNum) && pairwiseMinEnergyMatrix.pairs.pruned[a][s_AAtypes[a]][s_rots[a]][pos2][sAA][sRot] )
																isCompatible = false;
														}
														if(isCompatible){
															Kmaxmin = sKmin;
															s_AAtypes[pos2] = sAA;
															s_rots[pos2] = sRot;
														}
													}

											}
										}

										if(Kmaxmin == Double.NEGATIVE_INFINITY){
											checkSum = Double.NEGATIVE_INFINITY;
											break;//r will not be pruned using t
										}

										checkSum += Kmaxmin;//add max min K(j_u,j_s) to checkSum where j_u is compatible with i_r and j_s with i_t
									}
								}

								if ( checkSum > curEw)
									return true;//this rotamer can be pruned/eliminated
								else
									minDiff = Math.max(minDiff,checkSum);


							}
						}
					}
				}
			}
		}
		else //already pruned
			return true;

		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}



	//Check if a pair can be eliminated
	private boolean CanEliminate (int posNum1, int AANumAtPos1, int rotNumAtPos1, int posNum2,
			int AANumAtPos2, int rotNumAtPos2){


		if( pairwiseMinEnergyMatrix.getSinglePruned(posNum1, AANumAtPos1, rotNumAtPos1)
				|| pairwiseMinEnergyMatrix.getSinglePruned(posNum2, AANumAtPos2, rotNumAtPos2)
				|| pairwiseMinEnergyMatrix.pairs.pruned[posNum1][AANumAtPos1][rotNumAtPos1][posNum2][AANumAtPos2][rotNumAtPos2] )
			return true;//Already pruned


			for (int altAA1=0; altAA1<pairwiseMinEnergyMatrix.singles.E[posNum1].length; altAA1++){

				for (int altRot1=0; altRot1<pairwiseMinEnergyMatrix.singles.E[posNum1][altAA1].length; altRot1++){

					//if t and r are not actually the same rotamer of the same AA
					//if (!((altAA1==AANumAtPos1)&&(altRot1==rotNumAtPos1))){

					if (!pairwiseMinEnergyMatrix.getSinglePruned(posNum1, altAA1, altRot1)){ //not pruned

						for (int altAA2=0; altAA2<pairwiseMinEnergyMatrix.singles.E[posNum2].length; altAA2++){

							if( (!typeDependent) || ( (altAA1==AANumAtPos1) && (altAA2==AANumAtPos2) ) ){

								for (int altRot2=0; altRot2<pairwiseMinEnergyMatrix.singles.E[posNum2][altAA2].length; altRot2++){

									//if t and r are not actually the same rotamer of the same AA
									//if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){

										if ( (!pairwiseMinEnergyMatrix.getSinglePruned(posNum2 ,altAA2, altRot2))
												&& (!pairwiseMinEnergyMatrix.pairs.pruned[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2]) ){ //not pruned

											if( canEliminateUsing( posNum1, AANumAtPos1, rotNumAtPos1, altAA1, altRot1,
													posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2 ) )
												return true;
											else if (magicBullet) //magic bullet pairs, so no further checks
												return false;
										}
										//}
								}
							}
						}
					}
					//}
				}
			}

			//We have tried all of the other rotamers at the current position and none
			//of them is able to prune the given rotamer, so we return false
			return false;
	}



	//Check if a pair (indicated by r) can be eliminated using another pair (indicated by t)
	//(residues indicated by 1 and 2)
	private boolean canEliminateUsing (int pos1, int r1AA, int r1Rot, int t1AA, int t1Rot,
			int pos2, int r2AA, int r2Rot, int t2AA, int t2Rot) {

		double checkSum = 0;

		int s_AAtypes[] = new int[numMutable];//AA types and rotamers used to prune (in the Kmaxmin terms)
		int s_rots[] = new int[numMutable];

		//First add up K's at other positions (these are stored in pairK)
		for(int posj=0; posj<numMutable; posj++){

			if( inZ[posj] && (posj != pos1) && (posj != pos2) ){

				double Kmaxmin = Double.NEGATIVE_INFINITY;

				for(int sAA=0; sAA<pairwiseMinEnergyMatrix.singles.E[posj].length; sAA++){

					for(int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[posj][sAA].length; sRot++){

						if( pairwiseMinEnergyMatrix.getSinglePruned(posj, sAA, sRot)
								|| pairwiseMinEnergyMatrix.pairs.pruned[posj][sAA][sRot][pos1][t1AA][t1Rot]
										|| pairwiseMinEnergyMatrix.pairs.pruned[posj][sAA][sRot][pos2][t2AA][t2Rot] )
							continue;//Don't use this j_s

							if( isPrunedTriple(posj,sAA,sRot,pos1,t1AA,t1Rot,pos2,t2AA,t2Rot) )
								continue;
							//This is checked after lower-order pruning to make sure we don't call this function for pruned rotamers


							double sKmin = Double.POSITIVE_INFINITY;

							for(int uAA=0; uAA<pairwiseMinEnergyMatrix.singles.E[posj].length; uAA++){

								for(int uRot=0; uRot<pairwiseMinEnergyMatrix.singles.E[posj][uAA].length; uRot++){

									if( pairwiseMinEnergyMatrix.getSinglePruned(posj, uAA, uRot)
											|| pairwiseMinEnergyMatrix.pairs.pruned[posj][uAA][uRot][pos1][r1AA][r1Rot]
													|| pairwiseMinEnergyMatrix.pairs.pruned[posj][uAA][uRot][pos2][r2AA][r2Rot] )
										continue;

									if( isPrunedTriple(posj,uAA,uRot,pos1,r1AA,r1Rot,pos2,r2AA,r2Rot) )
										continue;

									if( pairK[posj][uAA][uRot][sAA][sRot] < sKmin )
										sKmin = pairK[posj][uAA][uRot][sAA][sRot];
								}
							}

							if ( sKmin == Double.POSITIVE_INFINITY )
								return true;
							//In this case, there is no rotamer u at posj that is compatible with both r1 at pos1 and r2 at pos1
							//So the rotamer pair (r1 at pos1, r2 at pos2) is impossible


							if(sKmin > Kmaxmin){
								boolean isCompatible = true;//is index_s compatible with the rest of the s-indices used?
										for(int a=0;a<posj;a++){
											if( inZ[a] && (a != pos1) && (a != pos2) && pairwiseMinEnergyMatrix.pairs.pruned[a][s_AAtypes[a]][s_rots[a]][posj][sAA][sRot] )
												isCompatible = false;
										}
										if(isCompatible){
											Kmaxmin = sKmin;
											s_AAtypes[posj] = sAA;
											s_rots[posj] = sRot;
										}
							}
					}
				}

				if(Kmaxmin == Double.NEGATIVE_INFINITY)
					return false;
				//This messes up the pruning of (r1 at pos1, r2 at pos2)
				//So go try to find some other competitor pair

				checkSum += Kmaxmin;//add max min K(j_u,j_s) to checkSum where j_u is compatible with i_r and j_s with i_t
			}
		}



		//Now add in K(i_r,i_t) (bold i,r,t)
		checkSum += pairwiseMinEnergyMatrix.getIntraAndShellE( pos1, r1AA, r1Rot ) 	//intra and shell energies: 1st residue
				- pairwiseMinEnergyMatrix.getIntraAndShellE( pos1, t1AA, t1Rot )
				+ pairwiseMinEnergyMatrix.getIntraAndShellE( pos2, r2AA, r2Rot ) //2nd residue
				- pairwiseMinEnergyMatrix.getIntraAndShellE( pos2, t2AA, t2Rot );


		//Pairwise energy within the pair
		checkSum +=  pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, pos2, r2AA, r2Rot )
				- pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, pos2, t2AA, t2Rot );


		//Now sum the other pairwise energy differences within the pruning zone
		for(int posj=0; posj<pos2; posj++){
			//We are assuming pos2 > pos1

			if( inZ[posj] && (posj != pos1) ){

				double minTerm = Double.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) over j_s in R_(j,i_r)
				double maxTerm = Double.NEGATIVE_INFINITY;//This will be the maximum E(i_t,j_s) over j_s in R(j,i_t)

				for(int sAA=0; sAA<pairwiseMinEnergyMatrix.singles.E[posj].length; sAA++){
					for (int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[posj][sAA].length; sRot++){

						if( pairwiseMinEnergyMatrix.getSinglePruned(posj, sAA, sRot) )//Don't consider pruned rotamers
						continue;

						if( ! ( pairwiseMinEnergyMatrix.pairs.pruned[pos1][r1AA][r1Rot][posj][sAA][sRot]
								|| pairwiseMinEnergyMatrix.pairs.pruned[pos2][r2AA][r2Rot][posj][sAA][sRot]
										|| isPrunedTriple(pos1,r1AA,r1Rot,pos2,r2AA,r2Rot,posj,sAA,sRot) ) ) {
							double checkMin = pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, posj, sAA, sRot );
							if(posj<pos1)
								checkMin += pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, posj, sAA, sRot );
							if ( checkMin < minTerm )
								minTerm = checkMin;
						}
						if( ! ( pairwiseMinEnergyMatrix.pairs.pruned[pos1][t1AA][t1Rot][posj][sAA][sRot] || pairwiseMinEnergyMatrix.pairs.pruned[pos2][t2AA][t2Rot][posj][sAA][sRot]
								|| isPrunedTriple(pos1,t1AA,t1Rot,pos2,t2AA,t2Rot,posj,sAA,sRot) ) ){
							double checkMax = pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, posj, sAA, sRot );
							if(posj<pos1)
								checkMax += pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, posj, sAA, sRot );
							if ( checkMax > maxTerm )
								maxTerm = checkMax;
						}
					}
				}


				if(minTerm == Double.POSITIVE_INFINITY)
					return true;//(r1, r2) can be pruned because it is incompatible with all rotamers at posj

				if(maxTerm == Double.NEGATIVE_INFINITY){
					pairwiseMinEnergyMatrix.pairs.pruned[pos1][t1AA][t1Rot][pos2][t2AA][t2Rot] = true;
					pairwiseMinEnergyMatrix.pairs.pruned[pos2][t2AA][t2Rot][pos1][t1AA][t1Rot] = true;
					prunedPairsCurRun++;
					return false;
				}
				//(t1, t2) can be pruned because it is incompatible with all rotamers at posj
				//So, as above, we cannot meaningfully continue trying to prune (r1,r2) with it



				checkSum += minTerm - maxTerm;
			}
		}


		//Finally add terms for residues outside the pruning zone
		for(int posj=0;posj<numMutable;posj++){

			if ( !inZ[posj] ) {//posj is not in the pruning zone

				double minTerm = Double.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) - E(i_t,j_s)  over j_s in R_j

				for(int sAA=0; sAA<pairwiseMinEnergyMatrix.singles.E[posj].length; sAA++){

					for (int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[posj][sAA].length; sRot++){

						if( ! ( pairwiseMinEnergyMatrix.getSinglePruned(posj, sAA, sRot)
								|| pairwiseMinEnergyMatrix.pairs.pruned[pos1][r1AA][r1Rot][posj][sAA][sRot]
										|| pairwiseMinEnergyMatrix.pairs.pruned[pos2][r2AA][r2Rot][posj][sAA][sRot] ) ){

							if( ! isPrunedTriple(pos1,r1AA,r1Rot,pos2,r2AA,r2Rot,posj,sAA,sRot) ){

								double checkMin = + pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, posj, sAA, sRot )
										+ pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, posj, sAA, sRot )
										- pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, posj, sAA, sRot )
										- pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, posj, sAA, sRot );

								if ( checkMin < minTerm )
									minTerm = checkMin;
							}
						}
					}
				}


				if(minTerm == Double.POSITIVE_INFINITY)
					return true;//(r1, r2) can be pruned because it is incompatible with all rotamers at posj

				checkSum += minTerm;
			}
		}


		if ( checkSum > curEw)
			return true;//this rotamer can be pruned/eliminated
		else
			minDiff = Math.max(minDiff,checkSum);

		return false;

	}


	private double precomputePairK(int curPos, int curAA, int curRot, int tAA, int tRot){
		//Computes K(i_r,i_t); i = curPos; curAA & curRot give r; AAt and tRot give t
		//This should only be called if i_r and i_t have not been pruned

		//First set K = E(i_r) - E(i_t)
		double K = pairwiseMinEnergyMatrix.getIntraAndShellE( curPos, curAA, curRot ) 	//intra energy for i_r
				- pairwiseMinEnergyMatrix.getIntraAndShellE( curPos, tAA, tRot );   // i_t now

		//Now sum pairwise energy differences within the pruning zone
		for(int pos2=0; pos2<curPos; pos2++){

			if( inZ[pos2] ){

				double minTerm = Double.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) over j_s in R_(j,i_r)
				double maxTerm = Double.NEGATIVE_INFINITY;//This will be the maximum E(i_t,j_s) over j_s in R(j,i_t)

				for(int AAs=0; AAs<pairwiseMinEnergyMatrix.singles.E[pos2].length; AAs++){

					for (int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[pos2][AAs].length; sRot++){

						if( pairwiseMinEnergyMatrix.getSinglePruned(pos2, AAs, sRot) )//Don't consider pruned rotamers
						continue;

						if( ! pairwiseMinEnergyMatrix.pairs.pruned[curPos][curAA][curRot][pos2][AAs][sRot] ){
							if ( pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, AAs, sRot ) < minTerm )
								minTerm = pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, AAs, sRot );
						}
						if( ! pairwiseMinEnergyMatrix.pairs.pruned[curPos][tAA][tRot][pos2][AAs][sRot] ){
							if ( pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, AAs, sRot ) > maxTerm )
								maxTerm = pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, AAs, sRot );
						}
					}
				}

				if( minTerm == Double.POSITIVE_INFINITY || maxTerm == Double.NEGATIVE_INFINITY ){

					if(minTerm == Double.POSITIVE_INFINITY){
						pairwiseMinEnergyMatrix.setSinglePruned(curPos,curAA,curRot,true);//ir can be pruned because it is incompatible with all rotamers at position j (pos2)
						prunedSinglesCurRun++;
					}

					if(maxTerm == Double.NEGATIVE_INFINITY){
						if(!pairwiseMinEnergyMatrix.getSinglePruned(curPos,tAA,tRot)){
							pairwiseMinEnergyMatrix.setSinglePruned(curPos,tAA,tRot,true);//it can be pruned (same argument)
							prunedSinglesCurRun++;
						}
					}

					return Double.NaN;//Either way, the pruning makes K(ir,it) useless
				}



				K += minTerm - maxTerm;
			}
		}

		//and finally terms from outside the pruning zone
		for(int pos2=0; pos2<curPos; pos2++){

			if( !inZ[pos2] ){//pos2 is not in the pruning zone

				double minTerm = Double.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) - E(i_t,j_s)  over j_s in R_j

				for(int AAs=0; AAs<pairwiseMinEnergyMatrix.singles.E[pos2].length; AAs++){

					for (int sRot=0; sRot<pairwiseMinEnergyMatrix.singles.E[pos2][AAs].length; sRot++){

						if( ! ( pairwiseMinEnergyMatrix.getSinglePruned(pos2, AAs, sRot)
								|| pairwiseMinEnergyMatrix.pairs.pruned[curPos][curAA][curRot][pos2][AAs][sRot] ) ) {

							double checkMin = pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, AAs, sRot )
									- pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, AAs, sRot );

							if ( checkMin < minTerm )
								minTerm = checkMin;
						}
					}
				}


				if(minTerm == Double.POSITIVE_INFINITY){
					pairwiseMinEnergyMatrix.setSinglePruned(curPos, curAA, curRot, true);//ir can be pruned because it is incompatible with all rotamers at position j (pos2)
					prunedSinglesCurRun++;
					return Double.NaN;
				}

				K += minTerm;
			}
		}

		return K;
	}

}