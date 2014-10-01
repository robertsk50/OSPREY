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
//	CETMatrix.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;


//Matrix of ContETerms for the various intra+shell and pairwise energies

public class CETMatrix implements Serializable {

    //these arrays are indexed by (residue #'d within active site, AA index, rotamer/RC #)
    //simlar to PairwiseEnergyMatrix

    ContETerm intraAndShellBounds[][][] = null;
    ContETerm pairwiseBounds[][][][][][] = null;
    
    int[][][][] supRot;
    
    int numRes;
    double ivalCutoff;//RCs unpruned at this ival are guaranteed to be included in this matrix
    //lets us know if we need to recompute the matrix 

    DegreeOfFreedom[] DOFList;
    
    AARotamerType resAATypes[][];//Names of AA types for each residue



//    public CETMatrix(int numMutable, MutableResParams strandMut, RotamerSearch rs){
//
//        numRes = numMutable;
//        DOFList = rs.m.DOFs;
//
//        intraAndShellBounds = new ContETerm[numRes][][];
//        pairwiseBounds = new ContETerm[numRes][][][][][];
//        
//        resAATypes = new AARotamerType[numRes][];
//
//        //First residue number within active site 
//        for(int p1=0; p1<strandMut.allMut.length;p1++){
//        	
//        	int str1 = strandMut.resStrand[p1];
//        	int strResNum1 = strandMut.resStrandNum[p1];
//
//                intraAndShellBounds[p1] = new ContETerm[rs.strandRot[str1].rl.getNumAAallowed()][];
//                pairwiseBounds[p1] = new ContETerm[rs.strandRot[str1].rl.getNumAAallowed()][][][][];
//                
//                resAATypes[p1] = rs.strandRot[str1].rl.getAAtypesAllowed();
//
//                for (int a1=0; a1<rs.strandRot[str1].getNumAllowable(strResNum1); a1++){
//                    int curAAind1 = rs.strandRot[str1].getIndexOfNthAllowable(strResNum1,a1);
//                    int numRot1 = rs.getNumRot(str1, strResNum1, curAAind1);
//
//                    intraAndShellBounds[p1][curAAind1] = new ContETerm[numRot1];
//                    pairwiseBounds[p1][curAAind1] = new ContETerm[numRot1][][][];
//
//                    for (int r1=0; r1<numRot1; r1++){
//                        pairwiseBounds[p1][curAAind1][r1] = new ContETerm[numRes][][];
//
//                        
//                        for (int p2=0;p2<strandMut.allMut.length;p2++){
//                        	int str2 = strandMut.resStrand[p2];
//                        	int strResNum2 = strandMut.resStrandNum[p2];
//
//                                        pairwiseBounds[p1][curAAind1][r1][p2] = new ContETerm[rs.strandRot[str2].rl.getNumAAallowed()][];
//                                        for (int a2=0; a2<rs.strandRot[str2].getNumAllowable(strResNum2); a2++){
//
//                                            int curAAind2 = rs.strandRot[str2].getIndexOfNthAllowable(strResNum2,a2);
//                                            int numRot2 = rs.getNumRot( str2, strResNum2, curAAind2 );
//                                            pairwiseBounds[p1][curAAind1][r1][p2][curAAind2] = new ContETerm[numRot2];
//                                        }
//                                
//                                
//                        }
//                    }
//                }
//
//            }
//    }

    public CETMatrix(Emat emat, Molecule m){

    	numRes = emat.resByPos.size();
        DOFList = m.DOFs;
        supRot = emat.singles.supRot;
    	
		intraAndShellBounds = new ContETerm[emat.singles.E.length][][];
        pairwiseBounds = new ContETerm[emat.singles.E.length][][][][][];

		//Runtime runtime = Runtime.getRuntime();
		for(int p1=0; p1<emat.singles.E.length;p1++){
			int[] p1ind = {p1};
			intraAndShellBounds[p1] = new ContETerm[emat.singles.E[p1].length][];
			pairwiseBounds[p1] = new ContETerm[emat.singles.E[p1].length][][][][];

			for(int a1=0;a1<emat.singles.E[p1].length;a1++){
				int[] p1a1ind = {p1,a1};
				intraAndShellBounds[p1][a1] = new ContETerm[emat.singles.E[p1][a1].length];
				pairwiseBounds[p1][a1] = new ContETerm[emat.singles.E[p1][a1].length][][][];
				for(int r1=0; r1<emat.singles.E[p1][a1].length;r1++){
					int[] p1a1r1ind = {p1,a1,r1};
					pairwiseBounds[p1][a1][r1] = new ContETerm[emat.singles.E.length][][];
					for(int p2=0; p2<emat.singles.E.length;p2++){
						if(p1!=p2 && emat.areNeighbors(p1,p2,m)  ){
							int[] p1a1r1p2ind = {p1,a1,r1,p2};
							pairwiseBounds[p1][a1][r1][p2] = new ContETerm[emat.singles.E[p2].length][];
							for(int a2=0;a2<emat.singles.E[p2].length;a2++){
								int[] p1a1r1p2a2ind = {p1,a1,r1,p2,a2};
								pairwiseBounds[p1][a1][r1][p2][a2] = new ContETerm[emat.singles.E[p2][a2].length];
							}
						}
					}
				}
			}
		}
	}



    

    public void mergeIn(CETMatrix M, int residueMutatable[]){
        //merge in bounds from a matrix M from a slave node,
        //which has the pairwise or intra+shell bounds for
        //the residue(s) indicated by 1s in residueMutatable


        int res1 = -1, res2 = -1;
        for(int i=numRes-1;i>=0;i--){//Go through the residues backwards to make sure res1>res2 (matching the matrix organization)
            //actually not needed now for LSBMatrix
                if (residueMutatable[i]==1) {
                        if (res1 == -1)
                                res1 = i;
                        else
                                res2 = i;
                }
        }


        if(res1 == -1)//Template run: nothing to merge in
            return;
        else if(res2 == -1) { //We only compute the CETs that need to be computed
        	
        		for(int a1=0; a1<intraAndShellBounds[res1].length;a1++){
        			for(int r1=0; r1<intraAndShellBounds[res1][a1].length;r1++){
        				if(M.intraAndShellBounds[res1][a1][r1] != null)
        					intraAndShellBounds[res1][a1][r1] = M.intraAndShellBounds[res1][a1][r1];
        			}
        		}
//            intraAndShellBounds[res1] = M.intraAndShellBounds[res1];
        }
        else{

        	for(int a1=0; a1<intraAndShellBounds[res1].length;a1++){
    			for(int r1=0; r1<intraAndShellBounds[res1][a1].length;r1++){
    				for(int a2=0; a2<intraAndShellBounds[res2].length;a2++){
    	    			for(int r2=0; r2<intraAndShellBounds[res2][a2].length;r2++){
    	    				if(M.pairwiseBounds[res1][a1][r1][res2][a2][r2] != null){
    	    					pairwiseBounds[res1][a1][r1][res2][a2][r2] = M.pairwiseBounds[res1][a1][r1][res2][a2][r2];
    	    					pairwiseBounds[res2][a2][r2][res1][a1][r1] = M.pairwiseBounds[res1][a1][r1][res2][a2][r2];
    	    				}
    	    			}
    	    		}
    			}
    		}
        	
//            for(int AAind=0; AAind<pairwiseBounds[res1].length; AAind++){
//                if( pairwiseBounds[res1][AAind] != null ){
//                    for(int rot=0; rot<pairwiseBounds[res1][AAind].length; rot++){
//                        pairwiseBounds[res1][AAind][rot][res2] = M.pairwiseBounds[res1][AAind][rot][res2];
//                    }
//
//                }
//            }
//
//            for(int AAind=0; AAind<pairwiseBounds[res2].length; AAind++){
//                if( pairwiseBounds[res2][AAind] != null ){
//                    for(int rot=0; rot<pairwiseBounds[res2][AAind].length; rot++){
//                        pairwiseBounds[res2][AAind][rot][res1] = M.pairwiseBounds[res2][AAind][rot][res1];
//                    }
//
//                }
//            }
        }
    }


    public void setShellRotE(int res, int AAind, int rot, ContETerm b){
        tryMolecCompression(res,AAind,rot,b);
        intraAndShellBounds[res][AAind][rot] = b;
    }
    
    public void setShellRotE(int[] i, ContETerm b) {
    	tryMolecCompression(i[0],i[1],i[2],b);
        intraAndShellBounds[i[0]][i[1]][i[2]] = b;
	}

    public void setPairwiseE(int res1, int AAind1, int rot1, int res2, int AAind2,
            int rot2, ContETerm b){
        
        tryMolecCompression(res1,AAind1,rot1,res2,AAind2,rot2,b);

        pairwiseBounds[res1][AAind1][rot1][res2][AAind2][rot2] = b;
        pairwiseBounds[res2][AAind2][rot2][res1][AAind1][rot1] = b;
    }
    
    public void setPairwiseE(int[] i, ContETerm b){
        
        tryMolecCompression(i[0],i[1],i[2],i[3],i[4],i[5],b);

        pairwiseBounds[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = b;
        pairwiseBounds[i[3]][i[4]][i[5]][i[0]][i[1]][i[2]] = b;
    }
    
    
    //tryMolecCompression: these methods are used when adding a new series b w/ SVE to the mtarix
    //we try to substitute in a molecule that is already in the matrix
    //to avoid storing a new one in b
    //we can do this iff there is another bound with SVE at the same residue with the same 
    //discreteDOFVals (including amino acid type)
    //since then setting continuous DOFs (during ContETerm evaluation) will make the relevant residues
    //of the molecules identical
    void tryMolecCompression(int res, int AAind, int rot, ContETerm b){
        if(b instanceof EPoly){
            EPoly ep = (EPoly)b;
            
            if(ep.sve != null){
                
                for(ContETerm b2 : intraAndShellBounds[res][AAind]){//can only share m with other residues at the same (res, AA type)
                    if(b2 instanceof EPoly){//and thus not null...
                        EPoly ep2 = (EPoly)b2;
                        
                        if(ep2.sve!=null){
                            
                            boolean match = true;
                            //since res and AA type match, b and b2 will have the same discrete DOFs
                            //check that the values match
                            for(int dd=0; dd<b.discreteDOFVals.size(); dd++){
                                
                                double dd1 = b.discreteDOFVals.get(dd);//want to compare doubles, not Doubles in the array list...need numerical comparison
                                double dd2 = b2.discreteDOFVals.get(dd);
                                
                                if(dd1!=dd2){
                                    match = false;
                                    break;
                                }
                            }
                            
                            if(match){
                                ep.sve.m = ep2.sve.m;
                                ep.sveOF.m = ep2.sve.m;
                                return;//succesfully compressed
                            }
                        }
                    }
                }
            }
        }
    }


    //pairwise version
    void tryMolecCompression(int res1, int AAind1, int rot1, int res2, int AAind2,
            int rot2, ContETerm b){
        
        if(b instanceof EPoly){
            EPoly ep = (EPoly)b;
            
            if(ep.sve != null){
                
                for(ContETerm c[][][] : pairwiseBounds[res1][AAind1]){//can only share m with other residues at the same (res, AA type)
                    if(c!=null){
                        for(ContETerm b2 : c[res2][AAind2]){
                            
                            if(b2 instanceof EPoly){//and thus not null...
                                EPoly ep2 = (EPoly)b2;

                                if(ep2.sve!=null){
                                    boolean match = true;
                                    //since res and AA type match, b and b2 will have the same discrete DOFs
                                    //check that the values match
                                    for(int dd=0; dd<b.discreteDOFVals.size(); dd++){
                                        double dd1 = b.discreteDOFVals.get(dd);//want to compare doubles, not Doubles in the array list...need numerical comparison
                                        double dd2 = b2.discreteDOFVals.get(dd);

                                        if(dd1!=dd2){
                                            match = false;
                                            break;
                                        }
                                    }
                                    if(match){
                                        ep.sve.m = ep2.sve.m;
                                        ep.sveOF.m = ep2.sve.m;
                                        return;//succesfully compressed
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    public ContETerm[] getCETList(Index3[] i){
        //given the AAs and rotamers assigned so far, get the list of LSBs
        //we'll use to compute the A* lower bound

        int numTerms = (i.length+1)*i.length/2;

        ContETerm terms[] = new ContETerm[numTerms];

        int termCount = 0;

        for(int res=0; res<i.length; res++){
            //Note: this works if all residues are assigned
            //or just the first several

            terms[termCount] = intraAndShellBounds[i[res].pos][i[res].aa][i[res].rot];
            termCount++;

            for(int res2=0; res2<res; res2++){
            	if(!areNeighbors(i[res].pos,i[res2].pos))
            		continue;
                terms[termCount] = pairwiseBounds[i[res].pos][i[res].aa][i[res].rot][i[res2].pos][i[res2].aa][i[res2].rot];
                termCount++;
            }
        }


        return terms;
    }
    
    
    
    /**
	 * Determines if two positions in the Emat are neighbors
	 * For now I just check if the pair matrix is defined for
	 * these two positions
	 * @param pos1
	 * @param pos2
	 * @return areNeighbors
	 */
	public boolean areNeighbors(int pos1, int pos2) {
		//Find a valid rotamer for pos1 and then check
		//if pos2 is defined for the rotamer
		for(int a1=0; a1<pairwiseBounds[pos1].length;a1++){
			for(int r1=0; r1<pairwiseBounds[pos1][a1].length;r1++){
				if(pairwiseBounds[pos1][a1][r1][pos2] == null)
					return false;
				else
					return true;
			}
		}

		return true;
	}
    
    
    
    public CETObjFunction getObjFunc( Index3[] indices, boolean includeMinE, boolean rcs, ContSCObjFunction sveOF){
        //Get the objective function for the given rotamers
        //if !rcs: given as (AA#, rot#) assignments for each residue (can be -1)
        //if rcs: given as RC set assignments, in AStarAxe.RCSets numbering, stored in rot (AAIndices can be null)
        //if sveOF isn't null we can use it to assign DOFs for any SVEs this objective function may call
        
        ContETerm terms[] = getCETList(indices);

        //we'll minimize over all DOFs that the terms depend on
        //constraints will be taken from DOFmin, DOFmax in the bounds
        //when getting constraints all bounds should agree (either using
        //single-voxel or all-unpruned voxels range of values for each DOF...)
        HashMap<Integer,double[]> DOFMap = new HashMap<Integer,double[]>();//map DOFNum-->constraints
        for(ContETerm b : terms){

            if(b!=null){

                for(int dof=0; dof<b.numDOFs; dof++){
                    double constr[] = new double[] {b.DOFmin.get(dof),b.DOFmax.get(dof)};
                    int DOFNum = b.DOFNums[dof];

                    //check against any previous bounds' constraints on same DOF...
                    if(DOFMap.containsKey(DOFNum)){
                        double prevConstr[] = DOFMap.get(DOFNum);
                        prevConstr[0] = Math.max(constr[0],prevConstr[0]);
                        prevConstr[1] = Math.min(constr[1],prevConstr[1]);
                        //Each assigned RC restricts the relevant DOFs to
                        //the range for that RC, and removes RCs at other positions
                        //in the case of pruned pairs.  So all the DOF upper and lower bounds
                        //in terms apply to the conformations set we're minimizing over--
                        //thus we use the intersection of the ranges
                        //that different terms may provide for our DOF
                        
                        if(prevConstr[0] > prevConstr[1]){
                            System.err.println("ERROR: inconsistent DOF constraints for objective function! DOFNum: "
                                    +DOFNum+" prevConstr: "+prevConstr[0]+" "+prevConstr[1]+" new constr: "+
                                    constr[0]+" "+constr[1]);
                        }
                    }
                    else
                        DOFMap.put(DOFNum, constr);
                }
            }
        }

        int numDOFs = DOFMap.size();

        int DOFNums[] = new int[numDOFs];
        double constraints[][] = new double[2][numDOFs];
        int count = 0;

        for(int DOFNum : DOFMap.keySet()){
            DOFNums[count] = DOFNum;
            double constr[] = DOFMap.get(DOFNum);
            constraints[0][count] = constr[0];
            constraints[1][count] = constr[1];
            count++;
        }

        DoubleMatrix1D DOFVals = DoubleFactory1D.dense.make(numDOFs);
        CETEnergyFunction ef = new CETEnergyFunction( DOFVals, DOFNums, terms, includeMinE );
        return new CETObjFunction(DOFNums,constraints,DOFList,ef,sveOF);
    }

    
    
    public void cleanupSVE(){
        //Whenever two ContETerms have the same amino-acid types,
        //AND any discreteDOFVals agree,
        //we can use the same molecule for sve (since the relevant res of the molecule will be exactly the same
        //once the continuous DOFs are set)
        //so let's link all the redundant molecules
        for(ContETerm[][] b : intraAndShellBounds){
            if(b!=null)
            for(ContETerm[] bb : b){
                if(bb!=null){
                    
                    //Molecule mutMolec = null;//molecule with given mutations
                    //should be the same for all bbb in bb
                    //WAIT have to account for discreteDOFVals
                    
                    ArrayList<Molecule> mutMolec = new ArrayList<>();
                    ArrayList<double[]> mutMolecDD = new ArrayList<>();//discrete DOFVals for each mutMolec
                    
                    for(ContETerm bbb : bb){
                        if(bbb instanceof EPoly){//not null and may have sve
                            EPoly ep = (EPoly)bbb;
                            if(ep.sve!=null){
                                
                                /*if(mutMolec==null)
                                    mutMolec = ep.sve.m;
                                else{
                                    ep.sve.m = mutMolec;
                                    ep.sveOF.m = mutMolec;
                                }*/
                                boolean alreadyStored = false;
                                for(int m=0; m<mutMolec.size(); m++){
                                    double[] ddVal = mutMolecDD.get(m);
                                    boolean match = true;
                                    for(int dd=0; dd<ddVal.length; dd++)//discreteDOFVals must match
                                        //if AA types, res match then there'll be the same number of discrete DOFs
                                        match = match && (ep.discreteDOFVals.get(dd)==ddVal[dd]);

                                    if(match){
                                        ep.sve.m = mutMolec.get(m);
                                        ep.sveOF.m = ep.sve.m;
                                        alreadyStored = true;
                                        break;
                                    }
                                }

                                if(!alreadyStored){
                                    mutMolec.add(ep.sve.m);
                                    double ddVals[] = new double[ep.discreteDOFVals.size()];
                                    for(int dd=0; dd<ddVals.length; dd++)
                                        ddVals[dd] = ep.discreteDOFVals.get(dd);
                                    mutMolecDD.add(ddVals);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        for(ContETerm[][][][][] b : pairwiseBounds){
            if(b!=null)
            for(ContETerm[][][][] bb : b){
                if(bb!=null){
                    
                    //to reuse molecule, first & second res & AA must match
                    //HashMap<Integer,Molecule> mutMolec = new HashMap<>();
                    //key for mutMolec is res2 + numRes*curAA2
                    //which is OK because res2 can't exceed numRes
                    
                    //We need a mutMolec for each (res2, curAA2, any discreteDOFVals)
                    ArrayList<Molecule> mutMolec = new ArrayList<>();
                    ArrayList<double[]> mutMolecInfo = new ArrayList<>();
                    //for each molecule in mutMolec, (res2, curAA2, any discreteDOFVals)
                    
                    for(ContETerm bbb[][][] : bb){
                        if(bbb!=null)
                        for(int res2=0; res2<bbb.length; res2++){
                            ContETerm bbbb[][] = bbb[res2];
                            if(bbbb!=null)
                            for(int curAA2=0; curAA2<bbbb.length; curAA2++){
                                ContETerm bbbbb[] = bbbb[curAA2];
                                if(bbbbb!=null)
                                for(ContETerm bbbbbb: bbbbb){
                                    if(bbbbbb instanceof EPoly){//not null and may have sve
                                        EPoly ep = (EPoly)bbbbbb;
                                        if(ep.sve!=null){
                                            
                                            /*if(mutMolec.get(res2 + numRes*curAA2)==null)
                                                mutMolec.put(res2 + numRes*curAA2, ep.sve.m);
                                            else{
                                                ep.sve.m = mutMolec.get(res2 + numRes*curAA2);
                                                ep.sveOF.m = ep.sve.m;
                                            }*/
                                            boolean alreadyStored = false;
                                            for(int m=0; m<mutMolec.size(); m++){
                                                double[] info = mutMolecInfo.get(m);
                                                boolean match = (info[0]==res2) && (info[1]==curAA2);
                                                for(int dd=2; dd<info.length; dd++)//discreteDOFVals must match
                                                    //if AA types, res match then there'll be the same number of discrete DOFs
                                                    match = match && (ep.discreteDOFVals.get(dd-2)==info[dd]);
                                                
                                                if(match){
                                                    ep.sve.m = mutMolec.get(m);
                                                    ep.sveOF.m = ep.sve.m;
                                                    alreadyStored = true;
                                                    break;
                                                }
                                            }
                                            
                                            if(!alreadyStored){
                                                mutMolec.add(ep.sve.m);
                                                double info[] = new double[ep.discreteDOFVals.size()+2];
                                                info[0] = res2;
                                                info[1] = curAA2;
                                                for(int dd=2; dd<info.length; dd++)
                                                    info[dd] = ep.discreteDOFVals.get(dd-2);
                                                mutMolecInfo.add(info);
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
    }
    

    
    static double squareDist(DoubleMatrix1D a, DoubleMatrix1D b){
        //square of distance between two vectors
        return a.zDotProduct(a) - 2*a.zDotProduct(b) + b.zDotProduct(b);
    }
    
    
    
    void analyzeFitTypes(){
        //Look at what types of fits are in this matrix
        System.out.println("Fit types in continuous energy term matrix (# of terms, description): ");
        
        HashMap<String,Integer> typeCounts = new HashMap<>();
        
        for(ContETerm c[][] : intraAndShellBounds){
            for(ContETerm cc[] : c){
                if(cc!=null){
                    for(ContETerm ccc : cc){
                        if(ccc!=null)
                            countFitType(ccc.fitDescription,typeCounts);
                    }
                }
            }
        }
        
        for(int res=0; res<numRes; res++){
            for(ContETerm cc[][][][] : pairwiseBounds[res]){
                if(cc!=null){
                    for(ContETerm ccc[][][] : cc){
                        if(ccc!=null){
                            for(int res2=0; res2<res; res2++){
                            	if(!areNeighbors(res, res2))
                            		continue;
                                for(ContETerm cccc[] : ccc[res2]){
                                    if(cccc!=null){
                                        for(ContETerm ccccc : cccc){
                                            if(ccccc!=null)
                                                countFitType(ccccc.fitDescription,typeCounts);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        for(String descr : typeCounts.keySet()){
            System.out.println(typeCounts.get(descr)+" "+descr);
        }
    }
    
    
    void countFitType(String descr, HashMap<String,Integer> typeCounts){
        //given the description of a fit, add it to the counts of descriptions in typeCounts
        if(typeCounts.containsKey(descr))
            typeCounts.put(descr, typeCounts.get(descr)+1);
        else
            typeCounts.put(descr, 1);
    }




    /**
     * Copy over the terms from an older CETMatrix 
     * @param cetm An old CET matrix that already has some CET entries computed 
     * @param emat Emat that was used to create the current CETMatrix
     */
	public void copyOverTerms(CETMatrix cetm) {
		//We copy over a term when we find two supRot entries that match between
		//the old cetm and the new emat
		
		//Loop through the new matrix since every term in the old matrix should be in
		//the new matrix
		for(int p1=0; p1<intraAndShellBounds.length;p1++ ){
			for(int a1=0; a1<intraAndShellBounds[p1].length;a1++ ){
				int oldR1 = 0;
				for(int r1=0; r1<intraAndShellBounds[p1][a1].length;r1++ ){
					if(oldR1 < cetm.supRot[p1][a1].length && supRot[p1][a1][r1][0] == cetm.supRot[p1][a1][oldR1][0]){ //0 index assumes no super rotamers;
						intraAndShellBounds[p1][a1][r1] = cetm.intraAndShellBounds[p1][a1][oldR1];
						
						for(int p2=p1+1; p2<intraAndShellBounds.length;p2++ ){
							if(!areNeighbors(p1, p2))
								continue;
							for(int a2=0; a2<intraAndShellBounds[p2].length;a2++ ){
								int oldR2 = 0;
								for(int r2=0; r2<intraAndShellBounds[p2][a2].length;r2++ ){
									if(oldR2 < cetm.supRot[p2][a2].length && supRot[p2][a2][r2][0] == cetm.supRot[p2][a2][oldR2][0]){
										pairwiseBounds[p1][a1][r1][p2][a2][r2] = cetm.pairwiseBounds[p1][a1][oldR1][p2][a2][oldR2];
										pairwiseBounds[p2][a2][r2][p1][a1][r1] = cetm.pairwiseBounds[p2][a2][oldR2][p1][a1][oldR1];
										oldR2++;
									}
								}
							}
						}
						oldR1++;
					}
				}
			}
		}
		
		
		
	}

	/**
     * Copy over all terms from an older CETMatrix 
     * @param cetm An larger CET matrix that has all CET entries computed for the current CET 
     * @param emat Emat that was used to create the current CETMatrix
     */
	public void copyAllTerms(CETMatrix cetm) {
		//We copy over a term when we find two supRot entries that match between
		//the old cetm and the new emat
		
		//Loop through the old (larger) matrix since every term in the new matrix should be in
		//the old matrix
		for(int p1=0; p1<cetm.intraAndShellBounds.length;p1++ ){
			for(int a1=0; a1<cetm.intraAndShellBounds[p1].length;a1++ ){
				int newR1 = 0;
				for(int r1=0; r1<cetm.intraAndShellBounds[p1][a1].length;r1++ ){
					if(newR1 < supRot[p1][a1].length && cetm.supRot[p1][a1][r1][0] == supRot[p1][a1][newR1][0]){ //0 index assumes no super rotamers;
						intraAndShellBounds[p1][a1][newR1] = cetm.intraAndShellBounds[p1][a1][r1];
						
						for(int p2=p1+1; p2<cetm.intraAndShellBounds.length;p2++ ){
							for(int a2=0; a2<cetm.intraAndShellBounds[p2].length;a2++ ){
								int newR2 = 0;
								for(int r2=0; r2<cetm.intraAndShellBounds[p2][a2].length;r2++ ){
									if(newR2 < supRot[p2][a2].length && cetm.supRot[p2][a2][r2][0] == supRot[p2][a2][newR2][0]){
										pairwiseBounds[p1][a1][newR1][p2][a2][newR2] = cetm.pairwiseBounds[p1][a1][r1][p2][a2][r2];
										pairwiseBounds[p2][a2][newR2][p1][a1][newR1] = cetm.pairwiseBounds[p2][a2][r2][p1][a1][r1];
										newR2++;
									}
								}
							}
						}
						newR1++;
					}
				}
			}
		}
		
		
		
	}


	/** 
	 * Set the runParams for each OneMutation to only compute Ematrix
	 * entries for RCs that don't have computed CETs
	 * @param mutArray OneMutation array that holds MPI information for Ematrix computation
	 */
	public void updateMutArray(OneMutation[] mutArray) {
		for(OneMutation mut: mutArray){
			if(mut.runParams != null){ //Will be null for TEMPL and INTRA runs
				mut.runParams.rotamers = new ArrayList<Index3>();
				for(int p1=0; p1<mut.resMut.length;p1++){
					if(mut.resMut[p1] > 0){
						for(int a1=0; a1<intraAndShellBounds[p1].length;a1++){
							for(int r1=0; r1<intraAndShellBounds[p1][a1].length;r1++){
								if(intraAndShellBounds[p1][a1][r1] == null)
									mut.runParams.rotamers.add(new Index3(p1,a1,r1));
							}
						}
					}
				}
			}
		}
		
	}






	


    
}