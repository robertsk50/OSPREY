/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.0
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
//	PGQueueNode.java
//
//	Version:           2.0
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;

import cern.colt.matrix.DoubleMatrix1D;

/**
 * Handles the data for a single node in the A* queue.
 *
 */
public class PGQueueNode implements Comparable,Serializable  {
	
	//the f(n) score associated with the current node
	public double fScore;
	
	//corresponding level
	public int level;
	
	//the numbers of the nodes in the considered conformation up to the current level
	public int confSoFar[];
	
	//the number of the corresponding node at that level
	public int nodeNum;
	
	//a pointer to the previous and next nodes in the expansion list
	//public PGQueueNode prevNode;
	//public PGQueueNode nextNode;
	
	// Residue positions that have not been assigned
	ArrayList<Integer> emptyLevels;
	ArrayList<Integer> nonEmptyLevels;
	
	LinkedList<EnergyTuple> curTuples = null;
	EMatrixEntryWIndex[] actualConf = null;
	
	 boolean FSTermIncluded;//There is a fit-series terms included in the fScore
     DoubleMatrix1D optFSPoint = null;//the values of the DOFs where the FS term is optimized
     //if the FS term hasn't been computed for this node yet, it is inherited from the parent

	
	PGQueueNode (int totalNumResidues, int curConf[], double fn, int newPos, int newRotAtPos) {

//		super(newRotAtPos,newPos,curConf,fn);
		
		nodeNum = newRotAtPos;
		level = newPos;
		
		
		confSoFar = new int[totalNumResidues];
		
		for (int i=0; i<totalNumResidues; i++){
			confSoFar[i] = curConf[i];
		}	
		
		confSoFar[newPos] = newRotAtPos;
		
		
		fScore = fn;
		
		//prevNode = null;
		//nextNode = null;
		
		emptyLevels = new ArrayList<Integer>();
		nonEmptyLevels = new ArrayList<Integer>();
		// Check which are the empty levels for this node.
		for(int i = 0; i < confSoFar.length; i++){
			if(confSoFar[i] == -1){
				emptyLevels.add(i);
			}
			else{
				nonEmptyLevels.add(i);
			}
		}		
	}
	
	
	public int compareTo(Object otherObject) throws ClassCastException{
        if(!(otherObject instanceof PGQueueNode)){
                throw new ClassCastException("A PGQueueNode object expected.");
        }
        PGQueueNode other = (PGQueueNode)otherObject;
        if(this.fScore > other.fScore){
                return 1;
        }
        else if(this.fScore < other.fScore){
                return -1;
        }
        else{ // Nodes have the same score
                if(this.level > other.level || this.nodeNum != other.nodeNum){
                        return -1; // but different levels, this is larger by default
                }
                else if(checkConf(this, other)){
                        return 0;
                }
                else{ // Two distinct nodes have the same fScore, say this one is larger.
                        return 1;
                }
        }
	}
	
	//Checks if the two given nodes have the same partially assigned conformation
	private boolean checkConf(PGQueueNode node1, PGQueueNode node2){
	
		if (node1.level!=node2.level) //different level
			return false;
		
		for (int l=0; l<node1.confSoFar.length; l++){
			if (node1.confSoFar[l]!=node2.confSoFar[l])
				return false;
		}
		
		//The partially assigned conformations are the same
		return true;		
	}
	
	public boolean empty(int level){
		if(confSoFar[level] == -1)
			return true;
		else
			return false;
	}
	
}
