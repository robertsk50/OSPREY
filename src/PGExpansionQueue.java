import java.util.concurrent.PriorityBlockingQueue;

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
//	PGExpansionQueue.java
//
//	Version:           2.0
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

/**
 * This queue is ordered in terms of increasing f(n) values of the nodes in the A* expansion tree;
 * 		only the visible nodes are contained in the queue.
 */
public class PGExpansionQueue {
	
	private PriorityBlockingQueue<PGQueueNode> thequeue;
	long totalNodes = 0;
	
		//constructor
	PGExpansionQueue () {
		 thequeue = new PriorityBlockingQueue<PGQueueNode>(50000);
	}
	
	PGExpansionQueue (int i) {
		 thequeue = new PriorityBlockingQueue<PGQueueNode>(i);
	}
	
	 public void insert(PGQueueNode newNode){
         thequeue.add(newNode);
         totalNodes++;
	 }

	 public PGQueueNode getMin(){
         return thequeue.poll();
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

	public int numNodes() {
		return thequeue.size();
	}
}
