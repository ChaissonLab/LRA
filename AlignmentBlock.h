#ifndef ALIGNMENT_BLOCK_H_
#define ALIGNMENT_BLOCK_H_

#include <iostream>
#include <fstream>
#include "Types.h"

using namespace std;
class Block {
 public:
	//
	// An alignment is a collection of blocks. The qPos and tPos in a block
	// is relative to the beginning of the alignment rather than the
	// target or query.
	//
	 
	GenomePos qPos, tPos, length;
	Block() {qPos=tPos=length=0;}
  Block(GenomePos q, GenomePos t, GenomePos l) : qPos(q), tPos(t), length(l) {}
	Block& Assign(Block &rhs) {
		qPos = rhs.qPos;
		tPos = rhs.tPos;
		length = rhs.length;
		return *this;
	}

	GenomePos QEnd() {
		return qPos + length;
	}

	GenomePos TEnd() {
		return tPos + length;
	}

  void Clear() {
    qPos = tPos =  length = 0;
  }
};



#endif
