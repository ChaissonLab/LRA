#ifndef _REFINE_Breakpoint_H_
#define _REFINE_Breakpoint_H_
#include <vector>
#include "AlignmentBlock.h"

void PrependBlocks(vector<Block> &src, vector<Block> &dest) {
  if (src.size() == 0) { return;}
  if (dest.size() == 0) { dest = src; return;}
  
  int last=src.size()-1;
  // If the last match is a gapless extension, that can mess up some cigar parsing.
  if (src[last].tPos + src[last].length == dest[0].tPos and
      src[last].qPos + src[last].length == dest[0].qPos) {
    dest[0].tPos -= src[last].length;
    dest[0].qPos -= src[last].length;
    dest[0].length += src[last].length;
    src.resize(last);
  }
  
  dest.resize(dest.size()+ src.size());
  for (int i=dest.size(); i > src.size(); i--) {
    dest[i-1] = dest[i-src.size()-1];
  }
  for (int i=0; i < src.size(); i++) {
    dest[i] = src[i];
  }
}

void AppendBlocks(vector<Block> &src, vector<Block> &dest) {
  if (src.size() == 0) { return;}
  if (dest.size() == 0) { dest = src; return;}
  
  int last=dest.size()-1;
  // If the last match is a gapless extension, that can mess up some cigar parsing.
  int srcStart=0;
  if (dest[last].tPos + dest[last].length == src[0].tPos and
      dest[last].qPos + dest[last].length == src[0].qPos) {
    dest[last].length += src[0].length;
    srcStart=1;
  }
  int destEnd=dest.size();
  dest.resize(dest.size() + src.size() - srcStart);
  for (int i=srcStart; i < src.size(); i++) {
    dest[destEnd] = src[i];
    destEnd++;
  }
}
  

void PathToBlocks(vector<int> &path, vector<Block> &blocks) {
  int i=0;
  int left=1,down=2,diag=3;
  int q=0, t=0;
  while (i < path.size() and path[i] != diag and (path[i] == left or path[i] == down)) {
    if (path[i] == left) { q++;}
    if (path[i] == down) { t++;}
    i++;
  }
  
  while (i < path.size()) {
    int ml=0;
    int qs=q;
    int ts=t;
    while (i < path.size() and path[i] == diag) {
      ml++;
      q++;
      t++;
      i++;
    }

    while (i < path.size() and (path[i] == left or path[i] == down)) {
      if (path[i] == left) { q++;}
      if (path[i] == down) { t++;}
      i++;
    }
    int match=min(q-qs, t-ts);
    if (match > 0) {
      blocks.push_back(Block(qs,ts,match));
    }
  }
}

void PrintMat(vector<int> &mat, int r) {
  for (int i=0; i < mat.size(); i++) {
    cout << "\t" << mat[i];
    if ((i+1)%r == 0) { cout << endl;}
  }
}

void TraceBack(vector<int> &path, int q, int t, int r, vector<int> &tb) {
  q++;
  t++;
  int i=t*r+q;
  int left=1;
  int down=2;
  int diag=3;
  while (q >0 or t > 0) {
    if (path[i] == diag) {
      q--;
      t--;
      tb.push_back(diag);
      
    }
    if (path[i] == left) {
      q--;
      tb.push_back(left);
    }
    if (path[i] == down) {
      t--;
      tb.push_back(down);
    }
    i=(t)*r+(q);
  }
  reverse(tb.begin(), tb.end());
}

void StoreQScoreVect(vector<int> &score, vector<int> &path, int q, int t, int r, vector<int> &qv, vector<int> &index) {
  qv.resize(r-1);
  index.resize(r-1);
  fill(qv.begin(), qv.end(), 0);
  int i=(t+1)*r+q+1;
  int left=1;
  int down=2;
  int diag=3;
  // index into matrix is up by 1
  q++;
  t++;
  while (i > 0) {
    if (path[i] == diag or path[i] == left) {
      assert(q > 0);
      qv[q-1] = score[i];
      index[q-1] = i;
    }
    if (path[i] == diag) {
      q--;
      t--;
    }
    if (path[i] == left) {
      q--;
    }
    if (path[i] == down) {
      t--;
    }
    i=(t)*r+(q);
  }
}
    
    

void RSdp(string &q, string &t, vector<int> &path, vector<int> &score,
	  int mat, int mis, int indel) {
  
  path.resize((q.size()+1)*(t.size()+1), -1);
  score.resize((q.size()+1)*(t.size()+1),0);
  int qs=q.size();
  int ts=t.size();
  int left=1;
  int down=2;
  int diag=3;
  int row=qs+1;
  for (int i=1; i < qs+1;i++) {
    path[i] = left;
    score[i] = score[i-1] + indel;
  }
  for (int i=1; i < ts+1; i++) {
    path[row*i] = down;
    score[i*row] = score[(i-1)*row] + indel;
  }
  for (int i=0; i < ts; i++) {
    for (int j=0; j < qs; j++) {
      int diagScore = score[i*row+j];
      if (q[j] == t[i]) {
	diagScore += mat;
      }
      else {
	diagScore += mis;
      }
      int leftScore = score[(i+1)*row + j] + indel;
      int downScore = score[i*row + (j+1)] + indel;
      int maxScore =max(diagScore, max(leftScore, downScore));
      score[(i+1)*row+(j+1)] = maxScore;
      if (maxScore == diagScore) {
	path[(i+1)*row+(j+1)] = diag;
      }
      else if (maxScore == leftScore) {
	path[(i+1)*row+(j+1)] = left;
      }
      else {
	path[(i+1)*row+(j+1)] = down;
      }
    }
  }
  //
  // No need for trace back here.
  //
}

int FindMax(vector<int> &score, int row, int &q, int &t) {
  if (score.size() == 0) {
    q=t=0;
    return 0;
  }
  vector<int>::iterator itr = max_element(score.begin(), score.end());
  int index=itr-score.begin();
  t=index/row-1;
  q=index%row-1;
  return score[index];
}

void RefineBreakpoint(Read &read,
		      Genome &genome,
		      Alignment &leftAln,
		      Alignment &rightAln,
		      const Options &opts) {
  int lqs=0, lqe=0, lts=0, lte=0;
  int rqs=0, rqe=0, rts=0, rte=0;

  lqs=leftAln.GetQStart();
  lqe=leftAln.GetQEnd();
  lts=leftAln.GetTStart();
  lte=leftAln.GetTEnd();
  int flqs=0, flqe=0, rlqs=0, rlqe=0;
  if (leftAln.strand == 0) {
    flqs=lqs;
    flqe=lqe;
  }
  else {
    flqs = read.length - lqe;
    flqe = read.length - lqs;
  }

  rqs=rightAln.GetQStart();
  rqe=rightAln.GetQEnd();
  rts=rightAln.GetTStart();
  rte=rightAln.GetTEnd();
  int frqs=0, frqe=0;
  if (rightAln.strand == 0) {
    frqs=rqs;
    frqe=rqe;
  }
  else {
    frqs = read.length - rqe;
    frqe = read.length - rqs;
  }
  /*  
  cerr << "LEFT\t" << (int) leftAln.strand  << "\t" << leftAln.blocks.size() << "\t" << flqs << "\t" << flqe << "\t"
       << "RIGHT\t" << (int) rightAln.strand << "\t" << rightAln.blocks.size() << "\t" << frqs << "\t" << frqe << endl;
  */
  int MAX_GAP=500;
  
  if (frqs > flqe and frqs - flqe < MAX_GAP) {
    //
    // The two endpoints are close enough to attempt to refine the breakpoint.
    //
    // For now just do standard dp, speed up later if need be.
    int span  = frqs - flqe;
    //
    // Determine spans that are refined for left.
    //
    string lqString, ltString;
    char *tChrom=genome.seqs[leftAln.chromIndex];
    bool lPrefixExtend=false;
    int lqExtStart, lqExtEnd;
    int ltExtStart, ltExtEnd;
    
    if (leftAln.strand == 0) {
      //
      // Missed segment is at right side of left aln
      //
      lqExtStart = lqe;
      lqExtEnd   = lqe+span;
      lqString   = string(leftAln.read + lqe, span);      

      //
      // Left align is forward strand. Refining will go from the end of alignment on
      //
      ltExtStart = leftAln.GetTEnd();
      int tSpan  = min(genome.lengths[leftAln.chromIndex]-leftAln.GetTEnd(), span);
      ltExtEnd   = ltExtStart+tSpan;
      ltString   = string(tChrom+ltExtStart, tSpan);
    }
    else {
      assert(lqs-span > 0);
      // missed segment is at left side of rev strand
      lqExtStart = lqs-span;
      lqExtEnd   = lqs;
      lqString   = string(leftAln.read + lqExtStart, span);

      //
      // Left align is reverse strand. The gap goes forward in the read, which means it exends back in target
      ltExtEnd   = leftAln.GetTStart();
      ltExtStart = max(0,ltExtEnd-span);
      ltString   = string(tChrom+ltExtStart, ltExtEnd-ltExtStart);
      lPrefixExtend=true;
      reverse(lqString.begin(), lqString.end());
      reverse(ltString.begin(), ltString.end());
    }
        
    vector<int> lPath, lScore, rPath, rScore;
    int mat=2;
    int mis=-2;
    int gap=-4;
    /*
    cerr << "query" << endl;
    cerr << lqString << endl;
    cerr << "target" << endl;
    cerr << ltString << endl;
    */
    RSdp(lqString, ltString, lPath, lScore, mat,mis,gap);

    //
    // Right-hand side logic is the reverse.
    //
    string rqString, rtString;
    string rqStringCopy, rtStringCopy;
    char *rtChrom=genome.seqs[rightAln.chromIndex];
    bool rPrefixExtend=false;
    int rqExtStart, rqExtEnd;
    int rtExtStart, rtExtEnd;
    if (rightAln.strand == 0) {
      //
      
      rqExtStart = rqs - span;
      rqExtEnd   = rqs;
      rqString=string(rightAln.read + rqExtStart, span);

      int rtSpan = min(rightAln.GetTStart(), span);
      rtExtStart = rightAln.GetTStart() - rtSpan;
      rtExtEnd   = rightAln.GetTStart();
      rtString   = string(rtChrom+rtExtStart, rtSpan);
      reverse(rqString.begin(), rqString.end());
      reverse(rtString.begin(), rtString.end());
      rPrefixExtend=true;
    }
    else {
      rqExtStart = rightAln.GetQEnd();
      rqExtEnd   = rqExtStart+span;
      assert(rqExtStart+span <= read.length);
      rqString   = string(rightAln.read + rqExtStart, span);

      rtExtStart = rightAln.GetTEnd();
      int tSpan  = span;
      if (rtExtStart+span >= genome.lengths[rightAln.chromIndex]) {
	tSpan = genome.lengths[rightAln.chromIndex] - rtExtStart;
      }
      rtExtEnd   = rtExtStart+tSpan;
      rtString   = string(rtChrom+rtExtStart, tSpan);
    }      
    /*
    cerr << "rquery" << endl;
    cerr << rqString << endl;
    cerr << "rtarget" << endl;
    cerr << rtString << endl;
    */

    RSdp(rqString, rtString, rPath, rScore, mat,mis,gap);

    
    int mls, mlq, mlt, mrs, mrq, mrt;
    mls=FindMax(lScore, span+1, mlq, mlt);
    mrs=FindMax(rScore, span+1, mrq, mrt);
    //    cerr << "left " << mlq << "\t" << mlt << "\t" << mls << "\tright\t" << mrq << "\t" << mrt << "\t" << mrs << endl;
    // Now to merge the two results.
    //
    // Case 1, the local alignments do not overlap
    int maxLIndex=0;
    int maxRIndex=0;

    if (mlq < span - mrq ) {
      //      cerr << "Alignments do not overlap " << endl;
    }
    else {
      //      cerr << "Alignments do overlap. Optimize" << endl;
      vector<int> lqScores, rqScores, lqIndex, rqIndex;
      StoreQScoreVect(lScore, lPath, mlq, mlt, span+1, lqScores, lqIndex);
      StoreQScoreVect(rScore, rPath, mrq, mrt, span+1, rqScores, rqIndex);
      assert(lqScores.size() == rqScores.size());
      int maxScore=0;
      for (int i=0; i < lqScores.size(); i++ ) {
	if (lqScores[i] + rqScores[lqScores.size()-i-1] > maxScore) {
	  maxScore = lqScores[i] + rqScores[lqScores.size()-i-1];
	  maxLIndex=i;
	  maxRIndex=lqScores.size()-i-1;
	}
      }
      int maxLI=lqIndex[maxLIndex];
      int maxRI=rqIndex[maxRIndex];
      mlq=maxLIndex;
      mlt=maxLI/(span+1)-1;
      mrq=maxRIndex;
      mrt=maxRI/(span+1)-1;
      
    }
    vector<int> ltb, rtb;
    TraceBack(lPath, mlq, mlt, span+1, ltb);
    TraceBack(rPath, mrq, mrt, span+1, rtb);
    /*
    for (int i=0;i<ltb.size(); i++){ cerr << ltb[i];  }
    cerr << endl;
    for (int i=0;i<rtb.size(); i++){ cerr << rtb[i];  }
    cerr << endl;
    */
    //
    // Now append these to the alignments.
    //
    vector<Block> lBlocks, rBlocks;
    int lqBlockStart=-1, ltBlockStart=-1;
    if (lPrefixExtend == true) {
      reverse(ltb.begin(), ltb.end());
      reverse(lqString.begin(), lqString.end());
      reverse(ltString.begin(), ltString.end());
      int qs=lqExtEnd-1;
      int ts=ltExtEnd-1;
      for (int i=0; i < ltb.size(); i++) {
	if (ltb[i] == 3) { qs--; ts--;}
	if (ltb[i] == 2) { ts--;}
	if (ltb[i] == 1) { qs--;}
      }
      lqBlockStart=leftAln.GetQStart()-mlq-1;
      ltBlockStart=leftAln.GetTStart()-mlt-1;
      //      mlq=lqString.size() - mlq - 1;
      //      mlt=ltString.size() - mlt - 1;      
    }
    else {
      lqBlockStart=leftAln.GetQEnd();
      ltBlockStart=leftAln.GetTEnd();
    }
    
    PathToBlocks(ltb, lBlocks);
    for (int i=0; i < lBlocks.size(); i++) { lBlocks[i].qPos+=lqBlockStart; lBlocks[i].tPos+=ltBlockStart;}

    if (lPrefixExtend) {
      PrependBlocks(lBlocks, leftAln.blocks);
    }
    else {
      AppendBlocks(lBlocks, leftAln.blocks);
    }
    

    int rqBlockStart=-1, rtBlockStart=-1;
    if (rPrefixExtend == true) {
      reverse(rtb.begin(), rtb.end());
      rqBlockStart=rightAln.GetQStart()-mrq-1;
      rtBlockStart=rightAln.GetTStart()-mrt-1;
      //      mrq=rqString.size() - mrq - 1;
      //      mrt=rtString.size() - mrt - 1;    
    }
    else {
      rqBlockStart = rightAln.GetQEnd();
      rtBlockStart = rightAln.GetTEnd();
    }
    PathToBlocks(rtb, rBlocks);
    for (int i=0; i < rBlocks.size(); i++) { rBlocks[i].qPos+=rqBlockStart; rBlocks[i].tPos+=rtBlockStart;}


    if (rPrefixExtend) {
      PrependBlocks(rBlocks, rightAln.blocks);
    }
    else {
      AppendBlocks(rBlocks, rightAln.blocks);
    }
  }
  
}




#endif
