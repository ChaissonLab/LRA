#ifndef TUPLE_OPS_H_
#define TUPLE_OPS_H_
#include <string>
#include "SeqUtils.h"
#include "Types.h"
#include <vector>
using namespace std;

Tuple mask=0;

template<typename Tup>
void TupleToString(Tup t, int k, string &s) {
	s="";
	for (int i = 0; i < k; i++) {
		s = s + binMap[(t & 3)];
		t >>=2;
	}
}

#define LOCAL_POS_BITS 12
class LocalTuple {
 public:
  uint32_t t: 32-LOCAL_POS_BITS;
  uint32_t pos: LOCAL_POS_BITS;
  static uint32_t for_mask_s;
  static uint32_t rev_mask_s;
	LocalTuple() {
		t=0;
		pos=0;
	}
	bool operator<(const LocalTuple &b) const {
		return t < b.t;
	}
	bool operator>(const LocalTuple &b) {
		return t > b.t;
	}
	friend int operator >(LocalTuple &a, LocalTuple &b) {
		return a.t > b.t;
	}
	friend int operator!=(LocalTuple &a, LocalTuple &b) {
		return a.t != b.t;
	}
	void ToString(int k, string &s) {
		TupleToString(t, k, s);
	}

};

uint32_t LocalTuple::for_mask_s=0;
uint32_t LocalTuple::rev_mask_s=0;

template<typename Tup> 
int DiagOffset(Tup &a, Tup &b) {
	unsigned long ad = a.first.pos - a.second.pos;
	unsigned long bd = b.first.pos - b.second.pos;
	int res = ad-bd;
	return res;
}

template<typename Tup> 
int DiagGap(Tup &a, Tup &b) {
	int firstGap  = b.first.pos - a.first.pos;
	int secondGap = b.second.pos - a.second.pos;
	int res= max(firstGap, secondGap) - min(firstGap, secondGap);
	return res;
}

class GenomeTuple {
public:
	Tuple t; // used to store kmers, this uses a 64 bit representation of k-mers (A=00, C=01, G=10, T=11)
	GenomePos pos; // unsigned 32 bit position in a genome (signed integers cannot index entire human genome!) 
	static Tuple for_mask_s;
	static Tuple rev_mask_s;
	GenomeTuple() {t = 0; pos = 0;} 
 	GenomeTuple(Tuple _t, GenomePos _p): t(_t), pos(_p) {}
	bool operator<(const GenomeTuple &b) const {
		return (t & GenomeTuple::for_mask_s) < (b.t & GenomeTuple::for_mask_s);
	}

	friend int operator > (GenomeTuple &a, GenomeTuple &b) {
		return (a.t & GenomeTuple::for_mask_s) > (b.t & GenomeTuple::for_mask_s);
	}
	friend int operator!= (GenomeTuple &a, GenomeTuple &b) {
		return (a.t & GenomeTuple::for_mask_s) != (b.t & GenomeTuple::for_mask_s);
	}
	void ToString(int k, string &s) {
		TupleToString(t, k, s);
	}
};

Tuple GenomeTuple::for_mask_s = 0;
Tuple GenomeTuple::rev_mask_s = 0;

template< typename Tup>
int InitMask(Tup &m, int k) {
	m.t=0;
	for (int i = 0; i < k; i++) {
		m.t<<=2;
		m.t+=3;
	}
	return m.t;
}

template<typename tup> void StoreTuple(char *seq, int pos, int k, tup &t) {
	t.t=0;
	unsigned int n;
	for (int p=pos; p <= pos+k-1 ; p++) {
		t.t <<=2;
		n = seqMap[seq[p]];
		t.t+=n;
	}
}

template<typename tup> void ShiftOne(char *seq, int pos, tup mask, tup &t) {
	t.t= (t.t << 2) & (Tuple) mask.t;
	t.t += (Tuple)seqMap[seq[pos]];	

}

template<typename tup> void ShiftOneRC(char *seq, int pos, int k, tup &t) {
	t.t >>= 2;
	unsigned int n=seqMap[seq[pos]];
	t.t += (~(n) & (Tuple)3) << (2*((Tuple)k-1));
}

template <typename tup> void TupleRC(tup a, tup &b, int k) {
	int i;
	tup nucMask;
	nucMask.t=3;
	tup least;
	b.t=0;
	for (i=0; i<k; i++) {
		least.t = ~ (a.t & nucMask.t) & nucMask.t;
		a.t >>=2;
		b.t <<=2;
		b.t += least.t;
	}
}

typedef pair<GenomeTuple, GenomeTuple> GenomePair;
typedef vector<GenomePair> GenomePairs;
typedef pair<LocalTuple, LocalTuple> LocalPair;
typedef vector<LocalPair> LocalPairs;

template<typename List>
void AppendValues(GenomePairs &dest, typename List::iterator sourceStart, typename List::iterator sourceEnd, GenomePos queryOffset, 
					GenomePos targetOffset) {
	int i=dest.size();
	dest.resize(dest.size() + sourceEnd-sourceStart);
	typename List::iterator sourceIt=sourceStart;
	for (; i < dest.size(); i++, ++sourceIt) {
		dest[i].first.pos = sourceIt->first.pos+ queryOffset;
		dest[i].second.pos = sourceIt->second.pos + targetOffset;
		dest[i].first.t = sourceIt->first.t;
		dest[i].second.t = sourceIt->second.t;
	}
}

template<typename List>
void AppendValues(GenomePairs &dest, typename List::iterator sourceStart, typename List::iterator sourceEnd, GenomePos &queryOffset,
					 GenomePos &targetOffset, int64_t &maxDiagNum, int64_t &minDiagNum, GenomePos qs, GenomePos qe, 
					 GenomePos ts, GenomePos te, GenomePos &prev_readStart, GenomePos &prev_readEnd) {
	//int i=dest.size();
	//dest.resize(dest.size() + sourceEnd-sourceStart);
	typename List::iterator sourceIt=sourceStart;
	int64_t diag;
	for (; sourceIt < sourceEnd; ++sourceIt) {
	  diag = (int64_t)(sourceIt->second.pos + targetOffset) - (int64_t)(sourceIt->first.pos + queryOffset);
	  if ( diag >= minDiagNum and diag <= maxDiagNum
	       and sourceIt->first.pos + queryOffset >= qs 
	       and sourceIt->first.pos + queryOffset < qe 
	       and sourceIt->second.pos + targetOffset >= ts 
	       and sourceIt->second.pos + targetOffset < te) {

	    dest.push_back(GenomePair(GenomeTuple(sourceIt->first.t, sourceIt->first.pos + queryOffset), 
				      GenomeTuple(sourceIt->second.t, sourceIt->second.pos + targetOffset)));
	    prev_readEnd = max(prev_readEnd, sourceIt->first.pos + queryOffset);
	    prev_readStart = min(prev_readStart, sourceIt->first.pos + queryOffset);
	    
	  }
	  else {
	    /*
	    cout << "Not adding " << "\t" << diag << "\t" << minDiagNum << "\t" << maxDiagNum << endl;
	    cout << "target: " << sourceIt->second.pos + targetOffset << "\t" << sourceIt->second.pos << "\t" << targetOffset << "\t" << ts << "\t" << te << endl;
	    cout << "query: " << sourceIt->first.pos + queryOffset << "\t" << sourceIt->first.pos << "\t" << queryOffset << "\t" << qs << "\t" << qe << endl;
	    cout <<  "diag " << (int) (diag >= minDiagNum and diag <= maxDiagNum) << "\t"
		 << "qs "<< (int) (sourceIt->first.pos + queryOffset >= qs ) << "\t"
		 << "qe " << (int) ( sourceIt->first.pos + queryOffset < qe ) << "\t"
		 << "ts " <<  (int)( sourceIt->second.pos + targetOffset >= ts ) << "\t"
		 << "te " << (int) (sourceIt->second.pos + targetOffset < te ) << endl;
	      }
	    */
	  }
	}
}
#endif
