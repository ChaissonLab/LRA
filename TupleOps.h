#ifndef TUPLE_OPS_H_
#define TUPLE_OPS_H_

typedef uint32_t GenomePos;
typedef uint64_t Tuple;
typedef uint16_t SmallTuple;
typedef uint16_t SmallPos;

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
  uint32_t t:32-LOCAL_POS_BITS;
  uint32_t pos: LOCAL_POS_BITS;

	bool operator<(const LocalTuple &b) {
		return t < b.t;
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


class GenomeTuple {
public:
	Tuple t;
	GenomePos pos;
	bool operator<(const GenomeTuple &b) {
		return t < b.t;
	}

	friend int operator >(GenomeTuple &a, GenomeTuple &b) {
		return a.t > b.t;
	}
	friend int operator!=(GenomeTuple &a, GenomeTuple &b) {
		return a.t != b.t;
	}
	void ToString(int k, string &s) {
		TupleToString(t, k, s);
	}
};

template< typename Tup>
int InitMask(Tup &m, int k) {
	m.t=0;
	for (int i = 0; i < k; i++) {
		m.t<<=2;
		m.t+=3;
	}
}

template<typename tup> void StoreTuple(char *seq, int pos, int k, tup &t) {
	t.t=0;
	for (int p=pos; p <= pos+k-1; p++) {
		t.t <<=2;
		t.t+=seqMap[seq[p]];
	}

}

template<typename tup> void ShiftOne(char *seq, int pos, tup mask, tup &t) {
	t.t= (t.t << 2) & (Tuple) mask.t;
	t.t += (Tuple)seqMap[seq[pos]];	
}

template<typename tup> void ShiftOneRC(char *seq, int pos, int k, tup &t) {
	t.t >>= 2;
	t.t += (~(seqMap[seq[pos]]) & (Tuple)3) << (2*((Tuple)k-1));
}


template <typename tup> void TupleRC(tup a, tup &b, int k) {
	int i;
	tup nucMask;
	nucMask.t=3;
	tup least;
	b.t=0;
	for (i=0; i <k; i++) {
		least.t = ~ (a.t & nucMask.t) & nucMask.t;
		a.t >>=2;
		b.t <<=2;
		b.t += least.t;
	}

}


#endif
