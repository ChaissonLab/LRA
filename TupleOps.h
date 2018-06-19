#ifndef TUPLE_OPS_H_
#define TUPLE_OPS_H_

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


class LocalTuple {
 public:
	SmallTuple t;
	SmallPos   pos;

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
	int pos;
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
	m=0;
	for (int i = 0; i < k; i++) {
		m<<=2;
		m+=3;
	}
}

template<typename tup> void StoreTuple(char *seq, int pos, int k, tup &t) {
	t = 0;
	for (int p=pos; p <= pos+k-1; p++) {
		t <<=2;
		t+=(tup)seqMap[seq[p]];
	}
}

template<typename tup> void ShiftOne(char *seq, int pos, tup mask, tup &t) {
	t = (t << 2) & (tup) mask;
	t += (tup)seqMap[seq[pos]];	
}

template<typename tup> void ShiftOneRC(char *seq, int pos, int k, tup &t) {
	t >>= 2;
	t += (~(seqMap[seq[pos]]) & (tup)3) << (2*((tup)k-1));
}


template <typename tup> void TupleRC(tup a, tup &b, int k) {
	int i;
	tup nucMask=3;
	tup least;
	b=0;
	for (i=0; i <k; i++) {
		least = ~ (a & nucMask) & nucMask;
		a >>=2;
		b <<=2;
		b += least;
	}
}


#endif
