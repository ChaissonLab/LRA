#ifndef TUPLE_OPS_H_
#define TUPLE_OPS_H_

typedef uint64_t Tuple;
Tuple mask=0;
class GenomeTuple {
public:
	Tuple tuple;
	int pos;
	bool operator<(const GenomeTuple &b) {
		return tuple < b.tuple;
	}

	friend int operator >(GenomeTuple &a, GenomeTuple &b) {
		return a.tuple > b.tuple;
	}
	friend int operator!=(GenomeTuple &a, GenomeTuple &b) {
		return a.tuple != b.tuple;
	}

};

int InitMask(int k) {
	mask=0;
	for (int i = 0; i < k; i++) {
		mask<<=2;
		mask+=3;
	}
}

template<typename tup> void StoreTuple(char *seq, int pos, int k, tup &tuple) {
	tuple = 0;
	for (int p=pos; p <= pos+k-1; p++) {
		tuple <<=2;
		tuple+=seqMap[seq[p]];
	}
}

template<typename tup> void ShiftOne(char *seq, int pos, tup &tuple) {
	tuple = (tuple << 2) & mask;
	tuple += seqMap[seq[pos]];	
}

template<typename tup> void ShiftOneRC(char *seq, int pos, int k, tup &tuple) {
	tuple >>= 2;
	tuple += (~(seqMap[seq[pos]]) & (uint64_t)3) << (2*((tup)k-1));
}


template <typename tup> void TupleRC(tup a, tup &b, int k) {
	int i;
	unsigned int nucMask=3;
	unsigned int least;
	b=0;
	for (i=0; i <k; i++) {
		least = ~ (a & nucMask) & nucMask;
		a >>=2;
		b <<=2;
		b += least;
	}
}


#endif
