#ifndef PFRAGMENT_INFO_H_
#define PFRAGMENT_INFO_H_

// TODO(Jingwen): combine Value Vector with Fragment_Info Pair
#include <vector>
#include <ostream>


class Fragment_Info
{
public:
	float val;
	int clusterNum; // store the index of the Cluster which the current anchor comes from;
	int matchstartNum;
	long int prev_sub; // the previous subproblem's number
	long int prev_ind; // the index in the Ev of the previous subproblem
	bool prev; //if prev == TRUE then the previous subproblem is row subproblem. Else it's col subproblem 
	bool inv; // if inv == TRUE, then the previous subproblem is dividing (s1, e1). Else it's dividing (s2, e2)
	bool orient; // if orient = 0 means reverse orientated anchor
	std::vector<unsigned int> SS_A_R1;
	std::vector<unsigned int> SS_B_R1;
	unsigned int counter_A_R1;
	unsigned int counter_B_R1;
	std::vector<unsigned int> SS_A_C1;
	std::vector<unsigned int> SS_B_C1;
	unsigned int counter_A_C1;
	unsigned int counter_B_C1;
	std::vector<unsigned int> SS_A_R2;
	std::vector<unsigned int> SS_B_R2;
	unsigned int counter_A_R2;
	unsigned int counter_B_R2;
	std::vector<unsigned int> SS_A_C2;
	std::vector<unsigned int> SS_B_C2;
	unsigned int counter_A_C2;
	unsigned int counter_B_C2;
	Fragment_Info();
	~Fragment_Info() {}; // deconstructor
	friend std::ostream & operator<<(std::ostream & os, const Fragment_Info & M);
};


Fragment_Info::Fragment_Info () {
	prev_sub = -1;
	prev_ind = -1;
	prev = 1;	
	inv = 1;
	orient = 1;
}


std::ostream & operator<<(std::ostream & os, const Fragment_Info & M) {
	os << "val: " << M.val << ", prev_sub: " << M.prev_sub << ", prev_ind: " << M.prev_ind << ", prev: " << M.prev << "\n";
	os << "SS_A_R1: " << M.SS_A_R1 << "\n";
	os << "SS_B_R1: " << M.SS_B_R1 << "\n";
	os << "counter_A_R1: " << M.counter_A_R1 << "\n";
	os << "counter_B_R1: " << M.counter_B_R1 << "\n";
	os << "SS_A_C1: " << M.SS_A_C1 << "\n";
	os << "SS_B_C1: " << M.SS_B_C1 << "\n";
	os << "counter_A_C1: " << M.counter_A_C1 << "\n";
	os << "counter_B_C1: " << M.counter_B_C1 << "\n";
	return os;
}


class Fragment_valueOrder {
 public:
	vector<float> fragments_value;
	vector<int> index;
	
 	Fragment_valueOrder(const vector<Fragment_Info> *c) {
 		vector<float> fv(c->size(), 0);
 		for (unsigned int i = 0; i < c->size(); i++) {
 			fv[i] = (*c)[i].val;
 		}
 		fragments_value = fv;
 		assert(fragments_value.size() == c->size());
		index.resize(c->size());
		for (unsigned int i = 0;i < index.size(); i++) { index[i] = i;}
		Sort();
	}

	int operator()(const int i, const int j) {
		assert(i < fragments_value.size());
		assert(j < fragments_value.size());
		return fragments_value[i] > fragments_value[j];			
	}

	void Sort() {
		sort(index.begin(), index.end(), *this);
	}

	float & operator[](int i) {
		return fragments_value[index[i]];
	}
	
	int size() {
		return index.size();
	}
};


#endif