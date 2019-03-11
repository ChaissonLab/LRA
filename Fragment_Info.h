#ifndef PFRAGMENT_INFO_H_
#define PFRAGMENT_INFO_H_

// TODO(Jingwen): combine Value Vector with Fragment_Info Pair
#include <vector>
#include <ostream>


class Fragment_Info
{
public:
	float val;
	long int prev_sub;
	long int prev_ind;
	bool prev; //if prev == TRUE then the previous subproblem is row subproblem. Else it's col subproblem 
	std::vector<unsigned int> SS_A_R;
	std::vector<unsigned int> SS_B_R;
	unsigned int counter_A_R;
	unsigned int counter_B_R;
	std::vector<unsigned int> SS_A_C;
	std::vector<unsigned int> SS_B_C;
	unsigned int counter_A_C;
	unsigned int counter_B_C;
	Fragment_Info();
	~Fragment_Info() {}; // deconstructor
	friend std::ostream & operator<<(std::ostream & os, const Fragment_Info & M);
};


Fragment_Info::Fragment_Info () {
	prev_sub = -1;
	prev_ind = -1;
	prev = 1;	
}


std::ostream & operator<<(std::ostream & os, const Fragment_Info & M) {
	os << "val: " << M.val << ", prev_sub: " << M.prev_sub << ", prev_ind: " << M.prev_ind << ", prev: " << M.prev << "\n";
	os << "SS_A_R: " << M.SS_A_R << "\n";
	os << "SS_B_R: " << M.SS_B_R << "\n";
	os << "counter_A_R: " << M.counter_A_R << "\n";
	os << "counter_B_R: " << M.counter_B_R << "\n";
	os << "SS_A_C: " << M.SS_A_C << "\n";
	os << "SS_B_C: " << M.SS_B_C << "\n";
	os << "counter_A_C: " << M.counter_A_C << "\n";
	os << "counter_B_C: " << M.counter_B_C << "\n";
	return os;
}


#endif