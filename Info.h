#ifndef INFO_H_
#define INFO_H_

#include <vector>
#include <ostream>

using std::vector;

typedef std::pair<unsigned int, unsigned int> Pair;


class info
 {
 public:
 	unsigned int pstart;
 	unsigned int pend;
 	unsigned int rc_num; // rc_num means the row/col number 
 	unsigned int num; // num means the subproblem Sub[num] which the current row/col belongs to 
 	vector<unsigned int> SS_A; // SS_A stores the subproblem number which end points on the current row are in Di
 	vector<unsigned int> SS_B; //  SS_B stores the subproblem number which start points the current row are in Ei
 	info(unsigned int s, unsigned int e, unsigned int n) : pstart(s), pend(e), rc_num(n), num(0) {}  // constructor
 	~info() {};
 	friend std::ostream & operator<<(std::ostream & os, const info & t); // overload of operator <<
 }; 

std::ostream & operator<<(std::ostream & os, const info & M) {
	os << "{pstart: " << M.pstart << ", pend: " << M.pend << ", rc_num: " << M.rc_num << ", num: " << M.num << endl;
	os << "SS_A: " << M.SS_A << endl;
	os << "SS_B: " << M.SS_B << "} "<< endl;

	return os;
}

#endif