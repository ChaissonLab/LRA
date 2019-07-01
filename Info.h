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
 	vector<unsigned int> SS_A1; // SS_A1 stores the subproblem number which end points (e1) on the current row are in Di
 	vector<unsigned int> SS_B1; //  SS_B1 stores the subproblem number which start points (s1) the current row are in Ei
 	vector<unsigned int> SS_A2; // SS_A2 stores the subproblem number which end points (e2) on the current row are in Di
 	vector<unsigned int> SS_B2; //  SS_B2 stores the subproblem number which start points (s2) the current row are in Ei
 	info(unsigned int s, unsigned int e, unsigned int n) : pstart(s), pend(e), rc_num(n), num(0) {}  // constructor
 	~info() {};
 	friend std::ostream & operator<<(std::ostream & os, const info & t); // overload of operator <<
 }; 

std::ostream & operator<<(std::ostream & os, const info & M) {
	os << "{pstart: " << M.pstart << ", pend: " << M.pend << ", rc_num: " << M.rc_num << ", num: " << M.num << endl;
	os << "SS_A1: " << M.SS_A1 << endl;
	os << "SS_B1: " << M.SS_B1 << "} "<< endl;

	return os;
}

#endif