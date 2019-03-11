#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_


#include <vector>	
#include <list>
#include <utility>
#include <limits>
#include <stack>
#include "overload.h"


typedef std::pair<unsigned int, unsigned int> Pair;
typedef std::pair<long int, long int> LPair;


class Subproblem
{
public:
	unsigned int num; // the number of subproblems 
	//unsigned int counter; // if (counter == 0) then the D array of this subproblem is fully filled out.
	unsigned int now;
	long int last;
	std::vector<long int> Di, Ei;
	std::vector<float> Dv, Ev;
	std::vector<unsigned int> Dp; // store the index of the point which gives the maximum value 
	std::vector<unsigned int> Ep; // store the Point which gives the maximum value 
	std::vector<long int> Eb; // Ei[j] gives the index of point in D array such that Ei[j] is the first element in Ei array that is >= Di[Ei[j]]
	std::vector<long int> Db;
	std::vector<unsigned int> D; // storing the index [0,1,2,...., n]
	std::vector<unsigned int> E;
	std::vector<std::pair<long int, long int>> Block; // store the block for the maximization structure
	std::stack<LPair> S_1; // S_1 will be used in the middle of computing maximization structure
	Subproblem(); // default constructor
	Subproblem(unsigned int num1);
	//Subproblem(unsigned int & num1, std::vector<long int> & Di1, std::vector<long int> & Ei1); // for leaf case 
	//Subproblem(unsigned int & num1, std::vector<long int> & Di1, std::vector<long int> & Ei1, std::vector<long int> & Eb1, std::vector<long int> & Db1); // for non-leaf case
	~Subproblem() {};
	friend std::ostream & operator<<(std::ostream & os, const Subproblem & M);
};

Subproblem::Subproblem() {}

//initialization 
Subproblem::Subproblem(unsigned int num1) {

	num = num1;
	last = -1;
	now = 0;
}




///////////////////////////////////////////////////////////////////////////////

/*
//initialization for leaf case
Subproblem::Subproblem(unsigned int & num1, std::vector<long int> & Di1, std::vector<long int> & Ei1) {
	unsigned int h = Ei1.size();

	num = num1;
	Ei = Ei1;
	E = Ei;
	std::iota(E.begin(), E.end(), 0);
	std::vector<float> p(h, 10);
	std::vector<long> z(h, 0);
	Ev = p;
	Ep = z;
	last = -1;
	now = 0;
}

// initializationfor non-leaf case
Subproblem::Subproblem(unsigned int & num1, std::vector<long int> & Di1, std::vector<long int> & Ei1, std::vector<long int> & Eb1, std::vector<long int> & Db1) {
	unsigned int h = Ei1.size();
	unsigned int l = Di1.size();

	now = 0;
	last = -1;
	Di = Di1;
	std::vector<float> v(l, 0);
	std::vector<unsigned int> w(l, 0);
	Db = Db1;
	Dv = v; 
	Dp = w;
	num = num1;
	Ei = Ei1;
	D = Di;
	E = Ei;
	std::iota(D.begin(), D.end(), 0);
	std::iota(E.begin(), E.end(), 0);
	std::vector<float> p(h, 10);
	std::vector<long> z(h, -1);
	Ev = p;
	Ep = z;
	Eb = Eb1;
	std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
	S_1.push(dummy_pair); 
}
*/



std::ostream & operator<<(std::ostream & os, const Subproblem & M) {
	os << "{num: " << M.num << "\n";
	os << "Di: " << M.Di << ", Dp: " << M.Dp << ", Dv: " << M.Dv <<  ", Db: " << M.Db << ", Ei: " << M.Ei << ", Ep: " << M.Ep << ", Ev: " << M.Ev << ", Eb: " << M.Eb << ", last: " << M.last << ", now: "
	 << M.now << ", Block: " << M.Block << ", S_1: " << M.S_1 << "\n"; 
	return os;
}



#endif