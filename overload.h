// C++ program to print vector objects 
// by overloading "<<" operator 

#ifndef OVERLOADING_H_
#define OVERLOADING_H_


#include <ostream>
#include <vector> 
#include <stack>

using std::pair;

// C++ template to print vector container elements 
template <typename T> 
std::ostream & operator<<(std::ostream & os, const std::vector<T> & v) 
{ 
    os << "["; 
    for (int i = 0; i < v.size(); ++i) { 
    	if (i != v.size() - 1) {os << v[i] << " "; }
        else {os << v[i]; }
    } 
    os << "]\n"; 
    return os; 
} 

// print out Pair
template <typename T>
std::ostream & operator<<(std::ostream & os, const std::pair<T,T> & M) {

	os << "(" << M.first << ", " << M.second << ")";
	return os;
}

// print out stack
template <typename T>
std::ostream & operator<<(std::ostream & os, const std::stack<T> & N) {

	std::stack<T> M;
	M = N;
	unsigned int m = M.size();
	os << "(";
	for (unsigned int i = 0; i < m; ++i) {
		if (i != M.size() - 1) {
			os << M.top() << ", ";
		}
		else {
			os << M.top();
		}
		M.pop();
	}
	os << ")\n";
	return os;
}



#endif

  