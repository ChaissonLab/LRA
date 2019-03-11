#ifndef POINT_H_
#define POINT_H_

#include "Types.h"

class Point
{
public:

	Pair se; //store the start and end of a point
	bool ind; // ind = 1 means this is a start
	unsigned int frag_num; // store the index of the fragment that contains this point
	Point(unsigned int & frag_num1);
	Point() {};
	~Point() {};
	
	friend std::ostream & operator<<(std::ostream & os, const Point & t); // overload of operator <<
};


Point::Point(unsigned int & frag_num1) {
	frag_num = frag_num1;
}

std::ostream & operator<<(std::ostream & os, const Point & M) {
	os << "Point: { Pair:" << M.se << ";  ind: " << M.ind << ";  frag_num: " << M.frag_num << "\n";  
	return os;
}


#endif