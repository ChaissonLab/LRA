#ifndef POINT_H_
#define POINT_H_

#include "Types.h"

class Point
{
public:

	Pair se; //store the coordinates of a point (q, t)
	bool orient; //  if orient = 0 means reverse orientated anchor
	bool ind; // ind = 1 means this is a start; ind = 0 means this is an end
	bool inv; // inv = 1 means this is forward directiion; inv = 0 means this is a backward direction
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