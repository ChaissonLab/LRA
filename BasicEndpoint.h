#ifndef BASIC_ENDPOINT_H_
#define BASIC_ENDPOINT_H_




/*
 An endpoint is one of the ends of a fragment, where
 a fragment is an exact match between two genomes. 
 So, a fragment is a 2D object that has a position and length,
 and an endpoint is 1D, where it just has a position.
 A fragment may be associated with a score that is the score
 of the fragment in a maximum scoring chain.  When finding a 
 maximum scoring chain using priority search trees, one must 
 be able to set the score of a fragment when indexing solely 
 a point. 
*/


class Coordinate {
 private:
	UInt x;
  UInt y; 
 public:
 UInt GetX() const { return x;}
 UInt GetY() const { return y;}
 UInt SetX(UInt _x) { return (x = _x);}
 UInt SetY(UInt _y) { return (y = _y);}
 int operator<(const Coordinate &rhs) const {
	 if (x == rhs.GetX()) return y < rhs.GetY();
	 else return x < rhs.GetX();
 }

 int operator<=(const Coordinate &rhs) const {
	 return (*this < rhs) or (x == rhs.x && y == rhs.y);
 }

 int Equals(const Coordinate &rhs) const {
	 return (x == rhs.GetX() and y == rhs.GetY());
 }

 //
 // Synonym for Equals.
 //
	int operator==(const Coordinate &rhs) const {
		return this->Equals(rhs);
	}

 Coordinate &operator=(const Coordinate &rhs) {
	 this->x = rhs.x;
	 this->y = rhs.y;
	 return *this;
 }
};
	

template<typename T_ScoredFragment>
 class BasicEndpoint {
	 T_ScoredFragment *fragmentPtr;

	public:
	 enum WhichEnd {Start, End};
	 //	typedef Coordinate KeyType;
	 typedef UInt KeyType;
	class LessThan {
	public:
		int operator()(const BasicEndpoint<T_ScoredFragment> &lhs, const BasicEndpoint<T_ScoredFragment> &rhs) const {
			return lhs.p < rhs.p;
		}
	};
	
 public:// private:
	Coordinate p;
	WhichEnd side;

	WhichEnd GetSide() { return side; }

	void FragmentPtrToStart(T_ScoredFragment *fragment) {
		p.SetX(fragment->GetX());
		p.SetY(fragment->GetY());
		side = Start;
		fragmentPtr = fragment;
	}

	void FragmentPtrToEnd(T_ScoredFragment *fragment) {
		p.SetX(fragment->GetX() + fragment->GetXLength());
		p.SetY(fragment->GetY() + fragment->GetYLength());
		side = End;
		fragmentPtr = fragment;
	}

	int GetScore() {
		return fragmentPtr->GetScore();
	}
	int SetScore(int score) {
		return (fragmentPtr->SetScore(score));
	}

	T_ScoredFragment* SetScoredReference(T_ScoredFragment *_fragmentPtr) {
		return (fragmentPtr = _fragmentPtr);
	}
	int operator<(const BasicEndpoint &rhs) const {
		return p < rhs.p;
	}		
	KeyType GetKey() {
		return p.GetY();
	}
	T_ScoredFragment* GetFragmentPtr() {
		return fragmentPtr;
	}
	void SetChainPrev(T_ScoredFragment *prevChainFragment) {
		fragmentPtr->SetChainPrev(prevChainFragment);
	}
};
 


#endif
