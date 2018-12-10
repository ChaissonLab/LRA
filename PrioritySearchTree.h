#ifndef PRIORITY_SEARCH_TREE_H_
#define PRIORITY_SEARCH_TREE_H_
#include <iostream>

/*
 * Define a priority search tree on a point that implements 
 * the following interface:
 *
 * int T_point::GetIndex()
 *    - Return the index of the point in a list of points.
 * int T_point::GetKey()
 *    - Return the key value that the points are sorted by (x-value in a 2D query)
 * int T_point::GetValue()
 *    - Return the value of a point.
 * int T_point::SetValue(int value)
 *    - sets the value of a point.
 *
 * This class implements a query FindMax(key), which returns
 * the index of the point with greatest value of all points with key [0...key).
 *
 * 
 */
#include <assert.h>
template<typename T_Point>
class PSTVertex {
 public: 
	unsigned int leftChildIndex;
	unsigned int rightChildIndex;
	unsigned int isALeaf;
	typename T_Point::KeyType medianKey;
	typename T_Point::KeyType maxKey;
	int pointIndex;
	int maxScoreNode;
	int fx,fy;
	PSTVertex() {
		isALeaf         = 0;
		leftChildIndex  = 0;
		rightChildIndex = 0;
		maxScoreNode = -1;
    maxKey = 0;
    medianKey = 0;
    pointIndex = -1;
		fx=-1,fy=-1;
	}
};

template<typename T_Point>
class PrioritySearchTree {
 public:
	void Print() {
		int i;
		cerr << "index\tleaf\tleft\tright\tmedian\tmax\tpoint\tnode\tx\ty" << endl;
		for (i=0; i < tree.size(); i++) {
			cerr << "pst: " << i << "\t" << (int) tree[i].isALeaf << "\t" << tree[i].leftChildIndex << "\t" << tree[i].rightChildIndex << "\t" << tree[i].medianKey << "\t" << tree[i].maxKey << "\t" << tree[i].pointIndex << "\t" << tree[i].maxScoreNode << "\t" << tree[i].fx << "\t" << tree[i].fy << endl;
		}
	}

 private:

	vector<PSTVertex<T_Point> > tree;
	vector<PSTVertex<T_Point> > *treePtr;
	int GetMedianIndex(int start, int end) {
		return (end + start) / 2;
	}
	typename T_Point::KeyType CreateTree(vector<T_Point> &points, int start, int end, unsigned int &iterativeIndex) {

		assert(iterativeIndex < (*treePtr).size());
		//
		// Look to see if this vertex is the parent of a leaf -- when there are only
		// two points below.
		//

		int medianIndex    = GetMedianIndex(start, end);
		int curVertexIndex = iterativeIndex;
		(*treePtr)[curVertexIndex].medianKey = points[medianIndex].GetKey();

		if (end == start) {
			// No children for this node, done.
			(*treePtr)[curVertexIndex].pointIndex = start;
			(*treePtr)[curVertexIndex].fx=points[medianIndex].x;
			(*treePtr)[curVertexIndex].fy=points[medianIndex].y;
			return (*treePtr)[curVertexIndex].medianKey;
		}
		// 
		// Check to see if the current
		// node is a leaf node.  No recursion on this node.
		//
		
		if (end - start == 1) {
			(*treePtr)[curVertexIndex].isALeaf    = 1;
			(*treePtr)[curVertexIndex].medianKey  = points[start].GetKey();
			(*treePtr)[curVertexIndex].pointIndex = start;
			(*treePtr)[curVertexIndex].fx=points[medianIndex].x;
			(*treePtr)[curVertexIndex].fy=points[medianIndex].y;

			//
			// Return the key of this vertex.  The parent
			// will know what to do with it.  If this is
			// a left child, the parent will use the key to
			// distinguish what is on the left side of the branches.
			// If it is the right side of a (*treePtr), it is ignored.
			//
			return (*treePtr)[curVertexIndex].medianKey;
		}
		else {
			//
			// This vertex contains at least two children, so it is not 
			// a leaf.  Recurse assigning leaves.
			//

			(*treePtr)[curVertexIndex].isALeaf = 0;
			(*treePtr)[curVertexIndex].leftChildIndex = ++iterativeIndex;
			typename T_Point::KeyType leftTreeKey, rightTreeKey;
			leftTreeKey  = CreateTree(points, start, medianIndex, iterativeIndex);
			
			//
			// The leftTreeKey separates the branches BELOW this vertex.
    	//
			(*treePtr)[curVertexIndex].medianKey = leftTreeKey;
			
			(*treePtr)[curVertexIndex].rightChildIndex = ++iterativeIndex;
			rightTreeKey = CreateTree(points, medianIndex, end, iterativeIndex);
			//
			// The rightTreeKey will separate the parent's left tree from the right.
			//
			(*treePtr)[curVertexIndex].maxKey = rightTreeKey;
			return rightTreeKey;
		}
	}
	
	int FindIndexOfMaxPoint(int curVertexIndex, 
													vector<T_Point> &points, typename T_Point::KeyType maxKey, int &maxPointValue, int &maxPointIndex) {

		//
		// Attempt to find the leaf vertex beneath this vertex that has the largest 
		// score, with a key less than max key.
		//
		// On return: 
		//   Return 1 if a value is assigned to maxPointValue, 0 otherwise.
		//   If a value is assigned to maxPointValue, this sets:
		//      maxPointValue is the score of the maximum point.
		//      maxPointIndex the index of the point in 'points' that has the maximum score.
		//   

		//
		// The vertex at curVertexIndex has a max score node beneath it, if it has been
		// initialized.  If the maxScoreNode has a key less than the current maxKey, then we 
		// know the maximum value is contained beneath this vertex, AND that its key
		// is within the range in the rage maximum query.
		// That means that there is no need to continue the search below here.
		//
		if ((*treePtr)[curVertexIndex].maxScoreNode == -1) {
			return 0;
		}
		if (points[(*treePtr)[curVertexIndex].maxScoreNode].GetKey() < maxKey) {
			if (points[(*treePtr)[curVertexIndex].maxScoreNode].score > maxPointValue) {
				maxPointValue = points[(*treePtr)[curVertexIndex].maxScoreNode].score;
				maxPointIndex = (*treePtr)[curVertexIndex].maxScoreNode;
				return 1;
		 	}
			else {
				return 0;
			}
		}
		//
		// Otherwise, the maximum scoring node beneath this node has a key greater than
		// the max key. That means that the search must continue for the maximum value
		// node with a key less than 'maxKey'.  The search has two cases:
		//
		// First, if the median key of this node is greater than the maxKey, all keys
		// on the right side of the tree are greater than maxKey, so do not search there.
		//
		// If the median key of this node si less than maxKey, there may be a node on the left
		// or right child of the current node with a maximum key.  Search both to the 
		// left and right.
		//
		else {
			if (!(*treePtr)[curVertexIndex].isALeaf) {
				if (maxKey <= (*treePtr)[curVertexIndex].medianKey) {
					return FindIndexOfMaxPoint((*treePtr)[curVertexIndex].leftChildIndex, points, maxKey, maxPointValue, maxPointIndex);
				}
				else {
					int foundValueLeft, foundValueRight;
					foundValueLeft = FindIndexOfMaxPoint((*treePtr)[curVertexIndex].leftChildIndex, points, maxKey, maxPointValue, maxPointIndex);
					foundValueRight = FindIndexOfMaxPoint((*treePtr)[curVertexIndex].rightChildIndex, points, maxKey, maxPointValue, maxPointIndex);
					return (foundValueLeft or foundValueRight);
				}
			}
			else {
				// 
				// The current node is a leaf node, but due to the condition from before, its key
				// is greater than or equal to the max key, therefore its score cannot 
				// be used for the maximum score.
				// Returning 0 here signifies that this search-branch did not turn up any candidates for
				// the maximum scoring node.
				return 0;
			}
		}
	}
	
 public:
    PrioritySearchTree() {
        treePtr = NULL;
    }

	void CreateTree(vector<T_Point> &points, vector<PSTVertex<T_Point> > *bufTreePtr=NULL) {
		/*
		 * Precondition: points is sorted according to key.
		 */

		// 
		// The tree is a binary tree containing all the points.  The 
		// perfectly balanced tree is of maximum size points.size()-1, so 
		// go ahead and preallocate that now.
		//
		if (bufTreePtr != NULL) {
			treePtr = bufTreePtr;
		}
		else {
			treePtr = &tree;
		}
 		treePtr->resize((points.size() * 2) - 1);
		unsigned int curVertexIndex = 0;
		CreateTree(points, 0, points.size(), curVertexIndex);
	}


	//
	// Implement the tree as an array of interior nodes.
	// Since there is already space allocated for the 


	int FindPoint(typename T_Point::KeyType pointKey, int curVertexIndex, int &pointVertexIndex) {
		if ((*treePtr)[curVertexIndex].isALeaf) {
			pointVertexIndex = curVertexIndex;
			return (*treePtr)[curVertexIndex].medianKey == pointKey;
		}
		else {
			if (pointKey <= (*treePtr)[curVertexIndex].medianKey) {
				return FindPoint(pointKey, (*treePtr)[curVertexIndex].leftChildIndex, pointVertexIndex);
			}
			else {
				return FindPoint(pointKey, (*treePtr)[curVertexIndex].rightChildIndex, pointVertexIndex);
			}
		}
	}

	void Activate(vector<T_Point> &points, int pointIndex) {
		int pointScore = points[pointIndex].score;
		// Now, update the pMax scores in the tree

		int curVertexIndex = 0;
		typename T_Point::KeyType pointKey = points[pointIndex].GetKey();
		unsigned int itIndex = 0;
		while (pointIndex != -1 and
					 (*treePtr)[curVertexIndex].isALeaf == 0) {
			assert(itIndex < (*treePtr).size());
			if ((*treePtr)[curVertexIndex].maxScoreNode == -1 or 
					points[(*treePtr)[curVertexIndex].maxScoreNode].score <= pointScore) {
				int tmpPMaxIndex = (*treePtr)[curVertexIndex].maxScoreNode;
				(*treePtr)[curVertexIndex].maxScoreNode = pointIndex;
				pointIndex = tmpPMaxIndex;
			}
			
			if (pointKey <= (*treePtr)[curVertexIndex].medianKey) {
				curVertexIndex = (*treePtr)[curVertexIndex].leftChildIndex;
			}
			else {
				curVertexIndex = (*treePtr)[curVertexIndex].rightChildIndex;
			}

			// Keep track of the number of times this loop is executed... an 
			// infinite loop will bomb.
			++itIndex;
		}
	}

	int FindIndexOfMaxPoint(vector<T_Point> &points, typename T_Point::KeyType maxPointKey, int &maxPointIndex) {
		// start at the root
		int curVertexIndex = 0;
		if ((*treePtr)[curVertexIndex].maxScoreNode == -1) {
			//
			// This case can only be hit if none of the points have been
			// activated. 
			//
			return 0;
		}
		int maxPointValue = -1;
		return FindIndexOfMaxPoint(0, points, maxPointKey, maxPointValue, maxPointIndex);
	}
};

#endif
