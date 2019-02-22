#ifndef FRAGMENT_H_
#define FRAGMENT_H_

class Fragment {
 public:
	int xl, yl, xh, yh;
	int score;
	int prev;
	int index;
	int GetScore() {
		return score;
	}
 Fragment(int _xl, int _yl, int _xh, int _yh, int _s, int _idx) : xl(_xl), yl(_yl), xh(_xh), yh(_yh), score(_s), index(_idx) {prev=-1;}
};
#endif
