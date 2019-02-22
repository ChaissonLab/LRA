using namespace std;
#include "GlobalChain.h"
#include "Fragment.h"
#include <vector>
#include <iostream>

typedef vector<vector<int> > v;

void TestGlobalChain(v ep) {
	int i;
	vector<Endpoint> endpoints;
	vector<Fragment> fragments;
	vector<int> opt;
	for (i=0;i<ep.size();i++) {
		fragments.push_back(Fragment(ep[i][0],ep[i][1],ep[i][2],ep[i][3],ep[i][2]-ep[i][0],0));
	}
	GlobalChain( fragments, opt, endpoints);
	cout << "Opt of size " << opt.size() << endl;
	for (i=0;i<opt.size(); i++) {
		cout << fragments[opt[i]].xl << "\t" 
				 << fragments[opt[i]].yl << "\t" 		 
				 << fragments[opt[i]].xh << "\t" 		 
				 << fragments[opt[i]].yh << "\t"
				 << fragments[opt[i]].score << endl;
	}
}

int main() {
	TestGlobalChain({ 
			{0,0,10,10}, 
		  {20,20,30,30},
			{40,40,50,50}, 
			{60,60,70,70}, 
			{80,80,90,90}, 
			{100,100,110,110}, 
			{120,120,130,130}, 
			{140,140,150,150}, 
				{81,31, 91,41}});
}
