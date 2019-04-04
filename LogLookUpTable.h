#ifndef LOG_LOOK_UP_TABLE_H_
#define LOG_LOOK_UP_TABLE_H_

#include<vector>
#include<cmath>

void 
CreateLookUpTable(std::vector<float> & LookUpTable){
	for (int i = 501; i <= 10001; i = i + 5) {
		//cerr << "i: " << i << endl;
		LookUpTable.push_back(std::log(i));
	}

}


#endif