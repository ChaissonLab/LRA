#ifndef LOG_LOOK_UP_TABLE_H_
#define LOG_LOOK_UP_TABLE_H_

#include<vector>
#include<cmath>

// static vector<float> LookUpTable;

void 
CreateLookUpTable(std::vector<float> & LookUpTable){
	for (int i = 1; i <= 10001; i = i + 5) {
		LookUpTable.push_back(logf(i));
	}

}


#endif