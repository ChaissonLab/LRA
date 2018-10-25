#ifndef READ_H_
#define READ_H_
using namespace std;
class Read {
 public:
	char *seq;
	char *qual;
	int  length;
	char *passthrough;
	string name;
	void Clear() {
		if (seq != NULL) {
			delete[] seq;
			seq=NULL;
		}
		if (qual != NULL) {
			delete[] qual;
			qual=NULL;
		}
		if (passthrough != NULL) {
			delete[] passthrough;
			passthrough=NULL;
		}
		length=0;
		name="";
	}
	Read() {
		seq=NULL;
		length=0;
		qual=NULL;
		name="";
		passthrough=NULL;
	}
	Read(char* _sq, int _len, string _name, char*_qual=NULL) {
		seq=_sq;
		length=_len;
		qual=_qual;
		passthrough=NULL;
		name=_name;
	}
};

#endif
