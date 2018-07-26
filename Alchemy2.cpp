#include <zlib.h>
#include <stdint.h>
#include "htslib/kseq.h"
#include <cstdlib>
KSEQ_INIT(gzFile, gzread)
#include "htslib/sam.h"
#include <numeric>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>
#include "Genome.h"
#include <random>

long ParseNumber(string n) {
	long mult=1;
	if (n.find("G") != n.npos or n.find("g") != n.npos) {
		mult=1000000000;
	}
	else if (n.find("M") != n.npos or n.find("m") != n.npos) {
		mult=1000000;
	}
	else if (n.find("k") != n.npos or n.find("k") != n.npos) {
		mult=1000;
	}
	return atoi(n.c_str())*mult;
}

using namespace std;
void HelpAlchemy2() {
	cout << "alchemy2: A model based read simulator. Models are built from a bam file and " << endl
			 << "          are aligner and genome dependent. Each model is a histogram of sequencer  " << endl
			 << "          values from a k-mer context. Simulated output is sampled according " << endl
			 << "          to the histogram." << endl;
	cout << "Usage: alchemy2 [options]" << endl
			 << " Storing a model: " << endl
			 << "   -b (string) Aligned reads bam file." << endl
			 << "   -k (int)    Context k-mer." << endl
			 << "   -L (int)    Maximum number of reads to try sampling (100,000)" << endl
			 << "   -m (string) Output model." << endl
			 << "   -g (int)    Maximum gap to model. This is intended to skip adding SVs to the model (30)" << endl
			 << "   -s (int)    Minimum number of output values to sample" << endl
			 << " Simulating reads: (default without a bam file)" << endl
			 << "   -m (string) Input model." << endl		
			 << "   -g (string) Input genome"  << endl
			 << "   -o (string) Output reads." << endl
			 << "   -r (string) Simulate from these reads not the genome." << endl
			 << "   -R (string) Simulate from regions (bed file)" << endl
			 << "   -l (int)    Minimum read length (1000)." << endl
			 << "   -p (string) Simulate from positions." << endl
			 << "   -E (string) Use emprical read length distribution (false=log-normal)." << endl
			 << "   -B (int)    Number of bases to simulated. Ignores -r" << endl
			 << "   -u (int)    Override model mean read length" << endl
			 << "   -V (int)    Override model read length variance" << endl << endl;
	cout << " Examples: " << endl
			 << "   alchemy2  -b reads.bam -g hg38.fa -s 1000 -k 5 -L 50000 -m hg38.alc" << endl
			 << "   alchemy2  -g hg38.fa  -m hg38.alc -B 20G  -o genome-wide.fasta" << endl
			 << "   alchemy2  -g hg38.fa  -m hg38.alc -B 300M -R regions.bed  -o regions.fasta" << endl;
}


class Output {
public:
	int op, len;
	string seq;
	int count;
};

float Rand() {
	return ((float)rand()) / RAND_MAX;
}

long RandInt(long max) {
	float f=Rand();
	return (long)(f*max);
}

int IsDel(string s) {
	return s.find("del-") == 0;
}

int ParseDel(string del) {
	string delLenStr=del.substr(4);
	return atoi(delLenStr.c_str());
}

class ContextOutput {
public:
	vector<string> seqs;
	vector<int> nobs, cumObs;


	ContextOutput() {
		size=0;
	}
	map<string, int> counts;
	int size;
	int totObs;
	string Sample() {
		float r = Rand();
		int samp=r*totObs;
		int i=max(0,(int)(std::lower_bound(cumObs.begin(), cumObs.end(), samp)-cumObs.begin()-1));
		assert(i <= cumObs.size());
		return seqs[i];
	}

	void Increment(int op, int length, string seq="") {
		if (size == 100000) { return;}
		if (op == 1) {
			stringstream strm;
			strm << "del-" << length;
			seq=strm.str();
		}
		if (counts.find(seq) == counts.end()) {
			counts[seq] = 0;
		}
		counts[seq]++;
		size++;
	}

	void Store(ostream &out) {
		int n=counts.size();
		out.write((char*)&n, sizeof(int));
		ToVector();
		out.write((char*)&nobs[0], sizeof(int)*nobs.size());
		vector<int> seqlen;
		string allseq="";
		for (int i=0; i < seqs.size(); i++) {
			seqlen.push_back(seqs[i].size());
			allseq+=seqs[i];
		}
		out.write((char*)&seqlen[0], sizeof(int) *seqlen.size());
		int alllen=allseq.size();
		out.write((char*)&alllen, sizeof(int));
		out.write((char*)&allseq[0], sizeof(char)*allseq.size());

		
	}

	void Read(istream &fin) {
		int n;
		fin.read((char*)&n, sizeof(int));
		nobs.resize(n);
		fin.read((char*)&nobs[0], sizeof(int)*n);
		vector<int> seqlen(n);
		fin.read((char*)&seqlen[0], sizeof(int)*n);
		string allseq;
		int lenAllseq;
		fin.read((char*)&lenAllseq, sizeof(int));
		allseq.resize(lenAllseq);
		fin.read((char*)&allseq[0], sizeof(char)*lenAllseq);
		int p=0;
		for (int i=0; i < seqlen.size(); i++) {
			seqs.push_back(string(&allseq[p], seqlen[i]));
			p+=seqlen[i];
		}
		totObs=0;
		cumObs.resize(nobs.size()+1);
		cumObs[0] = 0;
		for(int i=0; i < nobs.size(); i++) { totObs+=nobs[i]; cumObs[i+1] = totObs;}
	}


	void ToVector() {
		vector<pair< int, string> > op;
		for (map<string, int>::iterator it = counts.begin();
				 it != counts.end();
				 ++it) {
			op.push_back(pair<int, string>(it->second, it->first));
		}
		sort(op.begin(), op.end());
		seqs.resize(op.size());
		nobs.resize(op.size());
		for (int i=0;i< op.size();i++) {
			nobs[i] = op[i].first;
			seqs[i] = op[i].second;
		}
	}
	void Print() {
		ToVector();
		for (int i = 0; i < seqs.size(); i++) {
			cout << nobs[i] << "\t" << seqs[i] << endl;
		}
	}
};

typedef 	map<string, ContextOutput> Context;
class Model {
public:
	long nSamples;
	map<string, ContextOutput> model;
	int k;
	int maxGap;
	string desc;
	int avgReadLength;
	int readVar;
	double avgLogReadLength;
	double readLogVar;
	vector<int> readLengths;

	Model() {
		desc="Alchemy2-model-1";		
		nSamples=0;
	}

	string Sample(string kmer) {
		if ( model.find(kmer) == model.end()) {
			cout << "ERROR: " << kmer << " is not in the model." << endl;
			exit(1);
		}
		return model[kmer].Sample();
	}

	void Init(int _k, int _mg) {
		k=_k;
		maxGap=_mg;
	}
	void WriteModel(string &outName) {
		ofstream outfile(outName.c_str(), std::ios::out|std::ios::binary);
		outfile.write((char*)&desc[0], sizeof(char)*desc.size());
		outfile.write((char*)&k,sizeof(int));
		int n=model.size();
		outfile.write((char*)&n,sizeof(int));
		for (Context::iterator it=model.begin(); it!=model.end(); ++it) {
			outfile.write((char*) &it->first[0], sizeof(char)*k);
			it->second.Store(outfile);
		}
		long nbases=0;
		long ssq=0;
		for (int i=0; i< readLengths.size(); i++) {
			nbases+=readLengths[i];
			ssq+=readLengths[i]*readLengths[i];
		}
		double lnSum=0;
		double lnSumSq=0;
		for (int i=0; i < readLengths.size(); i++) {
			double l=log(readLengths[i]);
			lnSum+=l;
			lnSumSq+=l*l;
		}
			
		avgReadLength=0;
		readVar=0;
		if (readLengths.size() > 0) {
			avgReadLength=nbases/readLengths.size();
			readVar =(int) std::sqrt(ssq/readLengths.size() - avgReadLength*avgReadLength);
			
			avgLogReadLength=lnSum/readLengths.size();
			readLogVar=(double) std::sqrt(lnSumSq/readLengths.size()- avgLogReadLength*avgLogReadLength);
		}
		outfile.write((char*)&avgReadLength, sizeof(int));
		outfile.write((char*)&readVar, sizeof(int));
		outfile.write((char*)&avgLogReadLength, sizeof(double));
		outfile.write((char*)&readLogVar, sizeof(double));
		int nrl=readLengths.size();
		outfile.write((char*)&nrl, sizeof(int));
		outfile.write((char*)&readLengths[0], sizeof(int)*readLengths.size());
		outfile.close();
	}

	void ReadModel(string &inName) {
		ifstream infile(inName.c_str(), std::ios::in|std::ios::binary);
		string tmp;
		tmp.resize(desc.size());
		infile.read((char*)&tmp[0], sizeof(char)*desc.size());

		if (tmp != desc) {
			cout << "ERROR, model file does not appear to be an Alechemy model " << desc << endl;
			exit(1);
		}
		infile.read((char*) &k, sizeof(int));
		int n;
		infile.read((char*) &n, sizeof(int));
		cerr << "Alchemy2 reading " << n << " contexts" << endl;
		for (int i=0; i < n; i++) {
			string kmer;
			kmer.resize(k);
			infile.read((char*)&kmer[0], sizeof(char)*k);
			ContextOutput context;
			context.Read(infile);
			model[kmer] = context;
		}
		infile.read((char*)&avgReadLength, sizeof(int));
		infile.read((char*)&readVar, sizeof(int));			
		infile.read((char*)&avgLogReadLength, sizeof(double));
		infile.read((char*)&readLogVar, sizeof(double));
		int nrl;
		infile.read((char*)&nrl, sizeof(int));
		readLengths.resize(nrl);
		infile.read((char*)&readLengths[0], sizeof(int)*nrl);
		infile.close();
	}
	void Increment(string kmer, int op, int length, string seq="") {
		nSamples++;
		model[kmer].Increment(op, length, seq);		
	}

	void StoreRead(bam1_t *b, bam_hdr_t *h, Genome &genome) {
		if (b->core.l_qseq == 0) {
			return;
		}
		readLengths.push_back(b->core.l_qseq);
		int tStart=b->core.pos;
		int qPos=0;
		int tPos=tStart;
		int nCigar = b->core.n_cigar;
		uint32_t* cigar= bam_get_cigar(b);
		vector<int> posmap;
		bool started=false;
		int seqLen = b->core.l_qseq;
		if (seqLen < k) {
			return;
		}
		int i =0;
		for (i=0; i < nCigar; i++) {
			int op = cigar[i] & 0xF;
			int opLen = cigar[i] >> 4;
			if (op == 4 ) {
				qPos+= cigar[i] >> 4;
			}
			else if ( op == 0 or 
								op == 7 or
								op == 8 ) {				
				for (int j=0; j < opLen; j++) {
					posmap.push_back(qPos++);
					tPos++;
				}
			}
			else if ( op == 2 ) {
				for (int j = 0; j < opLen; j++) {
					posmap.push_back(qPos);
				}
				tPos += opLen;
			}
			else if ( op == 1 ) {
				qPos += opLen;
			}
		}
		
		//
		// Now sample k-mers from the reference
		//
		assert(k > 0);
		assert(k % 2 == 1);
		int chrom=genome.GetIndex(h->target_name[b->core.tid]);
		char *ref=&genome.seqs[chrom][tStart];
		string kmer = string(ref,k);
		for (int i=0; i < kmer.size(); i++) {
			kmer[i] = toupper(kmer[i]);
		}
		int seqlen=b->core.l_qseq;
		char *seq = new char[b->core.l_qseq];



		uint8_t *q = bam_get_seq(b);
		for (int i=0; i < seqlen; i++) {seq[i]=seq_nt16_str[bam_seqi(q,i)];	}
		for (int i = k/2+1; i < posmap.size() - (k/2+1); i++) {
			kmer=string(&ref[i-(k/2)], k);
			for (int j=0; j < k; j++) {
				kmer[j] = toupper(ref[i-(k/2)+j]);
			}
			
			if (posmap[i] == posmap[i-1]) {
				//
				// Found a deletion.
				//
				int j=i;
				while (j < seqLen and posmap[j] == posmap[i-1]) {
					j++;
				}
				if (j-i < maxGap) {
					if (kmer.find('N', 0) == kmer.npos) {
						Increment(kmer, 1, j-i, "");				
					}
				}
				i=j;
			}
			else {
				if (posmap[i+1]-posmap[i] < maxGap) {
					if (kmer.find('N', 0) == kmer.npos and posmap[i+1] > posmap[i]) {
						Increment(kmer, 0, posmap[i+1]-posmap[i], string(&seq[posmap[i]], posmap[i+1]-posmap[i]));

					}
				}
				i+= max(0, posmap[i+1] -posmap[i]-1);
			}				
		}
		delete[] seq;
	}
};



int main(int argc, char* argv[]) {

  opterr = 0;
	string bamFile="";
	int context=0;
	int samples=1000;
	string genomeFile="";
	string modelFile="", readsFile="", posFile="";
	string outFile="";
	string bedFile="";
	int numBases=0;
	int numReads=0;
	int maxGap=50;
	double avgReadLen=0;
	double readVar=0;
	char c;
	int maxSampledReads=50000;
	bool useEmpiricalReadLengths=false;
	int minReadLength=1000;
	if (argc == 1) {
		HelpAlchemy2();
		exit(1);
	}
  while ((c = getopt (argc, argv, "b:k:L:G:s:g:m:r:R:B:l:N:p:o:u:V:")) != -1) {
    switch (c)
      {
      case 'b':
				bamFile=optarg;
        break;
      case 'L':
        maxSampledReads=atoi(optarg);
        break;
      case 'k':
        context=atoi(optarg);
        break;
      case 'E':
				useEmpiricalReadLengths=true;
				break;
      case 'u':
        avgReadLen=atof(optarg);
				break;
			case 'V':
        readVar=atof(optarg);
				break;
      case 's':
        samples=atoi(optarg);
        break;
			case 'G':
				maxGap=atoi(optarg);
				break;
      case 'g':
				genomeFile=optarg;
				 break;
			 case 'm':
				 modelFile=optarg;
				 break;
			 case 'r':
				 readsFile=optarg;
				 break;
			case 'R':
				bedFile=optarg;
				break;
			case 'o':
				outFile=optarg;
				break;
			 case 'B':
				 numBases=ParseNumber(optarg);
				 break;
			 case '?':
				 HelpAlchemy2();
				 return 1;
			 default:
				 abort ();
			 }
	 }

	if (bamFile != "") {

		cerr << "Storing model" << endl;
		if (context == 0 or samples==0 or modelFile=="" or genomeFile=="") {
			HelpAlchemy2();
			cout << "Error. When storing a model, context k-mer (-k), number of samples (-s), output model (-m), and genome (-g) must be specified." << endl;
			exit(0);
		}
		Genome genome;
		genome.Read(genomeFile);

		htsFile *htsfp;
		htsfp = hts_open(bamFile.c_str(),"r");
		const htsFormat *fmt = hts_get_format(htsfp);
		if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
			cout << "Cannot determine format of input reads." << endl;
			exit(1);
		}


		bam_hdr_t *samHeader;			
		samHeader = sam_hdr_read(htsfp);

		Model model;
		model.Init(context, maxGap);
		bam1_t *b = bam_init1();
		int r=0;
		long nBase=0;
		while (true) {
			int res= sam_read1(htsfp, samHeader, b);
			if (res == 0) { break;}
			r++;
			if (r % 1000 == 0) {
				int nComplete=0;
				for (	map<string, ContextOutput>::iterator mit=model.model.begin(); 
							mit != model.model.end(); ++mit) {
					if (mit->second.size > samples) { nComplete++; }
				}
				cerr << "alchemy2 model: Parsed " << r << " alignments. " << model.model.size() 
						 << "/" << (1<<(model.k*2) )
						 << " kmers seen: " << std::setprecision(3) << 100*float(nComplete) / (1 << (model.k*2)) << "% sufficiently sampled." << endl;
				
			}
			model.StoreRead(b, samHeader, genome);
			nBase+= b->core.l_qseq;
			int large=9999999;
			int minCount=large;
			int avgCount=0;
			string minContext;
			map<string, ContextOutput>::iterator minIt = model.model.end();

			for (	map<string, ContextOutput>::iterator mit=model.model.begin(); 
						mit != model.model.end(); ++mit) {
				minCount = min(minCount, mit->second.size);				
				if (minCount == mit->second.size) {
					minContext = mit->first;
					minIt=mit;
				}
					avgCount+= mit->second.size;
			}
			if ((minCount != large and minCount >= samples) or r > maxSampledReads ) {
				cerr << "done " << minCount << endl;
				break;
			}
		}
		cerr << "Storing model " << model.nSamples << endl;
		model.WriteModel(modelFile);		
	}
	else {
		//
		// Running simulation, not storing model.
		//
		if (modelFile=="" or genomeFile == "") {
			HelpAlchemy2();
			cout << "Error. A model file (-m file) and genome file (-g genome) must be specified for simulations."<< endl;
			exit(1);
		}
		vector<int> testv;
		Model model;
		model.ReadModel(modelFile);
		Genome genome;
		genome.Read(genomeFile);

		vector<long> cumRegionLength;
		vector<string> chroms;
		vector<long> start, end;
		long c=0;
		long maxPos=0;
		if (bedFile == "") {
			for (int i=0; i < genome.lengths.size(); i++) {
				cumRegionLength.push_back(maxPos);
				maxPos += genome.lengths[i];
				chroms.push_back(genome.header.names[i]);
				start.push_back(0);
				end.push_back(genome.lengths[i]);
			}
		}
		else {
			ifstream bed(bedFile.c_str());
			string line;
			while (getline(bed, line)) {
				if (line=="") { break;}
				stringstream strm(line);
				string chrom;
				int s, e;				
				strm >> chrom >> s >> e;
				cumRegionLength.push_back(maxPos);
				maxPos+=e-s;
				chroms.push_back(chrom);
				start.push_back(s);
				end.push_back(e);
			}
		}

		ofstream readsFile(outFile.c_str());
		
		int nSimulatedBases=0;
		if (avgReadLen == 0) {
			avgReadLen = model.avgLogReadLength;
		}
		if (readVar == 0) {
			readVar = model.readLogVar;
		}
			
		std::default_random_engine generator;
		std::lognormal_distribution<double> distribution(avgReadLen, readVar);
		int readIndex=0;

		while (nSimulatedBases < numBases) {
			long i=RandInt(maxPos);
			int idx=max(0,(int)(lower_bound(cumRegionLength.begin(), cumRegionLength.end(), i)-cumRegionLength.begin() - 1));
			int offset=i-cumRegionLength[idx];

			int refIdx = idx;
			// Find offset into region
			string chrom=chroms[idx];
			if (bedFile != "") {
				// If simulating from regions, 
				offset+=start[idx];
				refIdx=genome.nameMap[chrom];
			}

			int maxLen=end[idx]-offset;

			int readLen;
			
			if (useEmpiricalReadLengths) {
				readLen = model.readLengths[RandInt(model.readLengths.size())];
			}
			else {
				readLen = distribution(generator);
			}
			readLen = min(maxLen, readLen);
			nSimulatedBases += readLen;
			if (readLen < minReadLength) {
				continue;
			}
			//			cout << nSimulatedBases << "\t" << chroms[idx] << "\t" << offset << "\t" << readLen << endl; 
			if (readLen < model.k) {
				continue;
			}
			char *ref=&genome.seqs[refIdx][offset];
			string sim="";
			string kmer;
			kmer.resize(model.k);
			for (int i=model.k/2; i < readLen-model.k/2; i++) {
				string simOut="";
				if (ref[i] == 'N') {
					simOut="N";
				}
				else {
					for (int j=0;j<model.k; j++) {
						kmer[j] = (char)toupper(ref[i-model.k/2+j]);
						if (kmer[j] == 'N') { kmer[j] = 'A';}
					}
					simOut=model.Sample(kmer);
				}
				if (IsDel(simOut)) {
					i+= ParseDel(simOut);
				}
				else {
					sim+=simOut;
				}
			}
			stringstream title;
			title << ">" << chrom << ":" << offset << "-" << offset+readLen << "/" << readIndex;
			readsFile << title.str() << endl;
			readsFile << sim << endl;			
			readIndex++;
			if (readIndex % 1000 == 0) {
				cerr << "Alchemy2: simulated " << nSimulatedBases << " bases." << endl;
			}
		}
		//
		// Done simulating
		//
		exit(0);
	}
	
	exit(0);
}
