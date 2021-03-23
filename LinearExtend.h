#ifndef LINEAR_EXTEND_H
#define LINEAR_EXTEND_H

#include <vector>
#include "TupleOps.h"
#include "SplitClusters.h"
#include "Clustering.h"

using namespace std;

class LongAnchors {
public:
	vector<int> anchorIndex;
	GenomePairs *matches;
	vector<int> *matchesLengths;
	int strand;

  	LongAnchors(GenomePairs *ms, vector<int> *mLs, int str) {
  		matches=ms;
  		strand=str;
  		matchesLengths=mLs;
	}	
	//
	// Cartesian sort of clusters. NOTE: if strand==1, sort on reverse
	//
	int operator()(const int i, const int j) {
		if (strand == 0) {
			if ((*matches)[i].first.pos != (*matches)[j].first.pos) {
				return (*matches)[i].first.pos < (*matches)[j].first.pos;
			}
			else {
				return (*matches)[i].second.pos < (*matches)[j].second.pos;
			}			
		}
		else { 
			if ((*matches)[i].first.pos + (*matchesLengths)[i] != (*matches)[j].first.pos + (*matchesLengths)[j]) {
				return (*matches)[i].first.pos + (*matchesLengths)[i] > (*matches)[j].first.pos + (*matchesLengths)[j];
			}
			else {
				return (*matches)[i].second.pos < (*matches)[j].second.pos;
			}				
		}
	}
	void Sort() {
		sort(anchorIndex.begin(), anchorIndex.end(), *this);
	}
};


void Checkbp(GenomePair &cur, GenomePair &next, Genome &genome, Read &read, int &ChromIndex, GenomePos& qe, GenomePos& te, 
				bool strand, Options &opts, int &mat, int K) {
	GenomePos curQ, curT, nextQ, nextT;
	if (strand == 0) {
		curQ = cur.first.pos + K;
		curT = cur.second.pos + K;
		nextQ = next.first.pos;
		nextT = next.second.pos;		
	}
	else {
		curQ = cur.first.pos + K;
		curT = cur.second.pos - 1;
		nextQ = next.first.pos;
		nextT = next.second.pos + K - 1;		
	}

	if (strand == 0) { //curT - genome.header.pos[ChromIndex]
		while (nextQ > curQ and nextT > curT and genome.seqs[ChromIndex][curT] == read.seq[curQ] ) {
			mat++;
			curQ++;
			curT++;
		}
	}
	else {
		while (nextQ > curQ and nextT < curT and genome.seqs[ChromIndex][curT] == read.seq[curQ] ) {
			mat++;
			curQ++;
			curT--;
		}
	}
	qe = curQ; te = curT;
}


int
CheckOverlap(GenomePair & match, Options & opts, vector<pair<GenomePos, bool>> & Set, bool & Ovp, int K) {

	for (int t = 0; t < Set.size(); t++) {
		if (Set[t].second == 0 and Set[t].first >= match.first.pos and Set[t].first < match.first.pos + K) {
			Ovp = 1;
			return 0;
		} 
		if (Set[t].second == 1 and Set[t].first >= match.second.pos and Set[t].first < match.second.pos + K) {
			Ovp = 1;
			return 0;
		}
	}
	return 0;
}


void 
DecideCoordinates (Cluster & cluster, bool strand, int chromIndex, float anchorfreq) {
	if (cluster.matches.size() == 0) return;
	GenomePos qStart, qEnd, tStart, tEnd;
	qStart = cluster.matches[0].first.pos;
	qEnd = qStart + cluster.matchesLengths[0];
	tStart = cluster.matches[0].second.pos;
	tEnd = tStart + cluster.matchesLengths[0];

	for (int n = 1; n < cluster.matches.size(); n++) {
		qStart = min(qStart, cluster.matches[n].first.pos);
		qEnd   = max(qEnd, cluster.matches[n].first.pos + cluster.matchesLengths[n]);
		tStart = min(tStart, cluster.matches[n].second.pos);
		tEnd   = max(tEnd, cluster.matches[n].second.pos + cluster.matchesLengths[n]);
	}
	cluster.qStart = qStart;
	cluster.qEnd = qEnd;
	cluster.tStart = tStart;
	cluster.tEnd = tEnd;
	cluster.strand = strand;
	cluster.chromIndex = chromIndex;
	cluster.anchorfreq = anchorfreq;
	return;
}


//
// This function implements linear extension for each Cluster on the chain;
// NOTE: The extension for anchors should avoid overlapping;
//
template<typename Tup>
void
LinearExtend(vector<Cluster*> clusters, vector<Cluster> & extCluster, vector<Tup> & chain, Options & opts, Genome & genome, Read & read, 
			int start, int &overlap, bool skiprepetitive, int K) {

	vector<pair<GenomePos, bool>> Set;
	int next, prev;
	for (int c = 0; c < chain.size(); c++) {
		int mat = 0;
		int cm = chain[c];
		if (clusters[cm]->matches.size() == 0) continue;

		Set.clear();
		prev = c - 1;
		next = c + 1;
		GenomePos qsb = clusters[cm]->qStart;
		GenomePos qeb = clusters[cm]->qEnd;
		GenomePos tsb = clusters[cm]->tStart;
		GenomePos teb = clusters[cm]->tEnd;
		//
		// Compute the ChromIndex for cluster
		//
		extCluster[c + start].chromIndex = clusters[cm]->chromIndex;
		extCluster[c + start].anchorfreq = clusters[cm]->anchorfreq;
		// cerr << "clusters[cm]->anchorfreq: " << clusters[cm]->anchorfreq << endl;
		
		if (skiprepetitive and extCluster[c + start].anchorfreq <= 1.1f ) { // Do not check overlap for repetitive clusters;
			//
			// Check the previous Cluster to get the overlapping points;
			//			
			GenomePos q1, q2, t1, t2;
			if (prev != -1) {
				q1 = clusters[chain[prev]]->qStart;
				q2 = clusters[chain[prev]]->qEnd;
				if (q1 > qsb and q1 < qeb) Set.push_back(make_pair(q1, 0)); 
				if (q2 > qsb and q2 < qeb) Set.push_back(make_pair(q2, 0));

				t1 = clusters[chain[prev]]->tStart;
				t2 = clusters[chain[prev]]->tEnd;
				if (t1 > tsb and t1 < teb) Set.push_back(make_pair(t1, 1)); 
				if (t2 > tsb and t2 < teb) Set.push_back(make_pair(t2, 1)); 
			} 

			//
			// Check the next Cluster to get the overlapping points;
			//
			if (next < chain.size()) {
				q1 = clusters[chain[next]]->qStart;
				q2 = clusters[chain[next]]->qEnd;
				if (q1 > qsb and q1 < qeb) Set.push_back(make_pair(q1, 0)); 
				if (q2 > qsb and q2 < qeb) Set.push_back(make_pair(q2, 0));

				t1 = clusters[chain[next]]->tStart;
				t2 = clusters[chain[next]]->tEnd;
				if (t1 > tsb and t1 < teb) Set.push_back(make_pair(t1, 1)); 
				if (t2 > tsb and t2 < teb) Set.push_back(make_pair(t2, 1)); 			
			}
		}
		

		//
		// Sort each Cluster
		// NOTICE: if the current Cluster get more anchors in the btwn space or end space, we need to sort them. Otherwise we do not sort them;
		//// TODO(Jingwen): Check whether do we need to sort for every Cluster!!!!!!!!!!!!!!!
		//
		
		if (clusters[cm]->strand == 0) { 
			//if (clusters[cm]->refinespace == 1) DiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
			DiagonalSort<GenomeTuple>(clusters[cm]->matches, 500);
			// DiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
		}
		else {
			//if (clusters[cm]->refinespace == 1) AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
			AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches, 500);
			// AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end(), read.length);
		}
		//
		// Linear Extension;
		// NOTICE: anchors in "clusters[cm]" have the same strand;
		//
		int n = 1;
		int m = 0;
		bool chm = 1; // chm == 1 means we need to check whether m overlaps with "Set";
		while (n < clusters[cm]->matches.size()) {
			//
			// Check whether clusters[cm]->matches[m] overlaps with "Set";
			//
			bool Ovp = 0;
			if (chm == 1) {
				CheckOverlap(clusters[cm]->matches[m], opts, Set, Ovp, K);

				if (Ovp == 1) {
					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
					extCluster[c + start].matchesLengths.push_back(K);	
					extCluster[c + start].overlap.push_back(1);
					overlap++;
					// cerr << "1: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[m].second.pos << endl;
					m = n;
					n++;		
					chm = 1;
					continue;	
				}
				else chm = 0;
			}

			//
			// Check whether clusters[cm]->matches[n] overlaps with "Set";
			//
			Ovp = 0;
			CheckOverlap(clusters[cm]->matches[n], opts, Set, Ovp, K);
			
			if (Ovp == 1) { // clusters[cm]->matches[n] overlaps with "Set";
				if (clusters[cm]->strand == 0) {
					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));					
				}
				else {
					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));	
					// cerr << "2: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
				}
				extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + K - clusters[cm]->matches[m].first.pos);
				extCluster[c + start].overlap.push_back(0);
				extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[n].first.pos), GenomeTuple(0, clusters[cm]->matches[n].second.pos)));
				extCluster[c + start].matchesLengths.push_back(K);
				extCluster[c + start].overlap.push_back(1);
				overlap++;
				// cerr << "3: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[n].first.pos << " " << clusters[cm]->matches[n].second.pos << endl;
				m = n + 1;
				n = m + 1;
				chm = 1;
				continue;
			}

			//
			// clusters[cm]->matches[n] does not overlap with "Set";
			//
			int64_t curDiag, nextDiag;
			if (clusters[cm]->strand == 0) {
				curDiag = (int64_t) clusters[cm]->matches[n-1].first.pos - (int64_t) clusters[cm]->matches[n-1].second.pos;
				nextDiag = (int64_t) clusters[cm]->matches[n].first.pos - (int64_t) clusters[cm]->matches[n].second.pos;
			}
			else {
				curDiag = (int64_t) clusters[cm]->matches[n-1].first.pos + (int64_t) clusters[cm]->matches[n-1].second.pos;
				nextDiag = (int64_t) clusters[cm]->matches[n].first.pos + (int64_t) clusters[cm]->matches[n].second.pos;			
			}

			if (curDiag == nextDiag) { // those two anchors are on the same diagonal

				if (clusters[cm]->matches[n].first.pos < clusters[cm]->matches[n-1].first.pos + K) { // anchor n-1 and anchor n are overlapped
					n++;
				}
				else { //  anchor n-1 and anchor n are not overlapped; Need to extend after anchor n-1
					GenomePos qe, te;
					Checkbp(clusters[cm]->matches[n-1], clusters[cm]->matches[n], genome, read, clusters[cm]->chromIndex, qe, te, clusters[cm]->strand, opts, mat, K);
					if (clusters[cm]->strand == 0 and qe == clusters[cm]->matches[n].first.pos and te == clusters[cm]->matches[n].second.pos) {
						n++;
					}
					else if (clusters[cm]->strand == 1 and qe == clusters[cm]->matches[n].first.pos and te == clusters[cm]->matches[n].second.pos + K - 1) {
						n++;
					}
					else {
						if (clusters[cm]->strand == 0) {
							extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
						}
						else {
							extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, te + 1)));
							// cerr << "4: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << te + 1 << endl;

						}
						extCluster[c + start].matchesLengths.push_back(qe - clusters[cm]->matches[m].first.pos);
						extCluster[c + start].overlap.push_back(0);
						m = n;
						n++;
					}
				}
			}
			else {
				if (clusters[cm]->strand == 0) {
					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));					
				}
				else {
					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));
					// cerr << "5: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
				}
				extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + K - clusters[cm]->matches[m].first.pos);
				extCluster[c + start].overlap.push_back(0);
				m = n;
				n++;
			}
			chm = 0;	
		}

		if (n == clusters[cm]->matches.size()) { // and m != n-1
			if (clusters[cm]->strand == 0) {
				extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
			}
			else {
				extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));
				// cerr << "6: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
			}
			extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + K - clusters[cm]->matches[m].first.pos);
			extCluster[c + start].overlap.push_back(0);
		}

		// extCluster[c + start].strand = clusters[cm]->strand;
		// extCluster[c + start].chromIndex = clusters[cm]->chromIndex;
		// extCluster[c + start].anchorfreq = clusters[cm]->anchorfreq;
		DecideCoordinates(extCluster[c + start], clusters[cm]->strand, clusters[cm]->chromIndex, clusters[cm]->anchorfreq);
		//cerr <<"c: " << c << "   mat:" << mat << endl;
		assert(clusters[cm]->matches.size() >= extCluster[c + start].matches.size());
		assert(extCluster[c + start].matches.size() == extCluster[c + start].overlap.size());
	}
}

// //
// // This function implements linear extension for each Cluster on the chain;
// // NOTE: The extension for anchors should avoid overlapping;
// //
// template<typename Tup>
// void
// LinearExtend_MergedCluster(vector<Merge_SplitChain> &mergeinfo, vector<Cluster*> clusters, vector<Cluster> & extCluster, Options & opts, Genome & genome, Read & read, 
// 			int start, int &overlap, bool skiprepetitive) {

// 	vector<pair<GenomePos, bool>> Set;
// 	int next, prev;
// 	for (int c = 0; c < mergeinfo.size(); c++) {
// 		int mat = 0;
// 		int cm = mergeinfo[c][m];
// 		int f = mergeinfo[c].back();
// 		int l = mergeinfo[c][0];
// 		// if (clusters[cm]->matches.size() == 0) continue;

// 		Set.clear();
// 		prev = c - 1;
// 		next = c + 1;
// 		GenomePos qsb = mergeinfo[c].qStart; 
// 		GenomePos qeb = mergeinfo[c].qEnd;
// 		GenomePos tsb = mergeinfo[c].tStart;
// 		GenomePos teb = mergeinfo[c].tEnd;
// 		//
// 		// Compute the ChromIndex for cluster
// 		//
// 		extCluster[c + start].chromIndex = clusters[l]->chromIndex;
		
// 		//
// 		// Check the previous Cluster to get the overlapping points;
// 		//			
// 		GenomePos q1, q2, t1, t2;
// 		if (prev != -1) {
// 			q1 = mergeinfo[prev].qStart;
// 			q2 = mergeinfo[prev].qEnd;
// 			if (q1 > qsb and q1 < qeb) Set.push_back(make_pair(q1, 0)); 
// 			if (q2 > qsb and q2 < qeb) Set.push_back(make_pair(q2, 0));

// 			t1 = mergeinfo[prev].tStart;
// 			t2 = mergeinfo[prev].tEnd;
// 			if (t1 > tsb and t1 < teb) Set.push_back(make_pair(t1, 1)); 
// 			if (t2 > tsb and t2 < teb) Set.push_back(make_pair(t2, 1)); 
// 		} 

// 		//
// 		// Check the next Cluster to get the overlapping points;
// 		//
// 		if (next < chain.size()) {
// 			q1 = mergeinfo[cur].qStart;
// 			q2 = clusters[chain[next]]->qEnd;
// 			if (q1 > qsb and q1 < qeb) Set.push_back(make_pair(q1, 0)); 
// 			if (q2 > qsb and q2 < qeb) Set.push_back(make_pair(q2, 0));

// 			t1 = clusters[chain[next]]->tStart;
// 			t2 = clusters[chain[next]]->tEnd;
// 			if (t1 > tsb and t1 < teb) Set.push_back(make_pair(t1, 1)); 
// 			if (t2 > tsb and t2 < teb) Set.push_back(make_pair(t2, 1)); 			
// 		}

// 		for (int m = 0; m < mergeinfo[c].size(); m++) {
// 			int cm = mergeinfo[c][m];
// 			extCluster[c + start].anchorfreq = clusters[cm]->anchorfreq;

// 			//
// 			// Sort each Cluster
// 			// NOTICE: if the current Cluster get more anchors in the btwn space or end space, we need to sort them. Otherwise we do not sort them;
// 			//// TODO(Jingwen): Check whether do we need to sort for every Cluster!!!!!!!!!!!!!!!
// 			//
			
// 			if (clusters[cm]->strand == 0) { 
// 				//if (clusters[cm]->refinespace == 1) DiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
// 				DiagonalSort<GenomeTuple>(clusters[cm]->matches, 500);
// 				// DiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
// 			}
// 			else {
// 				//if (clusters[cm]->refinespace == 1) AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end());
// 				AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches, 500);
// 				// AntiDiagonalSort<GenomeTuple>(clusters[cm]->matches.begin(), clusters[cm]->matches.end(), read.length);
// 			}
// 			//
// 			// Linear Extension;
// 			// NOTICE: anchors in "clusters[cm]" have the same strand;
// 			//
// 			int n = 1;
// 			int m = 0;
// 			bool chm = 1; // chm == 1 means we need to check whether m overlaps with "Set";
// 			while (n < clusters[cm]->matches.size()) {
// 				//
// 				// Check whether clusters[cm]->matches[m] overlaps with "Set";
// 				//
// 				bool Ovp = 0;
// 				if (chm == 1) {
// 					CheckOverlap(clusters[cm]->matches[m], opts, Set, Ovp);

// 					if (Ovp == 1) {
// 						extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
// 						extCluster[c + start].matchesLengths.push_back(opts.globalK);	
// 						extCluster[c + start].overlap.push_back(1);
// 						overlap++;
// 						// cerr << "1: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[m].second.pos << endl;
// 						m = n;
// 						n++;		
// 						chm = 1;
// 						continue;	
// 					}
// 					else chm = 0;
// 				}

// 				//
// 				// Check whether clusters[cm]->matches[n] overlaps with "Set";
// 				//
// 				Ovp = 0;
// 				CheckOverlap(clusters[cm]->matches[n], opts, Set, Ovp);
				
// 				if (Ovp == 1) { // clusters[cm]->matches[n] overlaps with "Set";
// 					if (clusters[cm]->strand == 0) {
// 						extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));					
// 					}
// 					else {
// 						extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));	
// 						// cerr << "2: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
// 					}
// 					extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + opts.globalK - clusters[cm]->matches[m].first.pos);
// 					extCluster[c + start].overlap.push_back(0);
// 					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[n].first.pos), GenomeTuple(0, clusters[cm]->matches[n].second.pos)));
// 					extCluster[c + start].matchesLengths.push_back(opts.globalK);
// 					extCluster[c + start].overlap.push_back(1);
// 					overlap++;
// 					// cerr << "3: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[n].first.pos << " " << clusters[cm]->matches[n].second.pos << endl;
// 					m = n + 1;
// 					n = m + 1;
// 					chm = 1;
// 					continue;
// 				}

// 				//
// 				// clusters[cm]->matches[n] does not overlap with "Set";
// 				//
// 				int64_t curDiag, nextDiag;
// 				if (clusters[cm]->strand == 0) {
// 					curDiag = (int64_t) clusters[cm]->matches[n-1].first.pos - (int64_t) clusters[cm]->matches[n-1].second.pos;
// 					nextDiag = (int64_t) clusters[cm]->matches[n].first.pos - (int64_t) clusters[cm]->matches[n].second.pos;
// 				}
// 				else {
// 					curDiag = (int64_t) clusters[cm]->matches[n-1].first.pos + (int64_t) clusters[cm]->matches[n-1].second.pos;
// 					nextDiag = (int64_t) clusters[cm]->matches[n].first.pos + (int64_t) clusters[cm]->matches[n].second.pos;			
// 				}

// 				if (curDiag == nextDiag) { // those two anchors are on the same diagonal

// 					if (clusters[cm]->matches[n].first.pos < clusters[cm]->matches[n-1].first.pos + opts.globalK) { // anchor n-1 and anchor n are overlapped
// 						n++;
// 					}
// 					else { //  anchor n-1 and anchor n are not overlapped; Need to extend after anchor n-1
// 						GenomePos qe, te;
// 						Checkbp(clusters[cm]->matches[n-1], clusters[cm]->matches[n], genome, read, clusters[cm]->chromIndex, qe, te, clusters[cm]->strand, opts, mat);
// 						if (clusters[cm]->strand == 0 and qe == clusters[cm]->matches[n].first.pos and te == clusters[cm]->matches[n].second.pos) {
// 							n++;
// 						}
// 						else if (clusters[cm]->strand == 1 and qe == clusters[cm]->matches[n].first.pos and te == clusters[cm]->matches[n].second.pos + opts.globalK - 1) {
// 							n++;
// 						}
// 						else {
// 							if (clusters[cm]->strand == 0) {
// 								extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
// 							}
// 							else {
// 								extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, te + 1)));
// 								// cerr << "4: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << te + 1 << endl;

// 							}
// 							extCluster[c + start].matchesLengths.push_back(qe - clusters[cm]->matches[m].first.pos);
// 							extCluster[c + start].overlap.push_back(0);
// 							m = n;
// 							n++;
// 						}
// 					}
// 				}
// 				else {
// 					if (clusters[cm]->strand == 0) {
// 						extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));					
// 					}
// 					else {
// 						extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));
// 						// cerr << "5: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
// 					}
// 					extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + opts.globalK - clusters[cm]->matches[m].first.pos);
// 					extCluster[c + start].overlap.push_back(0);
// 					m = n;
// 					n++;
// 				}
// 				chm = 0;	
// 			}

// 			if (n == clusters[cm]->matches.size()) { // and m != n-1
// 				if (clusters[cm]->strand == 0) {
// 					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[m].second.pos)));
// 				}
// 				else {
// 					extCluster[c + start].matches.push_back(GenomePair(GenomeTuple(0, clusters[cm]->matches[m].first.pos), GenomeTuple(0, clusters[cm]->matches[n-1].second.pos)));
// 					// cerr << "6: " << extCluster[c + start].matches.size() - 1 << " " << clusters[cm]->matches[m].first.pos << " " << clusters[cm]->matches[n-1].second.pos << endl;
// 				}
// 				extCluster[c + start].matchesLengths.push_back(clusters[cm]->matches[n-1].first.pos + opts.globalK - clusters[cm]->matches[m].first.pos);
// 				extCluster[c + start].overlap.push_back(0);
// 			}

// 			extCluster[c + start].strand = clusters[cm]->strand;
// 			DecideCoordinates(extCluster[c + start]);
// 			//cerr <<"c: " << c << "   mat:" << mat << endl;
// 			assert(clusters[cm]->matches.size() >= extCluster[c + start].matches.size());
// 			assert(extCluster[c + start].matches.size() == extCluster[c + start].overlap.size());
// 		}
// 	}
// }

//
// This function trim overlapped long anchors.
//
void 
TrimOverlappedAnchors(vector<Cluster> & extCluster, int start) {

	//
	// Get all long anchors
	//
	for (int c = 0 + start; c < extCluster.size(); c++) {
		LongAnchors longanchors(&(extCluster[c].matches), &(extCluster[c].matchesLengths), extCluster[c].strand);
		//LongAnchors longanchors(extCluster[c].matches&extCluster[c]);
		for (int ln = 0; ln < extCluster[c].matches.size(); ln++) {
			if (extCluster[c].matchesLengths[ln] >= 40) { // 50
				longanchors.anchorIndex.push_back(ln);
			}
		}

		//
		// Catersiansort anchors and trim overlapped adjacent two;
		//
		longanchors.Sort();
		for (int ln = 1; ln < longanchors.anchorIndex.size(); ln++) {

			int prev = longanchors.anchorIndex[ln-1];
			int cur = longanchors.anchorIndex[ln];

			//
			// If the adjacent anchors are overlapped less than 30b, trim the prev
			// 
			int overlap_r = 0, overlap_g = 0, overlap = 0;
			bool lr = 0, lg = 0;
			//
			// Check if two anchors are overlapped on read;
			//
			if (extCluster[c].strand == 0) {
				if (extCluster[c].matches[cur].first.pos < extCluster[c].matches[prev].first.pos + extCluster[c].matchesLengths[prev] 
					and 
					extCluster[c].matches[cur].first.pos >= extCluster[c].matches[prev].first.pos + extCluster[c].matchesLengths[prev] - 30) {

					overlap_r = extCluster[c].matches[prev].first.pos + extCluster[c].matchesLengths[prev] - extCluster[c].matches[cur].first.pos;
					assert(overlap_r <= 60); //assert(overlap_r > 0);//30
					lr = 1;
				}				
			}
			else {
				if (extCluster[c].matches[cur].first.pos + extCluster[c].matchesLengths[cur] > extCluster[c].matches[prev].first.pos
					and 
					extCluster[c].matches[cur].first.pos + extCluster[c].matchesLengths[cur] <= extCluster[c].matches[prev].first.pos + 30) {

					overlap_r = extCluster[c].matches[cur].first.pos + extCluster[c].matchesLengths[cur] - extCluster[c].matches[prev].first.pos;
					assert(overlap_r <= 60); //assert(overlap_r > 0);
					lr = 1;
				}				
			}
			//
			// Check if two anchors are overlapped on genome;
			//
			if (extCluster[c].matches[cur].second.pos < extCluster[c].matches[prev].second.pos + extCluster[c].matchesLengths[prev] 
				and 
				extCluster[c].matches[cur].second.pos >= extCluster[c].matches[prev].second.pos + extCluster[c].matchesLengths[prev] - 30) {

				overlap_g = extCluster[c].matches[prev].second.pos + extCluster[c].matchesLengths[prev] - extCluster[c].matches[cur].second.pos;
				assert(overlap_g <= 60); assert(overlap_g > 0);
				lg = 1;
			}
			
			if (overlap_r > 0 or overlap_g > 0) {
				overlap = max(overlap_r, overlap_g);
				if (extCluster[c].strand == 1) extCluster[c].matches[prev].first.pos += overlap+1;
				extCluster[c].matchesLengths[prev] -= overlap+1;
				if (lr == 1 and extCluster[c].strand == 0) assert(extCluster[c].matches[cur].first.pos >= extCluster[c].matches[prev].first.pos + extCluster[c].matchesLengths[prev]);
				if (lr == 1 and extCluster[c].strand == 1) assert(extCluster[c].matches[cur].first.pos + extCluster[c].matchesLengths[cur] <= extCluster[c].matches[prev].first.pos);
				if (lg == 1) assert(extCluster[c].matches[cur].second.pos >= extCluster[c].matches[prev].second.pos + extCluster[c].matchesLengths[prev]);
			
			}
		}
	}
}



//
// This function just do linear extension without avoiding overlapping points;
// Input: GenomePairs pairs; Output: GenomePairs ExtendPairs; vector<int> ExtendPairsMatchesLength;
// NOTE: Only forward strand;
//
void LinearExtend(GenomePairs * pairs, GenomePairs &Extendpairs, vector<int> &ExtendpairsMatchesLength, Options &opts, 
				Genome &genome, Read &read, int chromIndex, bool strand, bool skipsorting, int K) {
	//
	// Sort each Cluster
	//
	if (!skipsorting) DiagonalSort<GenomeTuple>((*pairs).begin(), (*pairs).end());
	//
	// Linear Extension;
	// NOTICE: anchors have the same strand;
	//
	int n = 1;
	int m = 0;
	// int strand = 0;
	int mat = 0;
	while (n < (*pairs).size()) {

		int64_t curDiag, nextDiag;
		if (strand == 0) {
			curDiag = (int64_t) (*pairs)[n-1].first.pos - (int64_t) (*pairs)[n-1].second.pos;
			nextDiag = (int64_t) (*pairs)[n].first.pos - (int64_t) (*pairs)[n].second.pos;			
		}
		else {
			curDiag = (int64_t) (*pairs)[n-1].first.pos + (int64_t) (*pairs)[n-1].second.pos;
			nextDiag = (int64_t) (*pairs)[n].first.pos + (int64_t) (*pairs)[n].second.pos;			
		}


		if (curDiag == nextDiag) { // those two anchors are on the same diagonal
			if ((*pairs)[n].first.pos < (*pairs)[n-1].first.pos + K) { // anchor n-1 and anchor n are overlapped
				n++;
			}
			else { //  anchor n-1 and anchor n are not overlapped; Need to extend after anchor n-1
				GenomePos qe, te;
				Checkbp((*pairs)[n-1], (*pairs)[n], genome, read, chromIndex, qe, te, strand, opts, mat, K);
				if (strand == 0 and qe == (*pairs)[n].first.pos and te == (*pairs)[n].second.pos) {n++;}
				else if (strand == 1 and qe == (*pairs)[n].first.pos and te == (*pairs)[n].second.pos + K - 1) {n++;}
				else {
					if (strand == 0) {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, (*pairs)[m].second.pos)));}
					else {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, te + 1)));}
					ExtendpairsMatchesLength.push_back(qe - (*pairs)[m].first.pos);
					m = n;
					n++;
				}
			}
		}
		else {
			if (strand == 0) {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, (*pairs)[m].second.pos)));}
			else {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, (*pairs)[n - 1].second.pos)));}
			ExtendpairsMatchesLength.push_back((*pairs)[n-1].first.pos + K - (*pairs)[m].first.pos);
			m = n;
			n++;
		}
	}
	if (n == (*pairs).size()) { // and m != n-1
		if (strand == 0) {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, (*pairs)[m].second.pos)));}
		else {Extendpairs.push_back(GenomePair(GenomeTuple(0, (*pairs)[m].first.pos), GenomeTuple(0, (*pairs)[n-1].second.pos)));}
		ExtendpairsMatchesLength.push_back((*pairs)[n-1].first.pos + K - (*pairs)[m].first.pos);		
	}
}

//
// This function takes GenomePairs ExtendPairs; vector<int> ExtendPairsMatchesLengths;
// Note: only forward strand;
//
void 
TrimOverlappedAnchors(GenomePairs &ExtendPairs, vector<int> &ExtendPairsMatchesLengths) {

	//
	// Get all long anchors
	//
	LongAnchors longanchors(&ExtendPairs, &ExtendPairsMatchesLengths, 0);
	for (int ln = 0; ln < ExtendPairs.size(); ln++) {
		if (ExtendPairsMatchesLengths[ln] >= 50) { 
			longanchors.anchorIndex.push_back(ln);
		}
	}

	//
	// Catersiansort anchors and trim overlapped adjacent two;
	//
	longanchors.Sort();
	for (int ln = 1; ln < longanchors.anchorIndex.size(); ln++) {

		int prev = longanchors.anchorIndex[ln-1];
		int cur = longanchors.anchorIndex[ln];

		//
		// If the adjacent anchors are overlapped less than 30b, trim the prev
		// 
		int overlap_r = 0, overlap_g = 0, overlap = 0;
		bool lr = 0, lg = 0;
		//
		// Check if two anchors are overlapped on read;
		//
		if (ExtendPairs[cur].first.pos < ExtendPairs[prev].first.pos + ExtendPairsMatchesLengths[prev] 
			and 
			ExtendPairs[cur].first.pos >= ExtendPairs[prev].first.pos + ExtendPairsMatchesLengths[prev] - 30) {

			overlap_r = ExtendPairs[prev].first.pos + ExtendPairsMatchesLengths[prev] - ExtendPairs[cur].first.pos;
			assert(overlap_r <= 30); assert(overlap_r > 0);
			lr = 1;
		}				
		//
		// Check if two anchors are overlapped on genome;
		//
		if (ExtendPairs[cur].second.pos < ExtendPairs[prev].second.pos + ExtendPairsMatchesLengths[prev] 
			and 
			ExtendPairs[cur].second.pos >= ExtendPairs[prev].second.pos + ExtendPairsMatchesLengths[prev] - 30) {

			overlap_g = ExtendPairs[prev].second.pos + ExtendPairsMatchesLengths[prev] - ExtendPairs[cur].second.pos;
			assert(overlap_g <= 30); assert(overlap_g > 0);
			lg = 1;
		}
		
		if (overlap_r > 0 or overlap_g > 0) {
			overlap = max(overlap_r, overlap_g);
			ExtendPairsMatchesLengths[prev] -= overlap+1;
			assert(ExtendPairsMatchesLengths[prev]>=0);
			if (lr == 1) assert(ExtendPairs[cur].first.pos >= ExtendPairs[prev].first.pos + ExtendPairsMatchesLengths[prev]);
			if (lg == 1) assert(ExtendPairs[cur].second.pos >= ExtendPairs[prev].second.pos + ExtendPairsMatchesLengths[prev]);
		}
	}	
}

template<typename Tup>
void
LinearExtend_chain (vector<Tup> &chain, vector<Cluster> &ExtendClusters, vector<Cluster*> &RefinedClusters, 
			Options &smallOpts, Genome &genome, Read &read, int start, int &overlap, bool skiprepetitive, int K) {
	//
	// Do linear extension for each anchors and avoid overlapping locations;
	// INPUT: RefinedClusters; OUTPUT: ExtendClusters;
	// NOTICE: ExtendClusters have members: strand, matches, matchesLengths, GenomePos, chromIndex;
	//
	LinearExtend(RefinedClusters, ExtendClusters, chain, smallOpts, genome, read, start, overlap, skiprepetitive, K);
	TrimOverlappedAnchors(ExtendClusters, start);
}

void 
MergeMatchesSameDiag (vector<Cluster> &ExtendClusters, vector<Cluster_SameDiag> &MergedClusters, Options &opts) {
	int original = 0, reduced = 0;
	for (int c = 0; c < ExtendClusters.size(); c++) {
		original += ExtendClusters[c].matches.size();
		MergedClusters.push_back(Cluster_SameDiag(&(ExtendClusters[c])));
		MergedClusters[c].start.push_back(0);
		MergedClusters[c].end.push_back(1);
		MergedClusters[c].coarse = c;
		long prev_diag = GetDiag(ExtendClusters[c].matches[0], ExtendClusters[c].matchesLengths[0], ExtendClusters[c].strand);
		GenomePos prev_qEnd = ExtendClusters[c].matches[0].first.pos + ExtendClusters[c].matchesLengths[0];
		int q = 1;
		while (q < ExtendClusters[c].matches.size()) {
			long cur_diag = GetDiag(ExtendClusters[c].matches[q],ExtendClusters[c].matchesLengths[q],  ExtendClusters[c].strand);
			if (ExtendClusters[c].overlap[q - 1] == 0 and ExtendClusters[c].overlap[q] == 0 
				and prev_diag == cur_diag and prev_qEnd < ExtendClusters[c].matches[q].first.pos
				and GapDifference(ExtendClusters[c].matches[q - 1], ExtendClusters[c].matches[q], ExtendClusters[c].matchesLengths[q - 1]) <= opts.merge_dist) { // 100 for contig
				MergedClusters[c].end.back() = q + 1;
			}
			else {
				MergedClusters[c].start.push_back(q);
				MergedClusters[c].end.push_back(q + 1);
			}
			prev_qEnd = ExtendClusters[c].matches[q].first.pos + ExtendClusters[c].matchesLengths[q];
			prev_diag = cur_diag;
			q++;
		}
		reduced += MergedClusters[c].size();
	}
	// cerr << "original: " << original << " reduced: " << reduced << endl;
}

#endif