#ifndef INDEXED_SEED_H_
#define INDEXED_SEED_H_
#include <thread>
#include "seqan/seeds.h"

struct IndexedSeed_;
typedef seqan::Tag<IndexedSeed_> IndexedSeed;

template<typename TConfig>
class seqan::Seed<IndexedSeed, TConfig> : public seqan::Seed<Simple, TConfig>{
public:
        int index;
        Seed(int i, int j, int k, int l) : seqan::Seed<Simple,TConfig>(i,j,k,l){index=-1;}
        Seed(int i, int j, int k) : seqan::Seed<Simple,TConfig>(i,j,k){index=-1;}
        Seed() : seqan::Seed<Simple,TConfig>(){index=-1;}
        Seed(int i, int j, int k, int l, int idx) : seqan::Seed<Simple,TConfig>(i,j,k,l), index(idx) {}
};

#endif
