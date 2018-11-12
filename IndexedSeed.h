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


typedef seqan::Seed<IndexedSeed> IndSeed;
typedef seqan::SeedSet<IndSeed> IndSeedSet;
typedef seqan::Iterator<IndSeedSet>::Type TIterator;


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
	beginPositionH(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionH;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
beginPositionV(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionV;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
endPositionH(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionH;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
endPositionV(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionV;
}


// Debugging code
template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return stream << "Seed<IndexedSeed, TConfig>(" << beginPositionH(seed)
                  << ", " << beginPositionV(seed) << ", "
                  << endPositionH(seed) << ", "
                  << endPositionV(seed) << ", lower diag = "
                  << lowerDiagonal(seed) << ", upper diag = "
                  << upperDiagonal(seed) << ")";
}

#endif
