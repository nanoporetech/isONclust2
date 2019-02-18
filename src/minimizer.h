#ifndef MINIMIZER_H_INCLUDED
#define MINIMIZER_H_INCLUDED

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "kmer_index.h"
#include "tbb/concurrent_vector.h"

#define MIN_DB_RESERVE 100000
#define MIN_SET_RESERVE 1000

class Seq;

typedef struct {
    unsigned Min;
    unsigned Pos;
    unsigned Index;
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(Min, Pos, Index);
    }
} Minimizer;

typedef std::vector<Minimizer> Minimizers;
bool operator==(const Minimizer& a, const Minimizer& b);

Minimizers GetKmerMinimizers(const KmerSeq& kmerSeq, int kmerSize,
			     int windowSize);

struct UnsignedHash {
    inline std::size_t operator()(const unsigned& u) const
    {
	return std::size_t(u);
    }
};

typedef std::pair<int, int> StrandedCluster;
struct StrandedClsHash {
    inline std::size_t operator()(const StrandedCluster& u) const
    {
	return std::size_t(int(u.first * u.second));
    }
};

typedef std::vector<unsigned> RepSet;
typedef std::unordered_map<unsigned, RepSet, UnsignedHash> MinimizerDB;
void AddMinimizers(const Minimizers& mins, unsigned cls, MinimizerDB& db);
typedef struct {
    unsigned Pos;
    unsigned Index;
} MinimizerHit;

typedef struct {
    unsigned Cls;
    MinimizerHit Hits;
} MinHitPair;

typedef std::vector<MinimizerHit> MinimizerHitVector;
typedef tbb::concurrent_vector<MinHitPair> RawMinimizerHits;
typedef std::unordered_map<StrandedCluster, MinimizerHitVector, StrandedClsHash>
    MinimizerHits;

MinimizerHits GetMinimizerHits(const Minimizers& mins,
			       const Minimizers& revMins,
			       const MinimizerDB& db);
void ConsolidateMinimizerHits(const RawMinimizerHits& hits, MinimizerHits& res,
			      int strand);
typedef std::vector<std::tuple<unsigned, unsigned, unsigned, std::string, int>>
    SortedHits;
void UpdateMinDB(int best, const Minimizers& oldMins, const Minimizers& newMins,
		 MinimizerDB& db);

#endif
