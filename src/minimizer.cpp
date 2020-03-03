#include "minimizer.h"
#include <algorithm>
#include <deque>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include "kmer_index.h"
#include "seq.h"
#include "tbb/parallel_for.h"

StrandedClsHash sch;

bool operator==(const Minimizer& a, const Minimizer& b)
{
    if (a.Pos != b.Pos) {
	return false;
    }
    if (a.Min != b.Min) {
	return false;
    }
    if (a.Index != b.Index) {
	return false;
    }
    return true;
}

void AddMinimizers(const Minimizers& mins, unsigned cls, MinimizerDB& db)
{
    for (const auto& m : mins) {
	auto v = db.find(m.Min);
	if (v == db.end()) {
	    db[m.Min] = std::move(RepSet{cls});
	}
	else if (v->second.size() == 0 || cls > v->second.back()) {
	    v->second.emplace_back(cls);
	}
    }
}

MinimizerHits GetMinimizerHits(const Minimizers& mins,
			       const Minimizers& revMins, const MinimizerDB& db)
{
    RawMinimizerHits hits;
    MinimizerHits res(20 * (mins.size() + revMins.size()), sch);

    hits.reserve(20 * mins.size());

    for (auto& m : mins) {
	auto it = db.find(m.Min);
	if (it != db.end()) {
	    for (auto& cls : it->second) {
		hits.emplace_back(
		    MinHitPair{cls, MinimizerHit{m.Pos, m.Index}});
	    }
	}
    }
    ConsolidateMinimizerHits(hits, res, 1);

    hits.clear();
    for (auto& rm : revMins) {
	auto it = db.find(rm.Min);
	if (it != db.end()) {
	    for (auto& cls : it->second) {
		hits.emplace_back(
		    MinHitPair{cls, MinimizerHit{rm.Pos, rm.Index}});
	    }
	}
    }
    ConsolidateMinimizerHits(hits, res, -1);

    return res;
}

Minimizers GetKmerMinimizers(const KmerSeq& kmerSeq, int kmerSize,
			     int windowSize)
{
    Minimizers minimizers;

    int initW = windowSize - kmerSize;
    minimizers.reserve(kmerSeq.size() - initW);

    unsigned index{0};

    std::deque<unsigned> windowKmers;
    for (int i = 0; i <= initW; i++) {
	windowKmers.push_back(kmerSeq[i]);
    }

    auto currIt = std::min_element(windowKmers.begin(), windowKmers.end());
    auto currMin = *currIt;
    unsigned pos = std::distance(windowKmers.begin(), currIt);
    auto min = Minimizer{currMin, pos, index};
    minimizers.push_back(min);
    index++;

    for (int i = initW + 1; i < (int)kmerSeq.size(); i++) {
	const auto& newKmer = kmerSeq[i];
	const auto& oldKmer = windowKmers.front();
	windowKmers.pop_front();
	windowKmers.push_back(newKmer);

	if (currMin == oldKmer) {
	    currIt = std::min_element(windowKmers.begin(), windowKmers.end());
	    pos = std::distance(windowKmers.begin(), currIt) + i - initW;
	    currMin = *currIt;
	    min = Minimizer{currMin, pos, index};
	    minimizers.push_back(min);
	    index++;
	}
	else if (newKmer < currMin) {
	    currMin = newKmer;
	    min = Minimizer{newKmer, (unsigned)i, index};
	    minimizers.push_back(min);
	    index++;
	}
    }

    return minimizers;
}
void UpdateMinDB(int best, const Minimizers& oldMins, const Minimizers& newMins,
		 MinimizerDB& db)
{
    std::set<unsigned> oldSet;
    std::set<unsigned> newSet;
    std::set<unsigned> toIns;
    std::set<unsigned> toDel;

    for (auto& m : oldMins) {
	oldSet.insert(m.Min);
    }
    for (auto& m : newMins) {
	newSet.insert(m.Min);
    }

    std::set_difference(oldSet.begin(), oldSet.end(), newSet.begin(),
			newSet.end(), std::inserter(toDel, toDel.begin()));
    std::set_difference(newSet.begin(), newSet.end(), oldSet.begin(),
			oldSet.end(), std::inserter(toIns, toIns.begin()));

    for (auto m : toDel) {
	auto& mins = db[m];
	std::set<unsigned> tmp(mins.begin(), mins.end());
	tmp.erase(best);
	mins.clear();
	mins.insert(mins.begin(), tmp.begin(), tmp.end());
	// if (tmp.size() == 0) {
	//    db.erase(m);
	//}
    }

    for (auto m : toIns) {
	auto& tv = db[m];
	tv.push_back(best);
	std::sort(tv.begin(), tv.end());
    }
}
