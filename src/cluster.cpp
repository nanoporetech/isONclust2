#include "cluster.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <numeric>

#include "args.h"
#include "cluster_data.h"
#include "consensus.h"
#include "kmer_index.h"
#include "minimizer.h"
#include "p_emp_prob.h"
#include "parasail.h"
#include "pbar.h"
#include "util.h"

using namespace std;
unsigned ALN_INVOKED{0};
unsigned CONS_INVOKED{0};
UnsignedHash uh;
extern std::unique_ptr<spoa::AlignmentEngine> SpoaEngine;

int ProcSeqWeight(ProcSeq& s) { return int(s.RawSeq->MeanQual()); }

void printSortedSizes(Clusters& cls)
{
    std::vector<unsigned> sizes;
    sizes.reserve(cls.size());
    for (auto c : cls) {
	auto s = c->size();
	if (s > 1) {
	    sizes.push_back(s);
	}
    }
    std::sort(sizes.rbegin(), sizes.rend());
    for (unsigned i = 0; i < sizes.size(); i++) {
	std::cerr << sizes[i];
	if (i != sizes.size() - 1) {
	    std::cerr << ",";
	}
    }
}

unsigned countNtClusters(Clusters& cls)
{
    unsigned count = 0;
    for (auto& c : cls) {
	if (c->size() > 1) {
	    count++;
	}
    }
    return count;
}

void dumpMinimizers(const Minimizers& mins, const std::string& readId,
		    unsigned kmerSize)
{
    for (auto& m : mins) {
	std::cerr << readId << "\t" << IndexToKmer(m.Min, kmerSize) << "\t"
		  << m.Pos << "\t" << m.Index << std::endl;
    }
}

void ClusterSortedReads(BatchP& leftBatch, BatchP& rightBatch, bool quiet,
			bool seqPurge)
{
    if (leftBatch->SortArgs != rightBatch->SortArgs) {
	std::cerr << "The left and right batches have been sorted with "
		     "different parameters! "
		  << std::endl;
	std::cerr << "Refusing to carry on with clustering as results would "
		     "not make sense! "
		  << std::endl;
	exit(1);
    }
    auto args = leftBatch->SortArgs;

    if (rightBatch->Depth > 0 &&
	rightBatch->BatchStart != (leftBatch->BatchEnd + 1)) {
	cerr << "Trying to merge non-consecutive batches! Giving up!" << endl;
	exit(1);
    }

    if (leftBatch->Depth > 0 && rightBatch->Depth > leftBatch->Depth) {
	cerr << "The left input batch must have higher depth!" << endl;
	exit(1);
    }

    if (leftBatch->MinDB.size() == 0) {
	leftBatch->MinDB = MinimizerDB(MIN_DB_RESERVE, uh);
    }

    rightBatch->MinDB = MinimizerDB(0, uh);

    auto& cls = leftBatch->Cls;
    leftBatch->ConsGs.reserve(cls.size());
    auto& reads = rightBatch->Cls;
    auto& minDB = leftBatch->MinDB;
    auto& consMaxSize = leftBatch->SortArgs.ConsMaxSize;

    auto sharedMinTab = InitMinSharedMap(args.KmerSize, args.WindowSize);

    if (args.Debug) {
	std::cerr
	    << "Iteration\tNrClusters\tMinDbSize\tCurrReadId\tClusterSizes"
	    << std::endl;
    }

    for (unsigned i = 0; i < reads.size(); i++) {
	if (reads[i]->size() == 0) {
	    continue;
	}
	auto& read = reads[i]->at(REP);
	auto& seq = read->RawSeq;
	const auto& hpcSeq = read->HpcSeq;
	const auto hpcErr = hpcSeq->ErrorRate();
	if (args.Debug) {
	    std::cerr << i << "\t";
	    std::cerr << countNtClusters(cls) << "\t";
	    std::cerr << minDB.size() << "\t";
	    std::cerr << seq->Name() << "\t";
	    printSortedSizes(cls);
	    std::cerr << std::endl;
	}

	int best = -2;

	if (VERBOSE && !quiet) {
	    Pbar((float)(i + 1) / float(reads.size()));
	}
	if (seq->Score() < 0) {
	    continue;
	}
	if (seq->Str().length() < unsigned(2 * args.KmerSize)) {
	    seq->SetScore(-1.0);
	    continue;
	}
	if (hpcSeq->Str().length() < unsigned(2 * args.KmerSize)) {
	    seq->SetScore(-1.0);
	    continue;
	}

	if ((-10 * log10(seq->ErrorRate())) <= args.MinQual) {
	    seq->SetScore(-1.0);
	    continue;
	}

	const auto& mins = read->Mins;
	const auto& revMins = read->RevMins;
	StrandedCluster stMatch;
	if (best != -1) {
	    stMatch = getBestCluster(i, leftBatch, rightBatch, sharedMinTab);
	    best = stMatch.first;
	}

	auto startIt = reads[i]->begin();
	auto readTmp = *startIt;

	auto readSeq = std::string(readTmp->RawSeq->Str());
	auto readRawErr = readTmp->RawSeq->ErrorRate();
	auto readHpcErr = readTmp->HpcSeq->ErrorRate();

	if (best == -1) {
	    auto newId = unsigned(cls.size());
	    auto nrReads = reads[i]->size();
	    AddMinimizers(mins, newId, minDB);
	    if (nrReads == 1) {
		auto nrep = new ProcSeq;
		auto rep = reads[i]->at(0);
		nrep->RawSeq = std::unique_ptr<Seq>(new Seq);
		*(nrep->RawSeq) = *(rep->RawSeq);
		nrep->HpcSeq = std::unique_ptr<Seq>(new Seq);
		*(nrep->HpcSeq) = *(rep->HpcSeq);
		nrep->Mins = rep->Mins;
		nrep->RevMins = rep->RevMins;
		nrep->MatchStrand = rep->MatchStrand;
		nrep->Id = rep->Id;
		auto repName = "rep_" + std::to_string(leftBatch->BatchNr) +
			       "_" + std::to_string(newId);
		nrep->RawSeq->SetName(repName);
		nrep->HpcSeq->SetName(repName);
		reads[i]->insert(reads[i]->begin(), 1,
				 std::unique_ptr<ProcSeq>(nrep));
	    }
	    leftBatch->ConsGs.emplace_back(spoa::createGraph());

	    AddSeqToGraph(reads[i]->at(0)->RawSeq->Str(),
			  leftBatch->ConsGs[newId].get(), SpoaEngine.get(), 1);

	    cls.emplace_back(reads[i]);
	    if (nrReads == 1 && cls[newId]->size() != 2) {
		std::cerr << "Inconsistent initial cluster size "
			  << cls[newId]->size()
			  << " for "
			     "read: "
			  << cls[newId]->at(0)->RawSeq->Name() << std::endl;
		std::cerr << "Aborting clustering!" << std::endl;
		exit(1);
	    }

	    leftBatch->NrCls++;
	    if (rightBatch->ConsGs.size() > 0 &&
		rightBatch->ConsGs[i] != nullptr) {
		rightBatch->ConsGs[i] = nullptr;
	    }
	}
	else {
	    cls[best]->reserve(cls[best]->size() + reads[i]->size());
	    // Free unecessary minimizers:
	    auto tmp = Minimizers(0);
	    for (unsigned j = 0; j < reads[i]->size(); j++) {
		auto& s = reads[i]->at(j);
		if (s == nullptr) {
		    std::cerr << "Null pointer at position " << j
			      << " in read array " << i
			      << " of size: " << reads[i]->size() << std::endl;
		    exit(1);
		}
		if (stMatch.second == -1) {
		    switch (s->MatchStrand) {
			case 1:
			    s->MatchStrand = -1;
			    break;
			case -1:
			    s->MatchStrand = 1;
			    break;
			default:
			    throw("Invalid match strand!");
		    }
		}
		s->Mins = tmp;
		s->RevMins = tmp;
		if (!seqPurge) {
		    s->RawSeq = nullptr;
		    s->HpcSeq = nullptr;
		}
	    }

	    if (reads[i]->size() > 1) {
		startIt++;
	    }

	    auto oldSize = reads[i]->size();
	    std::move(startIt, std::end(*(reads[i])),
		      std::back_inserter(*(cls[best])));

	    if (consMaxSize <= 0) {
		continue;
	    }

	    if ((leftBatch->Depth == -1) &&
		(leftBatch->SortArgs.ConsPeriod > 0) &&
		(int(cls[best]->size()) > leftBatch->SortArgs.ConsPeriod)) {
		continue;
	    }

	    auto consGraphLeft = leftBatch->ConsGs[best].get();
	    spoa::Graph* consGraphRight = nullptr;

	    if (rightBatch->ConsGs.size() > 0) {
		consGraphRight = rightBatch->ConsGs[i].get();
	    }

	    std::string consName = "cons_" +
				   std::to_string(leftBatch->BatchNr) + "_" +
				   std::to_string(i);

	    auto oldMins = cls.at(best)->at(REP)->Mins;
	    auto consMinSize = leftBatch->SortArgs.ConsMinSize;
	    if (leftBatch->Depth != -1) {
		consMinSize = 2;  // FIXME
	    }
	    auto ok = UpdateClusterConsensus(
		consName, *(cls[best]), consGraphLeft, consGraphRight, readSeq,
		readRawErr, readHpcErr, stMatch.second, consMinSize,
		consMaxSize, args.KmerSize, args.WindowSize);

	    if (ok) {
		CONS_INVOKED++;
		UpdateMinDB(best, oldMins, cls[best]->at(REP)->Mins, minDB);
	    }

	    if (ok && (int(consGraphLeft->num_sequences()) > consMaxSize)) {
		auto newGraph =
		    ConsPurge(consGraphLeft, SpoaEngine.get(), *(cls[best]));
		leftBatch->ConsGs[best].swap(newGraph);
	    }

	    if (rightBatch->ConsGs.size() > 0 &&
		rightBatch->ConsGs[i] != nullptr) {
		rightBatch->ConsGs[i] = nullptr;
	    }
	}
    }
    if (VERBOSE) {
	std::cout << std::endl;
    }
    leftBatch->Depth++;
    leftBatch->BatchEnd = rightBatch->BatchEnd;
    leftBatch->BatchBases = leftBatch->BatchBases + rightBatch->BatchBases;
}

double getMappedRatio(const Seq& hpcSeq, const Seq& clHpcSeq,
		      const Minimizers& mins, const MinimizerHitVector& hits,
		      const MinSharedMap& sharedMinTab, double minProbNoHits)
{
    auto clHpcErr = clHpcSeq.ErrorRate();
    double pError =
	1.0 - GetPMinShared(clHpcErr, hpcSeq.ErrorRate(), sharedMinTab);
    double totalMapped{0};

    if (pow(pError, hits[0].Index) >= minProbNoHits) {
	totalMapped += double(hits[0].Pos);
    }

    for (unsigned i = 0; i < hits.size() - 1; i++) {
	auto& h1 = hits[i];
	auto& h2 = hits[i + 1];
	double noMatchProb = pow(pError, double(h2.Index - (h1.Index + 1)));
	if (noMatchProb >= minProbNoHits) {
	    totalMapped += double(h2.Pos - h1.Pos);
	}
    }

    auto& h = hits[hits.size() - 1];
    if (pow(pError, double(mins.size() - (h.Index + 1))) >= minProbNoHits) {
	totalMapped += hpcSeq.Str().length() - h.Pos;
    }

    auto mr = totalMapped / (double)hpcSeq.Str().length();
    return mr;
}

StrandedCluster getBestClusterMapping(const ProcSeq& read,
				      const BatchP& leftBatch,
				      const MinimizerHits& hits,
				      const SortedHits& order,
				      const MinSharedMap& sharedMinTab)
{
    auto& hpcSeq = read.HpcSeq;
    auto hpcErr = read.HpcSeq->ErrorRate();
    auto& mins = read.Mins;
    auto& revMins = read.RevMins;
    auto& cls = leftBatch->Cls;
    auto minShared = leftBatch->SortArgs.MinShared;
    auto minFrac = leftBatch->SortArgs.MinFraction;
    auto minProbNoHits = leftBatch->SortArgs.MinProbNoHits;
    auto mappedTh = leftBatch->SortArgs.MappedThreshold;
    auto NEG = std::make_pair(int(-1), int(0));

    if (order.size() == 0) {
	return NEG;
    }

    auto nrTopHits = std::get<0>(order[0]);
    if (nrTopHits < (unsigned)minShared) {
	return NEG;
    }

    for (auto& c : order) {
	auto& nmHits = std::get<0>(c);
	auto& clId = std::get<2>(c);
	const Minimizers& m = mins;
	auto strand = int(std::get<3>(c));
	auto scl = std::make_pair(int(clId), int(strand));
	if (int(nmHits) < int((double)nrTopHits * minFrac)) {
	    return NEG;
	}
	float mr = 0.0;
	if (strand == 1) {
	    mr = getMappedRatio(*hpcSeq, *(cls.at(clId)->at(REP)->HpcSeq), mins,
				hits.at(scl), sharedMinTab, minProbNoHits);
	}
	else {
	    mr = getMappedRatio(*hpcSeq, *(cls.at(clId)->at(REP)->HpcSeq),
				revMins, hits.at(scl), sharedMinTab,
				minProbNoHits);
	}
	if (mr >= mappedTh) {
	    return scl;
	}
    }

    return NEG;
}

parasail_result_t* ParasailAlign(const std::string& read,
				 const std::string& ref, int gapOpen,
				 int gapExtend,
				 const parasail_matrix_t* user_matrix)
{
    auto alnRes = parasail_sg_trace_scan_16(read.c_str(), read.length(),
					    ref.c_str(), ref.length(), gapOpen,
					    gapExtend, user_matrix);
    if (parasail_result_is_saturated(alnRes)) {
	alnRes = parasail_sg_trace_scan_32(read.c_str(), read.length(),
					   ref.c_str(), ref.length(), gapOpen,
					   gapExtend, user_matrix);
    }

    return alnRes;
}

int setGapOpen(double e)
{
    if (e <= 0.01) {
	return 5;
    }
    else if ((0.01 < e) && (e <= 0.04)) {
	return 4;
    }
    else if ((0.04 < e) && (e <= 0.1)) {
	return 3;
    }
    else if (e > 0.1) {
	return 2;
    }
    throw;
}

double getAlnRatio(const std::string& comp, double e, unsigned slen,
		   unsigned kmerSize)
{
    double aligned{};
    auto limit = floor((1.0 - e) * kmerSize);
    auto i = comp.begin();
    auto j = std::next(i, kmerSize);

    while (j != comp.end()) {
	auto nm = std::count(i, j, '|');
	if (nm >= limit) {
	    aligned++;
	}
	i = std::next(i);
	j = std::next(j);
    }
    return aligned / slen;
}

StrandedCluster getBestClusterAln(const ProcSeq& read,
				  const SortedHits& hitOrder,
				  const BatchP& leftBatch)
{
    auto& clsLeft = leftBatch->Cls;
    auto alignedTh = leftBatch->SortArgs.AlignedThreshold;
    auto kmerSize = leftBatch->SortArgs.KmerSize;
    auto NEG = std::make_pair(int(-1), int(0));
    if (hitOrder.size() == 0) {
	return NEG;
    }
    auto topHit = std::get<0>(hitOrder[0]);
    auto& readSeq = read.RawSeq->Str();

    int match = 2;
    int mismatch = -2;
    int gapExtend = 1;

    auto user_matrix = parasail_matrix_create("ACGT", match, mismatch);

    for (auto& c : hitOrder) {
	if (std::get<0>(c) < topHit) {
	    break;
	}
	auto strand = std::get<3>(c);
	auto clId = unsigned(std::get<2>(c));
	auto& rep = clsLeft[clId]->at(REP)->RawSeq;
	auto repSeq = std::string(rep->Str());
	if (strand == -1) {
	    repSeq = RevComp(repSeq);
	}
	auto e2 = rep->ErrorRate();
	auto& clsQual = rep->Qual();

	auto e1 = read.RawSeq->ErrorRate();
	int gapOpen = setGapOpen(e1 + e2);

	auto alnRes =
	    ParasailAlign(readSeq, repSeq, gapOpen, gapExtend, user_matrix);
	auto tr = parasail_result_get_traceback(
	    alnRes, readSeq.c_str(), readSeq.length(), repSeq.c_str(),
	    repSeq.length(), user_matrix, '|', ' ', ' ');
	parasail_result_free(alnRes);

	std::string Comp(tr->comp);
	auto alnRatio = getAlnRatio(Comp, e1 + e2, readSeq.length(), kmerSize);
	parasail_traceback_free(tr);
	if (alnRatio >= alignedTh) {
	    return std::make_pair(std::get<2>(c), strand);
	}
    }

    parasail_matrix_free(user_matrix);
    return NEG;
}

void dumpSortedHits(const SortedHits& order, const std::string& readId,
		    const Clusters& cls)
{
    unsigned i = 0;
    std::cerr << readId << std::endl;
    for (auto& h : order) {
	std::cerr << "\t" << i << "\t"
		  << cls[std::get<2>(h)]->at(REP)->RawSeq->Name() << "\t"
		  << std::get<0>(h) << "\t" << std::get<1>(h) << "\t"
		  << std::get<2>(h) << std::endl;
	i++;
    }
}

StrandedCluster getBestCluster(const unsigned rightId, BatchP& leftBatch,
			       BatchP& rightBatch,
			       const MinSharedMap& sharedMinTab)
{
    auto mode = leftBatch->SortArgs.Mode;
    auto minShared = leftBatch->SortArgs.MinShared;
    auto minProbNoHits = leftBatch->SortArgs.MinProbNoHits;
    auto& read = rightBatch->Cls[rightId]->at(REP);
    auto&& hits = GetMinimizerHits(read->Mins, read->RevMins, leftBatch->MinDB);
    auto&& hitOrder = SortMinimizerHits(hits, leftBatch->Cls);
    auto NEG = std::make_pair(-1, 0);
    if (hitOrder.size() == 0) {
	return NEG;
    }

    if ((mode == Sahlin) || (mode == Fast)) {
	auto mapCluster = getBestClusterMapping(*read, leftBatch, hits,
						hitOrder, sharedMinTab);
	if (mapCluster.first > -1) {
	    return mapCluster;
	}
    }

    auto shared = std::get<0>(hitOrder[0]);
    if (shared < unsigned(minShared)) {
	return NEG;
    }

    if (mode == Fast) {
	return NEG;
    }

    if ((mode == Furious) || (mode == Sahlin)) {
	ALN_INVOKED++;
	auto alnCluster = getBestClusterAln(*read, hitOrder, leftBatch);
	return alnCluster;
    }
    return NEG;
}

void SortClustersBySize(Clusters& cls)
{
    std::stable_sort(
	cls.begin(), cls.end(),
	[](const shared_ptr<Cluster> a, const shared_ptr<Cluster> b) {
	    if (a->size() == b->size()) {
		return a->at(REP)->RawSeq->Score() >
		       b->at(REP)->RawSeq->Score();
	    }
	    return a->size() > b->size();
	});
}

unsigned AlnInvoked() { return ALN_INVOKED; }
unsigned ConsInvoked() { return CONS_INVOKED; }
double AlnInvokedPerc(int total)
{
    if (AlnInvoked() == 0) {
	return 0.0;
    }
    return double(AlnInvoked()) / double(total) * 100;
}

double ConsInvokedPerc(int total)
{
    if (ConsInvoked() == 0) {
	return 0.0;
    }
    return double(ConsInvoked()) / double(total) * 100;
}

unsigned sumPositions(const MinimizerHitVector& hits)
{
    unsigned sum = 0;
    for (auto& m : hits) {
	sum += m.Pos;
    }
    return sum;
}

void ConsolidateMinimizerHits(const RawMinimizerHits& hits, MinimizerHits& res,
			      int strand)
{
    for (auto& mp : hits) {
	res[std::make_pair(mp.Cls, strand)].emplace_back(std::move(mp.Hits));
    }
}

bool hitsSortUtil(const std::tuple<unsigned, unsigned, unsigned, int>& a,
		  const std::tuple<unsigned, unsigned, unsigned, int>& b)
{
    auto sa = std::get<0>(a);
    auto sb = std::get<0>(b);
    if (sa == sb) {
	auto pa = std::get<1>(a);
	auto pb = std::get<1>(b);
	return pa > pb;
    }
    return sa > sb;
}

SortedHits SortMinimizerHits(const MinimizerHits& hits, const Clusters& cls)
{
    SortedHits sorted;
    sorted.reserve(hits.size());
    for (auto& hit : hits) {
	sorted.push_back(
	    std::make_tuple(hit.second.size(), sumPositions(hit.second),
			    unsigned(hit.first.first), int(hit.first.second)));
    }
    std::stable_sort(sorted.begin(), sorted.end(), hitsSortUtil);
    return sorted;
}

