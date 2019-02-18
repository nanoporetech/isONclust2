#include "consensus.h"
#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>
#include "hpc.h"
#include "util.h"

bool UpdateClusterConsensus(std::string& consName, Cluster& cl, int minSize,
			    int maxSize, int lastSize, int kmerSize,
			    int windowSize)
{
    int period = 0;
    int nrReads = int(cl.size());
    int sampleSize = maxSize;
    if (nrReads <= minSize) {
	sampleSize = nrReads - 1;
	period = 1;
    }
    else if (nrReads >= maxSize) {
	sampleSize = maxSize - 1;
	period = nrReads / sampleSize;
    }
    else {
	period = 1;
    }
    std::vector<std::string> inSeq;

    inSeq.reserve(sampleSize + 1);
    inSeq.push_back(cl[0].RawSeq.Str());

    double hpcErr = cl[0].HpcSeq.ErrorRate();
    double rawErr = cl[0].RawSeq.ErrorRate();

    for (int i = 1; i < int(cl.size()); i++) {
	if (i % period != 0) {
	    continue;
	}
	auto pick = i - (rand() % (1 + period));
	if (pick < 1) {
	    pick = 1;
	}

	auto rs = cl[pick].RawSeq.Str();
	if (cl[pick].MatchStrand == -1) {
	    RevComp(rs);
	}

	inSeq.push_back(rs);
	hpcErr += cl[pick].HpcSeq.ErrorRate();
	rawErr += cl[pick].RawSeq.ErrorRate();
	if (int(inSeq.size()) >= maxSize) {
	    break;
	}
    }

    hpcErr = hpcErr / double(inSeq.size());
    ;
    rawErr = rawErr / double(inSeq.size());

    auto rng = std::default_random_engine{};
    auto sit = std::begin(inSeq);
    sit++;
    if (rand_double() < 0.33) {
	std::shuffle(sit, std::end(inSeq), rng);
    }

    std::int8_t m = 10;
    std::int8_t n = -5;
    std::int8_t g = -8;
    std::int8_t e = -2;
    std::int8_t q = -20;
    std::int8_t c = -1;

    std::uint8_t algorithm = 2;
    std::uint8_t result = 0;

    auto alignment_engine = spoa::createAlignmentEngine(
	static_cast<spoa::AlignmentType>(algorithm), m, n, g, e, q, c);

    auto graph = spoa::createGraph();

    for (const auto& it : inSeq) {
	auto alignment = alignment_engine->align(it, graph);
	graph->add_alignment(alignment, it);
    }

    std::string cons = graph->generate_consensus();
    auto consLen = cons.length();
    auto& rep = cl[0];

    if (consLen < rep.RawSeq.Str().length()) {
	return false;
    }

    rep.RawSeq.SetStr(cons);
    rep.RawSeq.SetName(consName);
    rep.RawSeq.SetErrorRate(rawErr);
    rep.RawSeq.SetScore(rawErr * double(cons.length()));
    auto fixedQualHpc = std::to_string(int(-10 * log10(hpcErr)) + 33)[0];
    auto fixedQualRaw = std::to_string(int(-10 * log10(rawErr)) + 33)[0];
    rep.RawSeq.SetQual(std::string(cons.length(), fixedQualRaw));

    Seq hpcSeq;

    if (cons.length() > unsigned(2 * kmerSize) ||
	cons.length() >= unsigned(windowSize)) {
	hpcSeq = HomopolymerCompressObj(rep.RawSeq);
	hpcSeq.SetErrorRate(hpcErr);
	hpcSeq.SetScore(hpcErr * double(hpcSeq.Str().length()));
	rep.HpcSeq.SetQual(std::string(hpcSeq.Str().length(), fixedQualHpc));
	if (hpcSeq.Str().length() < unsigned(2 * kmerSize) ||
	    hpcSeq.Str().length() < unsigned(windowSize)) {
	    hpcSeq.SetScore(-1.0);
	    rep.RawSeq.SetScore(-1.0);
	    rep.RawSeq.SetErrorRate(0.9999);
	    hpcSeq.SetErrorRate(0.9999);
	}
    }

    const auto& kmerSeq = KmerEncodeSeq(hpcSeq.Str(), kmerSize);
    const auto& revKmerSeq = KmerEncodeSeq(RevComp(hpcSeq.Str()), kmerSize);
    hpcSeq.SetErrorRate(hpcErr);
    rep.Mins = std::move(GetKmerMinimizers(kmerSeq, kmerSize, windowSize));
    rep.RevMins =
	std::move(GetKmerMinimizers(revKmerSeq, kmerSize, windowSize));
    return true;
}

