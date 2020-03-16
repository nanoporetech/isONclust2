#include "consensus.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>

#include "cluster.h"
#include "hpc.h"
#include "util.h"

std::unique_ptr<spoa::AlignmentEngine> SpoaEngine;

void AddSeqToGraph(const std::string& seq, spoa::Graph* graphPtr,
		   spoa::AlignmentEngine* ae, std::uint32_t weight)
{
    auto graph = std::unique_ptr<spoa::Graph>(graphPtr);
    auto alignment = ae->align(seq, graph);
    graph->add_alignment(alignment, seq, weight);
    graph.release();
}

void AddSeqToGraphWeight(const std::string& seq,
			 const std::vector<std::uint32_t>& w,
			 spoa::Graph* graphPtr, spoa::AlignmentEngine* ae)
{
    auto graph = std::unique_ptr<spoa::Graph>(graphPtr);
    auto alignment = ae->align(seq, graph);
    graph->add_alignment(alignment, seq, w);
    graph.release();
}

bool UpdateClusterConsensus(std::string& consName, Cluster& cl,
			    spoa::Graph* leftGraphPtr,
			    spoa::Graph* rightGraphPtr, std::string& readSeq,
			    double readRawErr, double readHpcErr,
			    int matchStrand, int consMinSize, int consMaxSize,
			    int kmerSize, int windowSize)
{
    auto leftSize = leftGraphPtr->num_sequences();
    auto rightSize = leftSize;
    rightSize = 1;

    auto rs = readSeq;

    if (matchStrand == -1) {
	RevComp(rs);
    }

    if (rightGraphPtr != nullptr) {
	rightSize = rightGraphPtr->num_sequences();
    }

    auto clsRep = cl.at(0).get();

    double hpcErr = (clsRep->HpcSeq->ErrorRate() * double(leftSize) +
		     readHpcErr * double(rightSize)) /
		    double(leftSize + rightSize);

    double rawErr = (clsRep->RawSeq->ErrorRate() * double(leftSize) +
		     readRawErr * double(rightSize)) /
		    double(leftSize + rightSize);

    /*
    double MAX_ERR = pow(10, -(30 / 10));
    if ((leftSize + rightSize) > 50) {
	if (rawErr > MAX_ERR) {
	    rawErr = MAX_ERR;
	}
	if (hpcErr > MAX_ERR) {
	    hpcErr = MAX_ERR;
	}
    }
    */

    if (rightGraphPtr == nullptr) {
	AddSeqToGraph(rs, leftGraphPtr, SpoaEngine.get(), 1);
    }
    else {
	AddSeqToGraph(rs, leftGraphPtr, SpoaEngine.get(), rightSize);
    }

    if (int(leftGraphPtr->num_sequences()) < consMinSize) {
	return false;
    }

    std::string cons = leftGraphPtr->generate_consensus();

    cons.reserve(cons.size());
    auto consLen = cons.length();
    auto& rep = cl[0];

    rep->RawSeq->SetStr(cons);
    rep->RawSeq->SetName(consName);
    rep->RawSeq->SetErrorRate(rawErr);
    rep->RawSeq->SetScore(rawErr * double(cons.length()));
    auto fixedQualHpc = std::to_string(int(-10 * log10(hpcErr)) + 33)[0];
    auto fixedQualRaw = std::to_string(int(-10 * log10(rawErr)) + 33)[0];
    rep->RawSeq->SetQual(std::string(cons.length(), fixedQualRaw));

    auto hpcSeq = std::unique_ptr<Seq>(new Seq);

    if (cons.length() > unsigned(2 * kmerSize) ||
	cons.length() >= unsigned(windowSize)) {
	*hpcSeq = HomopolymerCompressObj(*(rep->RawSeq));
	hpcSeq->SetErrorRate(hpcErr);
	hpcSeq->SetScore(hpcErr * double(hpcSeq->Str().length()));
	rep->HpcSeq->SetQual(std::string(hpcSeq->Str().length(), fixedQualHpc));
	if (hpcSeq->Str().length() < unsigned(2 * kmerSize) ||
	    hpcSeq->Str().length() < unsigned(windowSize)) {
	    hpcSeq->SetScore(-1.0);
	    rep->RawSeq->SetScore(-1.0);
	    rep->RawSeq->SetErrorRate(0.9999);
	    hpcSeq->SetErrorRate(0.9999);
	}
    }

    const auto& kmerSeq = KmerEncodeSeq(hpcSeq->Str(), kmerSize);
    const auto& revKmerSeq = KmerEncodeSeq(RevComp(hpcSeq->Str()), kmerSize);
    hpcSeq->SetErrorRate(hpcErr);
    rep->HpcSeq = std::move(hpcSeq);
    rep->Mins = GetKmerMinimizers(kmerSeq, kmerSize, windowSize);
    rep->RevMins = GetKmerMinimizers(revKmerSeq, kmerSize, windowSize);
    return true;
}

std::unique_ptr<spoa::Graph> ConsPurge(spoa::Graph* graphPtr,
				       spoa::AlignmentEngine* ae, Cluster& cl)
{
    auto& repSeq = cl[0]->RawSeq->Str();
    auto w = graphPtr->num_sequences();
    graphPtr->clear();
    auto newGraph = spoa::createGraph();
    AddSeqToGraph(repSeq, newGraph.get(), ae, w);
    return newGraph;
}

