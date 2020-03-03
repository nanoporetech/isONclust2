#include "qualscore.h"
#include <math.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "hpc.h"
#include "isONclust2_config.h"
#include "minimizer.h"
#include "tbb/parallel_for.h"
#include "util.h"

using namespace std;

void FillQualScores(SequencesP& sequences, int kmerSize, int windowSize,
		    const QualTab& qualTab, const QualTab& qualTabNomin)
{
    tbb::parallel_for(
	tbb::blocked_range<int>(0, sequences.size()),
	[&](tbb::blocked_range<int> r) {
	    for (int i = r.begin(); i < r.end(); ++i) {
		auto& s = sequences[i];
		if (s->Str().length() > unsigned(2 * kmerSize)) {
		    auto qs = CalcQualScore(*s, kmerSize, qualTab);
		    if (qs <= 0) {
			qs = -1.0;
		    }

		    s->SetScore(qs);
		    s->SetErrorRate(CalcErrorRate(s->Qual(), qualTabNomin));
		}
		else {
		    s->SetScore(-1.0);
		    s->SetErrorRate(1.0);
		}
	    }
	});
}

Batch* PrepareSortedBatch(SequencesP& sequences, int batchStart, int batchEnd,
			  int batchSize, int kmerSize, int windowSize,
			  double minQual, const QualTab& qualTab,
			  const QualTab& qualTabNomin)
{
    int size = 1 + batchEnd - batchStart;
    auto batch = new Batch;

    batch->Cls = Clusters(size);
    tbb::parallel_for(
	tbb::blocked_range<int>(0, size), [&](tbb::blocked_range<int> r) {
	    for (int i = r.begin(); i < r.end(); ++i) {
		auto j = batchStart + i;
		auto& s = sequences[j];
		if (batch->Cls[i] == nullptr) {
		    batch->Cls[i] = make_shared<Cluster>(Cluster());
		}
		if ((-10 * log10(s->ErrorRate())) <= minQual) {
		    batch->Cls[i]->emplace_back(std::make_shared<ProcSeq>(
			ProcSeq{*s, *s, Minimizers{}, Minimizers{}, 0}));
		    continue;
		}
		if (s->Str().length() > unsigned(2 * kmerSize) ||
		    s->Str().length() >= unsigned(windowSize)) {
		    auto hpcSeq = std::unique_ptr<Seq>(HomopolymerCompress(s));
		    if (hpcSeq->Str().length() < unsigned(2 * kmerSize) ||
			hpcSeq->Str().length() < unsigned(windowSize)) {
			s->SetScore(-1.0);
			hpcSeq->SetScore(-1.0);
			batch->Cls[i]->emplace_back(
			    make_shared<ProcSeq>(ProcSeq{
				*s, *hpcSeq, Minimizers{}, Minimizers{}, 0}));
			continue;
		    }
		    const auto& kmerSeq =
			KmerEncodeSeq(hpcSeq->Str(), kmerSize);
		    const auto& revKmerSeq =
			KmerEncodeSeq(RevComp(hpcSeq->Str()), kmerSize);
		    const auto& hpcErr =
			CalcErrorRate(hpcSeq->Qual(), qualTabNomin);
		    hpcSeq->SetErrorRate(hpcErr);
		    auto mins =
			GetKmerMinimizers(kmerSeq, kmerSize, windowSize);
		    auto revMins =
			GetKmerMinimizers(revKmerSeq, kmerSize, windowSize);
		    auto p = sequences[j].release();
		    batch->Cls[i]->emplace_back(make_shared<ProcSeq>(
			ProcSeq{std::move(*p), std::move(*hpcSeq),
				std::move(mins), std::move(revMins), 1}));
		}
		else {
		    s->SetScore(-1.0);
		    batch->Cls[i]->emplace_back(make_shared<ProcSeq>(
			ProcSeq{*s, *s, Minimizers{}, Minimizers{}, 0}));
		}
	    }
	});

    batch->NrCls = int(batch->Cls.size());
    batch->BatchStart = batchStart;
    batch->BatchEnd = batchEnd;
    batch->Depth = -1;
    batch->ConsGs = ConsGraphs(0);
    return batch;
}

double CalcQualScore(const Seq& s, int kmerSize, const QualTab& qualTab)
{
    auto quals = s.Qual();
    std::deque<double> win;
    if ((int)quals.length() <= kmerSize) {
	return -1.0;
    }

    for (int i = 0; i < kmerSize; i++) {
	win.push_back(1.0 - qualTab.at(int(quals[i])));
    }

    double currentNoError = 1.0;
    for (auto& p : win) {
	currentNoError *= p;
    }

    double sumOfExpectations = currentNoError;

    for (unsigned i = kmerSize; i < quals.size(); i++) {
	double p_enter = 1.0 - qualTab.at(int(quals[i]));
	double p_leave = win.front();
	win.pop_front();

	currentNoError *= (p_enter / p_leave);
	sumOfExpectations += currentNoError;
	win.push_back(p_enter);
    }
    return sumOfExpectations;
}

void SortByQualScores(SequencesP& sequences)
{
    std::stable_sort(
	sequences.begin(), sequences.end(),
	[](const std::unique_ptr<Seq>& a, const std::unique_ptr<Seq>& b) {
	    return a->Score() > b->Score();
	});
}

double CalcErrorRate(const std::string& quals, const QualTab& qualTab)
{
    double pSum{};
    for (auto q : quals) {
	pSum += qualTab.at(q);
    }
    return pSum / quals.length();
}

QualTab InitQualTab()
{
    QualTab res(129);
    double min = 0.79433;
    for (int i = 33; i <= 128; i++) {
	auto tmp = pow(10, -((i - 33) / 10.0));
	if (tmp > min) {
	    tmp = min;
	}
	res[i] = tmp;
    }

    return res;
}

QualTab InitQualTabNomin()
{
    QualTab res(129);
    for (int i = 33; i <= 128; i++) {
	auto tmp = pow(10, -((i - 33) / 10.0));
	res[i] = tmp;
    }

    return res;
}

