#ifndef QUALSCORE_H_INCLUDED
#define QUALSCORE_H_INCLUDED

#include <stdlib.h>
#include <deque>
#include <unordered_map>
#include <vector>

#include "seq.h"
#include "serialize.h"

typedef std::vector<double> QualTab;

void FillQualScores(SequencesP& sequences, int kmerSize, int windowSize,
		    const QualTab& qualTab, const QualTab& qualTabNomin);
double CalcQualScore(const Seq& s, int kmerSize, const QualTab& qualTab);
Batch* PrepareSortedBatch(SequencesP& sequences, int batchStart, int batchEnd,
			  int batchSize, int kmerSize, int windowSize,
			  double minQual, const QualTab& qualTab,
			  const QualTab& qualTabNomin);
QualTab InitQualTab();
QualTab InitQualTabNomin();
void SortByQualScores(SequencesP& sequences);
double CalcErrorRate(const std::string& quals, const QualTab& qualTab);

#endif
