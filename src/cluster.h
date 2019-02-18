#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#include <condition_variable>
#include <mutex>
#include <vector>
#include "cluster_data.h"
#include "minimizer.h"
#include "p_emp_prob.h"
#include "parasail.h"
#include "serialize.h"

#define REP 0

void ClusterSortedReads(BatchP& leftBatch, BatchP& rightBatch, bool quiet);

StrandedCluster getBestCluster(const unsigned rightId, BatchP& leftBatch,
			       BatchP& rightBatch,
			       const MinSharedMap& sharedMinTab);

double getMappedRatio(const Seq& hpcSeq, const Seq& clHpcSeq,
		      const Minimizers& mins, const MinimizerHitVector& hits,
		      const MinSharedMap& sharedMinTab, double minProbNoHits);
int setGapOpen(double e);
parasail_result_t* ParasailAlign(const std::string& ref,
				 const std::string& read, int gapOpen,
				 int gapExtend,
				 const parasail_matrix_t* user_matrix);
double getAlnRatio(const std::string& comp, double e, unsigned slen,
		   unsigned kmerSize);
void SortClustersBySize(Clusters& cls);
unsigned AlnInvoked();
double AlnInvokedPerc(int total);
SortedHits SortMinimizerHits(const MinimizerHits& hits, const Clusters& cls);

#endif

