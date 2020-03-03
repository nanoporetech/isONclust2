#ifndef CLUSTER_DATA_H_INCLUDED
#define CLUSTER_DATA_H_INCLUDED

#include <vector>
#include "minimizer.h"
#include "seq.h"

using namespace std;

extern bool VERBOSE;

typedef struct {
    Seq RawSeq;
    Seq HpcSeq;
    Minimizers Mins;
    Minimizers RevMins;
    int MatchStrand;
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(RawSeq, HpcSeq, Mins, RevMins, MatchStrand);
    };
} ProcSeq;

typedef std::vector<shared_ptr<ProcSeq>> Cluster;
typedef std::vector<shared_ptr<Cluster>> Clusters;

#endif

