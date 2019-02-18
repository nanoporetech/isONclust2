#ifndef CLUSTER_DATA_H_INCLUDED
#define CLUSTER_DATA_H_INCLUDED

#include <vector>
#include "minimizer.h"
#include "seq.h"

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

typedef std::vector<ProcSeq> Cluster;
typedef std::vector<Cluster> Clusters;

#endif

