#ifndef CLUSTER_DATA_H_INCLUDED
#define CLUSTER_DATA_H_INCLUDED

#include <vector>
#include "minimizer.h"
#include "seq.h"

using namespace std;

extern bool VERBOSE;

typedef std::unique_ptr<Seq> SeqUptr;

typedef struct {
    SeqUptr RawSeq;
    SeqUptr HpcSeq;
    Minimizers Mins;
    Minimizers RevMins;
    int MatchStrand;
    std::string Id;
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(RawSeq, HpcSeq, Mins, RevMins, MatchStrand, Id);
    };
} ProcSeq;

typedef std::vector<shared_ptr<ProcSeq>> Cluster;
typedef std::vector<shared_ptr<Cluster>> Clusters;

#endif

