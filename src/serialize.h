#ifndef SERIALIZE_H_INCLUDED
#define SERIALIZE_H_INCLUDED
#define CEREAL_THREAD_SAFE 1
#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>
#include <string>

#include "args.h"
#include "cluster_data.h"
#include "minimizer.h"
#include "spoa/spoa.hpp"

typedef std::vector<ProcSeq> SortedProcSeqs;
typedef std::unique_ptr<spoa::Graph> ConsGraph;
typedef std::vector<std::unique_ptr<spoa::Graph>> ConsGraphs;

class Batch {
public:
    int BatchNr;
    int BatchStart;
    int BatchEnd;
    unsigned long BatchBases;
    int TotalReads;
    int NrCls{};
    CmdArgs SortArgs;
    std::string LeftLeaf{""};
    std::string RightLeaf{""};
    int Depth{0};
    MinimizerDB MinDB;
    Clusters Cls;
    ConsGraphs ConsGs;
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(BatchNr, BatchStart, BatchEnd, BatchBases, TotalReads, NrCls,
		SortArgs, LeftLeaf, RightLeaf, Depth, MinDB, Cls, ConsGs);
    };

    int NrClusters()
    {
	if (NrCls == int(Cls.size())) {
	    int count = 0;
	    for (const auto& c : Cls) {
		if (c->at(0)->RawSeq.Score() > -1) {
		    count++;
		}
	    }
	    return count;
	}
	else {
	    std::cerr << "Inconsistent batch state: NrCluster " << NrCls
		      << " vs " << Cls.size() << std::endl;
	    exit(1);
	}
	return -1;
    };
    int NrNontrivialClusters()
    {
	if (NrCls == int(Cls.size())) {
	    int count = 0;
	    for (const auto& c : Cls) {
		if ((c->at(0)->RawSeq.Score() > -1) && (c->size() > 2)) {
		    count++;
		}
	    }
	    return count;
	}
	else {
	    std::cerr << "Inconsistent batch state: NrCluster " << NrCls
		      << " vs " << Cls.size() << std::endl;
	    exit(1);
	}
	return -1;
    };
    int NrFilteredReads()
    {
	if (NrCls == int(Cls.size())) {
	    int count = 0;
	    for (const auto& c : Cls) {
		if (c->at(0)->RawSeq.Score() < 0) {
		    count++;
		}
	    }
	    return count;
	}
	else {
	    std::cerr << "Inconsistent batch state: NrCluster " << NrCls
		      << " vs " << Cls.size() << std::endl;
	    exit(1);
	}
	return -1;
    };
    int MinDBSize() { return int(MinDB.size()); };
};

void SaveBatch(const std::unique_ptr<Batch>& b, std::string outf);
typedef std::unique_ptr<Batch> BatchP;
BatchP LoadBatch(std::string inf);
BatchP CreatePseudoBatch(std::unique_ptr<Batch>& inBatch);

#endif
