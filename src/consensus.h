#ifndef CONSENSUS_H_INCLUDED
#define CONSENSUS_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "qualscore.h"
#include "serialize.h"
#include "spoa/spoa.hpp"

bool UpdateClusterConsensus(std::string& consName, Cluster& cl,
			    spoa::Graph* leftGraphPtr,
			    spoa::Graph* rightGraphPtr, ProcSeq* readRep,
			    int matchStrand, int consMinSize, int consMaxSize,
			    int kmerSize, int windowSize);

void AddSeqToGraph(const std::string& seq, spoa::Graph* graphPtr,
		   spoa::AlignmentEngine* ae, std::uint32_t weight);

void AddSeqToGraphWeight(const std::string& seq,
			 const std::vector<std::uint32_t>& w,
			 spoa::Graph* graphPtr, spoa::AlignmentEngine* ae);

std::unique_ptr<spoa::Graph> ConsPurge(spoa::Graph* graphPtr,
				       spoa::AlignmentEngine* ae, Cluster& cl);

#endif
