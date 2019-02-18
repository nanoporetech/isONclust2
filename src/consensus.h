#ifndef CONSENSUS_H_INCLUDED
#define CONSENSUS_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "qualscore.h"
#include "spoa/spoa.hpp"

bool UpdateClusterConsensus(std::string& consName, Cluster& cl, int minSize,
			    int maxSize, int lastSize, int kmerSize,
			    int windowSize);

#endif
