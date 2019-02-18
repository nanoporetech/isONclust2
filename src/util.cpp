#include "util.h"
#include <algorithm>
#include <cstdlib>
#include <string>

double round(double number, int precision)
{
    int decimals = std::pow(10, precision);
    return (std::round(number * decimals)) / decimals;
}

// From: https://stackoverflow.com/a/42549427/5504770
std::string RevComp(const std::string& oseq)
{
    auto DNAseq = std::string(oseq);
    std::reverse(DNAseq.begin(), DNAseq.end());
    for (std::size_t i = 0; i < DNAseq.length(); ++i) {
	switch (DNAseq[i]) {
	    case 'A':
		DNAseq[i] = 'T';
		break;
	    case 'C':
		DNAseq[i] = 'G';
		break;
	    case 'G':
		DNAseq[i] = 'C';
		break;
	    case 'T':
		DNAseq[i] = 'A';
		break;
	    default:
		throw("Invalid base encountered: " + std::string(1, DNAseq[i]) +
		      "\n");
		break;
	}
    }
    return DNAseq;
}

double rand_double()
{
    double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    return r;
}

