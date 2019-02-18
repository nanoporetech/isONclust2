#include "kmer_index.h"
#include <string>
#include <vector>

KmerSeq KmerEncodeSeq(const std::string seq, unsigned kmerSize)
{
    std::vector<unsigned> res;
    if (seq.length() < kmerSize) {
	return res;
    }
    res.reserve(seq.length() - kmerSize);
    for (unsigned i = 0; i < seq.length() - kmerSize; i++) {
	auto kmer = seq.substr(i, kmerSize);
	res.emplace_back(KmerToIndex(kmer, kmer.end()));
    }
    return res;
}
