#ifndef KMER_INDEX_H_INCLUDED
#define KMER_INDEX_H_INCLUDED

#include<string>
#include<vector>
#include<iostream>

typedef std::vector<unsigned> KmerSeq;
KmerSeq KmerEncodeSeq(const std::string seq, unsigned kmerSize);

inline char NumberToBase(unsigned i) {
    switch (i){
        case 0:
            return 'A'; 
            break;
        case 1:
            return 'C';
            break;
        case 2:
            return 'G';
            break;
        case 3:
            return 'T';
            break;
    }
    return '!';
}

inline unsigned BaseToNumber(char c){
    switch(c) {
        case 'A':
            return 0;
            break;
        case 'C':
            return 1;
            break;
        case 'G':
            return 2;
            break;
        case 'T':
            return 3;
            break;
    }
    return -1;
}

inline std::string IndexToKmer(unsigned index, unsigned k) {
    if(k == 1){
        return std::string(1, NumberToBase(index));
    }
    auto prefixIndex = index / 4; 
    auto r = index % 4;
    auto symbol = NumberToBase(r);
    auto prefixPattern = IndexToKmer(prefixIndex, k-1);
    return prefixPattern + symbol;

}

inline unsigned KmerToIndex(const std::string& kmer, std::string::iterator end) {
	if (kmer.begin() == end) {
		return 0;
	}    
    --end;
	auto sym = *end;
	return 4 * KmerToIndex(kmer, end) + BaseToNumber(sym);
}

#endif
