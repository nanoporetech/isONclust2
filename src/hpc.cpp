#include "hpc.h"
#include <iostream>

Seq* HomopolymerCompress(const std::unique_ptr<Seq>& s)
{
    std::string compSeq;
    std::string compQual;
    auto& seq = s->Str();
    compSeq.reserve(seq.length());
    compQual.reserve(seq.length());
    auto& qual = s->Qual();
    char currBase = seq[0];
    char currQual = qual[0];

    compSeq += currBase;

    for (unsigned i = 1; i < seq.length(); i++) {
	if (seq[i] != currBase) {
	    currBase = seq[i];
	    compSeq += currBase;
	    compQual += currQual;
	    currQual = qual[i];
	}
	else if (currQual < qual[i]) {
	    currQual = qual[i];
	}
    }
    compQual += currQual;

    auto res = new Seq(s->Name(), compSeq, compQual, s->Score());
    return res;
}

Seq HomopolymerCompressObj(const Seq& s)
{
    std::string compSeq;
    std::string compQual;
    auto& seq = s.Str();
    compSeq.reserve(seq.length());
    compQual.reserve(seq.length());
    auto& qual = s.Qual();
    char currBase = seq[0];
    char currQual = qual[0];

    compSeq += currBase;

    for (unsigned i = 1; i < seq.length(); i++) {
	if (seq[i] != currBase) {
	    currBase = seq[i];
	    compSeq += currBase;
	    compQual += currQual;
	    currQual = qual[i];
	}
	else if (currQual < qual[i]) {
	    currQual = qual[i];
	}
    }
    compQual += currQual;

    auto res = Seq(s.Name(), compSeq, compQual, s.Score());
    return res;
}

