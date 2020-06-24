#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include <stdlib.h>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <unordered_map>
#include <vector>
#include "cluster.h"
#include "pbar.h"
#include "seq.h"

#include "serialize.h"

typedef struct {
public:
    std::string Fastq{};
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(Fastq);
    };
} SortedIdx;

void SequencesPToFastq(SequencesP& sequences, const std::string& outFastq,
		       const std::string& indexTab,
		       const std::string& indexCer);
void CreateOutdir(const std::string& outDir);
unsigned WriteFastqRecord(const Seq& s, std::ofstream& out);
void CreateFile(const std::string& outFile, std::ofstream& outfile);
void OpenFile(const std::string& inFile, std::ofstream& infile);
void WriteScores(SequencesP& sequences, const std::string& outFile);
typedef struct {
    unsigned Cls;
    int Strand;
} IdInfo;
typedef std::unordered_map<std::string, std::unique_ptr<IdInfo>> IdMap;

typedef struct {
    std::string Id;
    std::string Header;
    std::string Seq;
    std::string Plus;
    std::string Qual;
} FqRec;

typedef std::unique_ptr<FqRec> FqRecP;
typedef std::unordered_map<unsigned, std::vector<FqRecP>> SeqCache;
void WriteClusters(BatchP& cls, const std::string& outFile, SortedIdx* idx,
		   IdMap& idToCls);
std::unique_ptr<SortedIdx> LoadIndex(std::string inf);

#endif
