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
    std::unordered_map<std::string, unsigned long long int> IdxMap;
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(Fastq, IdxMap);
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
void WriteClusters(Clusters& cls, const std::string& outFile, SortedIdx* idx);
std::unique_ptr<SortedIdx> LoadIndex(std::string inf);

#endif
