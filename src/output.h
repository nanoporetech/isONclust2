#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include "cluster.h"
#include "pbar.h"
#include "seq.h"

void SequencesPToFastq(SequencesP& sequences, const std::string& outFastq);
void CreateOutdir(const std::string& outDir);
void WriteFastqRecord(const Seq& s, std::ofstream& out);
void CreateFile(const std::string& outFile, std::ofstream& outfile);
void WriteScores(SequencesP& sequences, const std::string& outFile);
void WriteClusters(Clusters& cls, const std::string& outFile);

#endif
