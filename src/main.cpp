#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include "args.h"
#include "bioparser/bioparser.hpp"
#include "cluster.h"
#include "minimizer.h"
#include "output.h"
#include "p_emp_prob.h"
#include "qualscore.h"
#include "serialize.h"
#include "util.h"

#define MAIN
int mainSort(int argc, char* argv[]);
int mainCluster(int argc, char* argv[]);
int mainDump(int argc, char* argv[]);
int mainInfo(int argc, char* argv[]);
void printBatchInfo(BatchP& b);
void dumpBatchInfo(BatchP& b, std::string outfile);
void dumpClusters(BatchP& b, std::string outdir);

UnsignedHash uh;

using namespace std;
int main(int argc, char* argv[])
{
    if (argc < 2) {
	std::cerr << "No subcommand specified!" << endl;
	std::cerr
	    << "Valid subcommands are: sort, cluster, dump, info, version, help"
	    << std::endl;
	std::cerr << "Inkove \"isONclust2 help\" for more details."
		  << std::endl;
	exit(1);
    }
    auto subCmd = std::string(argv[1]);

    if (subCmd == "help" || subCmd == "-h") {
	print_help();
	exit(0);
    }
    else if (subCmd == "sort") {
	mainSort(argc, argv);
    }
    else if (subCmd == "cluster") {
	mainCluster(argc, argv);
    }
    else if (subCmd == "dump") {
	mainDump(argc, argv);
    }
    else if (subCmd == "info") {
	mainInfo(argc, argv);
    }
    else if (subCmd == "version") {
	print_version();
	exit(0);
    }
    else {
	std::cerr << "Invalid subcommand: " << subCmd << std::endl;
	std::cerr
	    << "Valid subcommands are: sort, cluster, dump, info, version, help"
	    << std::endl;
	std::cerr << "Inkove \"isONclust2 help\" for more details."
		  << std::endl;
	exit(1);
    }

    return 0;
}

int mainSort(int argc, char* argv[])
{
    auto cmdArgs = ParseArgsSort(argc, argv);
    VERBOSE = cmdArgs->Verbose;

    if (VERBOSE) {
	cerr << "isONclust2 version: " << isONclust2_VERSION << endl;
	cerr << "Batches output directory: " << cmdArgs->BatchOutFolder << endl;
	cerr << "Minimum batch size: " << cmdArgs->BatchSize << " kilobases"
	     << endl;
	cerr << "Kmer size: " << cmdArgs->KmerSize << endl;
	cerr << "Window size: " << cmdArgs->WindowSize << endl;
	cerr << "Minimum average quality: " << cmdArgs->MinQual << endl;
	cerr << "Minimum shared minimizers: " << cmdArgs->MinShared << endl;
	cerr << "Minimum fraction of top minimizer hit: "
	     << cmdArgs->MinFraction << endl;
	cerr << "Mapping threshold: " << cmdArgs->MappedThreshold << endl;
	cerr << "Alignment threshold: " << cmdArgs->AlignedThreshold << endl;
	cerr << "Minimum probability no hit: " << cmdArgs->MinProbNoHits
	     << endl;
	cerr << "Debug output: " << (cmdArgs->Debug ? "on" : "off") << endl;
    }

    auto batchDir = cmdArgs->BatchOutFolder + "/batches";
    CreateOutdir(cmdArgs->BatchOutFolder);
    CreateOutdir(batchDir);

    auto fqParser =
	bioparser::createParser<bioparser::FastqParser, Seq>(cmdArgs->InFastq);

    SequencesP sequences;
    sequences.reserve(10000000);
    fqParser->parse(sequences, -1, false);

    if (VERBOSE) {
	cerr << "Parsed " << sequences.size() << " sequences." << endl;
    }

    auto qualTab = InitQualTab();
    auto qualTabNomin = InitQualTabNomin();
    FillQualScores(sequences, cmdArgs->KmerSize, cmdArgs->WindowSize, qualTab,
		   qualTabNomin);
    SortByQualScores(sequences);

    if (VERBOSE) {
	cerr << "Finished sorting sequences." << endl;
    }
    string sortedFastq = cmdArgs->BatchOutFolder + "/sorted_reads.fastq";
    SequencesPToFastq(sequences, sortedFastq);
    if (VERBOSE) {
	cerr << "Sorted sequences written to: " << sortedFastq << endl;
    }
    string scoresTsv = cmdArgs->BatchOutFolder + "/scores.tsv";
    WriteScores(sequences, scoresTsv);
    if (VERBOSE) {
	cerr << "Scores written to: " << scoresTsv << endl;
    }

    if (VERBOSE) {
	cerr << "Preparing batches:" << endl;
    }
    unsigned long batchBases{0};
    int batchSeqs{0};
    int nrBatches{0};
    int batchStart{0};
    unsigned i = 0;

    for (i = 0; i < sequences.size(); i++) {
	batchBases += sequences[i]->Str().length();
	batchSeqs++;

	if ((cmdArgs->BatchSize > 0) &&
	    ((batchBases > (unsigned long)(cmdArgs->BatchSize * 1000)) ||
	     ((cmdArgs->BatchMaxSeq > 0) &&
	      batchSeqs >= cmdArgs->BatchMaxSeq))) {
	    const auto batch = std::unique_ptr<Batch>(PrepareSortedBatch(
		sequences, batchStart, i, batchBases, cmdArgs->KmerSize,
		cmdArgs->WindowSize, cmdArgs->MinQual, qualTab, qualTabNomin));
	    batch->BatchNr = nrBatches;
	    batch->BatchBases = batchBases;
	    batch->SortArgs = *cmdArgs;

	    auto outFile =
		batchDir + "/isONbatch_" + std::to_string(nrBatches) + ".cer";
	    SaveBatch(batch, outFile);

	    if (VERBOSE) {
		cerr << "\tWritten batch " << nrBatches << " with "
		     << (i - batchStart + 1) << " sequences and "
		     << int((double(batchBases) / 1000.0)) << " kilobases."
		     << endl;
	    }

	    batchBases = 0;
	    batchSeqs = 0;
	    batchStart = i + 1;
	    nrBatches++;
	}
    }

    if (batchStart < int(sequences.size())) {
	auto batch = std::unique_ptr<Batch>(PrepareSortedBatch(
	    sequences, batchStart, sequences.size() - 1, batchBases,
	    cmdArgs->KmerSize, cmdArgs->WindowSize, cmdArgs->MinQual, qualTab,
	    qualTabNomin));
	batch->BatchNr = nrBatches;
	batch->BatchBases = batchBases;
	batch->SortArgs = *cmdArgs;
	auto outFile =
	    batchDir + "/isONbatch_" + std::to_string(nrBatches) + ".cer";
	SaveBatch(batch, outFile);

	if (VERBOSE) {
	    cerr << "\tWritten batch " << nrBatches << " with "
		 << (i - batchStart + 1) << " sequences and "
		 << int((double(batchBases) / 1000.0)) << " kilobases." << endl;
	}
    }

    return 0;
}

int mainDump(int argc, char* argv[])
{
    auto cmdArgs = ParseArgsDump(argc, argv);
    VERBOSE = cmdArgs->Verbose;
    auto b = LoadBatch(cmdArgs->InCereal);
    if (VERBOSE) {
	cerr << "Loaded batch from " << cmdArgs->InCereal << ":" << std::endl;
	printBatchInfo(b);
    }
    b->MinDB = MinimizerDB(0, uh);
    CreateOutdir(cmdArgs->OutDir);
    dumpBatchInfo(b, cmdArgs->OutDir + "/batch_info.tsv");
    if (VERBOSE) {
	cerr << "Sorting clusters by size." << std::endl;
    }
    if (b->Cls.size() > 0) {
	SortClustersBySize(b->Cls);
	dumpClusters(b, cmdArgs->OutDir);
    }
    if (VERBOSE) {
	cerr << "Dump complete." << std::endl;
    }

    return 0;
}

int mainCluster(int argc, char* argv[])
{
    auto cmdArgs = ParseArgsCluster(argc, argv);
    VERBOSE = cmdArgs->Verbose;
    bool SINGLE = false;
    if (cmdArgs->RightCereal == "") {
	SINGLE = true;
    }

    auto leftBatch = LoadBatch(cmdArgs->LeftCereal);
    if (VERBOSE) {
	cerr << "Loaded input batch from " << cmdArgs->LeftCereal << ":"
	     << std::endl;
	printBatchInfo(leftBatch);
    }
    std::unique_ptr<Batch> rightBatch;
    if (!SINGLE) {
	rightBatch = LoadBatch(cmdArgs->RightCereal);
	cerr << "Loaded input batch from " << cmdArgs->RightCereal << ":"
	     << std::endl;
	rightBatch->MinDB = MinimizerDB(MIN_DB_RESERVE);
	;
	printBatchInfo(rightBatch);
    }
    else {
	rightBatch = CreatePseudoBatch(leftBatch);
	if (VERBOSE) {
	    cerr << "Created pseudo-batch for single clustering:" << std::endl;
	    printBatchInfo(rightBatch);
	    cerr << "Resetting input clusters." << endl;
	}
	leftBatch->Cls.clear();
	leftBatch->NrCls = 0;
    }

    if (leftBatch->SortArgs.Mode != cmdArgs->Mode) {
	leftBatch->SortArgs.Mode = cmdArgs->Mode;
    }

    if (rightBatch->SortArgs.Mode != cmdArgs->Mode) {
	rightBatch->SortArgs.Mode = cmdArgs->Mode;
    }

    ClusterSortedReads(leftBatch, rightBatch, cmdArgs->Quiet);

    if (VERBOSE) {
	cerr << "Finished clustering!" << endl;
	cerr << "Alignment invocation count: " << AlnInvoked() << " (";
	cerr << AlnInvokedPerc(rightBatch->Cls.size()) << "%)" << endl;

	unsigned count{};
	for (auto& c : leftBatch->Cls) {
	    if (c.size() > 1) {
		count++;
	    }
	}
	cerr << "Number of clusters larger than 1: " << count << endl;
	cerr << "Output batch statistics:" << endl;
	printBatchInfo(leftBatch);
    }
    leftBatch->LeftLeaf = cmdArgs->LeftCereal;
    leftBatch->RightLeaf = cmdArgs->RightCereal;
    SaveBatch(leftBatch, cmdArgs->OutCereal);
    if (VERBOSE) {
	cerr << "Output batch written to: " << cmdArgs->OutCereal << endl;
    }
    return 0;
}

int mainInfo(int argc, char* argv[])
{
    auto target = std::string(argv[2]);
    if (target == "-h") {
	cerr << "help message" << endl;
	exit(0);
    }
    auto b = LoadBatch(target);
    cerr << "Loaded batch from " << target << ":" << std::endl;
    printBatchInfo(b);
    return 0;
}

void printBatchInfo(BatchP& b)
{
    cerr << "\tBatch number: " << b->BatchNr << endl;
    cerr << "\tBatch range: [" << b->BatchStart << "," << b->BatchEnd << "]"
	 << endl;
    cerr << "\tDepth: " << b->Depth << endl;
    cerr << "\tNr sequences: " << b->BatchEnd - b->BatchStart + 1 << endl;
    cerr << "\tNr bases: " << b->BatchBases << endl;
    cerr << "\tNr clusters: " << b->NrClusters() << endl;
    cerr << "\tMinimizers in database: " << b->MinDBSize() << endl;
}

void dumpBatchInfo(BatchP& b, std::string outfile)
{
    std::ofstream out;
    CreateFile(outfile, out);
    out << "Name\tValue" << endl;
    out << "BatchNumber\t" << b->BatchNr << endl;
    out << "BatchStart\t" << b->BatchStart << endl;
    out << "BatchEnd\t" << b->BatchEnd << endl;
    out << "Depth\t" << b->Depth << endl;
    out << "NrBases\t" << b->BatchBases << endl;
    out << "NrClusters\t" << b->NrClusters() << endl;
    out << "MinDBsize\t" << b->MinDBSize() << endl;
    out.close();
}

void dumpClusters(BatchP& b, std::string outdir)
{
    std::ofstream outInfo;
    std::string outfile = outdir + "/clusters_info.tsv";
    std::string clsdir = outdir + "/cluster_fastq";
    CreateFile(outfile, outInfo);
    CreateOutdir(clsdir);
    outInfo << "ClusterId\tSize" << endl;
    unsigned i = 0;
    for (auto& c : b->Cls) {
	outInfo << i << "\t" << c.size() - 1 << endl;
	i++;
    }
    outInfo.close();
    WriteClusters(b->Cls, outdir);
}
