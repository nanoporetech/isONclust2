#include "output.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "cluster_data.h"
#include "qualscore.h"
#include "seq.h"
#include "util.h"

// From:
// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
int dirExists(const std::string path)
{
    struct stat info;

    if (stat(path.c_str(), &info) != 0)
	return 0;
    else if (info.st_mode & S_IFDIR)
	return 1;
    else
	return 0;
}

void SequencesPToFastq(SequencesP& sequences, const std::string& outFastq)
{
    std::ofstream outfile;
    CreateFile(outFastq, outfile);

    for (auto& s : sequences) {
	WriteFastqRecord(*s, outfile);
    }
    outfile.close();
}

void CreateFile(const std::string& outFile, std::ofstream& outfile)
{
    outfile.open(outFile);
    if (!outfile.is_open()) {
	std::cerr << "Failed to open " + outFile + "!" << std::endl;
	exit(1);
    }
}

void WriteFastqRecord(const Seq& s, std::ofstream& out)
{
    out << "@" << s.Name() << std::endl;
    out << s.Str() << std::endl;
    out << "+" << std::endl;
    out << s.Qual() << std::endl;
}

void WriteScores(SequencesP& sequences, const std::string& outFile)
{
    std::ofstream outfile;
    CreateFile(outFile, outfile);

    for (auto& s : sequences) {
	outfile << s->Name() << "\t" << s->Score() << std::endl;
    }
    outfile.close();
}

void CreateOutdir(const std::string& outDir)
{
    if (dirExists(outDir)) {
	std::cerr << "Warning: reusing existing output directory: " << outDir
		  << std::endl;
	return;
    }
    int result = mkdir(outDir.c_str(), 0755);
    if (result != 0) {
	std::cerr << "Failed to create output directory!"
		  << "Error: " << result << std::endl;
	exit(1);
    }
}

void sortMembers(Cluster& v)
{
    auto it = v.begin();
    it++;
    std::stable_sort(
	it, v.end(),
	[](const shared_ptr<ProcSeq> a, const shared_ptr<ProcSeq> b) {
	    return a->RawSeq.Score() > b->RawSeq.Score();
	});
}

void WriteClusters(Clusters& cls, const std::string& outDir)
{
    std::ofstream outfile;
    std::ofstream outfq;
    std::ofstream outcons;
    std::string outFile = outDir + "/clusters.tsv";
    std::string outCons = outDir + "/cluster_cons.fq";
    CreateFile(outFile, outfile);
    CreateFile(outCons, outcons);
    outfile << "ClusterId\tStrand\tRead\tLength\tScore\tMeanQual" << std::endl;
    if (VERBOSE) {
	std::cerr << "Writing out clusters:" << std::endl;
    }
    for (unsigned i = 0; i < cls.size(); i++) {
	sortMembers(*cls[i]);
	CreateFile(
	    outDir + "/cluster_fastq/isONcluster_" + std::to_string(i) + ".fq",
	    outfq);
	if (VERBOSE) {
	    Pbar((float)(i + 1) / float(cls.size()));
	}
	int nrReads = 0;
	for (auto& read : *cls[i]) {
	    if (nrReads == 0) {
		auto& s = read->RawSeq;
		if (s.Score() < 0) {
		    continue;
		}
		auto seq = std::string(read->RawSeq.Str());
		if (read->MatchStrand == -1) {
		    seq = RevComp(seq);
		}
		outcons << "@cluster_" << i << " origin=" << s.Name() << ":"
			<< read->MatchStrand << " length=" << seq.length()
			<< " size=" << cls[i]->size() - 1 << std::endl;
		outcons << seq << std::endl;
		outcons << "+" << std::endl;
		outcons << s.Qual() << std::endl;
		nrReads++;
		continue;
	    }
	    outfile << i << "\t" << read->MatchStrand << "\t"
		    << read->RawSeq.Name();
	    outfile << "\t" << read->RawSeq.Str().length() << "\t"
		    << read->RawSeq.Score() << "\t"
		    << -10 * log10(read->RawSeq.ErrorRate()) << std::endl;
	    auto& s = read->RawSeq;
	    auto seq = std::string(read->RawSeq.Str());
	    if (read->MatchStrand == -1) {
		seq = RevComp(seq);
	    }
	    outfq << "@" << s.Name() << " cluster=" << i << ":"
		  << read->MatchStrand << std::endl;
	    outfq << seq << std::endl;
	    outfq << "+" << std::endl;
	    outfq << s.Qual() << std::endl;
	    nrReads++;
	}
	outfq.close();
    }
    std::cerr << std::endl;
    outfile.close();
}

