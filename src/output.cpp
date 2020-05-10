#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include "cluster_data.h"
#include "qualscore.h"
#include "seq.h"
#include "util.h"

#include "output.h"

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

void SequencesPToFastq(SequencesP& sequences, const std::string& outFastq,
		       const std::string& indexTab, const std::string& indexCer)
{
    std::ofstream outfile;
    std::ofstream outTsv;
    CreateFile(outFastq, outfile);
    CreateFile(indexTab, outTsv);
    outTsv << "Id\tPos" << std::endl;
    unsigned long long int seeker = 0;
    SortedIdx idx;
    idx.Fastq = outFastq;

    for (auto& s : sequences) {
	idx.IdxMap[s->Name()] = seeker;
	outTsv << s->Name() << "\t" << seeker << std::endl;
	seeker += WriteFastqRecord(*s, outfile);
    }
    outfile.close();
    outTsv.close();
    std::ofstream os(indexCer, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(idx);
}

void CreateFile(const std::string& outFile, std::ofstream& outfile)
{
    outfile.open(outFile);
    if (!outfile.is_open()) {
	std::cerr << "Failed to open " + outFile + "!" << std::endl;
	exit(1);
    }
}

void OpenFile(const std::string& inFile, std::ifstream& infile)
{
    infile.open(inFile);
    if (!infile.is_open()) {
	std::cerr << "Failed to open " + inFile + "!" << std::endl;
	exit(1);
    }
}

unsigned WriteFastqRecord(const Seq& s, std::ofstream& out)
{
    out << "@" << s.Name() << std::endl;
    out << s.Str() << std::endl;
    out << "+" << std::endl;
    out << s.Qual() << std::endl;
    return (s.Name().length() + s.Str().length() + s.Qual().length() + 6);
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

std::unique_ptr<SortedIdx> LoadIndex(std::string inf)
{
    std::ifstream instream(inf, std::ios::binary);

    auto p = std::unique_ptr<SortedIdx>(new SortedIdx);
    cereal::BinaryInputArchive iarchive(instream);
    try {
	iarchive(*p);
    }
    catch (std::runtime_error e) {
	std::cerr << "Failed to load index " << inf << ":" << e.what()
		  << std::endl;
	exit(1);
    }
    return p;
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

std::vector<std::string> RecordAtPos(unsigned pos, std::ifstream& infq)
{
    std::vector<std::string> res(4);
    infq.seekg(pos, infq.beg);
    std::getline(infq, res[0]);
    std::getline(infq, res[1]);
    std::getline(infq, res[2]);
    std::getline(infq, res[3]);
    return res;
}

void WriteClusters(Clusters& cls, const std::string& outDir, SortedIdx* idx)
{
    std::ofstream outfile;
    std::ofstream outfq;
    std::ofstream outcons;
    std::ifstream infq;
    OpenFile(idx->Fastq, infq);
    std::string outFile = outDir + "/clusters.tsv";
    std::string outCons = outDir + "/cluster_cons.fq";
    CreateFile(outFile, outfile);
    CreateFile(outCons, outcons);
    outfile << "ClusterId\tStrand\tRead\tLength\tScore\tMeanQual" << std::endl;
    if (VERBOSE) {
	std::cerr << "Writing out clusters:" << std::endl;
    }
    for (unsigned i = 0; i < cls.size(); i++) {
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
		if (s->Score() < 0) {
		    continue;
		}
		auto seq = std::string(read->RawSeq->Str());
		if (read->MatchStrand == -1) {
		    seq = RevComp(seq);
		}
		outcons << "@cluster_" << i << " origin=" << s->Name() << ":"
			<< read->MatchStrand << " length=" << seq.length()
			<< " size=" << cls[i]->size() - 1 << std::endl;
		outcons << seq << std::endl;
		outcons << "+" << std::endl;
		outcons << s->Qual() << std::endl;  // FIXME
		nrReads++;
		continue;
	    }
	    outfile << i << "\t" << read->MatchStrand << "\t" << read->Id
		    << std::endl;

	    auto rec = RecordAtPos(idx->IdxMap[read->Id], infq);
	    if (read->MatchStrand == -1) {
		// std::cerr << read->Id << std::endl;
		rec[1] = RevComp(rec[1]);
		reverse(rec[3].begin(), rec[3].end());
	    }
	    stringstream ss;
	    ss << rec[0];
	    ss << " cluster=" << i << ":" << read->MatchStrand << std::endl;
	    outfq << ss.str();
	    outfq << rec[1] << std::endl;
	    outfq << rec[2] << std::endl;
	    outfq << rec[3] << std::endl;
	}
	outfq.close();
    }
    std::cerr << std::endl;
    outfile.close();
}

