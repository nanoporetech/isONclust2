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
	if (s->Score() < 0) {
	    continue;
	}
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

typedef std::unique_ptr<FqRec> FqRecP;

FqRecP GetRecords(std::ifstream& infq)
{
    FqRecP rec(new FqRec);
    std::getline(infq, rec->Header);
    if (infq.eof()) {
	return nullptr;
    }
    std::getline(infq, rec->Seq);
    std::getline(infq, rec->Plus);
    std::getline(infq, rec->Qual);
    auto ididx = rec->Header.find(" ");
    rec->Id = rec->Header.substr(1, ididx);
    return rec;
}

void WriteFqRec(FqRecP& r, std::ofstream& fh)
{
    fh << r->Header << std::endl;
    fh << r->Seq << std::endl;
    fh << r->Plus << std::endl;
    fh << r->Qual << std::endl;
}

void WriteClusters(BatchP& b, const std::string& outDir, SortedIdx* idx,
		   IdMap& idToCls)
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

    outfile << "ClusterId\tStrand\tRead" << std::endl;
    if (VERBOSE) {
	std::cerr << "Writing out cluster information:" << std::endl;
    }

    auto& cls = b->Cls;
    for (unsigned i = 0; i < cls.size(); i++) {
	if (VERBOSE) {
	    Pbar((float)(i + 1) / float(cls.size()));
	}
	auto reads = cls[i];
	if (reads == nullptr) {
	    std::cerr << "Null pointer instead of cluster rep at index: " << i
		      << std::endl;
	    exit(1);
	}
	auto read = reads->at(0);
	if (read->RawSeq == nullptr) {
	    std::cerr << "Null pointer instead of cluster rep sequence "
			 "at index: "
		      << i << std::endl;
	    exit(1);
	}
	auto& s = read->RawSeq;
	if (s->Score() < 0) {
	    continue;
	}
	auto seq = std::string(read->RawSeq->Str());
	auto qual = std::string(read->RawSeq->Qual());
	if (read->MatchStrand == -1) {
	    seq = RevComp(seq);
	    std::reverse(qual.begin(), qual.end());
	}
	outcons << "@cluster_" << i << " origin=" << s->Name() << ":"
		<< read->MatchStrand << " length=" << seq.length()
		<< " size=" << cls[i]->size() - 1 << std::endl;
	outcons << seq << std::endl;
	outcons << "+" << std::endl;
	outcons << s->Qual() << std::endl;  // FIXME
    }
    outcons.close();
    if (VERBOSE) {
	std::cerr << std::endl;
    }
    b->Cls = Clusters();

    SeqCache seqCache;

    if (VERBOSE) {
	std::cerr << "Loading fastq records:" << std::endl;
    }
    unsigned j = 0;
    unsigned jj = 0;
    for (auto rec = GetRecords(infq); rec != nullptr; rec = GetRecords(infq)) {
	if (VERBOSE) {
	    auto now = ((float)(j + 1) / float(idToCls.size()));
	    if (unsigned(now) > jj || j == 0) {
		Pbar(now);
		jj = unsigned(now);
	    }
	}
	auto readId = rec->Header.substr(1, rec->Header.size());
	auto v = idToCls.find(readId);
	if (v == idToCls.end()) {
	    continue;
	}

	if (v->second->Strand == -1) {
	    rec->Seq = RevComp(rec->Seq);
	    std::reverse(rec->Qual.begin(), rec->Qual.end());
	}

	outfile << v->second->Cls << "\t" << v->second->Strand << "\t" << readId
		<< std::endl;
	seqCache[v->second->Cls].push_back(std::move(rec));
	j++;
    }
    outfile.close();
    if (VERBOSE) {
	std::cerr << std::endl;
    }
    idToCls.clear();

    if (VERBOSE) {
	std::cerr << "Writing out cluster fastqs:" << std::endl;
    }

    unsigned k = 0;
    unsigned kk = 0;
    for (auto& c : seqCache) {
	if (VERBOSE) {
	    auto now = ((float)(k + 1) / float(seqCache.size()));
	    if (unsigned(now) > kk || k == 0 || kk == 0) {
		Pbar(now);
		kk = unsigned(now);
	    }
	}
	std::ofstream outfq;
	CreateFile(outDir + "/cluster_fastq/" + std::to_string(c.first) + ".fq",
		   outfq);
	for (auto& r : c.second) {
	    WriteFqRec(r, outfq);
	}
	outfq.flush();
	outfq.close();
	k++;
    }

    if (VERBOSE) {
	std::cerr << std::endl;
    }
}

