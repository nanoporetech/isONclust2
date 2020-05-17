#ifndef ARGS_H_INCLUDED
#define ARGS_H_INCLUDED
#include "isONclust2_config.h"

#include <string>

typedef enum { Sahlin, Fast, Furious, None } ClsMode;

struct CmdArgs {
    bool Verbose{};
    bool Debug{};
    std::string InFastq;
    int KmerSize{11};
    int BatchSize{50000};
    int BatchMaxSeq{30000};
    int WindowSize{15};
    int MinShared{5};
    int ConsMinSize{50};
    int ConsMaxSize{-150};
    int ConsPeriod{500};
    int MinClsSize{3};
    double MinQual{7.0};
    double MappedThreshold{0.65};
    double AlignedThreshold{0.2};
    double MinFraction{0.8};
    double MinProbNoHits{0.1};
    std::string BatchOutFolder{"isONclust2_batches"};
    ClsMode Mode{Sahlin};
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(Verbose, Debug, InFastq, KmerSize, BatchSize, BatchMaxSeq,
		WindowSize, MinShared, ConsMinSize, ConsMaxSize, ConsPeriod,
		MinClsSize, MinQual, MappedThreshold, AlignedThreshold,
		MinFraction, MinProbNoHits, BatchOutFolder, Mode);
    }
};

bool operator==(const CmdArgs& lhs, const CmdArgs& rgh);
bool operator!=(const CmdArgs& lhs, const CmdArgs& rgh);

struct CmdArgsCluster {
    bool Verbose{};
    bool Quiet{};
    bool Debug{};
    bool MinPurge{};
    bool SeqPurge{};
    int MinClsSize{-1};
    std::string LeftCereal{""};
    std::string RightCereal{""};
    std::string OutCereal{""};
    ClsMode Mode{None};
    int SpoaAlgo{2};
};

struct CmdArgsDump {
    bool Verbose{};
    bool Debug{};
    std::string InCereal{""};
    std::string Index{""};
    std::string OutDir{""};
};

struct CmdArgs* ParseArgsSort(int argc, char** argv);
struct CmdArgsCluster* ParseArgsCluster(int argc, char** argv);
struct CmdArgsDump* ParseArgsDump(int argc, char** argv);
void print_help();
void print_version();
void print_help_sort();
void print_help_cluster();
void print_help_dump();
void print_help_info();

#endif
