#ifndef ARGS_H_INCLUDED
#define ARGS_H_INCLUDED
#include "isONclust2_config.h"

#include <string>

typedef enum { Sahlin, Fast, Furious } ClsMode;

struct CmdArgs {
    bool Verbose{};
    bool Debug{};
    std::string InFastq;
    int KmerSize{13};
    int BatchSize{50000};
    int BatchMaxSeq{30000};
    int WindowSize{20};
    int MinShared{5};
    int ConsMinSize{50};
    int ConsMaxSize{150};
    int ConsPeriod{-1};
    double MinQual{7.0};
    double MappedThreshold{0.7};
    double AlignedThreshold{0.4};
    double MinFraction{0.8};
    double MinProbNoHits{0.1};
    std::string BatchOutFolder{"isONclust2_batches"};
    ClsMode Mode{Sahlin};
    template <class Archive>
    void serialize(Archive& archive)
    {
	archive(Verbose, Debug, InFastq, KmerSize, BatchSize, WindowSize,
		MinShared, ConsMinSize, ConsMaxSize, MinQual, MappedThreshold,
		AlignedThreshold, MinFraction, MinProbNoHits, BatchOutFolder);
    }
};

bool operator==(const CmdArgs& lhs, const CmdArgs& rgh);
bool operator!=(const CmdArgs& lhs, const CmdArgs& rgh);

struct CmdArgsCluster {
    bool Verbose{};
    bool Quiet{};
    bool Debug{};
    std::string LeftCereal{""};
    std::string RightCereal{""};
    std::string OutCereal{""};
    ClsMode Mode{Sahlin};
};

struct CmdArgsDump {
    bool Verbose{};
    bool Debug{};
    std::string InCereal{""};
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
