#include "args.h"
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;
bool VERBOSE = false;

/// Parse command line arguments.
struct CmdArgs* ParseArgsSort(int argc, char* argv[])
{
    const struct option longopts[] = {
	{"version", no_argument, 0, 'V'},
	{"verbose", no_argument, 0, 'v'},
	{"debug", no_argument, 0, 'd'},
	{"mode", optional_argument, 0, 'x'},
	{"help", no_argument, 0, 'h'},
	{"kmer-size", required_argument, 0, 'k'},
	{"window-size", required_argument, 0, 'w'},
	{"min-shared", required_argument, 0, 'm'},
	{"mapped-threshold", required_argument, 0, 'r'},
	{"aligned-threshold", required_argument, 0, 'a'},
	{"min-fraction", required_argument, 0, 'f'},
	{"min-prob-no-hits", required_argument, 0, 'p'},
	{"min-qual", required_argument, 0, 'q'},
	{"low-cons-size", required_argument, 0, 'g'},
	{"max-cons-size", required_argument, 0, 'c'},
	{"cons-period", required_argument, 0, 'P'},
	{"outfolder", required_argument, 0, 'o'},
	{"batch-size", required_argument, 0, 'B'},
	{"batch-max-seq", required_argument, 0, 'M'},
	{0, 0, 0, 0},
    };

    int index;
    int iarg = 0;
    const string mSahlin = std::string("sahlin");
    const string mFast = std::string("fast");
    const string mFurious = std::string("furious");

    auto res = new struct CmdArgs;
    argc--;
    auto sargv = new char*[argc];
    std::copy(argv + 1, argv + 1 + argc, sargv);

    while (iarg != -1) {
	iarg = getopt_long(
	    argc, sargv, "k:w:dhvo:m:r:a:f:p:q:B:x:g:c:M:P:", longopts, &index);

	switch (iarg) {
	    case 'h':
		print_help_sort();
		exit(0);
		break;
	    case 'k':
		res->KmerSize = atoi(optarg);
		break;
	    case 'm':
		res->MinShared = atoi(optarg);
		break;
	    case 'r':
		res->MappedThreshold = atof(optarg);
		break;
	    case 'a':
		res->AlignedThreshold = atof(optarg);
		break;
	    case 'f':
		res->MinFraction = atof(optarg);
		break;
	    case 'p':
		res->MinProbNoHits = atof(optarg);
		break;
	    case 'q':
		res->MinQual = atof(optarg);
		break;
	    case 'w':
		res->WindowSize = atoi(optarg);
		break;
	    case 'B':
		res->BatchSize = atoi(optarg);
		break;
	    case 'M':
		res->BatchMaxSeq = atoi(optarg);
		break;
	    case 'P':
		res->ConsPeriod = atoi(optarg);
		break;
	    case 'g':
		res->ConsMinSize = atoi(optarg);
		break;
	    case 'c':
		res->ConsMaxSize = atoi(optarg);
		break;
	    case 'o':
		res->BatchOutFolder = optarg;
		break;
	    case 'v':
		res->Verbose = true;
		break;
	    case 'd':
		res->Debug = true;
		break;
	    case 'x':
		string m = string(optarg);
		if (m == mSahlin) {
		    res->Mode = Sahlin;
		}
		else if (m == mFast) {
		    res->Mode = Fast;
		}
		else if (m == mFurious) {
		    res->Mode = Furious;
		}
		else {
		    cerr << "Illegal clustering mode: " << m << endl;
		    print_help_sort();
		    exit(1);
		}

		break;
	}
    }

    if (argc - optind != 1) {
	cerr << "Please specify one input fastq file!" << endl;
	exit(1);
    }

    if (res->KmerSize > 31) {
	cerr << "Maximum supported kmer size is 31!" << endl;
	exit(1);
    }

    if (res->KmerSize < 10) {
	cerr << "Minimum supported kmer size is 10!" << endl;
	exit(1);
    }

    if (res->KmerSize > res->WindowSize) {
	cerr << "Kmer size cannot be larger than the window size!" << endl;
	exit(1);
    }

    res->InFastq = sargv[optind];
    return res;
}

/// Parse command line arguments.
struct CmdArgsCluster* ParseArgsCluster(int argc, char* argv[])
{
    const struct option longopts[] = {
	{"quiet", no_argument, 0, 'Q'},
	{"version", no_argument, 0, 'V'},
	{"verbose", no_argument, 0, 'v'},
	{"debug", no_argument, 0, 'd'},
	{"spoa-algo", optional_argument, 0, 'A'},
	{"mode", optional_argument, 0, 'x'},
	{"help", no_argument, 0, 'h'},
	{"outfile", required_argument, 0, 'o'},
	{"left-batch", required_argument, 0, 'l'},
	{"right-batch", optional_argument, 0, 'r'},
	{0, 0, 0, 0},
    };

    int index;
    int iarg = 0;
    const string mSahlin = std::string("sahlin");
    const string mFast = std::string("fast");
    const string mFurious = std::string("furious");

    auto res = new struct CmdArgsCluster;
    argc--;
    auto sargv = new char*[argc];
    std::copy(argv + 1, argv + 1 + argc, sargv);

    while (iarg != -1) {
	iarg = getopt_long(argc, sargv, "Vdhvo:l:r:Qx:A:", longopts, &index);

	switch (iarg) {
	    case 'h':
		print_help_cluster();
		exit(0);
		break;
	    case 'o':
		res->OutCereal = optarg;
		break;
	    case 'l':
		res->LeftCereal = optarg;
		break;
	    case 'r':
		res->RightCereal = optarg;
		break;
	    case 'v':
		res->Verbose = true;
		break;
	    case 'Q':
		res->Quiet = true;
		break;
	    case 'd':
		res->Debug = true;
		break;
	    case 'A':
		res->SpoaAlgo = atoi(optarg);
		break;
	    case 'x':
		string m = string(optarg);
		if (m == mSahlin) {
		    res->Mode = Sahlin;
		}
		else if (m == mFast) {
		    res->Mode = Fast;
		}
		else if (m == mFurious) {
		    res->Mode = Furious;
		}
		else {
		    cerr << "Illegal clustering mode: " << m << endl;
		    print_help_sort();
		    exit(1);
		}

		break;
	}
    }

    if (res->LeftCereal == "") {
	cerr << "Specifying left input batch is mandatory!" << endl;
	print_help();
	exit(1);
    }

    if (res->OutCereal == "") {
	cerr << "Specifying output batch file is mandatory!" << endl;
	print_help();
	exit(1);
    }

    return res;
}

// Parse command line arguments.
struct CmdArgsDump* ParseArgsDump(int argc, char* argv[])
{
    const struct option longopts[] = {
	{"verbose", no_argument, 0, 'v'},
	{"debug", no_argument, 0, 'd'},
	{"help", no_argument, 0, 'h'},
	{"outdir", required_argument, 0, 'o'},
	{0, 0, 0, 0},
    };

    int index;
    int iarg = 0;

    auto res = new struct CmdArgsDump;
    argc--;
    auto sargv = new char*[argc];
    std::copy(argv + 1, argv + 1 + argc, sargv);

    while (iarg != -1) {
	iarg = getopt_long(argc, sargv, "dhvo:", longopts, &index);

	switch (iarg) {
	    case 'h':
		print_help_dump();
		exit(0);
		break;
	    case 'o':
		res->OutDir = optarg;
		break;
	    case 'v':
		res->Verbose = true;
		break;
	    case 'd':
		res->Debug = true;
		break;
	}
    }

    if (res->OutDir == "") {
	cerr << "Specifying output directory is mandatory!" << endl;
	print_help();
	exit(1);
    }

    res->InCereal = sargv[optind];
    return res;
}

void print_help()
{
    print_version();
    cout << "Available subcommands: sort, cluster, dump, info, help, version"
	 << endl;
    print_help_sort();
    print_help_cluster();
    print_help_dump();
    print_help_info();
    cout << "\nhelp - print help message\n";
    cout << "\nversion - print version\n";
}

int ConsMinSize{20};
int ConsMaxSize{100};
int ConsPeriod{400};
void print_help_sort()
{
    cout << "\nsort - sort reads and write out batches:\n"
	    "\t-B --batch-size        Batch size in kilobases (default: "
	    "\t-M --batch-max-seq     Maximum number of sequences per batch "
	    "(default: "
	    "3000).\n"
	    "\t-k --kmer-size         Kmer size (default: 11).\n"
	    "\t-w --window-size       Window size (default: 15).\n"
	    "\t-m --min-shared        Minimum number of minimizers shared "
	    "between "
	    "read and cluster (default: 5).\n"
	    "\t-q --min-qual          Minimum average quality value "
	    "(default: "
	    "7.0).\n"
	    "\t-x --mode  Clustering mode:\n"
	    "\t           * sahlin (default): use minimizers first, alignment "
	    "second \n"
	    "\t           * fast: use minimizers only\n"
	    "\t           * furious: always use alignment\n"
	    "\t-g --low-cons-size     Use all sequences for consensus below "
	    "this size (default: "
	    "20).\n"
	    "\t-c --max-cons-size     Maximum number of sequences used for "
	    "consensus (default: "
	    "150).\n"
	    "\t-P --cons-period       Do not recalculate consensus after this "
	    "many seuqences added (default: 500).\n"
	    "\t-r --mapped-threshold  Minmum mapped fraction of read to be "
	    "\tincluded in cluster (default: 0.65).\n"
	    "\t-a --aligned-threshold Minimum aligned fraction of read to "
	    "be included in cluster (default: 0.2).\n"
	    "\t-f --min-fraction      Minimum fraction of minimizers shared "
	    "compared to best hit, in order to continue mapping (default: "
	    "0.8).\n"
	    "\t-p --min-prob-no-hits  Minimum probability for i consecutive "
	    "\tminimizers to be different between read and representative "
	    "(default: 0.1)\n"
	    "\t-o --outfolder         Output folder (default: "
	    "\t./isONclust2_batches).\n"
	    "\t-h --help              Print help.\n"
	    "\t-v --verbose           Verbose output.\n"
	    "\t-d --debug             Print debug info.\n"
	    "\t[positional argument]  Input fastq file (required).\n";
}

void print_help_cluster()
{
    cout << "\ncluster - cluster and/or merge batches:\n"
	    "\t-l --left-batch        Left input batch (mandatory).\n"
	    "\t-r --right-batch       Right input batch (optional).\n"
	    "\t-o --outfile           Output batch.\n"
	    "\t-x --mode  Clustering mode:\n"
	    "\t           * sahlin (default): use minimizers first, alignment "
	    "second \n"
	    "\t           * fast: use minimizers only\n"
	    "\t           * furious: use alignment only\n"
	    "\t-A --spoa-algo  spoa alignment algorithm:\n"
	    "\t           * 0 (default): local\n"
	    "\t           * 1 : global\n"
	    "\t           * 1 : semi-global\n"
	    "\t-v --verbose           Verbose output.\n"
	    "\t-Q --quiet             Supress progress bar.\n"
	    "\t-d --debug             Print debug info.\n"
	    "\t-h --help              Print help.\n";
};
void print_help_dump()
{
    cout << "\ndump - dump clustered batch:\n"
	    "\t-o --outdir            Output directory.\n"
	    "\t-v --verbose           Verbose output.\n"
	    "\t-d --debug             Print debug info.\n"
	    "\t-h --help              Print help.\n";
};
void print_help_info()
{
    cout << "\ninfo:\n"
	    "\t-h --help              Print help.\n"
	    "\t[positional argument]  Input serialized batch (required).\n";
};

void print_version()
{
    std::cout << "isONclust2 version: ";
    std::cout << std::string(isONclust2_VERSION) << endl;
};

bool operator==(const CmdArgs& lhs, const CmdArgs& rgh)
{
    if (lhs.Verbose != rgh.Verbose) {
	return false;
    }
    if (lhs.Debug != rgh.Debug) {
	return false;
    }
    if (lhs.KmerSize != rgh.KmerSize) {
	return false;
    }
    if (lhs.BatchSize != rgh.BatchSize) {
	return false;
    }
    if (lhs.WindowSize != rgh.WindowSize) {
	return false;
    }
    if (lhs.MinShared != rgh.MinShared) {
	return false;
    }
    if (lhs.MinQual != rgh.MinQual) {
	return false;
    }
    if (lhs.MappedThreshold != rgh.MappedThreshold) {
	return false;
    }
    if (lhs.MinFraction != rgh.MinFraction) {
	return false;
    }

    return true;
}

bool operator!=(const CmdArgs& lhs, const CmdArgs& rgh)
{
    return !(lhs == rgh);
}
