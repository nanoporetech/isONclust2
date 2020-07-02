![ONT_logo](/ONT_logo.png)

-----------------------------

isONclust2 - a tool for de novo clustering of long transcriptomic reads
=======================================================================

`isONclust2` is a tool for clustering long transcriptomic reads in to gene families.
The tool is based on the approach pioneered by [isONclust](https://github.com/ksahlin/isONclust), using minimizers and occasional pairwise alignment.

`isONclust2` is implemented in C++, which makes it fast enough to cluster large transcriptomic datasets produced on PromethION sequencers. The tool is not a mere re-implementation of the original `isONclust` approach, as it deals with the strandedness of the reads and provides further optional features. 

*WARNING: In order to be able to handle large datasets, `isONclust2` splits the input data into batches which have to be processed in a specified order to obtain the results. Hence the use of `isONclust2` as a standalone tool is highly discouraged and one should always use it through the de novo transcriptomics pipeline at [https://github.com/nanoporetech/pipeline-nanopore-denovo-isoforms](https://github.com/nanoporetech/pipeline-nanopore-denovo-isoforms).*

Getting Started
===============

## Installation

The best way to install `isONclust2` is from bioconda:

- Make sure you have [miniconda3](https://docs.conda.io/en/latest/miniconda.html) installed.
- Install the tool by issuing `conda install -c bioconda isonclust2`

## Compiling from source

- Clone the repository: `git clone --recursive https://github.com/nanoporetech/isONclust2.git`
- Make sure you have cmake v3.1 or later.
- Issue `cd isONclust2; mkdir build; cd build; cmake ..; make -j`
- The produced binary is static with no library dependencies. Link this under your path to use the tool.

## Usage

```
isONclust2 version: v2.3-a0e5b32
Available subcommands: sort, cluster, dump, info, help, version

sort - sort reads and write out batches:
        -B --batch-size        Batch size in kilobases (default: 50000)
        -M --batch-max-seq     Maximum number of sequences per batch (default: 3000).
        -k --kmer-size         Kmer size (default: 11).
        -w --window-size       Window size (default: 15).
        -m --min-shared        Minimum number of minimizers shared between read and cluster (default: 5).
        -q --min-qual          Minimum average quality value (default: 7.0).
        -x --mode  Clustering mode:
                   * sahlin (default): use minimizers first, alignment second
                   * fast: use minimizers only
                   * furious: always use alignment
        -g --low-cons-size     Use all sequences for consensus below this size (default: 20).
        -c --max-cons-size     Maximum number of sequences used for consensus (default: 150).
        -P --cons-period       Do not recalculate consensus after this many seuqences added (default: 500).
        -r --mapped-threshold  Minmum mapped fraction of read to be     included in cluster (default: 0.65).
        -a --aligned-threshold Minimum aligned fraction of read to be included in cluster (default: 0.2).
        -f --min-fraction      Minimum fraction of minimizers shared compared to best hit, in order to continue mapping (default: 0.8).
        -p --min-prob-no-hits  Minimum probability for i consecutive    minimizers to be different between read and representative (default: 0.1)
        -F --min-cls-size      Skip clusters smaller than this in the left batch (default: 3).
        -o --outfolder         Output folder (default:  ./isONclust2_batches).
        -h --help              Print help.
        -v --verbose           Verbose output.
        -d --debug             Print debug info.
        [positional argument]  Input fastq file (required).

cluster - cluster and/or merge batches:
        -l --left-batch        Left input batch (mandatory).
        -r --right-batch       Right input batch (optional).
        -o --outfile           Output batch.
        -x --mode  Clustering mode:
                   * sahlin (default): use minimizers first, alignment second
                   * fast: use minimizers only
                   * furious: use alignment only
        -A --spoa-algo  spoa alignment algorithm:
                   * 0 (default): local
                   * 1 : global
                   * 1 : semi-global
        -z --min-purge         Purge minimizer database from output batch.
        -j --keep-seq          Do not purge non-representative sequences from output batches.
        -F --min-cls-size      Skip clusters smaller than this in the left batch.
        -v --verbose           Verbose output.
        -Q --quiet             Supress progress bar.
        -d --debug             Print debug info.
        -h --help              Print help.

dump - dump clustered batch:
        -o --outdir            Output directory.
        -i --index             Index of sorted reads.
        -v --verbose           Verbose output.
        -d --debug             Print debug info.
        -h --help              Print help.

info:
        -h --help              Print help.
        [positional argument]  Input serialized batch (required).

help - print help message

version - print version
```

Help
====

## Licence and Copyright

(c) 2020 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips


## References and Supporting Information

See the post announcing the transcriptomics tools at the Oxford Nanopore Technologies community [here](https://community.nanoporetech.com/posts/new-transcriptomics-analys).

## Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.

