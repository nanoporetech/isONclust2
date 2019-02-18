#ifndef HPC_H_INCLUDED
#define HPC_H_INCLUDED

#include "seq.h"

Seq* HomopolymerCompress(const std::unique_ptr<Seq>& s);
Seq HomopolymerCompressObj(const Seq& s);

#endif
