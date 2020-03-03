#include "serialize.h"
#include "cluster.h"

extern UnsignedHash uh;
void SaveBatch(const std::unique_ptr<Batch>& b, std::string outf)
{
    std::ofstream os(outf, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(*b);
}

BatchP LoadBatch(std::string inf)
{
    std::ifstream instream(inf, std::ios::binary);

    auto p = std::unique_ptr<Batch>(new Batch);
    cereal::BinaryInputArchive iarchive(instream);
    try {
	iarchive(*p);
    }
    catch (std::runtime_error e) {
	std::cerr << "Failed to load batch " << inf << ":" << e.what()
		  << std::endl;
	exit(1);
    }
    return p;
}

BatchP CreatePseudoBatch(std::unique_ptr<Batch>& inBatch)
{
    auto nb = std::unique_ptr<Batch>(new Batch);
    nb->BatchNr = -inBatch->BatchNr;
    nb->BatchStart = inBatch->BatchStart;
    nb->BatchEnd = inBatch->BatchEnd;
    nb->BatchBases = 0;
    nb->TotalReads = 0;
    nb->SortArgs = inBatch->SortArgs;
    nb->Depth = -1;
    nb->Cls = Clusters(inBatch->Cls);
    nb->NrCls = int(nb->Cls.size());

    return nb;
}
