
#include "seq.h"

Seq::Seq(const char* name, uint32_t name_size, const char* data,
	 uint32_t data_size, const char* quality, uint32_t quality_size)
    : name(name, name_size), seq(data, data_size), qual(quality, quality_size)
{
}

Seq::Seq(const std::string& name, const std::string& seq,
	 const std::string& qual, double score)
    : name(name), seq(seq), qual(qual), score(score)
{
}
