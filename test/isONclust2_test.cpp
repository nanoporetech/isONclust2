#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "cluster.h"
#include "gtest/gtest.h"
#include "hpc.h"
#include "kmer_index.h"
#include "minimizer.h"
#include "output.h"
#include "p_emp_prob.h"
#include "parasail.h"
#include "qualscore.h"
#include "seq.h"

// Test sequence sorting.
TEST(SortingTest, SortingTest)
{
    std::vector<std::unique_ptr<Seq>> input;

    std::unique_ptr<Seq> s0(
	new Seq("s0", "ATGCGCATATGCGC", "@IIIIIIIIIIIII", 0.0));
    input.push_back(std::move(s0));

    std::unique_ptr<Seq> s1(
	new Seq("s1", "ATGCTGACATGCATGC", "@IIIIIIIIIIIIIII", 0.0));
    input.push_back(std::move(s1));

    std::unique_ptr<Seq> s2(new Seq("s2", "ATGCATGCCGATGTACATGCATGCATCGACGT",
				    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 0.0));
    input.push_back(std::move(s2));

    auto qualTab = InitQualTab();
    FillQualScores(input, 3, 5, qualTab, InitQualTabNomin());
    SortByQualScores(input);

    std::vector<std::string> valid_order{"s2", "s1", "s0"};

    for (unsigned i = 0; i < valid_order.size(); i++) {
	// std::cout << input[i]->Score() << std::endl;
	EXPECT_EQ(valid_order[i], input[i]->Name());
    }
}

// Test minimizer calculation.
TEST(MinimizerTest, MinimizerTest)
{
    auto min =
	GetKmerMinimizers(KmerEncodeSeq(std::string("ACGCCGATC"), 2), 2, 4);
    std::string ac = "AC";
    std::string cc = "CC";
    std::string at = "AT";
    Minimizers valid_result = {Minimizer{KmerToIndex(ac, ac.end()), 0, 0},
			       Minimizer{KmerToIndex(cc, cc.end()), 3, 1},
			       Minimizer{KmerToIndex(at, at.end()), 6, 2}};
    EXPECT_EQ(min, valid_result);
}

// Test homopolymer compression.
TEST(HpcTest, HpcTest)
{
    std::unique_ptr<Seq> stmp(
	new Seq("Foo", "AAATTTGCGTTAA", "++:+?++++++@+", 0.0));
    auto so = HomopolymerCompress(stmp);
    EXPECT_EQ(so->Str(), "ATGCGTA");
    EXPECT_EQ(so->Qual(), ":?++++@");
}

// Test error rate calculation.
TEST(ErrorRateTest, ErrorRateTest)
{
    auto res = CalcErrorRate("IIIIIIIIIIIIIIIIII", InitQualTab());
    EXPECT_DOUBLE_EQ(res, 0.0001);
}

// Test empirical probability lookup.
TEST(EmpProbLookupTest, EmpProbLookupTest)
{
    auto msm = InitMinSharedMap(13, 20);
    auto res = GetPMinShared(0.111, 0.131, msm);
    EXPECT_DOUBLE_EQ(res, 0.11736487779693013);
}

// Test minimizer matching.
TEST(MinMatchTest, MinMatchTest)
{
    std::string ref =
	"GGTAGTGGTGGCGGGTCTCCTTGAGAGCACTCGTCGAGTATGCCGAAAATATGTTAATGG"
	"CAGGAAGTTTGATTATAGCCATTAGCGTGTCATAATGTAGAAAGTCTCGATAATAAAGCT"
	"CAGGACGCGCCTCCGTTAAAGGAAGGCGGGATCCTGCGCGATGGCTATCTATAGTATGTA"
	"GTTACCTCTGATTGTCATGTGAACAGGAGGCCAGTACCACCTGATACGGCCTTGTAAACC"
	"TACCACTACTTCGCTTAAGACGGTGCTCCCCTCCCCATTTGCGGCCGTTCGTCGTGTCCC";

    std::unique_ptr<Seq> refSeq(
	new Seq("ref", ref, std::string(ref.length(), 'I'), 0.0));
    auto refHpc = HomopolymerCompress(refSeq);

    std::string read =
	"AGATATTATAGCCATACGTGTCATAATGTAGAAGTCTCGATAATAAAGCTCAGGACGCGC"
	"CTCCGTTAAGGAAGGCGGATCCGCGCGATGGGCTATCTATAGTATGTGGTTACCCTGATA"
	"GTCATGTGAGACAGGAGGCCAGTCCACCTGATACGGCTTGTAAACTACCACTACTTCGCT";

    std::unique_ptr<Seq> readSeq(
	new Seq("ref", read, std::string(read.length(), 'I'), 0.0));
    auto readHpc = HomopolymerCompress(readSeq);

    int kmerSize = 13;
    int windowSize = 20;

    MinimizerDB minDB;
    auto refMins = GetKmerMinimizers(KmerEncodeSeq(refHpc->Str(), kmerSize),
				     kmerSize, windowSize);
    AddMinimizers(refMins, 1, minDB);

    auto readMins = GetKmerMinimizers(KmerEncodeSeq(readHpc->Str(), kmerSize),
				      kmerSize, windowSize);
    auto hits = GetMinimizerHits(readMins, Minimizers(), minDB);
    Clusters cls;
    auto hitOrder = SortMinimizerHits(hits, cls);

    EXPECT_EQ(std::get<0>(hitOrder[0]), 14);

    auto msm = InitMinSharedMap(kmerSize, windowSize);
    auto qt = InitQualTab();
    readHpc->SetErrorRate(CalcErrorRate(readHpc->Qual(), qt));
    refHpc->SetErrorRate(CalcErrorRate(refHpc->Qual(), qt));

    double pError =
	1.0 - GetPMinShared(refHpc->ErrorRate(), readHpc->ErrorRate(), msm);
    EXPECT_DOUBLE_EQ(pError, 0.17140336964776648);

    auto mr = getMappedRatio(*readHpc, *refHpc, readMins,
			     hits[std::make_pair(1, 1)], msm, 0.1);
    EXPECT_DOUBLE_EQ(mr, 0.3835616438356164);
}

TEST(AlnRatioTest, AlnRatioTest)
{
    std::string ref =
	"GGTAGTGGTGGCGGGTCTCCTTGAGAGCACTCGTCGAGTATGCCGAAAATATGTTAATGG"
	"CAGGAAGTTTGATTATAGCCATTAGCGTGTCATAATGTAGAAAGTCTCGATAATAAAGCT"
	"CAGGACGCGCCTCCGTTAAAGGAAGGCGGGATCCTGCGCGATGGCTATCTATAGTATGTA"
	"GTTACCTCTGATTGTCATGTGAACAGGAGGCCAGTACCACCTGATACGGCCTTGTAAACC"
	"TACCACTACTTCGCTTAAGACGGTGCTCCCCTCCCCATTTGCGGCCGTTCGTCGTGTCCC";

    std::unique_ptr<Seq> refSeq(
	new Seq("ref", ref, std::string(ref.length(), 'I'), 0.0));

    std::string read =
	"AGATATTATAGCCATACGTGTCATAATGTAGAAGTCTCGATAATAAAGCTCAGGACGCGC"
	"CTCCGTTAAGGAAGGCGGATCCGCGCGATGGGCTATCTATAGTATGTGGTTACCCTGATA"
	"GTCATGTGAGACAGGAGGCCAGTCCACCTGATACGGCTTGTAAACTACCACTACTTCGCT";

    std::unique_ptr<Seq> readSeq(
	new Seq("ref", read, std::string(read.length(), 'I'), 0.0));

    int kmerSize = 13;
    int match = 2;
    int mismatch = -2;
    int gapExtend = 1;
    auto user_matrix = parasail_matrix_create("ACGT", match, mismatch);

    auto qualTab = InitQualTab();
    auto e1 = CalcErrorRate(refSeq->Qual(), qualTab);
    auto e2 = CalcErrorRate(readSeq->Qual(), qualTab);
    int gapOpen = setGapOpen(e1 + e2);
    auto alnRes = ParasailAlign(refSeq->Str(), readSeq->Str(), gapOpen,
				gapExtend, user_matrix);
    auto tr = parasail_result_get_traceback(
	alnRes, refSeq->Str().c_str(), refSeq->Str().length(),
	readSeq->Str().c_str(), readSeq->Str().length(), user_matrix, '|', ' ',
	' ');
    parasail_result_free(alnRes);

    std::string Comp(tr->comp);
    auto alnRatio =
	getAlnRatio(Comp, e1 + e2, readSeq->Str().length(), kmerSize);
    parasail_traceback_free(tr);
    parasail_matrix_free(user_matrix);
    EXPECT_DOUBLE_EQ(alnRatio, 0.7111111111111111);
}

// Test kmer transformation.
TEST(TestKmerTransform, TestKmerTransform)
{
    unsigned k = 4;
    std::vector<std::string> res;
    std::vector<unsigned> indices;
    for (unsigned i = 0; i < (unsigned)std::pow(4, k); i++) {
	auto kmer = IndexToKmer(i, k);
	res.push_back(kmer);
	indices.push_back(i);
    }
    std::vector<std::string> sorted_res(res);
    std::sort(sorted_res.begin(), sorted_res.end());
    EXPECT_EQ(res, sorted_res);

    std::vector<unsigned> res_indices;
    for (auto& kmer : res) {
	res_indices.push_back(KmerToIndex(kmer, kmer.end()));
    }
    EXPECT_EQ(indices, res_indices);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
