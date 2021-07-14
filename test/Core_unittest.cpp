#include <gtest/gtest.h>
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "../Core.cpp"
#include "CoreTest.h"

//Tests for function const bool chkQrm(UnitigColorMap<CoreInfo>&, const uint32_t&)//
//	1. The maximum color id is (not) already smaller than the quorum DONE
//	2. We are (not) dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors DONE
//	3. A color is (not) found for the given unitig DONE
//	4. The quorum is (not) reached DONE

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is already smaller than the quorum
TEST_F(ChkQrmTest, FewClrs){
	cdbg = ColoredCDBG<CoreInfo>(DEFAULT_TEST_K, DEFAULT_TEST_G);
	cdbgOpt.k = DEFAULT_TEST_K;
	cdbgOpt.g = DEFAULT_TEST_G;
	cdbgOpt.filename_seq_in.push_back("Test.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	EXPECT_FALSE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	2. We are dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	3. A color is found for the given unitig
//	4. The quorum is reached
TEST_F(ChkQrmTest, IterAllClrs){
	cdbg = ColoredCDBG<CoreInfo>(DEFAULT_TEST_K, DEFAULT_TEST_G);
	cdbgOpt.k = DEFAULT_TEST_K;
	cdbgOpt.g = DEFAULT_TEST_G;
	cdbgOpt.filename_seq_in.push_back("Test.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	qrm = 1;

	EXPECT_TRUE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	2. We are not dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	3. A color is found for the given unitig
//	4. The quorum is not reached
TEST_F(ChkQrmTest, ClrNtFnd){
	cdbg = ColoredCDBG<CoreInfo>(DEFAULT_TEST_K, DEFAULT_TEST_G);
	cdbgOpt.k = DEFAULT_TEST_K;
	cdbgOpt.g = DEFAULT_TEST_G;
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbgOpt.filename_seq_in.push_back("Test_color1.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	EXPECT_FALSE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	2. We are dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	3. A color is found for the given unitig
//	4. The quorum is not reached
TEST_F(ChkQrmTest, NtFstKmer){
	cdbg = ColoredCDBG<CoreInfo>(DEFAULT_TEST_K, DEFAULT_TEST_G);
	cdbgOpt.k = DEFAULT_TEST_K;
	cdbgOpt.g = DEFAULT_TEST_G;
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbgOpt.filename_seq_in.push_back("Test_color1.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.dist = 2;
	u.len = 1;

	EXPECT_FALSE(chkQrm(u, qrm));
	EXPECT_EQ("GGCAAACAC", u.mappedSequenceToString());
	i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(2, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	2. We are (not) dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	3. A color is (not) found for the given unitig
//	4. The quorum is reached
TEST_F(ChkQrmTest, SkpdCol){	
	cdbg = ColoredCDBG<CoreInfo>(11, 8);
	cdbgOpt.k = 11;
	cdbgOpt.g = 8;
	cdbgOpt.filename_seq_in.push_back("Test21.fa");
	cdbgOpt.filename_seq_in.push_back("Test21_color1.fa");
	cdbgOpt.filename_seq_in.push_back("Test21_color2.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.dist = 1;
	u.len = 1;

	EXPECT_TRUE(chkQrm(u, qrm));
	EXPECT_EQ("TTAAAACCCAC", u.mappedSequenceToString());
	i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(1, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(1, (*i).first);
	EXPECT_EQ(2, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests for function void markCore(ColoredCDBG<CoreInfo>&, const uint32_t&, const uint32_t&)//
//	1. The quorum is (not) fulfilled for a k-mer DONE
//	2. A new interval has (not) to be started DONE
//	3. The bridging path's length has (not) to be reseted DONE
//	4. Delta is (not) exceeded DONE

//Tests the function markCore under the following conditions
//	1. The quorum is (not) fulfilled for a k-mer
//	2. A new interval has to be started
//	3. The bridging path's length has (not) to be reseted
//	4. Delta is (not) exceeded
TEST_F(MarkCoreTest, twoClrs){
	cdbgOpt.filename_seq_in.push_back("Test_color1.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	qrm = 2;

	markCore(cdbg, qrm, dlt);

	i = cdbg.begin();

	EXPECT_EQ("AAGGCAAACAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(1, inter->first);
	EXPECT_EQ(1, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAGGCAAAGAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ(i, cdbg.end());
}

//Tests the function markCore under the following conditions
//	1. The quorum is fulfilled for a k-mer
//	2. A new interval has (not) to be started
//	3. The bridging path's length has not to be reseted
TEST_F(MarkCoreTest, snglClr){
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	markCore(cdbg, qrm, dlt);

	i = cdbg.begin();

	EXPECT_EQ("AAGGCAAACAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(2, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAGGCAAAGAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(2, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ(i, cdbg.end());
}

//Tests the function markCore under the following conditions
//	1. The quorum is (not) fulfilled for a k-mer
//	2. A new interval has to be started
//	3. The bridging path's length has (not) to be reseted
//	4. Delta is exceeded
TEST_F(MarkCoreTest, TwoInts){
	cdbgOpt.filename_seq_in.push_back("Test_color8.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	qrm = 2;

	markCore(cdbg, qrm, dlt);

	i = cdbg.begin();

	EXPECT_EQ("AAGGCAAACAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(2, inter->first);
	EXPECT_EQ(2, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAGGCAAAGAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ("GCAAACACA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ("GCAAACACC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ(i, cdbg.end());
}

//Tests for function TravTrackQueue detectCore(ColoredCDBG<CoreInfo>&, const uint32_t&, const uint32_t&)//
//	1. The quorum is (not) fulfilled for a k-mer DONE
//	2. A new interval has (not) to be started 1/0
//	3. The bridging path's length has (not) to be reseted DONE
//	4. Delta is (not) exceeded DONE
//	5. No/An interval has been found DONE

//Tests the function detectCore under the following conditions
//	1. The quorum is (not) fulfilled for a k-mer
//	2. A new interval has to be started
//	3. The bridging path's length has (not) to be reseted
//	4. Delta is (not) exceeded
//	5. No/An interval has been found
TEST_F(DetectCoreTest, twoClrs){
	cdbgOpt.filename_seq_in.push_back("Test_color1.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	qrm = 2;

	queue = detectCore(cdbg, qrm, dlt);

	i = cdbg.begin();

	EXPECT_EQ("AAGGCAAACAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(1, inter->first);
	EXPECT_EQ(1, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAGGCAAAGAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(1, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_EQ(i, cdbg.end());
	ASSERT_EQ(queue.size(), 8);
	//1. Element (1, Kmer("GGCAAAGAC"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAAGAC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//2. Element (1, Kmer("AAAGGCAAA"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("AAAGGCAAA", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//3. Element (1, Kmer("GCAAACACA"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACA", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//4. Element (1, Kmer("GCAAACACA"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACA", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//5. Element (1, Kmer("AAAGGCAAA"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("AAAGGCAAA", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//6. Element (2, Kmer("GGCAAACAC"), false)
	EXPECT_EQ(2, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//7. Element (2, Kmer("GGCAAACAC"), true)
	EXPECT_EQ(2, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//8. Element (3, Kmer("GGCAAAGAC"), true)
	EXPECT_EQ(3, queue.top().cDist);
	EXPECT_EQ("GGCAAAGAC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
}

//Tests the function detectCore under the following conditions
//	1. The quorum is fulfilled for a k-mer
//	2. A new interval has (not) to be started
//	3. The bridging path's length has not to be reseted
//	5. An interval has been found
TEST_F(DetectCoreTest, snglClr){
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	queue = detectCore(cdbg, qrm, dlt);

	i = cdbg.begin();

	EXPECT_EQ("AAGGCAAACAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(2, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAGGCAAAGAC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(1, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(2, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(2, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACA", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ("GCAAACACC", i->mappedSequenceToString());
	col = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, (*col).first);
	EXPECT_EQ(0, (*col).second);
	++col;
	EXPECT_EQ(col, i->getData()->getUnitigColors(*i)->end());
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	inter = i->getData()->getData(*i)->coreList.begin();
	EXPECT_EQ(0, inter->first);
	EXPECT_EQ(0, inter->second);
	++inter;
	EXPECT_EQ(inter, i->getData()->getData(*i)->coreList.end());
	++i;
	EXPECT_EQ(i, cdbg.end());
	ASSERT_EQ(queue.size(), 10);

	
	//6. Element (1, Kmer("GGCAAACAC"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//7. Element (1, Kmer("GGCAAACAC"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//1. Element (1, Kmer("GGCAAAGAC"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAAGAC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//2. Element (1, Kmer("GGCAAAGAC"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAAGAC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//2. Element (1, Kmer("AAAGGCAAA"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("AAAGGCAAA", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//3. Element (1, Kmer("GCAAACACA"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACA", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//4. Element (1, Kmer("GCAAACACA"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACA", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//5. Element (1, Kmer("AAAGGCAAA"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("AAAGGCAAA", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//3. Element (1, Kmer("GCAAACACC"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
	//4. Element (1, Kmer("GCAAACACC"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GCAAACACC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	exit(0);
}