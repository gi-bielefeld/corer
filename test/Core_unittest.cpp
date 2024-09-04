#include <gtest/gtest.h>
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "../src/Core.cpp"
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
//	2. A new interval has (not) to be started DONE
//	3. The bridging path's length has (not) to be reseted DONE
//	4. Delta is (not) exceeded DONE
//	5. No/An interval has been found DONE
//	6. A unitig has (a/no) successor(s) DONE
//	7. A unitig has (a/no) predecessor(s) DONE

//Tests the function detectCore under the following conditions
//	1. The quorum is (not) fulfilled for a k-mer
//	2. A new interval has to be started
//	3. The bridging path's length has (not) to be reseted
//	4. Delta is not exceeded
//	5. No/An interval has been found
//	6. A unitig has (a/no) successor(s)
//	7. A unitig has (a/no) predecessor(s)
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
	ASSERT_EQ(queue.size(), 5);

	for(c = 3; c > 0; --c){
		EXPECT_EQ(1, queue.top().cDist);

		//Element (1, Kmer("GGCAAAGAC"), false)
		if(!queue.top().track.toString().compare("GGCAAAGAC")){
			EXPECT_FALSE(queue.top().isSucTrav);
			EXPECT_FALSE(u2SeenPred);
			u2SeenPred = true;
		//Element (1, Kmer("AAAGGCAAA"), true)
		} else if(!queue.top().track.toString().compare("AAAGGCAAA") && queue.top().isSucTrav){
			if(!u3SeenSuc){
				u3SeenSuc = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u3SeenSuc);
			}
		//Element (1, Kmer("GCAAACACA"), false)
		} else if(!queue.top().track.toString().compare("GCAAACACA") && !queue.top().isSucTrav){
			if(!u4SeenPred){
				u4SeenPred = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u4SeenPred);
			}
		} else{
			//This should not happen...
			EXPECT_EQ("", queue.top().track.toString());
		}

		queue.pop();
	}

	for(c = 2; c > 0; --c){
		EXPECT_EQ(2, queue.top().cDist);
		EXPECT_EQ("GGCAAACAC", queue.top().track.toString());

		if(queue.top().isSucTrav){
			//Element (2, Kmer("GGCAAACAC"), true)
			EXPECT_FALSE(u1SeenSuc);
			u1SeenSuc = true;
		} else{
			//Element (2, Kmer("GGCAAACAC"), false)
			EXPECT_FALSE(u1SeenPred);
			u1SeenPred = true;
		}

		queue.pop();
	}
}

//Tests the function detectCore under the following conditions
//	1. The quorum is fulfilled for a k-mer
//	2. A new interval has (not) to be started
//	3. The bridging path's length has not to be reseted
//	5. An interval has been found
//	6. A unitig has (a/no) successor(s)
//	7. A unitig has (a/no) predecessor(s)
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
	ASSERT_EQ(queue.size(), 6);

	for(c = 6; c > 0; --c){
		EXPECT_EQ(1, queue.top().cDist);

		//Element (1, Kmer("GGCAAACAC"), true)
		if(!queue.top().track.toString().compare("GGCAAACAC") && queue.top().isSucTrav){
			if(!u1SeenSuc){
				u1SeenSuc = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u1SeenSuc);
			}
		//Element (1, Kmer("GGCAAACAC"), false)
		} else if(!queue.top().track.toString().compare("GGCAAACAC") && !queue.top().isSucTrav){
			if(!u1SeenPred){
				u1SeenPred = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u1SeenPred);
			}
		//Element (1, Kmer("GGCAAAGAC"), false)
		} else if(!queue.top().track.toString().compare("GGCAAAGAC") && !queue.top().isSucTrav){
			if(!u2SeenPred){
				u2SeenPred = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u2SeenPred);
			}
		//Element (1, Kmer("AAAGGCAAA"), true)
		} else if(!queue.top().track.toString().compare("AAAGGCAAA") && queue.top().isSucTrav){
			if(!u3SeenSuc){
				u3SeenSuc = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u3SeenSuc);
			}
		//Element (1, Kmer("GCAAACACA"), false)
		} else if(!queue.top().track.toString().compare("GCAAACACA") && !queue.top().isSucTrav){
			if(!u4SeenPred){
				u4SeenPred = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u4SeenPred);
			}
		//Element (1, Kmer("GCAAACACC"), false)
		} else if(!queue.top().track.toString().compare("GCAAACACC") && !queue.top().isSucTrav){
			if(!u5SeenPred){
				u5SeenPred = true;
			} else{
				//This should not happen...
				EXPECT_FALSE(u5SeenPred);
			}
		} else{
			//This should not happen...
			EXPECT_EQ("", queue.top().track.toString());
		}

		queue.pop();
	}
}

//Tests the function detectCore under the following conditions
//	1. The quorum is (not) fulfilled for a k-mer
//	2. A new interval has to be started
//	3. The bridging path's length has not to be reseted
//	4. Delta is exceeded
//	5. No/An interval has been found
//	6. A unitig has (a) successor(s)
//	7. A unitig has (a) predecessor(s)
TEST_F(DetectCoreTest, TwoInts){
	cdbgOpt.filename_seq_in.push_back("Test_color8.fa");

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
	ASSERT_EQ(queue.size(), 2);
	//1. Element (1, Kmer("GGCAAACAC"), true)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_TRUE(queue.top().isSucTrav);
	queue.pop();
	//2. Element (1, Kmer("GGCAAACAC"), false)
	EXPECT_EQ(1, queue.top().cDist);
	EXPECT_EQ("GGCAAACAC", queue.top().track.toString());
	EXPECT_FALSE(queue.top().isSucTrav);
	queue.pop();
}

//Tests for function void markKmers(ColoredCDBG<CoreInfo>&, vector<string>&, const uint32_t&)//
//	1. A sequence is (not) shorter than k DONE
//	2. A sequence contains only one/multiple k-mers DONE
//	3. A k-mer can(not) be found in the graph DONE
//	4. A k-mer is (not) found on the reference strand DONE
//	5. The first not matched character of the sequence can(not) be skipped DONE

//Tests the function markKmers under the following conditions
//	1. A sequence is (not) shorter than k
//	2. A sequence contains only one/multiple k-mers
//	3. A k-mer can(not) be found in the graph
//	4. A k-mer is (not) found on the reference strand
TEST_F(MarkKmersTest, UltTst){
	sl = {"", "ACTTTTCAA", "AAGGCAAACAC", "GTCTTTGCCT"};
	markKmers(cdbg, sl, false, 42);

	i = cdbg.begin();
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 2);
	++i;
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 2);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
}

//Tests the function markKmers under the following conditions
//	1. A sequence is not shorter than k
//	2. A sequence contains only multiple k-mers
//	3. A k-mer can(not) be found in the graph
//	4. A k-mer is found on the reference strand
//	5. The first not matched character of the sequence can(not) be skipped
TEST_F(MarkKmersTest, SkpNxt){
	sl = {"AAGGCAAAGGCAAA"};
	markKmers(cdbg, sl, false, 42);

	i = cdbg.begin();
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	++i;
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->coreList.empty());
}