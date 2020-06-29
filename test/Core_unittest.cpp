#include <gtest/gtest.h>
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "../Core.h"
#include "ChkQrmTest.h"

//Tests for function const bool chkQrm(UnitigColorMap<CoreInfo>&, const uint32_t&)//
//The maximum color id is (not) already smaller than the quorum DONE
//We are (not) dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors DONE
//A color is (not) found for the given unitig DONE
//The quorum is (not) reached DONE

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is already smaller than the quorum
TEST_F(ChkQrmTest, FewClrs){
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	EXPECT_FALSE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	UnitigColors::const_iterator i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	3. We are dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	4. A color is found for the given unitig
//	5. The quorum is reached
TEST_F(ChkQrmTest, IterAllClrs){
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	qrm = 1;

	EXPECT_TRUE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	UnitigColors::const_iterator i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests the function chkQrm under the following conditions
//	1. The maximum color id is not already smaller than the quorum
//	2. We are not dealing with the first color when the number of colors allowed to miss has to be decremented by the number of skipped colors
//	3. A color is (not) found for the given unitig
//	4. The quorum is not reached
TEST_F(ChkQrmTest, ClrNtFnd){
	cdbgOpt.filename_seq_in.push_back("Test_color1.fa");

	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);

	u = *cdbg.begin();
	u.len = 1;

	EXPECT_FALSE(chkQrm(u, qrm));
	EXPECT_EQ("AAGGCAAAC", u.mappedSequenceToString());
	UnitigColors::const_iterator i = u.getData()->getUnitigColors(u)->begin(u);
	EXPECT_EQ(0, (*i).first);
	EXPECT_EQ(0, (*i).second);
	++i;
	EXPECT_EQ(i, u.getData()->getUnitigColors(u)->end());
	EXPECT_TRUE(u.getData()->getData(u)->coreList.empty());
}

//Tests for function void markCore(ColoredCDBG<CoreInfo>&, const uint32_t&, const uint32_t&)//
//The quorum is (not) fulfilled for a k-mer 0/0
//A new interval has (not) to be started 0/0
//The bridging path's length has (not) to be reseted 0/0
//Delta is (not) exceeded 0/0

//Tests the function markCore under the following conditions
//	1.
TEST(markCoreTest, smplTst){

}