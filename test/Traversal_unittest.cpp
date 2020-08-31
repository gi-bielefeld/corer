#include <gtest/gtest.h>
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "../Traversal.h"
#include "TraversalTest.h"

//Tests for function const bool prioShrtst(const Path&, const Path&)//
//The first path is (not) longer than the second one DONE
//The second path is (not) longer than the first one DONE

//Tests the function prioShrtst under the following conditions
//	1. The first path is longer than the second
//	2. The second path is not longer than the first
TEST_F(PrioShrtstTest, FstLng){ EXPECT_TRUE(prioShrtst(p, q)); }

//Tests the function prioShrtst under the following conditions
//	1. The first path is not longer than the second
//	2. The second path is longer than the first
TEST_F(PrioShrtstTest, SndLng){ EXPECT_FALSE(prioShrtst(q, p)); }

//Tests for function const uint32_t getCoreDist(const neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false>&, const bool&)//
//The function is (not) called in the context of a successive traversal DONE
//The function is called for a unitig's reference/reverse complementary strand DONE

//Tests the function getCoreDist under the following conditions
//	1. The function is called in the context of a successive traversal
//	2. The function is called for a unitig's reference strand
TEST_F(GetCoreDistTest, SucRef){
	uni = *cdbg.begin()->getSuccessors().begin();
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));

	EXPECT_EQ(1 ,getCoreDist(cdbg.begin()->getSuccessors().begin(), true));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is called in the context of a successive traversal
//	2. The function is called for a unitig's reverse complementary strand
TEST_F(GetCoreDistTest, SucRev){
	i = cdbg.begin();
	++i;
	++i;
	uni = *i;
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	uni = *cdbg.begin();
	uni.strand = false;

	EXPECT_EQ(1, getCoreDist(uni.getSuccessors().begin(), true));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is not called in the context of a successive traversal
//	2. The function is called for a unitig's reference strand
TEST_F(GetCoreDistTest, PredRef){
	i = cdbg.begin();
	++i;
	++i;
	uni = *i;
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));

	EXPECT_EQ(1, getCoreDist(cdbg.begin()->getPredecessors().begin(), true));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is not called in the context of a successive traversal
//	2. The function is called for a unitig's reverse complementary strand
TEST_F(GetCoreDistTest, PredRev){
	uni = *cdbg.begin()->getSuccessors().begin();
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	uni = *cdbg.begin();
	uni.strand = false;

	EXPECT_EQ(1, getCoreDist(uni.getPredecessors().begin(), false));
}