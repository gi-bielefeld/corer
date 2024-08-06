#include <gtest/gtest.h>

#include "CoreInfoTest.h"
#include "CoreInfo.h"

//Tests for function void CoreInfo::updateCoreList(uint32_t, uint32_t, const uint32_t&)//
//	1. The core list is (not) empty DONE
//	2. The new interval's end is (not) too small compared to the current interval's start DONE
//	3. The interval is (not) too far behind the current one in the list DONE

//Tests the function updateCoreList under the following conditions
//	1. The core list is (not) empty
//	2. The new interval's end is (not) too small compared to the current interval's start
//	3. The interval is not too far behind the current one in the list
TEST_F(UpdateCoreListTest, EmpLst){
	i->getData()->getData(*i)->coreList.push_back(make_pair(2, 2));
	i->getData()->getData(*i)->updateCoreList(0, 0, 0);
	++i;
	i->getData()->getData(*i)->coreList.push_back(make_pair(0, 0));
	i->getData()->getData(*i)->coreList.push_back(make_pair(2, 2));
	i->getData()->getData(*i)->updateCoreList(1, 1, 0);
	++i;
	i->getData()->getData(*i)->updateCoreList(0, 0, 0);

	i = cdbg.begin();
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 2);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().first, 2);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().second, 2);
	++i;
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 2);
	++i;
	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 1);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
}

//Tests the function updateCoreList under the following conditions
//	1. The core list is not empty
//	2. The new interval's end is not too small compared to the current interval's start
//	3. The interval is too far behind the current one in the list
TEST_F(UpdateCoreListTest, FrBhnd){
	i->getData()->getData(*i)->coreList.push_back(make_pair(0, 0));
	i->getData()->getData(*i)->updateCoreList(2, 2, 0);

	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 2);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().first, 2);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().second, 2);
}

//Tests the function updateCoreList under the following conditions
//	1. The core list is not empty
//	2. The new interval's end is (not) too small compared to the current interval's start
//	3. The interval is (not) too far behind the current one in the list
TEST_F(UpdateCoreListTest, MultInt1){
	i->getData()->getData(*i)->coreList.push_back(make_pair(0, 0));
	i->getData()->getData(*i)->coreList.push_back(make_pair(3, 3));
	i->getData()->getData(*i)->coreList.push_back(make_pair(6, 7));
	i->getData()->getData(*i)->coreList.push_back(make_pair(10, 10));
	i->getData()->getData(*i)->coreList.push_back(make_pair(13, 14));
	i->getData()->getData(*i)->updateCoreList(5, 5, 1);
	i->getData()->getData(*i)->updateCoreList(12, 14, 1);

	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 3);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	in = ++(i->getData()->getData(*i)->coreList.begin());
	EXPECT_EQ(in->first, 3);
	EXPECT_EQ(in->second, 7);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().first, 10);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().second, 14);
}

//Tests the function updateCoreList under the following conditions
//	1. The core list is not empty
//	2. The new interval's end is (not) too small compared to the current interval's start
//	3. The interval is (not) too far behind the current one in the list
TEST_F(UpdateCoreListTest, MultInt2){
	i->getData()->getData(*i)->coreList.push_back(make_pair(0, 0));
	i->getData()->getData(*i)->coreList.push_back(make_pair(3, 3));
	i->getData()->getData(*i)->coreList.push_back(make_pair(6, 7));
	i->getData()->getData(*i)->coreList.push_back(make_pair(10, 10));
	i->getData()->getData(*i)->coreList.push_back(make_pair(13, 14));
	i->getData()->getData(*i)->updateCoreList(9, 9, 1);

	ASSERT_FALSE(i->getData()->getData(*i)->coreList.empty());
	EXPECT_EQ(i->getData()->getData(*i)->coreList.size(), 4);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().first, 0);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.front().second, 0);
	in = ++(i->getData()->getData(*i)->coreList.begin());
	EXPECT_EQ(in->first, 3);
	EXPECT_EQ(in->second, 3);
	++in;
	EXPECT_EQ(in->first, 6);
	EXPECT_EQ(in->second, 10);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().first, 13);
	EXPECT_EQ(i->getData()->getData(*i)->coreList.back().second, 14);
}