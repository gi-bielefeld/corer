#include <gtest/gtest.h>

#include "BridgingTest.h"
#include "../Bridging.cpp"

//Tests for function const bool lCrTooFar(const size_t&, const list<pair<uint32_t, uint32_t>>&, const uint32_t&)//
//The core list is empty and the unitig is (not) too long DONE
//The core list is not empty and the rightmost core k-mer (i.e. the k-mer with the highest start offset) is (not) too far DONE

//Tests the function lCrTooFar under the following conditions
//	1.The core list is empty
//	2.The unitig is too long
TEST_F(CrTooFarTest, LUni){
	d = 41;

	EXPECT_TRUE(lCrTooFar(len, ints, d));
}

//Tests the function lCrTooFar under the following conditions
//	1.The core list is empty
//	2.The unitig is not too long
TEST_F(CrTooFarTest, SUni){
	EXPECT_FALSE(lCrTooFar(len, ints, d));
}

//Tests the function lCrTooFar under the following conditions
//	1.The core list is not empty
//	2.The rightmost core k-mer is too far
TEST_F(CrTooFarTest, FarKmer){
	d = 41;
	ints.push_back(pair<uint32_t, uint32_t>(0, 0));

	EXPECT_TRUE(lCrTooFar(len, ints, d));
}

//Tests the function lCrTooFar under the following conditions
//	1.The core list is not empty
//	2.The rightmost core k-mer is not too far
TEST_F(CrTooFarTest, NearKmer){
	ints.push_back(pair<uint32_t, uint32_t>(0, 0));

	EXPECT_FALSE(lCrTooFar(len, ints, d));
}

//Tests for function const bool rCrTooFar(const size_t&, const list<pair<uint32_t, uint32_t>>&, const uint32_t&)
//The core list is empty and the unitig is (not) too long DONE
//The core list is not empty and the leftmost core k-mer (i.e. the k-mer with lowest start offset) is (not) too far DONE

//Tests the function rCrTooFar under the following conditions
//	1.The core list is empty
//	2.The unitig is too long
TEST_F(CrTooFarTest, LongU){
	d = 41;

	EXPECT_TRUE(rCrTooFar(len, ints, d));
}

//Tests the function rCrTooFar under the following conditions
//	1.The core list is empty
//	2.The unitig is not too long
TEST_F(CrTooFarTest, ShrtU){
	EXPECT_FALSE(rCrTooFar(len, ints, d));
}

//Tests the function rCrTooFar under the following conditions
//	1.The core list is not empty
//	2.The leftmost core k-mer is too far
TEST_F(CrTooFarTest, FKmer){
	d = 41;
	ints.push_back(pair<uint32_t, uint32_t>(41, 41));

	EXPECT_TRUE(rCrTooFar(len, ints, d));
}

//Tests the function rCrTooFar under the following conditions
//	1.The core list is not empty
//	2.The leftmost core k-mer is not too far
TEST_F(CrTooFarTest, NKmer){
	ints.push_back(pair<uint32_t, uint32_t>(41, 41));

	EXPECT_FALSE(rCrTooFar(len, ints, d));
}