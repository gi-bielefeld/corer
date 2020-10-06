#include <gtest/gtest.h>

#include "BridgingTest.h"
#include "../Bridging.cpp"

//Tests for function const bool lCrTooFar(const size_t&, const list<pair<uint32_t, uint32_t>>&, const uint32_t&)//
//	1. The core list is empty and the unitig is (not) too long DONE
//	2. The core list is not empty and the rightmost core k-mer (i.e. the k-mer with the highest start offset) is (not) too far DONE

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
//	1. The core list is empty and the unitig is (not) too long DONE
//	2. The core list is not empty and the leftmost core k-mer (i.e. the k-mer with lowest start offset) is (not) too far DONE

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

//Tests for function void markBrdg(const list<Path>&, const bool&)//
//	1. Our list of paths is (not) empty 1/0
//	2. Our list consists of one/many path(s) 0/0
//	3. We are (not) dealing with a list of successive paths 1/0
//	4. We are dealing with a list of successive paths and are (not) on the unitigs' reference strand 0/0
//	5. We are dealing with a list of predecessive paths and are (not) on the unitigs' reference strand 0/0

//Tests the function markBrdg under the following conditions
//	1. Our list of paths is empty
//	3. We are dealing with a list of successive paths
TEST_F(MarkBrdgTest, NoPths){
	markBrdg(l, true);
	EXPECT_TRUE(l.empty());
}

//Tests for function markBrdg under the following conditions
//	1. Our list of paths is not empty
//	2. Our list consists of one path
//	3. We are not dealing with a list of successive paths
//	5. We are dealing with a list of predecessive paths and are on the unitigs' reference strand
TEST_F(MarkBrdgTest, SnglPth){
	
	CCDBG_Build_opt cdbgOpt;
	cdbgOpt.k = DEFAULT_TEST_K;
	cdbgOpt.g = DEFAULT_TEST_G;
	cdbgOpt.filename_seq_in.push_back("Test.fa");
}