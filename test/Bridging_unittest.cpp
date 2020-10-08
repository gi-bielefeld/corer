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
//	1. Our list of paths is (not) empty DONE
//	2. Our list consists of one/many path(s) DONE
//	3. We are (not) dealing with a list of successive paths DONE
//	4. We are dealing with a list of successive paths and are (not) on the unitigs' reference strand DONE
//	5. We are dealing with a list of predecessive paths and are (not) on the unitigs' reference strand DONE

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
	l.push_back(Path(4, list<UnitigColorMap<CoreInfo>>()));
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->predCoreDist = 4;
	l.front().second.push_back(*i);
	++i;
	++i;
	i->getData()->getData(*i)->predCoreDist = 1;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	l.front().second.push_back(*i);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	l.front().second.push_front(*i);

	markBrdg(l, false);
	EXPECT_EQ(1, l.size());

	for(list<UnitigColorMap<CoreInfo>>::const_iterator u = l.front().second.begin(); u != l.front().second.end(); ++u) EXPECT_TRUE(u->getData()->getData(*u)->preBrdg);
}

//Tests for function markBrdg under the following conditions
//	1. Our list of paths is not empty
//	2. Our list consists of many paths
//	3. We are dealing with a list of successive paths
//	4. We are dealing with a list of successive paths and are on the unitigs' reference strand
TEST_F(MarkBrdgTest, MnyPths){
	l.push_back(Path(4, list<UnitigColorMap<CoreInfo>>()));
	l.push_back(Path(4, list<UnitigColorMap<CoreInfo>>()));
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,2));
	i->getData()->getData(*i)->sucCoreDist = 2;
	l.front().second.push_back(*i);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,2));
	i->getData()->getData(*i)->sucCoreDist = 2;
	l.back().second.push_back(*i);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	l.front().second.push_front(*i);
	l.back().second.push_front(*i);

	markBrdg(l, true);
	EXPECT_EQ(2, l.size());
	EXPECT_EQ(2, l.front().second.size());
	EXPECT_EQ(2, l.back().second.size());

	for(list<Path>::const_iterator p = l.begin(); p != l.end(); ++p){
		for(list<UnitigColorMap<CoreInfo>>::const_iterator u = l.front().second.begin(); u != l.front().second.end(); ++u) EXPECT_TRUE(u->getData()->getData(*u)->sufBrdg);
	}
}

//Tests for function markBrdg under the following conditions
//	1. Our list of paths is not empty
//	2. Our list consists of one path
//	3. We are dealing with a list of successive paths
//	4. We are dealing with a list of successive paths and are (not) on the unitigs' reference strand
TEST_F(MarkBrdgTest, SucRev){
	l.push_back(Path(4, list<UnitigColorMap<CoreInfo>>()));
	cdbgOpt.filename_seq_in.push_back("Test8.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	l.front().second.push_back(*i);
	++i;
	++i;
	uni = *i;
	uni.strand = false;
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	uni.getData()->getData(uni)->predCoreDist = 1;
	l.front().second.push_back(uni);

	markBrdg(l, true);
	EXPECT_EQ(1, l.size());
	EXPECT_EQ(2, l.front().second.size());
	uni = l.front().second.front();
	EXPECT_TRUE(uni.getData()->getData(uni)->sufBrdg);
	uni = l.front().second.back();
	EXPECT_TRUE(uni.getData()->getData(uni)->preBrdg);
}

//Tests for function markBrdg under the following conditions
//	1. Our list of paths is not empty
//	2. Our list consists of one path
//	3. We are not dealing with a list of successive paths
//	5. We are dealing with a list of predecessive paths and are (not) on the unitigs' reference strand
TEST_F(MarkBrdgTest, PredRev){
	l.push_back(Path(4, list<UnitigColorMap<CoreInfo>>()));
	cdbgOpt.filename_seq_in.push_back("Test11_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	l.front().second.push_back(*i);
	++i;
	uni = *i;
	uni.strand = false;
	uni.getData()->getData(uni)->coreList.push_back(pair<uint32_t, uint32_t>(2,3));
	uni.getData()->getData(uni)->sucCoreDist = 3;
	l.front().second.push_back(uni);

	markBrdg(l, false);
	EXPECT_EQ(1, l.size());
	EXPECT_EQ(2, l.front().second.size());
	uni = l.front().second.front();
	EXPECT_TRUE(uni.getData()->getData(uni)->preBrdg);
	uni = l.front().second.back();
	EXPECT_TRUE(uni.getData()->getData(uni)->sufBrdg);
}

//Tests for function void detectBrdg(ColoredCDBG<CoreInfo>&, const uint32_t&)//
//	1. A unitig's core list is (not) empty 0/0
//	2. The last k-mer on a unitig is (not) marked as core 0/0
//	3. A unitig's core list is empty, the unitig's suffix is already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	4. A unitig's core list is not empty, the unitig's suffix is already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	5. The last k-mer on a unitig is marked as core, the unitig's suffix is already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	6. The last k-mer on a unitig is not marked as core, the unitig's suffix is already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	7. A BFS on successive unitigs has to be done, is successful and there is a/no core k-mer on this unitig 0/0
//	8. A BFS on successive unitigs has to be done, fails and there is a/no core k-mer on this unitig 0/0
//	9. A unitig's suffix is (not) already marked as bridging 0/0
//	10. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	11. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	12. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	13. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large 0/0
//	14. A unitig's prefix is (not) already marked as bridging 0/0
//	15. A unitig's core list is empty, the unitig's prefix is already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	16. A unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	17. A unitig's core list is not empty, the unitig's prefix is already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	18. A unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	19. The first k-mer on a unitig is marked as core, the unitig's prefix is already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	20. The first k-mer on a unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	21. The first k-mer on a unitig is not marked as core, the unitig's prefix is already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	22. The first k-mer on a unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large 0/0
//	23. A BFS on predecessive unitigs has to be done and there is (no) core k-mer on the current unitig 0/0
//	24. A BFS on predecessive unitigs has to be done, is successful and there is a/no core k-mer on this unitig 0/0
//	25. A BFS on predecessive unitigs has to be done, fails and there is a/no core k-mer on this unitig 0/0

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
TEST_F(DetectBrdgTest, NoCore){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
}