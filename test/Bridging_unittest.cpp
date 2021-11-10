#include <gtest/gtest.h>

#include "BridgingTest.h"
#include "../src/Bridging.cpp"

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
	len = 41;
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
	len = 41;
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

//Tests for function void void markBrdg(const list<Path>&, const bool&, const uint32_t&)//
//	1. Our list of paths is (not) empty DONE
//	2. Our list consists of one/many path(s) DONE
//	3. We are (not) dealing with a list of successive paths DONE
//	4. We are dealing with a list of successive paths and are (not) on the unitigs' reference strand DONE
//	5. We are dealing with a list of predecessive paths and are (not) on the unitigs' reference strand DONE

//Tests the function markBrdg under the following conditions
//	1. Our list of paths is empty
//	3. We are dealing with a list of successive paths
TEST_F(MarkBrdgTest, NoPths){
	markBrdg(l, true, dlt);
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

	markBrdg(l, false, dlt);
	EXPECT_EQ(1, l.size());
	list<UnitigColorMap<CoreInfo>>::const_iterator u = l.front().second.begin();
	EXPECT_TRUE(u->getData()->getData(*u)->preBrdg);
	++u;

	for(; u != l.front().second.end(); ++u) EXPECT_TRUE(u->getData()->getData(*u)->sufBrdg);
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

	markBrdg(l, true, dlt);
	EXPECT_EQ(2, l.size());
	EXPECT_EQ(2, l.front().second.size());
	EXPECT_EQ(2, l.back().second.size());

	for(list<Path>::const_iterator p = l.begin(); p != l.end(); ++p){
		list<UnitigColorMap<CoreInfo>>::const_iterator u = l.front().second.begin();
		EXPECT_TRUE(u->getData()->getData(*u)->sufBrdg);
		++u;

		for(; u != l.front().second.end(); ++u) EXPECT_TRUE(u->getData()->getData(*u)->preBrdg);
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

	markBrdg(l, true, dlt);
	EXPECT_EQ(1, l.size());
	EXPECT_EQ(2, l.front().second.size());
	uni = l.front().second.front();
	EXPECT_TRUE(uni.getData()->getData(uni)->sufBrdg);
	uni = l.front().second.back();
	EXPECT_TRUE(uni.getData()->getData(uni)->sufBrdg);
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

	markBrdg(l, false, dlt);
	EXPECT_EQ(1, l.size());
	EXPECT_EQ(2, l.front().second.size());
	uni = l.front().second.front();
	EXPECT_TRUE(uni.getData()->getData(uni)->preBrdg);
	uni = l.front().second.back();
	EXPECT_TRUE(uni.getData()->getData(uni)->preBrdg);
}

//Tests for function void markBrdg(ColoredCDBG<CoreInfo>&, const uint32_t&)//
//	1. There is a/no core k-mer on the current unitig DONE
//	2. There is no core k-mer on the current unitig and predCoreDist was (never) updated DONE
//	3. There is no core k-mer on the current unitig and sucCoreDist was (never) updated DONE
//	4. There is no core k-mer on the current unitig, predCoreDist and sucCoreDist were both updated and this unitig is (not) bridging DONE
//	5. There is a core k-mer on the current unitig and predCoreDist was (never) updated DONE
//	6. There is a core k-mer on the current unitig and sucCoreDist was (never) updated DONE
//	7. There is a core k-mer on the current unitig, predCoreDist was updated and the unitig's prefix is (not) bridging DONE
//	8. There is a core k-mer on the current unitig, sucCoreDist was updated and the unitig's suffix is (not) bridging DONE

//Tests the function markBrdg under the following conditions
//	1. There is a/no core k-mer on the current unitig
//	2. There is no core k-mer on the current unitig and predCoreDist was updated
//	3. There is no core k-mer on the current unitig and sucCoreDist was (never) updated
//	4. There is no core k-mer on the current unitig, predCoreDist and sucCoreDist were both updated and this unitig is bridging
//	5. There is a core k-mer on the current unitig and predCoreDist was (never) updated
//	6. There is a core k-mer on the current unitig and sucCoreDist was (never) updated
//	7. There is a core k-mer on the current unitig, predCoreDist was updated and the unitig's prefix is bridging
//	8. There is a core k-mer on the current unitig, sucCoreDist was updated and the unitig's suffix is bridging
TEST_F(MarkBrdgTest1, DltFul){
	queue = detectCore(cdbg, 2, dlt);
	annotateDists(cdbg, queue, dlt);
	markBrdg(cdbg, dlt);
	i = cdbg.begin();

	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
}

//Tests the function markBrdg under the following conditions
//	1. There is a/no core k-mer on the current unitig
//	2. There is no core k-mer on the current unitig and predCoreDist was updated
//	3. There is no core k-mer on the current unitig and sucCoreDist was (never) updated
//	4. There is no core k-mer on the current unitig, predCoreDist and sucCoreDist were both updated and this unitig is not bridging
//	5. There is a core k-mer on the current unitig and predCoreDist was (never) updated
//	6. There is a core k-mer on the current unitig and sucCoreDist was (never) updated
//	7. There is a core k-mer on the current unitig, predCoreDist was updated and the unitig's prefix is not bridging
//	8. There is a core k-mer on the current unitig, sucCoreDist was updated and the unitig's suffix is not bridging
TEST_F(MarkBrdgTest1, NotBrd){
	dlt = 3;
	queue = detectCore(cdbg, 2, dlt);
	annotateDists(cdbg, queue, dlt);
	markBrdg(cdbg, dlt);
	i = cdbg.begin();

	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
}

//Tests the function markBrdg under the following conditions
//	1. There is no core k-mer on the current unitig
//	2. There is no core k-mer on the current unitig and predCoreDist was never updated
//	3. There is no core k-mer on the current unitig and sucCoreDist was never updated
TEST_F(MarkBrdgTest1, NoCr){
	queue = detectCore(cdbg, 3, dlt);
	annotateDists(cdbg, queue, dlt);
	markBrdg(cdbg, dlt);
	i = cdbg.begin();

	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
}

//Tests for function void detectBrdg(ColoredCDBG<CoreInfo>&, const uint32_t&)//
//	1. A unitig's core list is (not) empty DONE
//	2. The last k-mer on a unitig is (not) marked as core DONE
//	3. A BFS on successive unitigs has to be done, is successful and there is a/no core k-mer on this unitig DONE
//	4. A BFS on successive unitigs has to be done, fails and there is a/no core k-mer on this unitig DONE
//	5. A unitig's suffix is (not) already marked as bridging DONE
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large DONE
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large DONE
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large DONE
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large DONE
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is (not) already marked as bridging DONE
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large DONE
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large DONE
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large DONE
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large DONE
//	15. A BFS on predecessive unitigs has to be done and there is a/no core k-mer on the current unitig DONE
//	16. A BFS on predecessive unitigs has to be done, is successful and there is a/no core k-mer on this unitig DONE
//	17. A BFS on predecessive unitigs has to be done, fails and there is a/no core k-mer on this unitig DONE
//	18. A successive path list of an earlier processed unitig still exists when a BFS on predecessive unitigs is successful for the next unitig DONE

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	4. A BFS on successive unitigs has to be done, fails and there is no core k-mer on this unitig
//	5. A unitig's suffix is not already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is (not) already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large
TEST_F(DetectBrdgTest, NoCore){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(2,2));
	++i;
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	detectBrdg(cdbg, 2);

	for(i = cdbg.begin(); i != cdbg.end(); ++i){
		EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
		EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	}
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	3. A BFS on successive unitigs has to be done, is successful and there is no core k-mer on this unitig
//	4. A BFS on successive unitigs has to be done, fails and there is a/no core k-mer on this unitig
//	5. A unitig's suffix is (not) already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is (not) already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	15. A BFS on predecessive unitigs has to be done and there is no core k-mer on the current unitig
//	16. A BFS on predecessive unitigs has to be done, is successful and there is no core k-mer on this unitig
TEST_F(DetectBrdgTest, SucBFS){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	detectBrdg(cdbg, 4);

	i = cdbg.begin();
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(1, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(1, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	3. A BFS on successive unitigs has to be done, is successful and there is a core k-mer on this unitig
//	4. A BFS on successive unitigs has to be done, fails and there is no core k-mer on this unitig
//	5. A unitig's suffix is (not) already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is (not) already marked as bridging
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	15. A BFS on predecessive unitigs has to be done and there is a core k-mer on the current unitig
//	16. A BFS on predecessive unitigs has to be done, is successful and there is a core k-mer on this unitig
TEST_F(DetectBrdgTest, SucSucBFSCore){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	detectBrdg(cdbg, 2);

	i = cdbg.begin();
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(1, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(1, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	5. A unitig's suffix is not already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is not already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is already too large
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is already too large
TEST_F(DetectBrdgTest, NoDlt){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,2));
	detectBrdg(cdbg, 0);

	for(i = cdbg.begin(); i != cdbg.end(); ++i){
		EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
		EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	}
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is not marked as core
//	4. A BFS on successive unitigs has to be done, fails and there is a/no core k-mer on this unitig
//	5. A unitig's suffix is not already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is (not) already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is not already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is (not) already too large
//	15. A BFS on predecessive unitigs has to be done and there is a core k-mer on the current unitig
//	17. A BFS on predecessive unitigs has to be done, fails and there is a core k-mer on this unitig
TEST_F(DetectBrdgTest, PredBFS){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	detectBrdg(cdbg, 2);

	for(i = cdbg.begin(); i != cdbg.end(); ++i){
		EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
		EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	}
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	3. A BFS on successive unitigs has to be done, is successful and there is no core k-mer on this unitig
//	4. A BFS on successive unitigs has to be done, fails and there is no core k-mer on this unitig
//	5. A unitig's suffix is not already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is not already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	13. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	15. A BFS on predecessive unitigs has to be done and there is no core k-mer on the current unitig
//	17. A BFS on predecessive unitigs has to be done, fails and there is no core k-mer on this unitig
TEST_F(DetectBrdgTest, PredBFSNoCore){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	detectBrdg(cdbg, 4);

	for(i = cdbg.begin(); i != cdbg.end(); ++i){
		EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
		EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);

		if(i->mappedSequenceToString() == "GCAAACACA"){
			EXPECT_EQ(1, i->getData()->getData(*i)->sucCoreDist);
		} else{
			EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
		}
		
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	}
}

//Tests for function detectBrdg under the following conditions
//	1. A unitig's core list is (not) empty
//	2. The last k-mer on a unitig is (not) marked as core
//	3. A BFS on successive unitigs has to be done, is successful and there is no core k-mer on this unitig
//	5. A unitig's suffix is (not) already marked as bridging
//	6. A unitig's core list is empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	7. A unitig's core list is not empty, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	8. The last k-mer on a unitig is marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	9. The last k-mer on a unitig is not marked as core, the unitig's suffix is not already marked as bridging and the distance to bridge to the left side is not already too large
//	10. We have to continue with the current unitig after a potential BFS on successive unitigs and the current unitig's prefix is (not) already marked as bridging
//	11. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	12. We have to continue with the current unitig after a potential BFS on successive unitigs, the current unitig's core list is not empty, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	14. We have to continue with the current unitig after a potential BFS on successive unitigs, the first k-mer on the unitig is not marked as core, the unitig's prefix is not already marked as bridging and the distance to bridge to the right side is not already too large
//	15. A BFS on predecessive unitigs has to be done and there is a/no core k-mer on the current unitig
//	16. A BFS on predecessive unitigs has to be done, is successful and there is a/no core k-mer on this unitig
//	17. A BFS on predecessive unitigs has to be done, fails and there is a/no core k-mer on this unitig
//	18. A successive path list of an earlier processed unitig still exists when a BFS on predecessive unitigs is successful for the next unitig
TEST_F(DetectBrdgTest, OldSucPth){
	cdbgOpt.filename_ref_in.push_back("Test20_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	detectBrdg(cdbg, 3);

	i = cdbg.begin();
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(2, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(1, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_FALSE(i->getData()->getData(*i)->preBrdg);
	EXPECT_FALSE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_TRUE(i->getData()->getData(*i)->preBrdg);
	EXPECT_TRUE(i->getData()->getData(*i)->sufBrdg);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(2, i->getData()->getData(*i)->predCoreDist);
}