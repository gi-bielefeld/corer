#include <gtest/gtest.h>
#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "../Traversal.h"
#include "TraversalTest.h"

//Tests for function const bool prioShrtst(const Path&, const Path&)//
//	1. The first path is (not) longer than the second one DONE
//	2. The second path is (not) longer than the first one DONE

//Tests the function prioShrtst under the following conditions
//	1. The first path is longer than the second
//	2. The second path is not longer than the first
TEST_F(PrioShrtstTest, FstLng){ EXPECT_TRUE(prioShrtst(p, q)); }

//Tests the function prioShrtst under the following conditions
//	1. The first path is not longer than the second
//	2. The second path is longer than the first
TEST_F(PrioShrtstTest, SndLng){ EXPECT_FALSE(prioShrtst(q, p)); }

//Tests for function const bool prioShrtst(const TravTrack&, const TravTrack&)//
//	1. The first TravTrack does (not) belong to a longer path than the second one DONE
//	2. The second TravTrack does (not) belong to a longer path than the first one DONE

//Tests the function prioShrtst under the following conditions
//	1. The first TravTrack does belong to a longer path than the second one
//	2. The second TravTrack does not belong to a longer path than the first one
TEST_F(PrioShrtstTest1, FstLng){ EXPECT_TRUE(prioShrtst(t, s)); }

//Tests the function prioShrtst under the following conditions
//	1. The first TravTrack does not belong to a longer path than the second one
//	2. The second TravTrack does belong to a longer path than the first one
TEST_F(PrioShrtstTest1, SndLng){ EXPECT_FALSE(prioShrtst(s, t)); }

//Tests for function const uint32_t getCoreDist(const neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false>&, const bool&)//
//	1. The function is (not) called in the context of a successive traversal DONE
//	2. The function is called for a unitig's reference/reverse complementary strand DONE

//Tests the function getCoreDist under the following conditions
//	1. The function is called in the context of a successive traversal
//	2. The function is called for a unitig's reference strand
TEST_F(GetCoreDistTest, SucRef){
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;

	EXPECT_EQ(1 ,getCoreDist(i->getSuccessors().begin(), true));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is called in the context of a successive traversal
//	2. The function is called for a unitig's reverse complementary strand
TEST_F(GetCoreDistTest, SucRev){
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;
	++i;
	uni = *i;
	uni.strand = false;

	EXPECT_EQ(3, getCoreDist(uni.getSuccessors().begin(), true));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is not called in the context of a successive traversal
//	2. The function is called for a unitig's reference strand
TEST_F(GetCoreDistTest, PredRef){
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	++i;
	++i;
	++i;

	EXPECT_EQ(2, getCoreDist(i->getPredecessors().begin(), false));
}

//Tests the function getCoreDist under the following conditions
//	1. The function is not called in the context of a successive traversal
//	2. The function is called for a unitig's reverse complementary strand
TEST_F(GetCoreDistTest, PredRev){
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	++i;
	uni = *i;
	uni.strand = false;

	EXPECT_EQ(1, getCoreDist(uni.getPredecessors().begin(), false));
}

//Tests for function void expSucPths(priority_queue<Path>&, const uint32_t&, list<Path>&)//
//	1. The last unitig in the top priority path does (not) have successors DONE
//	2. The last unitig in the top priority path has one/many successors DONE
//	3. The last unitig in the top priority path has a successor at which we are (not) on the reference strand and for which the distance to the next core k-mer is already known to be too large DONE
//	4. The last unitig in the top priority path has a successor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large DONE
//	5. The last unitig in the top priority path has a successor on which a/no core k-mer lies DONE
//	6. The last unitig in the top priority path has a successor on which a core k-mer lies that is (not) close enough DONE
//	7. The last unitig in the top priority path has a successor for which adding it to the path does (not) make the path too long DONE
//	8. After processing the top priority path the queue is (not) empty DONE

//Tests the function expSucPths under the following conditions
//	1. The last unitig in the top priority path does have successors
//	2. The last unitig in the top priority path has many successors
//	5. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	6. The last unitig in the top priority path has a successor on which a core k-mer lies that is close enough
//	7. The last unitig in the top priority path has a successor for which adding it to the path does make the path too long
//	8. After processing the top priority path the queue is empty
TEST_F(ExpSucPthsTest, PthTooLng){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(2,2));
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));

	expSucPths(queue, 3, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_EQ(1, res.size());

	for(list<Path>::const_iterator p = res.begin(); p != res.end(); ++p){
		EXPECT_EQ(3, p->first);
		EXPECT_EQ(2, p->second.size());
		uint32_t nbUni = 0;

		for(list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin(); u != p->second.end(); ++u){
			++nbUni;

			if(nbUni == 1) EXPECT_EQ("AAAGGCAAA", u->mappedSequenceToString());

			if(nbUni == 2) EXPECT_EQ("AAGGCAAAGAC", u->mappedSequenceToString());
		}
	}
}

//Tests the function expSucPths under the following conditions
//	1. The last unitig in the top priority path does (not) have successors
//	2. The last unitig in the top priority path has many successors
//	4. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a successor on which no core k-mer lies
//	7. The last unitig in the top priority path has a successor for which adding it to the path does not make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpSucPthsTest, NoSuc){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));

	expSucPths(queue, 5, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function expSucPths under the following conditions
//	1. The last unitig in the top priority path does have successors
//	2. The last unitig in the top priority path has one successor
//	4. The last unitig in the top priority path has a successor at which we are not on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a successor on which a core k-mer lies
//	6. The last unitig in the top priority path has a successor on which a core k-mer lies that is not close enough
//	7. The last unitig in the top priority path has a successor for which adding it to the path does make the path too long
//	8. After processing the top priority path the queue is empty
TEST_F(ExpSucPthsTest, OneSuc){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;
	++i;
	l.push_back(*i);
	l.front().strand = false;
	queue.push(Path(0, l));

	expSucPths(queue, 2, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function expSucPths under the following conditions
//	1. The last unitig in the top priority path does have successors
//	2. The last unitig in the top priority path has one successor
//	3. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	4. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	6. The last unitig in the top priority path has a successor on which a core k-mer lies that is not close enough
//	7. The last unitig in the top priority path has a successor for which adding it to the path does neighborIterator make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpSucPthsTest, SndVis){
	cdbgOpt.filename_ref_in.push_back("Test17_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	i->getData()->getData(*i)->sucCoreDist = 2;
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));

	expSucPths(queue, 3, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function expSucPths under the following conditions
//	1. The last unitig in the top priority path does (not) have successors
//	2. The last unitig in the top priority path has one successor
//	3. The last unitig in the top priority path has a successor at which we are not on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	4. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a successor for which adding it to the path does not make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpSucPthsTest, RevVis){
	cdbgOpt.filename_ref_in.push_back("Test18_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->predCoreDist = 2;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));

	expSucPths(queue, 3, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests for function inline const uint32_t calcOff(const uint32_t&, const size_t&, const bool&)//
//	1. An offset on the reference/reverse complementary strand has to be calculated DONE

//Tests the function calcOff under the following conditions
//	1. An offset on the reference strand has to be calculated
TEST(CalcOffTest, refStnd){
	EXPECT_EQ(42, calcOff(42, 50, true));
}

//Tests the function calcOff under the following conditions
//	1. An offset on the reverse complementary strand has to be calculated
TEST(CalcOffTest, revStnd){
	EXPECT_EQ(4, calcOff(2, 7, false));
}

//Tests for function void addDists(const list<Path>&, const bool&)//
//	1. The path list is (not) empty DONE
//	2. A path's length is (not) larger than 2 DONE
//	3. We are dealing with a successive path that does (not) end on the reverse complementary strand DONE
//	4. We are dealing with a predecessive path that does (not) end on the reverse complementary strand DONE
//	5. We are dealing with a successive path and on some unitig not at the path's end we are (not) on the reverse complementary strand DONE
//	6. We are dealing with a predecessive path and on some unitig not at the path's end we are (not) on the reverse complementary strand DONE
//	7. The path list contains more than one successive path in which the same unitig occurs but with differing distance to the core k-mer at which the path ends and the unitig is (not) visited on the reverse complementary strand DONE
//	8. The path list contains more than one predecessive path in which the same unitig occurs but with differing distance to the core k-mer at which the path ends and the unitig is (not) visited on the reverse complementary strand DONE

//Tests the function addDists under the following conditions
//	1. The path list is empty
TEST_F(AddDistsTest, NoPth){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	addDists(pList, true);

	for( ; i != cdbg.end(); ++i){
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
		EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	}
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is (not) larger than 2
//	3. We are dealing with a successive path that does not end on the reverse complementary strand
//	5. We are dealing with a successive path and on some unitig not at the path's end we are not on the reverse complementary strand
TEST_F(AddDistsTest, SucPths){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uList.push_back(*i);
	pList.push_back(Path(4, uList));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,2));
	uList = list<UnitigColorMap<CoreInfo>>();
	uList.push_back(*i);
	pList.push_back(Path(2, uList));
	++i;
	pList.front().second.push_front(*i);
	pList.back().second.push_front(*i);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->sucCoreDist = 1;
	pList.front().second.push_back(*i);
	addDists(pList, true);

	i = cdbg.begin();
	EXPECT_EQ(4, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(2, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(1, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is not larger than 2
//	4. We are dealing with a predecessive path that does not end on the reverse complementary strand
TEST_F(AddDistsTest, PredPth){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uList.push_back(*i);
	pList.push_back(Path(1, uList));
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	pList.front().second.push_back(*i);
	addDists(pList, false);

	i = cdbg.begin();
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(1, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is larger than 2
//	4. We are dealing with a predecessive path that does not end on the reverse complementary strand
//	6. We are dealing with a predecessive path and on some unitig not at the path's end we are not on the reverse complementary strand
TEST_F(AddDistsTest, LpredP){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uList.push_back(*i);
	pList.push_back(Path(4, uList));
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->predCoreDist = 1;
	pList.front().second.push_back(*i);
	++i;
	pList.front().second.push_front(*i);
	addDists(pList, false);

	i = cdbg.begin();
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(4, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(1, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is not larger than 2
//	3. We are dealing with a successive path that does end on the reverse complementary strand
TEST_F(AddDistsTest, SucRev){
	cdbgOpt.filename_seq_in.push_back("Test8.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uList.push_back(*i);
	pList.push_back(Path(1, uList));
	++i;
	++i;
	u = *i;
	u.strand = false;
	u.getData()->getData(u)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	pList.front().second.push_back(u);
	addDists(pList, true);

	i = cdbg.begin();
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(1, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is larger than 2
//	4. We are dealing with a predecessive path that does end on the reverse complementary strand
//	6. We are dealing with a predecessive path and on some unitig not at the path's end we are (not) on the reverse complementary strand
TEST_F(AddDistsTest, LpRev){
	cdbgOpt.filename_seq_in.push_back("Test11_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->sucCoreDist = 6;
	i->getData()->getData(*i)->predCoreDist = 4;
	u = *i;
	u.strand = false;
	uList.push_back(u);
	uList.push_back(*i);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,2));
	i->getData()->getData(*i)->sucCoreDist = 2;
	uList.push_front(*i);
	u = *i;
	u.strand = false;
	uList.push_back(u);
	pList.push_back(Path(6, uList));
	addDists(pList, false);

	i = cdbg.begin();
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(6, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(4, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(2, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is larger than 2
//	3. We are dealing with a successive path that does not end on the reverse complementary strand
//	5. We are dealing with a successive path and on some unitig not at the path's end we are (not) on the reverse complementary strand
TEST_F(AddDistsTest, InterRev){
	cdbgOpt.filename_seq_in.push_back("Test11_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uList.push_back(*i);
	++i;
	uList.push_back(*i);
	u = *i;
	u.strand = false;
	uList.push_back(u);
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,2));
	uList.push_back(*i);
	pList.push_back(Path(6, uList));
	addDists(pList, true);

	i = cdbg.begin();
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(6, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(4, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(2, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(UINT32_MAX, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is larger than 2
//	3. We are dealing with a successive path that does (not) end on the reverse complementary strand
//	5. We are dealing with a successive path and on some unitig not at the path's end we are (not) on the reverse complementary strand
//	7. The path list contains more than one successive path in which the same unitig occurs but with differing distance to the core k-mer at which the path ends and the unitig is (not) visited on the reverse complementary strand
TEST_F(AddDistsTest, VarLenSuc){
	cdbgOpt.filename_seq_in.push_back("Test15_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();//1
	uList.push_back(*i);//[1]
	pList.push_back(Path(3, uList));//[(3, [1])]
	pList.push_back(Path(5, uList));//[(3, [1]), (5, [1])]
	++i;//2
	pList.front().second.push_front(*i);//[(3, [2, 1]), (5, [1])]
	pList.front().second.push_back(*i);//[(3, [2, 1, 2]), (5, [1])]
	pList.back().second.push_front(*i);//[(3, [2, 1, 2]), (5, [2, 1])]
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(4,5));
	++i;//3
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3])]
	i = cdbg.begin();//1
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1])]
	++i;//2
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2])]
	uList.back().strand = false;//[1c]
	pList.push_back(Path(3, uList));//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [1c])]
	u = *i;//2
	u.strand = false;//2c
	pList.back().second.push_front(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c])]
	pList.back().second.push_back(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c])]
	pList.push_back(Path(5, uList));//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [1c])]
	pList.back().second.push_front(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c])]
	++i;//3
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3])]
	pList.back().second.back().strand = false;//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c])]
	pList.back().second.push_back(*cdbg.begin());//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1])]
	pList.back().second.back().strand = false;//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1c])]
	pList.back().second.push_back(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1c, 2c])]
	addDists(pList, true);

	i = cdbg.begin();
	EXPECT_EQ(21, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(20, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(5, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(4, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(31, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(30, i->getData()->getData(*i)->predCoreDist);
}

//Tests the function addDists under the following conditions
//	1. The path list is not empty
//	2. A path's length is larger than 2
//	4. We are dealing with a predecessive path that does (not) end on the reverse complementary strand
//	6. We are dealing with a predecessive path and on some unitig not at the path's end we are (not) on the reverse complementary strand
//	8. The path list contains more than one predecessive path in which the same unitig occurs but with differing distance to the core k-mer at which the path ends and the unitig is (not) visited on the reverse complementary strand
TEST_F(AddDistsTest, VarLenPred){
	cdbgOpt.filename_seq_in.push_back("Test15_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();//1
	uList.push_back(*i);//[1]
	pList.push_back(Path(3, uList));//[(3, [1])]
	pList.push_back(Path(5, uList));//[(3, [1]), (5, [1])]
	++i;//2
	pList.front().second.push_front(*i);//[(3, [2, 1]), (5, [1])]
	pList.front().second.push_back(*i);//[(3, [2, 1, 2]), (5, [1])]
	pList.back().second.push_front(*i);//[(3, [2, 1, 2]), (5, [2, 1])]
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(4,5));
	++i;//3
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3])]
	i = cdbg.begin();//1
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1])]
	++i;//2
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2])]
	uList.back().strand = false;//[1c]
	pList.push_back(Path(3, uList));//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [1c])]
	u = *i;//2
	u.strand = false;//2c
	pList.back().second.push_front(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c])]
	pList.back().second.push_back(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c])]
	pList.push_back(Path(5, uList));//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [1c])]
	pList.back().second.push_front(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c])]
	++i;//3
	pList.back().second.push_back(*i);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3])]
	pList.back().second.back().strand = false;//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c])]
	pList.back().second.push_back(*cdbg.begin());//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1])]
	pList.back().second.back().strand = false;//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1c])]
	pList.back().second.push_back(u);//[(3, [2, 1, 2]), (5, [2, 1, 3, 1, 2]), (3, [2c, 1c, 2c]), (5, [2c, 1c, 3c, 1c, 2c])]
	addDists(pList, false);

	i = cdbg.begin();
	EXPECT_EQ(21, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(20, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(5, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(4, i->getData()->getData(*i)->predCoreDist);
	++i;
	EXPECT_EQ(31, i->getData()->getData(*i)->sucCoreDist);
	EXPECT_EQ(30, i->getData()->getData(*i)->predCoreDist);
}

//Tests for function const bool doSucBFS(const UnitigColorMap<CoreInfo>, const uint32_t, list<Path>&)//
//	1. The result path is (not) empty DONE
//	2. The last unitig in the top priority path does (not) have successors DONE
//	3. The last unitig in the top priority path has one/many successors DONE
//	4. The last unitig in the top priority path has a successor at which we are (not) on the reference strand and for which the distance to the next core k-mer is already known to be too large DONE
//	5. The last unitig in the top priority path has a successor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large DONE
//	6. The last unitig in the top priority path has a successor on which a/no core k-mer lies DONE
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is (not) close enough DONE
//	8. The last unitig in the top priority path has a successor for which adding it to the path does (not) make the path too long DONE
//	9. After processing the top priority path the queue is (not) empty DONE


//Tests the function doSucBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has many successors
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which no core k-mer lies
//	8. The last unitig in the top priority path has a successor for which adding it to the path does make the path too long
//	9. After processing the top priority path the queue is empty
TEST_F(DoSucBFStest, NoRes){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->sucCoreDist = 1;
	i = cdbg.begin();
	++i;
	++i;

	EXPECT_FALSE(doSucBFS(*i, 3, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is not empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has many successors
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is close enough
//	8. The last unitig in the top priority path has a successor for which adding it to the path does (not) make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoSucBFStest, SomeRes){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	uint32_t c = 0;

	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(2,2));
	++i;
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->sucCoreDist = 1;
	i = cdbg.begin();
	++i;
	++i;

	EXPECT_TRUE(doSucBFS(*i, 4, res));
	EXPECT_EQ(3, res.front().first);
	EXPECT_EQ("AAAGGCAAA", res.front().second.front().mappedSequenceToString());
	EXPECT_EQ("AAGGCAAAGAC", res.front().second.back().mappedSequenceToString());
	EXPECT_EQ(4, res.back().first);
	EXPECT_EQ(3, res.back().second.size());

	for(list<UnitigColorMap<CoreInfo>>::const_iterator u = res.back().second.begin(); u != res.back().second.end(); ++u, ++c){
		if(c == 0) EXPECT_EQ("AAAGGCAAA", u->mappedSequenceToString());

		if(c == 1) EXPECT_EQ("AAGGCAAACAC", u->mappedSequenceToString());

		if(c == 2) EXPECT_EQ("GCAAACACC", u->mappedSequenceToString());
	}
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is not empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has many successors
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is close enough
//	8. The last unitig in the top priority path has a successor for which adding it to the path does make the path too long
//	9. After processing the top priority path the queue is empty
TEST_F(DoSucBFStest, PthTooLng){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(2,2));
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i = cdbg.begin();
	++i;
	++i;

	EXPECT_TRUE(doSucBFS(*i, 3, res));
	EXPECT_EQ(1, res.size());
	EXPECT_EQ(3, res.front().first);
	EXPECT_EQ(2, res.front().second.size());
	EXPECT_EQ("AAAGGCAAA", res.front().second.front().mappedSequenceToString());
	EXPECT_EQ("AAGGCAAAGAC", res.front().second.back().mappedSequenceToString());
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does (not) have successors
//	3. The last unitig in the top priority path has many successors
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which no core k-mer lies
//	8. The last unitig in the top priority path has a successor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoSucBFStest, NoSuc){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;

	EXPECT_FALSE(doSucBFS(*i, 5, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has one successor
//	5. The last unitig in the top priority path has a successor at which we are not on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which a core k-mer lies
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is not close enough
//	8. The last unitig in the top priority path has a successor for which adding it to the path does make the path too long
//	9. After processing the top priority path the queue is empty
TEST_F(DoSucBFStest, OneSuc){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;
	++i;
	u = *i;
	u.strand = false;

	EXPECT_FALSE(doSucBFS(u, 2, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has one successor
//	4. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is not close enough
//	8. The last unitig in the top priority path has a successor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoSucBFStest, SndVis){
	cdbgOpt.filename_ref_in.push_back("Test17_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(1,1));
	i->getData()->getData(*i)->sucCoreDist = 2;
	++i;
	++i;

	EXPECT_FALSE(doSucBFS(*i, 3, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doSucBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have successors
//	3. The last unitig in the top priority path has one successor
//	4. The last unitig in the top priority path has a successor at which we are not on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	5. The last unitig in the top priority path has a successor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a successor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a successor on which a core k-mer lies that is not close enough
//	8. The last unitig in the top priority path has a successor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoSucBFStest, RevVis){
	cdbgOpt.filename_ref_in.push_back("Test18_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	i->getData()->getData(*i)->predCoreDist = 2;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	++i;

	EXPECT_FALSE(doSucBFS(*i, 3, res));
	EXPECT_TRUE(res.empty());
}

//Tests for function inline const uint32_t findMinPthLen(const list<Path>&)//
//	1. The list of paths is (not) empty DONE
//	2. The temporary minimum has (not) to be updated DONE

//Tests the function findMinPthLen under the following conditions
//	1. The list of paths is empty
TEST_F(FindMinPthLenTest, NoPth){
	EXPECT_EQ(UINT32_MAX, findMinPthLen(l));
}

//Tests the function findMinPthLen under the following conditions
//	1. The list of paths is not empty
//	2. The temporary minimum has (not) to be updated
TEST_F(FindMinPthLenTest, SomePths){
	l.push_back(Path(33, list<UnitigColorMap<CoreInfo>>()));
	l.push_back(Path(21, list<UnitigColorMap<CoreInfo>>()));
	l.push_back(Path(42, list<UnitigColorMap<CoreInfo>>()));

	EXPECT_EQ(21, findMinPthLen(l));
}

//Tests for function void expPredPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>&, const uint32_t&, list<Path>&)//
//	1. The last unitig in the top priority path does (not) have predecessors DONE
//	2. The last unitig in the top priority path has one/many predecessor(s) DONE
//	3. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is already known to be too large DONE
//	4. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large DONE
//	5. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies DONE
//	6. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is (not) close enough DONE
//	7. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long DONE
//	8. After processing the top priority path the queue is (not) empty DONE

//Tests the function expPredPths under the following conditions
//	1. The last unitig in the top priority path does have predecessors
//	2. The last unitig in the top priority path has one predecessor
//	4. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a predecessor on which a core k-mer lies
//	6. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is close enough
//	8. After processing the top priority path the queue is empty
TEST_F(ExpPredPthsTest, HasPred){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	l.push_back(*i);
	queue.push(Path(0, l));
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));

	expPredPths(queue, 42, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_EQ(1, res.size());

	for(list<Path>::const_iterator p = res.begin(); p != res.end(); ++p){
		EXPECT_EQ(1, p->first);
		EXPECT_EQ(2, p->second.size());
		
		EXPECT_EQ("AAGGCAAACAC", p->second.front().mappedSequenceToString());

		EXPECT_EQ("AAAGGCAAA", p->second.back().mappedSequenceToString());
	}
}

//Tests the function expPredPths under the following conditions
//	1. The last unitig in the top priority path does not have predecessors
//	8. After processing the top priority path the queue is empty
TEST_F(ExpPredPthsTest, NoPred){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));

	expPredPths(queue, 42, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function expPredPths under the following conditions
//	1. The last unitig in the top priority path does have predecessors
//	2. The last unitig in the top priority path has many predecessors
//	4. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies
//	6. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is not close enough
//	7. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpPredPthsTest, MnyPred){
	cdbgOpt.filename_seq_in.push_back("Test15_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(7,8));
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(13,14));
	l.push_back(*i);
	queue.push(Path(0, l));

	expPredPths(queue, 10, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function expPredPths under the following conditions
//	1. The last unitig in the top priority path does have predecessors
//	2. The last unitig in the top priority path has one/many predecessor(s)
//	3. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	4. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies
//	6. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is not close enough
//	7. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpPredPthsTest, FltdRef){
	cdbgOpt.filename_seq_in.push_back("Test11_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	i->getData()->getData(*i)->predCoreDist = 2;
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));

	expPredPths(queue, 5, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests the function  expPredPths under the following conditions
//	1. The last unitig in the top priority path does have predecessors
//	2. The last unitig in the top priority path has one predecessor
//	3. The last unitig in the top priority path has a predecessor at which we are not on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	4. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	5. The last unitig in the top priority path has a predecessor on which no core k-mer lies
//	7. The last unitig in the top priority path has a predecessor for which adding it to the path does not make the path too long
//	8. After processing the top priority path the queue is (not) empty
TEST_F(ExpPredPthsTest, FltdRev){
	cdbgOpt.filename_ref_in.push_back("Test19_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	++i;
	++i;
	++i;
	l.push_back(*i);
	queue.push(Path(0, l));
	++i;
	++i;
	i->getData()->getData(*i)->sucCoreDist = 2;

	expPredPths(queue, 3, res);
	EXPECT_TRUE(queue.empty());
	EXPECT_TRUE(res.empty());
}

//Tests for function const bool doPredBFS(const UnitigColorMap<CoreInfo>, const uint32_t, list<Path>&)//
//	1. The result path is (not) empty DONE
//	2. The last unitig in the top priority path does (not) have predecessors DONE
//	3. The last unitig in the top priority path has one/many predecessor(s) DONE
//	4. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is already known to be too large DONE
//	5. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large DONE
//	6. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies DONE
//	7. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is (not) close enough DONE
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long DONE
//	9. After processing the top priority path the queue is (not) empty DONE

//Tests the function doPredBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does (not) have predecessors
//	3. The last unitig in the top priority path has one predecessor
//	5. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which no core k-mer lies
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoPredBFStest, NoRes){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	++i;

	EXPECT_FALSE(doPredBFS(*i, 42, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is not empty
//	2. The last unitig in the top priority path does have predecessors
//	3. The last unitig in the top priority path has one predecessor
//	5. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is close enough
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoPredBFStest, SomeRes){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;

	EXPECT_TRUE(doPredBFS(*i, 4, res));
	EXPECT_EQ(1, res.size());
	EXPECT_EQ(4, res.front().first);
	EXPECT_EQ(3, res.front().second.size());
	list<UnitigColorMap<CoreInfo>>::const_iterator ui = res.front().second.begin();
	EXPECT_EQ("GCAAACACA", ui->mappedSequenceToString());
	++ui;
	EXPECT_EQ("AAGGCAAACAC", ui->mappedSequenceToString());
	EXPECT_EQ(4, ui->getData()->getData(*ui)->predCoreDist);
	++ui;
	EXPECT_EQ("AAAGGCAAA", ui->mappedSequenceToString());
	EXPECT_EQ(1, ui->getData()->getData(*ui)->predCoreDist);
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is not empty
//	2. The last unitig in the top priority path does have predecessors
//	3. The last unitig in the top priority path has one predecessor
//	5. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which a core k-mer lies
//	7. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is close enough
//	9. After processing the top priority path the queue is empty
TEST_F(DoPredBFStest, HasPred){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));

	EXPECT_TRUE(doPredBFS(*cdbg.begin(), 42, res));
	EXPECT_EQ(1, res.size());

	for(list<Path>::const_iterator p = res.begin(); p != res.end(); ++p){
		EXPECT_EQ(1, p->first);
		EXPECT_EQ(2, p->second.size());
		
		EXPECT_EQ("AAGGCAAACAC", p->second.front().mappedSequenceToString());

		EXPECT_EQ("AAAGGCAAA", p->second.back().mappedSequenceToString());
	}
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does not have predecessors
//	9. After processing the top priority path the queue is empty
TEST_F(DoPredBFStest, NoPred){
	cdbgOpt.filename_seq_in.push_back("Test.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	++i;
	++i;

	EXPECT_FALSE(doPredBFS(*i, 42, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have predecessors
//	3. The last unitig in the top priority path has one/many predecessors
//	5. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is not close enough
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoPredBFStest, MnyPred){
	cdbgOpt.filename_seq_in.push_back("Test15_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(7,8));
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(13,14));

	EXPECT_FALSE(doPredBFS(*i, 10, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have predecessors
//	3. The last unitig in the top priority path has one/many predecessor(s)
//	4. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	5. The last unitig in the top priority path has a predecessor at which we are (not) on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which a/no core k-mer lies
//	7. The last unitig in the top priority path has a predecessor on which a core k-mer lies that is not close enough
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does (not) make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoPredBFStest, FltdRef){
	cdbgOpt.filename_seq_in.push_back("Test11_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	i->getData()->getData(*i)->predCoreDist = 2;
	++i;
	++i;

	EXPECT_FALSE(doPredBFS(*i, 5, res));
	EXPECT_TRUE(res.empty());
}

//Tests the function doPredBFS under the following conditions
//	1. The result path is empty
//	2. The last unitig in the top priority path does have predecessors
//	3. The last unitig in the top priority path has one predecessor
//	4. The last unitig in the top priority path has a predecessor at which we are not on the reference strand and for which the distance to the next core k-mer is already known to be too large
//	5. The last unitig in the top priority path has a predecessor at which we are on the reference strand and for which the distance to the next core k-mer is not already known to be too large
//	6. The last unitig in the top priority path has a predecessor on which no core k-mer lies
//	8. The last unitig in the top priority path has a predecessor for which adding it to the path does not make the path too long
//	9. After processing the top priority path the queue is (not) empty
TEST_F(DoPredBFStest, FltdRev){
	cdbgOpt.filename_ref_in.push_back("Test19_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	i = cdbg.begin();
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,0));
	++i;
	i->getData()->getData(*i)->coreList.push_back(pair<uint32_t, uint32_t>(0,1));
	++i;
	++i;
	++i;
	++i;
	++i;
	i->getData()->getData(*i)->sucCoreDist = 2;
	i = cdbg.begin();
	++i;
	++i;
	++i;
	++i;

	EXPECT_FALSE(doPredBFS(*i, 3, res));
	EXPECT_TRUE(res.empty());
}