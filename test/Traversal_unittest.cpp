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
//	3. The last unitig in the top priority path at which we are (not) on the reference strand has a successor for which the distance to the next core k-mer is already known to be too large 1/0
//	4. The last unitig in the top priority path at which we are (not) on the reference strand has a successor for which the distance to the next core k-mer is not already known to be too large DONE
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
//	4. The last unitig in the top priority path at which we are on the reference strand has a successor for which the distance to the next core k-mer is not already known to be too large
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
//	4. The last unitig in the top priority path at which we are not on the reference strand has a successor for which the distance to the next core k-mer is not already known to be too large
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
//	3. The last unitig in the top priority path at which we are on the reference strand has a successor for which the distance to the next core k-mer is already known to be too large
//	4. The last unitig in the top priority path at which we are on the reference strand has a successor for which the distance to the next core k-mer is not already known to be too large
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
//	

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