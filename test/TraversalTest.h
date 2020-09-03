#ifndef TRAVERSAL_TEST_HPP
#define TRAVERSAL_TEST_HPP

#include "CoreTest.h"

class PrioShrtstTest : public ::testing::Test {

	protected:

		PrioShrtstTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
			cdbgOpt.filename_seq_in.push_back("Test.fa");
			cdbg.build(cdbgOpt);
			cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
			cdbg.buildColors(cdbgOpt);
			i = cdbg.begin();
			p.first = 4;
			p.second.push_back(*i);
			++i;
			++i;
			++i;
			p.second.push_back(*i);
			q.first = 3;
			q.second.push_back(*cdbg.begin());
		}

		//Some unitig iterator
		ColoredCDBG<CoreInfo>::iterator i;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//Paths to be compared
		Path p, q;
		
};

class GetCoreDistTest : public ::testing::Test {

protected:

	GetCoreDistTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
		cdbgOpt.k = DEFAULT_TEST_K;
		cdbgOpt.g = DEFAULT_TEST_G;
		cdbgOpt.filename_seq_in.push_back("Test.fa");
		cdbg.build(cdbgOpt);
		cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
		cdbg.buildColors(cdbgOpt);
	}

	//Some unitig iterator
	ColoredCDBG<CoreInfo>::iterator i;
	//Unitig with linked CoreInfo object
	UnitigColorMap<CoreInfo> uni;
	//Colored de Bruijn graph build options
	CCDBG_Build_opt cdbgOpt;
	//Compacted, colored de Bruijn graph with linked CoreInfo objects
	ColoredCDBG<CoreInfo> cdbg;

};

class ExpSucPthsTest : public ::testing::Test {

protected:

	ExpSucPthsTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
		cdbgOpt.k = DEFAULT_TEST_K;
		cdbgOpt.g = DEFAULT_TEST_G;
		cdbgOpt.filename_seq_in.push_back("Test.fa");
		cdbg.build(cdbgOpt);
		cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
		cdbg.buildColors(cdbgOpt);
		i = cdbg.begin();
	};

	//Some unitig iterator
	ColoredCDBG<CoreInfo>::iterator i;
	//Colored de Bruijn graph build options
	CCDBG_Build_opt cdbgOpt;
	//Compacted, colored de Bruijn graph with linked CoreInfo objects
	ColoredCDBG<CoreInfo> cdbg;
	//A unitig list
	list<UnitigColorMap<CoreInfo>> l;
	//A path list to save all results
	list<Path> res;
	//Priority queue of paths to explore
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue;
	
};

#endif