#ifndef COREINFO_TEST_HPP
#define COREINFO_TEST_HPP

#include "CoreTest.h"

class UpdateCoreListTest : public ::testing::Test {

	protected:

		UpdateCoreListTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
			cdbgOpt.filename_seq_in.push_back("Test.fa");
			cdbg.build(cdbgOpt);
			cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
			cdbg.buildColors(cdbgOpt);
			i = cdbg.begin();
		}

		//An iterator to iterate over core intervals
		list<pair<uint32_t, uint32_t>>::const_iterator in;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//Some unitig iterator
		ColoredCDBG<CoreInfo>::iterator i;
};

#endif