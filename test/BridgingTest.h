#ifndef BRIDGING_TEST_HPP
#define BRIDGING_TEST_HPP

#include <list>

#include "Traversal.h"
#include "CoreTest.h"
#include "Bridging.h"

using namespace std;

class CrTooFarTest : public ::testing::Test {

	protected:

		CrTooFarTest(): d(42), len(42) {}

		//Maximum distance from unitig border
		uint32_t d;
		//Length of the hypothetical unitig
		size_t len;
		//List of core k-mer intervals
		list<pair<uint32_t, uint32_t>> ints;
};

class MarkBrdgTest : public ::testing::Test {

	protected:

		MarkBrdgTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
		}

		//Some unitig iterator
		ColoredCDBG<CoreInfo>::iterator i;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//A unitig
		UnitigColorMap<CoreInfo> uni;
		//A list of paths
		list<Path> l;
};

class MarkBrdgTest1 : public ::testing::Test {

	protected:

		MarkBrdgTest1(): dlt(4), cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
			cdbgOpt.filename_ref_in.push_back("Test.fa");
			cdbgOpt.filename_ref_in.push_back("Test_color7.fa");
			cdbg.build(cdbgOpt);
			cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
			cdbg.buildColors(cdbgOpt);
		}

		//The used value of delta
		uint32_t dlt;
		//Some unitig iterator
		ColoredCDBG<CoreInfo>::iterator i;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//A priority queue for bridging k-mer detection
		TravTrackQueue queue;
};

class DetectBrdgTest : public ::testing::Test {

	protected:

		DetectBrdgTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
		}

		//Some unitig iterator
		ColoredCDBG<CoreInfo>::iterator i;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
};

#endif