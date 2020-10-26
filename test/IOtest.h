#ifndef IO_TEST_HPP
#define IO_TEST_HPP

#include "../Traversal.h"
#include "CoreTest.h"

class PrsArgsTest : public ::testing::Test {

	protected:

		PrsArgsTest(): nbArgs(0), qrm(0), dlt(DEFAULT_DELTA), thrds(DEFAULT_NB_THREADS), oSnps(OUTPUT_CORE_SNIPPETS_DEFAULT) {}

		// void TearDown() override {
		// 	for(uint16_t i = 0; i < nbArgs; ++i) free(&argv[i]);
		// }

		//Unitig snippet output flag
		bool oSnps;
		//Number of command line arguments
		int nbArgs;
		//Arrays with command line arguments
		char** argv;
		//Graph file prefix
		string filePref;
		//Output file prefix
		string oFilePref;
		//The core quorum
		uint32_t qrm;
		//Maximum path length
		uint32_t dlt;
		//Number of threads
		size_t thrds;
};

class GenCoreGraphTest : public ::testing::Test {

	protected:

		GenCoreGraphTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G), crGrph(DEFAULT_TEST_K, DEFAULT_TEST_G), qrm(2), oName("testCoreGraph"), thrds(DEFAULT_NB_THREADS) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
			cdbgOpt.filename_ref_in.push_back("Test.fa");
		}

		//Some unitig iterator
		ColoredCDBG<>::iterator i;
		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph representing the pangenome's core
		ColoredCDBG<> crGrph;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//Some UnitigColors iterator
		UnitigColors::const_iterator c;
		//The core quorum
		uint32_t qrm;
		//Name of output graph
		string oName;
		//A string vector iterator
		vector<string>::const_iterator n;
		//String vector to store loaded color names
		vector<string> cNms;
		//Number of threads
		size_t thrds;
};

#endif