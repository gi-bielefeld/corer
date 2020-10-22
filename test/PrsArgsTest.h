#ifndef PRS_ARGS_TEST_HPP
#define PRS_ARGS_TEST_HPP

#include "../Traversal.h"

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

#endif