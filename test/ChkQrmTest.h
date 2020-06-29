#ifndef CHK_QRM_TEST_HPP
#define CHK_QRM_TEST_HPP

#define DEFAULT_TEST_K 9
#define DEFAULT_TEST_G 4

class ChkQrmTest : public ::testing::Test {

	protected:

		ChkQrmTest(): cdbg(DEFAULT_TEST_K, DEFAULT_TEST_G), qrm(2) {
			cdbgOpt.k = DEFAULT_TEST_K;
			cdbgOpt.g = DEFAULT_TEST_G;
			cdbgOpt.filename_seq_in.push_back("Test.fa");
		}

		//Colored de Bruijn graph build options
		CCDBG_Build_opt cdbgOpt;
		//Compacted, colored de Bruijn graph with linked CoreInfo objects
		ColoredCDBG<CoreInfo> cdbg;
		//The unitig that shall be checked
		UnitigColorMap<CoreInfo> u;
		//The core quorum
		uint32_t qrm;
};

#endif