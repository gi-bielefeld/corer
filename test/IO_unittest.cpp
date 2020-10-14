#include <gtest/gtest.h>

#include "../IO.cpp"
#include "PrsArgsTest.h"

//Tests for function const bool prsArgs(int&, char**, string&, uint32_t&, uint32_t&)//
//	1. Graph prefix is (not) given DONE
//	2. Quorum is (not) given DONE
//	3. Given quorum is (not) positive DONE
//	4. Given quorum is (not) larger than INT32_MAX DONE
//	5. Delta is (not) given DONE
//	6. Given delta is (not) negative DONE 
//	7. Given delta is (not) larger than INT32_MAX DONE
//	8. Number of threads is (not) given DONE
//	9. Number of threads is (not) positive DONE
//	10. Help flag is (not) set DONE

//Tests the function prsArgs with no parameters
TEST_F(PrsArgsTest, NoParams){
	nbArgs = 1;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 1);
	EXPECT_EQ(filePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is given
//	3. Given quorum is negative
//	4. Given quorum is not larger than INT32_MAX
//	5. Delta is given
//	6. Given delta is positive
//	7. Given delta is not larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
TEST_F(PrsArgsTest, NegQrm){
	nbArgs = 7;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");
	argv[1] = strdup("-g");
	argv[2] = strdup("G");
	argv[3] = strdup("-q");
	argv[4] = strdup("-1");
	argv[5] = strdup("-d");
	argv[6] = strdup("1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 7);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is given
//	3. Given quorum is not negative
//	4. Given quorum is larger than INT32_MAX
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
TEST_F(PrsArgsTest, TooLrgQrm){
	nbArgs = 9;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[5] = strdup("-g");
	argv[6] = strdup("G");
	argv[7] = strdup("-q");
	argv[8] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 9);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	5. Delta is given
//	6. Given delta is not positive
//	8. Number of threads is not given
//	10. Help flag is not set
TEST_F(PrsArgsTest, NonPosDlt){
	nbArgs = 13;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[9] = strdup("-g");
	argv[10] = strdup("G");
	argv[11] = strdup("-d");
	argv[12] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 13);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	5. Delta is given
//	7. Given delta is larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
TEST_F(PrsArgsTest, TooLrgDlt){
	nbArgs = 17;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[13] = strdup("-g");
	argv[14] = strdup("G");
	argv[15] = strdup("-d");
	argv[16] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 17);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is set
TEST_F(PrsArgsTest, HlpFlgSet){
	nbArgs = 20;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[17] = strdup("-g");
	argv[18] = strdup("G");
	argv[19] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 20);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is positive
//	10. Help flag is not set
TEST_F(PrsArgsTest, PosNbThrds){
	nbArgs = 24;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[20] = strdup("-g");
	argv[21] = strdup("G");
	argv[22] = strdup("-t");
	argv[23] = strdup("2");

	EXPECT_TRUE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 24);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, 2);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is not positive
//	10. Help flag is not set
TEST_F(PrsArgsTest, NonPosNbThrds){
	nbArgs = 28;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[24] = strdup("-g");
	argv[25] = strdup("G");
	argv[26] = strdup("-t");
	argv[27] = strdup("0");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt, thrds));
	EXPECT_EQ(nbArgs, 28);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
}

//Tests for function void outputSnippets(const ColoredCDBG<CoreInfo>&)//
//	1. A unitig in the graph has (no) core k-mers DONE
//	2. A unitig has no core k-mers and is (not) marked as bridging DONE
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging DONE
//	4. A unitig has one/many core k-mer(s) 1/0
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging DONE

//Tests the function outputSnippets under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is (not) marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging
//Program calls:
//	Bifrost:
//		>Bifrost build -r OutputSnippetsTestGraphSources1.txt -o OutputSnippetsTestGraph1 -k 9 -m 4 -c -v
//	Corer:
//		>../Corer -q 2 -g OutputSnippetsTestGraph1
//Expected output:
//AAGGCAAACAC
//AAAGGCAAA
//GCAAACACA

//Tests the function outputSnippets under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging
//Program calls:
//	Bifrost:
//		>Bifrost build -r OutputSnippetsTestGraphSources2.txt -o OutputSnippetsTestGraph2 -k 9 -m 4 -c -v
//	Corer:
//		>../Corer -q 2 -g OutputSnippetsTestGraph2
//Expected output:
//AAGGCAAACAC
//AAGGCAAAG
//AAAGGCAAA
//GCAAACACA

//Tests the function outputSnippets under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is not marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is not marked as bridging
//Program calls:
//	Bifrost:
//		>Bifrost build -r OutputSnippetsTestGraphSources3.txt -o OutputSnippetsTestGraph3 -k 9 -m 4 -c -v
//	Corer:
//		>../Corer -q 2 -g OutputSnippetsTestGraph3
//Expected output:
//GGCAAACAC
//GCAAACACA

//Tests the function outputSnippets under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is not marked as bridging
//	4. A unitig has many core k-mers
//	5. A unitig has at least one core k-mer and its suffix is not marked as bridging
//Program calls:
//	Bifrost:
//		>Bifrost build -r OutputSnippetsTestGraphSources4.txt -o OutputSnippetsTestGraph4 -k 9 -m 4 -c -v
//	Corer:
//		>../Corer -q 2 -d 1 -g OutputSnippetsTestGraph4
//Expected output:
//AAGGCAAAC
//GGCAAACAC