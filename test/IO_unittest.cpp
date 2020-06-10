#include <gtest/gtest.h>

#include "../IO.h"
#include "PrsArgsTest.h"

//Tests for function const bool prsArgs(int&, char**, string&, uint32_t&, uint32_t&){//
//Graph prefix is (not) given DONE
//Quorum is (not) given DONE
//Given quorum is (not) positive DONE
//Given quorum is (not) larger than INT32_MAX DONE
//Delta is (not) given DONE
//Given delta is (not) negative DONE 
//Given delta is (not) larger than INT32_MAX DONE
//Help flag is (not) set DONE

//Tests the function prsArgs with no parameters
TEST_F(PrsArgsTest, NoParams){
	nbArgs = 1;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 1);
	EXPECT_EQ(filePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is given
//	3. Given quorum is negative
//	4. Given quorum is not larger than INT32_MAX
//	5. Delta is given
//	6. Given delta is positive
//	7. Given delta is not larger than INT32_MAX
//	8. Help flag is not set
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

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 7);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is given
//	3. Given quorum is not negative
//	4. Given quorum is larger than INT32_MAX
//	5. Delta is not given
//	6. Help flag is not set
TEST_F(PrsArgsTest, TooLrgQrm){
	nbArgs = 11;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");
	argv[7] = strdup("-g");
	argv[8] = strdup("G");
	argv[9] = strdup("-q");
	argv[10] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 11);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	3. Delta is given
//	4. Given delta is not positive
//	5. Help flag is not set
TEST_F(PrsArgsTest, NonPosDlt){
	nbArgs = 15;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[11] = strdup("-g");
	argv[12] = strdup("G");
	argv[13] = strdup("-d");
	argv[14] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 15);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	3. Delta is given
//	4. Given delta is larger than INT32_MAX
//	5. Help flag is not set
TEST_F(PrsArgsTest, TooLrgDlt){
	nbArgs = 19;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[15] = strdup("-g");
	argv[16] = strdup("G");
	argv[17] = strdup("-d");
	argv[18] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 19);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}

//Tests the function prsArgs under the following conditions
//	1. Graph prefix is given
//	2. Quorum is not given
//	3. Delta is not given
//	4. Help flag is set
TEST_F(PrsArgsTest, HlpFlgSet){
	nbArgs = 22;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[19] = strdup("-g");
	argv[20] = strdup("G");
	argv[21] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, qrm, dlt));
	EXPECT_EQ(nbArgs, 22);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
}