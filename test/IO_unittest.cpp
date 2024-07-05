#include <gtest/gtest.h>

#include "../src/IO.cpp"
#include "IOtest.h"
#include "../src/Bridging.h"

//Tests for function const bool prsArgs(int& nArgs, char** argList, string& inGfl, string& inCfl, string& outPref, uint32_t& qrm, 
//uint32_t& dlt, size_t& nThrds, bool& oSnps)//
//	1. Input graph sequence file is (not) given DONE
//	2. Quorum is (not) given DONE
//	3. Given quorum is (not) positive DONE
//	4. Given quorum is (not) larger than INT32_MAX DONE
//	5. Delta is (not) given DONE
//	6. Given delta is (not) negative DONE 
//	7. Given delta is (not) larger than INT32_MAX DONE
//	8. Number of threads is (not) given DONE
//	9. Number of threads is (not) positive DONE
//	10. Help flag is (not) set DONE
//	11. Output graph prefix is (not) given DONE
//	12.	Unitig snippet output is (not) requested DONE
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is (not) set as well DONE
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is (not) set as well DONE
//	15. Input graph color file is (not) given DONE
//	16. Input graph sequence file and output graph prefix are not set, and input graph color file is (not) set (as well) DONE
//	17. Input graph color file and output graph prefix are set, and input graph sequence file is (not) set (as well) DONE
//	18. Input graph files are given and output graph prefix is not set 0

//Tests the function prsArgs with no parameters
TEST_F(PrsArgsTest, NoParams){
	nbArgs = 1;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 1);
	EXPECT_EQ(gSeqFile, "");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is given
//	3. Given quorum is not positive
//	4. Given quorum is not larger than INT32_MAX
//	5. Delta is given
//	6. Given delta is positive
//	7. Given delta is not larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, NegQrm){
	nbArgs = 9;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");
	argv[1] = strdup("-i");
	argv[2] = strdup("G");
	argv[3] = strdup("-d");
	argv[4] = strdup("1");
	argv[5] = strdup("-o");
	argv[6] = strdup("O");
	argv[7] = strdup("-q");
	argv[8] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 9);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, 1);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is given
//	3. Given quorum is positive
//	4. Given quorum is larger than INT32_MAX
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, TooLrgQrm){
	nbArgs = 15;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[9] = strdup("-i");
	argv[10] = strdup("G");
	argv[11] = strdup("-o");
	argv[12] = strdup("O");
	argv[13] = strdup("-q");
	argv[14] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 15);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is given
//	6. Given delta is not positive
//  7. Given delta is not larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, NonPosDlt){
	nbArgs = 21;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[15] = strdup("-i");
	argv[16] = strdup("G");
	argv[17] = strdup("-o");
	argv[18] = strdup("O");
	argv[19] = strdup("-d");
	argv[20] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 21);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is given
//  6. Given delta is not negative
//	7. Given delta is larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, TooLrgDlt){
	nbArgs = 27;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[21] = strdup("-i");
	argv[22] = strdup("G");
	argv[23] = strdup("-o");
	argv[24] = strdup("O");
	argv[25] = strdup("-d");
	argv[26] = strdup("2147483648");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 27);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, HlpFlgSet){
	nbArgs = 32;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[27] = strdup("-i");
	argv[28] = strdup("G");
	argv[29] = strdup("-o");
	argv[30] = strdup("O");
	argv[31] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 32);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is positive
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, PosNbThrds){
	nbArgs = 38;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[32] = strdup("-i");
	argv[33] = strdup("G");
	argv[34] = strdup("-t");
	argv[35] = strdup("2");
	argv[36] = strdup("-o");
	argv[37] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 38);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, 2);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is not positive
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is set as well
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, NonPosNbThrds){
	nbArgs = 44;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[38] = strdup("-i");
	argv[39] = strdup("G");
	argv[40] = strdup("-o");
	argv[41] = strdup("O");
	argv[42] = strdup("-t");
	argv[43] = strdup("0");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 44);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is not given
//	12.	Unitig snippet output is requested
//	13. Input graph sequence file is set, input graph color file is not set and output graph prefix is not set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, NoOut){
	nbArgs = 47;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[44] = strdup("-i");
	argv[45] = strdup("G");
	argv[46] = strdup("-s");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 47);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(!OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is not given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	14. Output graph prefix is set, input graph color file is not set and input graph sequence file is not set as well
//	15. Input graph color file is not given
TEST_F(PrsArgsTest, NoGph){
	nbArgs = 49;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[47] = strdup("-o");
	argv[48] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 49);
	EXPECT_EQ(gSeqFile, "");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//  1. Input graph sequence file is not given
//  2. Quorum is not given
//  5. Delta is not given
//  8. Number of threads is not given
//  10. Help flag is not set
//  11. Output graph prefix is not given
//  12. Unitig snippet output is not requested
//  15. Input graph color file is given
//  16. Input graph sequence file and output graph prefix are not set, and input graph color file is set
TEST_F(PrsArgsTest, CflGvn){
	nbArgs = 51;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[49] = strdup("-c");
	argv[50] = strdup("testColorFile.color.bfg");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 51);
	EXPECT_EQ(gSeqFile, "");
	EXPECT_EQ(gColFile, "testColorFile.color.bfg");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//  1. Input graph sequence file is not given
//  2. Quorum is given
//	3. Given quorum is positive
//	4. Given quorum is not larger than INT32_MAX
//  5. Delta is not given
//  8. Number of threads is not given
//  10. Help flag is not set
//  11. Output graph prefix is not given
//  12. Unitig snippet output is not requested
//  15. Input graph color file is not given
//	16. Input graph sequence file and output graph prefix are not set, and input graph color file is not set as well
TEST_F(PrsArgsTest, NoInOut){
	nbArgs = 53;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[51] = strdup("-q");
	argv[52] = strdup("1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 53);
	EXPECT_EQ(gSeqFile, "");
	EXPECT_EQ(gColFile, "");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 1);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	15. Input graph color file is given
//	17. Input graph color file and output graph prefix are set, and input graph sequence file is set as well
TEST_F(PrsArgsTest, AllInOut){
	nbArgs = 59;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[53] = strdup("-i");
	argv[54] = strdup("G");
	argv[55] = strdup("-c");
	argv[56] = strdup("testColorFile.color.bfg");
	argv[57] = strdup("-o");
	argv[58] = strdup("O");

	EXPECT_TRUE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 59);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "testColorFile.color.bfg");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is not given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	15. Input graph color file is given
//	17. Input graph color file and output graph prefix are set, and input graph sequence file is not set
TEST_F(PrsArgsTest, NoSeqG){
	nbArgs = 63;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[59] = strdup("-c");
	argv[60] = strdup("testColorFile.color.bfg");
	argv[61] = strdup("-o");
	argv[62] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 63);
	EXPECT_EQ(gSeqFile, "");
	EXPECT_EQ(gColFile, "testColorFile.color.bfg");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph sequence file is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is not given
//	12.	Unitig snippet output is not requested
//	15. Input graph color file is given
//	18. Input graph files are given and output graph prefix is not set
TEST_F(PrsArgsTest, ClrInNoOut){
	nbArgs = 67;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[63] = strdup("-c");
	argv[64] = strdup("testColorFile.color.bfg");
	argv[65] = strdup("-i");
	argv[66] = strdup("G");

	EXPECT_FALSE(prsArgs(nbArgs, argv, gSeqFile, gColFile, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 67);
	EXPECT_EQ(gSeqFile, "G");
	EXPECT_EQ(gColFile, "testColorFile.color.bfg");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests for function void genCoreGraph(ColoredCDBG<CoreInfo>&, const string&, const size_t&)//
//	1. A unitig in the graph has (no) core k-mers DONE
//	2. A unitig has no core k-mers and is (not) marked as bridging DONE
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging DONE
//	4. A unitig has one/many core k-mer(s) DONE
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging DONE

//Tests the function genCoreGraph under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is (not) marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging
TEST_F(GenCoreGraphTest, HsCr){
	cdbgOpt.filename_ref_in.push_back("Test_color7.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	markCore(cdbg, qrm, DEFAULT_DELTA);
	detectBrdg(cdbg, DEFAULT_DELTA);
	genCoreGraph(cdbg, oName, thrds);

	ASSERT_TRUE(crGrph.read(oName + GFA_FILE_ENDING, oName + COLOR_FILE_ENDING, thrds));
	EXPECT_EQ(2, crGrph.getNbColors());
	cNms = crGrph.getColorNames();
	n = cNms.begin();
	EXPECT_EQ("Test.fa", *n);
	++n;
	EXPECT_EQ("Test_color7.fa", *n);
	EXPECT_EQ(1, crGrph.size());
	i = crGrph.begin();
	EXPECT_EQ("TGTGTTTGCCTTT", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(1, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(2, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(3, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(4, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(4, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
}

//Tests the function genCoreGraph under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is (not) marked as bridging
TEST_F(GenCoreGraphTest, TwPths){
	cdbgOpt.filename_ref_in.push_back("Test_color1.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	markCore(cdbg, qrm, DEFAULT_DELTA);
	detectBrdg(cdbg, DEFAULT_DELTA);
	genCoreGraph(cdbg, oName, thrds);

	ASSERT_TRUE(crGrph.read(oName + GFA_FILE_ENDING, oName + COLOR_FILE_ENDING, thrds));
	EXPECT_EQ(2, crGrph.getNbColors());
	cNms = crGrph.getColorNames();
	n = cNms.begin();
	EXPECT_EQ("Test.fa", *n);
	++n;
	EXPECT_EQ("Test_color1.fa", *n);
	EXPECT_EQ(3, crGrph.size());
	i = crGrph.begin();
	EXPECT_EQ("TGTGTTTGCCTT", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(1, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(2, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(3, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(2, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
	++i;
	EXPECT_EQ("AAGGCAAAG", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
	++i;
	EXPECT_EQ("AAAGGCAAA", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
}

//Tests the function genCoreGraph under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is not marked as bridging
//	4. A unitig has one core k-mer
//	5. A unitig has at least one core k-mer and its suffix is not marked as bridging
TEST_F(GenCoreGraphTest, NoPrefSufBrdg){
	cdbgOpt.filename_ref_in.push_back("Test_color9.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	markCore(cdbg, qrm, DEFAULT_DELTA);
	detectBrdg(cdbg, DEFAULT_DELTA);
	genCoreGraph(cdbg, oName, thrds);

	ASSERT_TRUE(crGrph.read(oName + GFA_FILE_ENDING, oName + COLOR_FILE_ENDING, thrds));
	EXPECT_EQ(2, crGrph.getNbColors());
	cNms = crGrph.getColorNames();
	n = cNms.begin();
	EXPECT_EQ("Test.fa", *n);
	++n;
	EXPECT_EQ("Test_color9.fa", *n);
	EXPECT_EQ(1, crGrph.size());
	i = crGrph.begin();
	EXPECT_EQ("TGTGTTTGCCT", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(1, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(2, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(2, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
}

//Tests the function genCoreGraph under the following conditions
//	1. A unitig in the graph has (no) core k-mers
//	2. A unitig has no core k-mers and is not marked as bridging
//	3. A unitig has at least one core k-mer and its prefix is not marked as bridging
//	4. A unitig has many core k-mers
//	5. A unitig has at least one core k-mer and its suffix is not marked as bridging
TEST_F(GenCoreGraphTest, SmlDlt){
	cdbgOpt.filename_ref_in.push_back("Test_color8.fa");
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	cdbg.buildColors(cdbgOpt);
	markCore(cdbg, qrm, 1);
	detectBrdg(cdbg, 1);
	genCoreGraph(cdbg, oName, thrds);

	ASSERT_TRUE(crGrph.read(oName + GFA_FILE_ENDING, oName + COLOR_FILE_ENDING, thrds));
	EXPECT_EQ(2, crGrph.getNbColors());
	cNms = crGrph.getColorNames();
	n = cNms.begin();
	EXPECT_EQ("Test.fa", *n);
	++n;
	EXPECT_EQ("Test_color8.fa", *n);
	EXPECT_EQ(2, crGrph.size());
	i = crGrph.begin();
	EXPECT_EQ("AAGGCAAAC", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
	++i;
	EXPECT_EQ("GGCAAACAC", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
	++c;
	EXPECT_EQ(i->getData()->getUnitigColors(*i)->end(), c);
}

//Tests for function void outputSnippets(const ColoredCDBG<CoreInfo>&)//
//	1. A unitig in the graph has (no) core k-mers DONE
//	2. A unitig has no core k-mers and is (not) marked as bridging DONE
//	3. A unitig has at least one core k-mer and its prefix is (not) marked as bridging DONE
//	4. A unitig has one/many core k-mer(s) DONE
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