#include <gtest/gtest.h>

#include "../IO.cpp"
#include "IOtest.h"
#include "../Bridging.h"

//Tests for function const bool prsArgs(int&, char**, string&, uint32_t&, uint32_t&)//
//	1. Input graph prefix is (not) given DONE
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
//	13. Input graph prefix is set and output graph prefix is (not) set (as well) DONE
//	14. Output graph prefix is set and input graph prefix is (not) set (as well) DONE

//Tests the function prsArgs with no parameters
TEST_F(PrsArgsTest, NoParams){
	nbArgs = 1;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 1);
	EXPECT_EQ(filePref, "");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is given
//	3. Given quorum is negative
//	4. Given quorum is not larger than INT32_MAX
//	5. Delta is given
//	6. Given delta is positive
//	7. Given delta is not larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, NegQrm){
	nbArgs = 9;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("Corer");
	argv[1] = strdup("-i");
	argv[2] = strdup("G");
	argv[3] = strdup("-q");
	argv[4] = strdup("-1");
	argv[5] = strdup("-d");
	argv[6] = strdup("1");
	argv[7] = strdup("-o");
	argv[8] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 9);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is given
//	3. Given quorum is not negative
//	4. Given quorum is larger than INT32_MAX
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, TooLrgQrm){
	nbArgs = 11;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[5] = strdup("-i");
	argv[6] = strdup("G");
	argv[7] = strdup("-q");
	argv[8] = strdup("2147483648");
	argv[9] = strdup("-o");
	argv[10] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 11);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is not given
//	5. Delta is given
//	6. Given delta is not positive
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, NonPosDlt){
	nbArgs = 15;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[9] = strdup("-i");
	argv[10] = strdup("G");
	argv[11] = strdup("-d");
	argv[12] = strdup("-1");
	argv[13] = strdup("-o");
	argv[14] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 15);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is not given
//	5. Delta is given
//	7. Given delta is larger than INT32_MAX
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, TooLrgDlt){
	nbArgs = 19;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[13] = strdup("-i");
	argv[14] = strdup("G");
	argv[15] = strdup("-d");
	argv[16] = strdup("2147483648");
	argv[17] = strdup("-o");
	argv[18] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 19);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, HlpFlgSet){
	nbArgs = 22;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[17] = strdup("-i");
	argv[18] = strdup("G");
	argv[19] = strdup("-h");
	argv[20] = strdup("-o");
	argv[21] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 22);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is positive
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, PosNbThrds){
	nbArgs = 26;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[20] = strdup("-i");
	argv[21] = strdup("G");
	argv[22] = strdup("-t");
	argv[23] = strdup("2");
	argv[24] = strdup("-o");
	argv[25] = strdup("O");

	EXPECT_TRUE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 26);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "O");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, 2);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is given
//	9. Number of threads is not positive
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is not requested
//	13. Input graph prefix is set and output graph prefix is set as well
//	14. Output graph prefix is set and input graph prefix is set as well
TEST_F(PrsArgsTest, NonPosNbThrds){
	nbArgs = 32;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[26] = strdup("-i");
	argv[27] = strdup("G");
	argv[28] = strdup("-t");
	argv[29] = strdup("0");
	argv[30] = strdup("-o");
	argv[31] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 32);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is (not) given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is not given
//	12.	Unitig snippet output is requested
//	13. Input graph prefix is set and output graph prefix is not set
TEST_F(PrsArgsTest, NoOut){
	nbArgs = 33;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[30] = strdup("-i");
	argv[31] = strdup("G");
	argv[32] = strdup("-s");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 33);
	EXPECT_EQ(filePref, "G");
	EXPECT_EQ(oFilePref, "");
	EXPECT_EQ(qrm, 0);
	EXPECT_EQ(dlt, DEFAULT_DELTA);
	EXPECT_EQ(thrds, DEFAULT_NB_THREADS);
	EXPECT_EQ(!OUTPUT_CORE_SNIPPETS_DEFAULT, oSnps);
}

//Tests the function prsArgs under the following conditions
//	1. Input graph prefix is not given
//	2. Quorum is not given
//	5. Delta is not given
//	8. Number of threads is not given
//	10. Help flag is not set
//	11. Output graph prefix is given
//	12.	Unitig snippet output is requested
//	14. Output graph prefix is set and input graph prefix is not set
TEST_F(PrsArgsTest, NoGph){
	nbArgs = 35;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[33] = strdup("-o");
	argv[34] = strdup("O");

	EXPECT_FALSE(prsArgs(nbArgs, argv, filePref, oFilePref, qrm, dlt, thrds, oSnps));
	EXPECT_EQ(nbArgs, 35);
	EXPECT_EQ(filePref, "");
	EXPECT_EQ(oFilePref, "O");
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

	//Testing
	cout << "Input graph is:" << endl;
	for(ColoredCDBG<CoreInfo>::iterator j = cdbg.begin(); j != cdbg.end(); ++j) cout << j->mappedSequenceToString() << endl;
	UnitigColorMap<CoreInfo> uni = cdbg.find(Kmer("TGTGTTTGC"));
	cout << "K-mer: TGTGTTTGC could " << (uni.isEmpty ? "" : "not ") << "be found" << endl;
	uni = cdbg.find(Kmer("GCAAACACA"));
	cout << "K-mer: GCAAACACA could " << (uni.isEmpty ? "" : "not ") << "be found" << endl;

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
	EXPECT_EQ("AAAGGCAAACACA", i->mappedSequenceToString());
	c = i->getData()->getUnitigColors(*i)->begin(*i);
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(0, c.getColorID());
	++c;
	EXPECT_EQ(0, c.getKmerPosition());
	EXPECT_EQ(1, c.getColorID());
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
	EXPECT_EQ(4, c.getKmerPosition());
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