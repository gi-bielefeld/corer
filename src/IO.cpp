#include <getopt.h>

#include "IO.h"

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& inGfl, string& inCfl, string& outPref, uint32_t& qrm, uint32_t& dlt, size_t& 
	nThrds, bool& oSnps){
	bool iFlGvn = false, cFlGvn = false, oFlGvn = false;
	int option_index = 0, a;

	//Check wheather arguments are given for anything at all
	if(nArgs < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"igraph",   required_argument,  0, 'i'},
        {"cgraph",   required_argument,  0, 'c'},
        {"ograph",   required_argument,  0, 'o'},
        {"quorum",   required_argument,  0, 'q'},
        {"delta",    required_argument,  0, 'd'},
        {"threads",  required_argument,  0, 't'},
        {"snippets", no_argument,        0, 's'},
        {"help",     no_argument,        0, 'h'},
        {0,          0,                  0,  0 }
    };

    //Testing
    // bool qGvn = false, dGvn = false, nTgvn = false, hFlgSt = false;

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'i':
				//Save input graph sequences file
				inGfl = optarg;
				//Note that we have the name of a graph sequences file for reading the input
				iFlGvn = true;
				break;
			case 'c':
				//Save input graph colors file
				inCfl = optarg;
				//Note that we have the name of a graph colors file for reading the input
				cFlGvn = true;
				break;
			case 'o':
				//Save output file prefix
				outPref = optarg;
				//Note that we have a file prefix for writing the output
				oFlGvn = true;
				break;
			case 'q':
				//Testing
				// qGvn = true;

				//A quorum has to be positive
				if(atol(optarg) <= 0 || atol(optarg) > INT32_MAX){
					//Testing
					// cout << "1 Option " << (iFlGvn? "1" : "2") << endl;
					// cout << "2 Option " << (qGvn? "1" : "2") << endl;
					// // cout << "optarg: " << optarg << endl;
					// // cout << "atol(optarg): " << atol(optarg) << endl;
					// // cout << "INT32_MAX: " << INT32_MAX << endl;
					// if(atol(optarg) <= 0){
					// 	cout << "3 Option 2" << endl;
					// } else{
					// 	cout << "3 Option 1" << endl;
					// }
					// if(atol(optarg) > INT32_MAX){
					// 	cout << "4 Option 1" << endl;
					// } else{
					// 	cout << "4 Option 2" << endl;
					// }
					// cout << "5 Option " << (dGvn? "1" : "2") << endl;
					// cout << "8 Option " << (nTgvn? "1" : "2") << endl;
					// cout << "10 Option " << (hFlgSt? "1" : "2") << endl;
					// cout << "11 Option " << (oFlGvn? "1" : "2") << endl;
					// cout << "12 Option " << (oSnps? "1" : "2") << endl;
					// if(iFlGvn){
					// 	cout << "13 Option " << (oFlGvn? "1" : "2") << endl;
					// 	if(cFlGvn){
					// 		if(oFlGvn){
					// 			cout << "17 Option 1" << endl;
					// 		} else{
					// 			cout << "18" << endl;
					// 		}
					// 	}
					// } else if(!oFlGvn){
					// 	cout << "16 Option " << (cFlGvn? "1" : "2") << endl;
					// }
					// if(oFlGvn){
					// 	cout << "14 Option " << (oFlGvn? "1" : "2") << endl;

					// 	if(!iFlGvn && cFlGvn) cout << "17 Option 2" << endl;
					// }
					// cout << "15 Option " << (cFlGvn? "1" : "2") << endl;

					cerr << "ERROR: Quorum value not applicable" << endl;
					return false;
				}

				//Testing
				// cout << "3 Option 1" << endl;
				// cout << "4 Option 2" << endl;

				qrm = atoi(optarg);
				break;
			case 'd':
				//Testing
				// dGvn = true;
				// cout << "optarg: " << optarg << endl;
				// cout << "atol(optarg): " << atol(optarg) << endl;

				//Distance has to be non-negative
				if(atol(optarg) < 0 || atol(optarg) > INT32_MAX){
					//Testing
					// cout << "1 Option " << (iFlGvn? "1" : "2") << endl;
					// cout << "2 Option " << (qGvn? "1" : "2") << endl;
					// cout << "5 Option " << (dGvn? "1" : "2") << endl;
					// if(atol(optarg) < 0){
					// 	cout << "6 Option 1" << endl;
					// } else{
					// 	cout << "6 Option 2" << endl;
					// }
					// if(atol(optarg) > INT32_MAX){
					// 	cout << "7 Option 1" << endl;
					// } else{
					// 	cout << "7 Option 2" << endl;
					// }
					// cout << "8 Option " << (nTgvn? "1" : "2") << endl;
					// cout << "10 Option 2" << endl;
					// cout << "11 Option " << (oFlGvn? "1" : "2") << endl;
					// cout << "12 Option " << (oSnps? "1" : "2") << endl;
					// if(iFlGvn){
					// 	cout << "13 Option " << (oFlGvn? "1" : "2") << endl;
					// 	if(cFlGvn){
					// 		if(oFlGvn){
					// 			cout << "17 Option 1" << endl;
					// 		} else{
					// 			cout << "18" << endl;
					// 		}
					// 	}
					// } else if(!oFlGvn){
					// 	cout << "16 Option " << (cFlGvn? "1" : "2") << endl;
					// }
					// if(oFlGvn){
					// 	cout << "14 Option " << (oFlGvn? "1" : "2") << endl;
					// 	if(!iFlGvn && cFlGvn) cout << "17 Option 2" << endl;
					// }
					// cout << "15 Option " << (cFlGvn? "1" : "2") << endl;

					cerr << "ERROR: Distance value not applicable" << endl;
					return false;
				}

				//Testing
				// cout << "6 Option 1" << endl;
				// cout << "7 Option 2" << endl;

				dlt = atoi(optarg);
				break;
			case 't':
				//Testing
				// nTgvn = true;

				//Number of threads needs to be a positive number
				if(atoi(optarg) < 1){
					//Testing
					// cout << "9 Option 2" << endl;

					cerr << "ERROR: Number of threads not applicable" << endl;
					return false;
				}

				//Testing
				// cout << "9 Option 1" << endl;

				nThrds = atoi(optarg);
				break;
			case 's':
				//Note that we will have to output the core as snippets
				oSnps = true;
				break;
			case 'h':
				//Testing
				// cout << "10 Option 1" << endl;
				// hFlgSt = true;

				return false;
			default:
				break;
		}
	}

	//Testing
	// cout << "1 Option " << (iFlGvn? "1" : "2") << endl;
	// cout << "2 Option " << (qGvn? "1" : "2") << endl;
	// cout << "5 Option " << (dGvn? "1" : "2") << endl;
	// cout << "8 Option " << (nTgvn? "1" : "2") << endl;
	// cout << "10 Option 2" << endl;
	// cout << "11 Option " << (oFlGvn? "1" : "2") << endl;
	// cout << "12 Option " << (oSnps? "1" : "2") << endl;
	// if(iFlGvn){
	// 	cout << "13 Option " << (oFlGvn? "1" : "2") << endl;
	// 	if(cFlGvn){
	// 		if(oFlGvn){
	// 			cout << "17 Option 1" << endl;
	// 		} else{
	// 			cout << "18" << endl;
	// 		}
	// 	}
	// } else if(!oFlGvn){
	// 	cout << "16 Option " << (cFlGvn? "1" : "2") << endl;
	// }
	// if(oFlGvn){
	// 	cout << "14 Option " << (oFlGvn? "1" : "2") << endl;
	// 	if(!iFlGvn && cFlGvn) cout << "17 Option 2" << endl;
	// }
	// cout << "15 Option " << (cFlGvn? "1" : "2") << endl;

	return iFlGvn && cFlGvn && oFlGvn;
}

//This function iterates over the given graph and outputs all core and bridging parts as snippets
void outputSnippets(ColoredCDBG<CoreInfo>& cdbg){
	//Start and end positions of intervals to extract
	size_t start, end;
	//Core interval list iterator
	list<pair<uint32_t, uint32_t>>::const_iterator intvl;

	//Iterate over all unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Check if unitig has no core k-mers
		if(i->getData()->getData(*i)->coreList.empty()){
			//Check if unitig's sequence is marked as bridging
			if(i->getData()->getData(*i)->preBrdg || i->getData()->getData(*i)->sufBrdg){
				//Output the complete sequence
				cout << i->mappedSequenceToString() << endl;
			}
		} else{
			//Get interval list iterator
			intvl = i->getData()->getData(*i)->coreList.begin();

			//Check if unitig's beginning is marked as bridging
			if(i->getData()->getData(*i)->preBrdg){
				//The first substring we have to output starts at the sequence's beginning
				start = 0;
			} else{
				//The first substring we have to output starts at the first core interval's beginning
				start = intvl->first;
			}

			//We assume that our substring will end at the first core interval's end
			end = intvl->second;
			//Move to next interval
			++intvl;

			//Keep outputting substings as long as intervals are left
			while(intvl != i->getData()->getData(*i)->coreList.end()){
				//Output last substring
				cout << i->mappedSequenceToString().substr(start, end - start + cdbg.getK()) << endl;
				//The next substring starts at the current interval
				start = intvl->first;
				//We assume it ends with the current interval
				end = intvl->second;
				//Move to next interval
				++intvl;
			}

			//Check if unitig's suffix is marked as bridging
			if(i->getData()->getData(*i)->sufBrdg){
				//Output last substring reaching to sequence's end
				cout << i->mappedSequenceToString().substr(start) << endl;
			} else{
				//Output last substring reaching to interval's end
				cout << i->mappedSequenceToString().substr(start, end - start + cdbg.getK()) << endl;
			}
		}
	}
}

//This function constructs a graph only consisting of a detected core and writes it to the specified output file
void genCoreGraph(ColoredCDBG<CoreInfo>& cdbg, const string& oName, const size_t& thrds){
	size_t start, end;
	UnitigColorMap<void> ogUni;
	UnitigColorMap<CoreInfo> igUni;
	CCDBG_Build_opt oGBO;
	ColoredCDBG<> oGrph(cdbg.getK(), cdbg.getG());
	list<pair<uint32_t, uint32_t>>::const_iterator intvl;

	//Add sequences to the graph//
	
	//Iterate over all unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Check if unitig has no core k-mers
		if(i->getData()->getData(*i)->coreList.empty()){
			//Check if unitig's sequence is marked as bridging
			if(i->getData()->getData(*i)->preBrdg || i->getData()->getData(*i)->sufBrdg){
				//Add complete sequence to output graph
				oGrph.add(i->referenceUnitigToString());
			}
		} else{
			//Get interval list iterator
			intvl = i->getData()->getData(*i)->coreList.begin();

			//Check if unitig's beginning is marked as bridging
			if(i->getData()->getData(*i)->preBrdg){
				//The first substring we have to output starts at the sequence's beginning
				start = 0;
			} else{
				//The first substring we have to output starts at the first core interval's beginning
				start = intvl->first;
			}

			//We assume that our substring will end at the first core interval's end
			end = intvl->second;
			//Move to next interval
			++intvl;

			//Keep outputting substings as long as intervals are left
			while(intvl != i->getData()->getData(*i)->coreList.end()){
				//Add substring to output graph
				oGrph.add(i->referenceUnitigToString().substr(start, end - start + cdbg.getK()));
				//The next substring starts at the current interval
				start = intvl->first;
				//We assume it ends with the current interval
				end = intvl->second;
				//Move to next interval
				++intvl;
			}

			//Check if unitig's suffix is marked as bridging
			if(i->getData()->getData(*i)->sufBrdg){
				//Add last substring reaching to sequence's end
				oGrph.add(i->referenceUnitigToString().substr(start));
			} else{
				//Add last substring reaching to interval's end
				oGrph.add(i->referenceUnitigToString().substr(start, end - start + cdbg.getK()));
			}
		}
	}

	//Copy colors from old graph//

	//Copy color names
	oGBO.filename_seq_in = cdbg.getColorNames();

	//Try to initialize color matrices and throw an error if neccessary
	if(!oGrph.buildColors(oGBO)){
		cerr << "ERROR: Color matrices of output graph could not be initialized!" << endl;
		exit(EXIT_FAILURE);
	}

	//Iterate over unitigs of output graph
	for(ColoredCDBG<>::iterator i = oGrph.begin(); i != oGrph.end(); ++i){
		//Get current unitig
		ogUni = *i;
		//Set length to 1 k-mer
		ogUni.len = 1;

		//Iterate over unitig's k-mers
		for(uint32_t j = 0; j <= ogUni.size - oGrph.getK(); ++j){
			//Update k-mer's position
			ogUni.dist = j;
			//Search for corresponding k-mer in input graph
			igUni = cdbg.find(Kmer(ogUni.mappedSequenceToString().c_str()));

			//Check if we could not find the k-mer
			if(igUni.isEmpty){
				cerr << "ERROR: K-mer " << ogUni.mappedSequenceToString() << " not found in input graph!" << endl;
				exit(EXIT_FAILURE);
			}

			//Iterate over k-mer's colors
			for(UnitigColors::const_iterator k = igUni.getData()->getUnitigColors(igUni)->begin(igUni); k != igUni.getData()->getUnitigColors(igUni)->end(); ++k)
				//Add color to output graph's k-mer
				ogUni.getData()->getUnitigColors(ogUni)->add(ogUni, k.getColorID());
		}
	}

	//Write graph to file//

	//Try to write graph and throw an error if neccessary
	if(!oGrph.write(oName, thrds, BIFROST_VERBOSE_MODE)){
		cerr << "ERROR: Output graph could not be written to file!" << endl;
		exit(EXIT_FAILURE);
	}
}