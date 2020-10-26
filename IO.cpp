#include <getopt.h>

#include "IO.h"

const bool prsArgs(int& nArgs, char** argList, string& inPref, string& outPref, uint32_t& qrm, uint32_t& dlt, size_t& nThrds, bool& oSnps){
	bool iFlGvn = false, oFlGvn = false;
	int option_index = 0, a;

	//Check wheather arguments are given for anything at all
	if(nArgs < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"igraph",   required_argument,  0, 'i'},
        {"ograph",   required_argument,  0, 'o'},
        {"quorum",   required_argument,  0, 'q'},
        {"delta",    required_argument,  0, 'd'},
        {"threads",  required_argument,  0, 't'},
        {"snippets", no_argument,        0, 's'},
        {"help",     no_argument,        0, 'h'},
        {0,          0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'i':
				//Save input file prefix
				inPref = optarg;
				//Note that we have a file prefix for reading the input
				iFlGvn = true;
				break;
			case 'o':
				//Save output file prefix
				outPref = optarg;
				//Note that we have a file prefix for writing the output
				oFlGvn = true;
				break;
			case 'q':
				//A quorum has to be positive
				if(atoi(optarg) <= 0 || atoi(optarg) > INT32_MAX){
					cerr << "ERROR: Quorum value not applicable" << endl;
					return false;
				}

				qrm = atoi(optarg);
				break;
			case 'd':
				//Distance has to be non-negative
				if(atoi(optarg) < 0 || atoi(optarg) > INT32_MAX){
					cerr << "ERROR: Distance value not applicable" << endl;
					return false;
				}

				dlt = atoi(optarg);
				break;
			case 't':
				//Number of threads needs to be a positive number
				if(atoi(optarg) < 1){
					cerr << "ERROR: Number of threads not applicable" << endl;
					return false;
				}

				nThrds = atoi(optarg);
				break;
			case 's':
				//Note that we will have to output the core as snippets
				oSnps = true;
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return iFlGvn && oFlGvn;
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

//This function
void genCoreGraph(ColoredCDBG<CoreInfo>& cdbg, const string& oName){
	size_t start, end;
	UnitigColorMap<CoreInfo> ogUni, igUni;
	ColoredCDBG<> oGrph(cdbg.getK());
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
				oGrph.add(i->mappedSequenceToString().substr(start));
			} else{
				//Add last substring reaching to interval's end
				oGrph.add(i->mappedSequenceToString().substr(start, end - start + cdbg.getK()));
			}
		}
	}

	//Copy colors from old graph//

	//Iterate over unitigs of output graph
	for(ColoredCDBG<CoreInfo>::iterator i = oGrph.begin(); i != oGrph.end(); ++i){
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

			//Make sure we could find the demanded k-mer
			if(!igUni.empty()){
				//Copy k-mer's colors
				igUni->getData()<-Is that possible?!
			}
		}
	}

	//Write graph to file//
}