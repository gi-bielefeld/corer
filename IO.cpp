#include <getopt.h>

#include "IO.h"

const bool prsArgs(int& nArgs, char** argList, string& filePref, uint32_t& qrm, uint32_t& dlt, size_t& nThrds){
	bool success = false;
	int option_index = 0, a;

	//Check wheather arguments are given for anything at all
	if(nArgs < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"graph",   required_argument,  0, 'g'},
        {"quorum",  required_argument,  0, 'q'},
        {"delta",   required_argument,  0, 'd'},
        {"threads", required_argument,  0, 't'},
        {"help",    no_argument,        0, 'h'},
        {0,         0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'g':
				//Save file prefix
				filePref = optarg;
				//A file prefix is all we need to continue
				success = true;
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
			case 'h':
				return false;
			default:
				break;
		}
	}

	return success;
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