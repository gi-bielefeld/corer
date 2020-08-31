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