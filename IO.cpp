#include <getopt.h>

#include "IO.h"

const bool prsArgs(int& nb_args, char** argList, string& filePref, uint32_t& qrm, uint32_t& dlt){
	bool success = false;
	int option_index = 0, a;

	//Check wheather arguments are given for anything at all
	if(nb_args < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"graph",   required_argument,  0, 'g'},
        {"quorum",  required_argument,  0, 'q'},
        {"delta",   required_argument,  0, 'd'},
        {"help",    no_argument,        0, 'h'},
        {0,         0,                  0,  0 }
    };

   	int b;

    //Parse all parameters given
	while ((a = getopt_long(nb_args, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'g':
				//Save file prefix
				filePref = optarg;
				//A file prefix is all we need to continue
				success = true;
				break;
			case 'q':
				//Testing
				// b = atoi(optarg);
				// cout << "optarg: " << *optarg << endl;

				//A quorum has to be positive
				if(atoi(optarg) <= 0 || atoi(optarg) > UINT32_MAX){
					cerr << "ERROR: Quorum value not applicable" << endl;
					return false;
				}

				//Testing
				// cout << "Do we get here as well?" << endl;
				
				qrm = atoi(optarg);
				break;
			case 'd':
				//Distance has to be non-negative
				if(atoi(optarg) < 0 || atoi(optarg) > UINT32_MAX){
					cerr << "ERROR: Distance value not applicable" << endl;
					return false;
				}

				dlt = atoi(optarg);
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return success;
}