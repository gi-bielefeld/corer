#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Core.h"
#include "Traversal.h"
#include "IO.cpp"

int main(int argc, char **argv){
	uint32_t qrm = 0;
	uint32_t dlt = DEFAULT_DELTA;
	string gFilePref;

	//Parse arguments
	if(!prsArgs(argc, argv, gFilePref, qrm, dlt)){
		//Display help message
		dspHlp();
		return EXIT_FAILURE;
	}

	cout << "Loaded parameters" << endl;
	cout << "Parameters are gFilePref: " << gFilePref << " qrm: " << qrm << " dlt: " << dlt << endl;

	return EXIT_SUCCESS;
}