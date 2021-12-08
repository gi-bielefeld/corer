#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#define COLOR_FILE_ENDING ".bfg_colors"
#define GFA_FILE_ENDING ".gfa"

//This function checks if the given unitig fulfills the given quorum and returns true in this case; false otherwise
//ATTENTION: This function only works correctly if the given unitig only consists of 1 k-mer!
const bool chkQrm(UnitigColorMap<void> &u, const uint32_t& q){
	uint32_t cnt;
	size_t curID, lstID;
	int32_t allwdToMs;

	//If the maximum color id (+1, since ids start at 0) is already smaller than our quorum it can never be fulfilled
	if(q > u.getData()->getUnitigColors(u)->colorMax(u) + 1) return false;

	//Calculate how many colors we are allowed to miss before it is clear that we cannot fulfill the quorum anymore
	allwdToMs = u.getData()->getUnitigColors(u)->colorMax(u) + 1 - q;
	//Get a color iterator
	UnitigColors::const_iterator i = u.getData()->getUnitigColors(u)->begin(u);
	//Set color counter
	cnt = 0;
	//Set current id to maximum
	curID = SIZE_MAX;
	//Initialize last id
	lstID = 0;

	//Iterate over unitig's colors until too many colors have already been missed
	while(i != u.getData()->getUnitigColors(u)->end() && allwdToMs >= 0){
		//Update current id
		curID = i.getColorID();

		//Last id is always smaller than the current id except we are dealing with the first color
		if(lstID < curID){
			//Decrement by the number of colors we have skipped between last and current color if they are not consecutive
			allwdToMs -= curID - lstID - 1;
		} else{
			//Decrement by the number of colors we have skipped starting from 0
			allwdToMs -= curID;
		}

		//Check if the current color is present at this unitig
		if(u.getData()->getUnitigColors(u)->contains(u, curID)){
			//Increment counter and check if we are done
			if(++cnt == q) return true;
		} else{
			//Decrement number of colors we are still allowed to miss
			--allwdToMs;
		}

		//Update last id
		lstID = curID;
		//Move to the next color
		++i;
	}

	return false;
}

//This program loads a Bifrost graph and counts the number of core k-mers
int main(int argc, char **argv){
	uint32_t qrm;
	size_t nbCoreKmers = 0, nbKmers = 0;
	UnitigColorMap<void> uni;

	//Check if there are enough input parameters
	if(argc < 3){
		cerr << "ERROR: Not enough input parameters given\nWe need the graph file prefix and a quorum" << endl;
		return 1;
	}

	//Load graph
	string graphFilePref = argv[1];
	ColoredCDBG<> cdbg = ColoredCDBG<>();
	if(!cdbg.read((graphFilePref + GFA_FILE_ENDING).c_str(), (graphFilePref + COLOR_FILE_ENDING).c_str())){
		cout << "ERROR: Graph could not be loaded" << endl;
		return 1;
	}

	//Load quorum
	//A quorum has to be positive
	if(atoi(argv[2]) <= 0 || atoi(argv[2]) > INT32_MAX){
		cerr << "ERROR: Quorum value not applicable" << endl;
		return 1;
	}
	qrm = atoi(argv[2]);

	//Iterate over unitigs
	for(ColoredCDBG<>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		uni = *i;
		//Set length to 1 k-mer
		uni.len = 1;

		//Iterate over unitig's k-mers
		for(uint32_t j = 0; j <= uni.size - cdbg.getK(); ++j){
			//Update k-mer's position
			uni.dist = j;
			++nbKmers;

			//Check if quorum is fulfilled
			if(chkQrm(uni, qrm)) ++nbCoreKmers;
		}
	}

	cout << "Number of core " << cdbg.getK() << "-mers in this graph is " << nbCoreKmers << " of " << nbKmers << " " << cdbg.getK() << "-mers in total" << endl;

	return 0;
}