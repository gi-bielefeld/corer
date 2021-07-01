#include "Core.h"

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig
void markCore(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& qrm, const uint32_t& dlt){
	uint32_t nBrd;
	int32_t l, r = -1;
	UnitigColorMap<CoreInfo> uni;

	//Iterate over unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get current unitig
		uni = *i;
		//Set length to 1 k-mer
		uni.len = 1;
		//Reset left interval border
		l = -1;
		//Reset bridging path length
		nBrd = 0;

		//Iterate over unitig's k-mers
		for(uint32_t j = 0; j <= uni.size - cdbg.getK(); ++j){
			//Update k-mer's position
			uni.dist = j;

			//Testing
			// if(!uni.mappedSequenceToString().compare("GTGGGTTTTAA") || !uni.mappedSequenceToString().compare("TGGGTTTTAAG") || !uni.mappedSequenceToString().compare("CTTAAAACCCA") || !uni.mappedSequenceToString().compare("TTAAAACCCAC")) cerr << "markCore: Checking k-mer " << uni.mappedSequenceToString() << endl;

			//Check if quorum is fulfilled
			if(chkQrm(uni, qrm)){
				//Testing
				// if(!uni.mappedSequenceToString().compare("GTGGGTTTTAA") || !uni.mappedSequenceToString().compare("TGGGTTTTAAG") || !uni.mappedSequenceToString().compare("CTTAAAACCCA") || !uni.mappedSequenceToString().compare("TTAAAACCCAC")){
				// 	cerr << "markCore: Quorum is fulfilled" << endl;
				// 	cerr << "markCore: Color set is" << endl;
				// 	for(UnitigColors::const_iterator c = uni.getData()->getUnitigColors(uni)->begin(uni); c != uni.getData()->getUnitigColors(uni)->end(); ++c) cerr << "Position:" << (*c).first << " Color:" << (*c).second << endl;
				// }

				//Check if there is no current interval yet and set left border
				if(l < 0) l = j;

				//Update right border
				r = j;

				//Reset bridging path length if necessary
				if(nBrd > 0) nBrd = 0;
			} else{
				//Testing
				// if(!uni.mappedSequenceToString().compare("GTGGGTTTTAA") || !uni.mappedSequenceToString().compare("TGGGTTTTAAG") || !uni.mappedSequenceToString().compare("CTTAAAACCCA") || !uni.mappedSequenceToString().compare("TTAAAACCCAC")){
				// 	cerr << "markCore: Quorum is not fulfilled" << endl;
				// 	cerr << "markCore: Color set is" << endl;
				// 	for(UnitigColors::const_iterator c = uni.getData()->getUnitigColors(uni)->begin(uni); c != uni.getData()->getUnitigColors(uni)->end(); ++c) cerr << "Position:" << (*c).first << " Color:" << (*c).second << endl;
				// }

				//Check if increased path length exceeds delta and an interval has already started
				if(++nBrd >= dlt && l > -1){
					//Testing
					// cout << r << endl;

					//Add interval
					uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));
					//Reset left interval borders
					l = -1;
					//Reset path length
					nBrd = 0;
				}
			}
		}

		//Check if there exists an open interval which was not yet added and add it
		if(l > -1){
			//Testing
			// cout << r << endl;

			uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));
		}

		//Testing
		// if(!i->mappedSequenceToString().compare("CTTAAAACCCAC") || !i->mappedSequenceToString().compare("GTGGGTTTTAAG")){
		// 	cerr << "markCore: Found unitig " << i->mappedSequenceToString() << endl;
		// 	cerr << "markCore: Core list is " << (i->getData()->getData(*i)->coreList.empty() ? "" : "not ") << "empty" << endl;
		// 	cerr << "markCore: Core list is" << endl;
		// 	for(list<pair<uint32_t, uint32_t>>::const_iterator intvl = i->getData()->getData(*i)->coreList.begin(); intvl != i->getData()->getData(*i)->coreList.end(); ++intvl) cerr << '[' << intvl->first << ',' << intvl->second << ']' << endl;
		// 	exit(EXIT_SUCCESS);	
		// }
	}
}

//This function checks if the given unitig fulfills the given quorum and returns true in this case; false otherwise
//ATTENTION: This function only works correctly if the given unitig only consists of 1 k-mer!
const bool chkQrm(UnitigColorMap<CoreInfo> &u, const uint32_t& q){
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