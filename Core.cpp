#include "Core.h"
#include "Traversal.h"

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

			//Check if quorum is fulfilled
			if(chkQrm(uni, qrm)){
				//Check if there is no current interval yet and set left border
				if(l < 0) l = j;

				//Update right border
				r = j;

				//Reset bridging path length if necessary
				if(nBrd > 0) nBrd = 0;
			} else{
				//Check if increased path length exceeds delta and an interval has already started
				if(++nBrd >= dlt && l > -1){
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
		if(l > -1) uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));
	}
}

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig. It outputs a priority queue containing all unitigs
//with core parts as TravTracks for a graph traversal.
TravTrackQueue detectCore(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& qrm, const uint32_t& dlt){
	uint32_t nBrd;
	int32_t l, r = -1;
	UnitigColorMap<CoreInfo> uni;
	TravTrackQueue queue(prioShrtst);

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

			//Check if quorum is fulfilled
			if(chkQrm(uni, qrm)){
				//Testing
				cout << "1 Option 1" << endl;

				//Check if there is no current interval yet and set left border
				if(l < 0){
					//Testing
					cout << "2 Option 1" << endl;

					l = j;
				} else {
					//Testing
					cout << "2 Option 2" << endl;
				}

				//Update right border
				r = j;

				//Reset bridging path length if necessary
				if(nBrd > 0){
					//Testing
					cout << "3 Option 1" << endl;

					nBrd = 0;
				} else{
					//Testing
					cout << "3 Option 2" << endl;
				}
			} else{
				//Testing
				cout << "1 Option 2" << endl;

				//Check if increased path length exceeds delta and an interval has already started
				if(++nBrd >= dlt && l > -1){
					//Testing
					cout << "4 Option 1" << endl;

					//Add interval
					uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));
					//Reset left interval borders
					l = -1;
					//Reset path length
					nBrd = 0;
				} else{
					//Testing
					cout << "4 Option 2" << endl;
				}
			}
		}

		//If there remains an interval not added yet add it
		if(l > -1) uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));

		if(!uni.getData()->getData(uni)->coreList.empty()){
			//Testing
			cout << "5 Option 1" << endl;
			cout << "i->len - r: " << i->len - r << endl;
			cout << "l + 1: " << l + 1 << endl;

			//Push a TravTrack to the queue (forward and backward traversal)
			queue.push(TravTrack(i->len - r, Kmer(uni.mappedSequenceToString().c_str()), true));
			queue.push(TravTrack(uni.getData()->getData(uni)->coreList.front().first + 1, Kmer(uni.mappedSequenceToString().c_str()), false));
		} else{
			//Testing
			cout << "5 Option 2" << endl;
		}
	}

	return queue;
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