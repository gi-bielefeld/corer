#include "Core.h"
#include "Traversal.h"
#include "CoreInfo.cpp"

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

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig. 
//It outputs a priority queue containing all unitigs with core parts as TravTracks for a graph traversal.
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

		//If there remains an interval not added yet add it
		if(l > -1) uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));

		if(!uni.getData()->getData(uni)->coreList.empty()){
			//If a unitig has no successors we do not need a traversal on them
			if(uni.getSuccessors().hasSuccessors()) queue.push(TravTrack(i->len - r, Kmer(uni.mappedSequenceToString().c_str()), 
				true));

			//If a unitig has no predecessors we do not need a traversal on them
			if(uni.getPredecessors().hasPredecessors()) queue.push(TravTrack(uni.getData()->getData(uni)->coreList.front().first + 1
				, Kmer(uni.mappedSequenceToString().c_str()), false));
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

//This function iterates over the given list of sequences, searches for all k-mer matches in the graph and marks them as core if 
//present. Close core k-mers on the same unitig are also connected by bridging k-mers if possible.
void markKmers(ColoredCDBG<CoreInfo>& cdbg, vector<string>& seqList, const bool& inexact, const uint32_t& dlt){
	vector<pair<size_t, UnitigColorMap<CoreInfo>>> foundKmers;

	//Iterate over all sequences
	for(vector<string>::const_iterator s = seqList.begin(); s != seqList.end(); ++s){
		//Search for sequence matches
		foundKmers = cdbg.searchSequence(*s, true, inexact, inexact, inexact, true);

		//Iterate over findings
		for(vector<pair<size_t, UnitigColorMap<CoreInfo>>>::const_iterator m = foundKmers.begin(); m != foundKmers.end(); ++m){
			//Mark matched k-mer and potentially bridging ones as core
			m->second.getData()->getData(m->second)->updateCoreList(m->second.dist, m->second.dist + m->second.len - 1, dlt);
		}
	}
}
