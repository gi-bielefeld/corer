#include "Core.h"

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig
void markCore(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& qrm, const uint32_t& dlt){
	uint32_t nBrd;
	int32_t l, r;
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
				//Testing
				cout << "Quorum is fulfilled for a k-mer" << endl;

				//Check if there is no current interval yet
				if(l < 0){
					//Testing
					cout << "A new interval has to be started" << endl;

					//Set left border
					l = j;
				} else{
					//Testing
					cout << "A new interval has not to be started" << endl;
				}

				//Update right border
				r = j;

				//Reset bridging path length if necessary
				if(nBrd > 0){
					//Testing
					cout << "Bridging path's length has to be reseted" << endl;

					nBrd = 0;
				} else{
					//Testing
					cout << "Bridging path's length has not to be reseted" << endl;
				}
			} else{
				//Testing
				cout << "Quorum is not fulfilled for a k-mer" << endl;

				//Increase path length and check if delta is exceeded
				if(++nBrd > dlt){
					//Testing
					cout << "Delta is exceeded" << endl;

					//Update right interval border
					r = j - nBrd;
					//Add interval
					uni.getData()->getData(uni)->coreList.push_back(make_pair(l, r));
					//Reset interval borders
					l = -1;
					r = -1;
					//Reset path length
					nBrd = 0;
				} else{
					//Testing
					cout << "Delta is not exceeded" << endl;
				}
			}
		}
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
	while(allwdToMs >= 0){
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
		i.nextColor();
	}

	return false;
}