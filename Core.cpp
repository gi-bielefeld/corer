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
			if(chkQrm(uni, qrm)){//TODO Implement this function!
				//Check if there is no current interval yet
				if(l < 0){
					//Set left border
					l = j;
				}

				//Update right border
				r = j;

				//Reset bridging path length if necessary
				if(nBrd > 0) nBrd = 0;
			} else{
				//Increase path length and check if delta is exceeded
				if(++nBrd > dlt){
					//Update right interval border
					r = j - nBrd;
					//Add interval
					uni.getData()->getData(uni)->addInt(make_pair(l, r));//TODO Implement this function!
					//Reset interval borders
					l = -1;
					r = -1;
					//Reset path length
					nBrd = 0;
				}
			}
		}
	}
}