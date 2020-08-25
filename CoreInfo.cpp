#include "CoreInfo.h"

//This function iterates over all paths of the given path list and adds distance information to all unitigs involved representing the distance to the core k-mer at the end of each corresponding path
void addDists(const list<Path>& pthLst, const bool& isSucPth){
	uint32_t lstDst;
	Path::const_iterator j;

	//Iterate over all paths
	for(list<Path>::const_iterator i = pthLst.second.begin(); i != pthLst.second.end(); ++i){
		//Get reverse iterator for path
		j = i->rbegin();

		//Check which distance we have to set
		if(isSucPth ^ j->strand){
			//Set distance
			j->getData()->getData(*j)->predCoreDist = j->len - j->getData()->getData(*j)->coreList.back().second;
			//Save distance for next unitig
			lstDst = j->getData()->getData(*j)->predCoreDist;
		} else{
			//Set distance
			j->getData()->getData(*j)->sucCoreDist = j->getData()->getData(*j)->coreList.front().first + 1;
			//Save distance for next unitig
			lstDst = j->getData()->getData(*j)->sucCoreDist;
		}

		//Go to next unitig
		++j;

		//Iterate over all remaining unitigs in current path except for the last one
		while(j->getData()->getData(*j)->coreList.empty()){
			//Update saved distance for next unitig
			lstDst += j->len;

			//Check which distance we have to set
			if(isSucPth ^ j->strand){
				//Set distance
				j->getData()->getData(*j)->predCoreDist = lstDst;
			} else{
				//Set distance
				j->getData()->getData(*j)->sucCoreDist = lstDst;
			}

			//Go to next unitig
			++j;
		}
	}
}