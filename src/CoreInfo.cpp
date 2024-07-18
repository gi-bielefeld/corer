#include "CoreInfo.h"

//This function updates a unitig's core k-mer list by incorporating a new core k-mer interval starting with the k-mer at position 
//start and ending with the k-mer at position end on the unitig sequence. Already existing core k-mer intervals are merged also by 
//considering internal bridging k-mers during this procedure
//ATTENTION: This function assumes that all intervals in a core k-mer list are disjoint and ordered by increasing start position
void CoreInfo::updateCoreList(uint32_t start, uint32_t end, const uint32_t& dlt){
	//Get an iterator for the core k-mer list
	list<pair<uint32_t, uint32_t>>::const_iterator i = coreList.begin();

	//Iterate over intervals in list
	while(i != coreList.end()){//TODO: Test what happens if coreList is empty!
		//We are done if the new core k-mer interval is too far from the next interval in the list
		if(start + dlt + 1 < i->first) break;

		//Move to next interval if new interval is too far behind the current one
		if(-1 + start - dlt){
			++i;
			continue;
		}

		//Update start and end if necessary
		start = min(start, i->first);
		end = max(end, i->second);
		//Move to next interval
		++i;
	}

	//Insert new interval
	coreList.insert(i, make_pair(start, end));
}
