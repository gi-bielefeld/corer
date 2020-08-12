#include "Bridging.h"
#include "Traversal.h"

//This function checks if the distance from the end of a unitig to either its nearest core k-mer or its beginning is already delta
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt){//TODO: Test this function!
	if(!coreList.empty()) return ulen - (coreList.back().second + 1) >= dlt;

	return ulen > dlt;
}

//This function detects and marks all bridging k-mers in the graph
void markBrd(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& dlt){
	//Pointer to current unitig's CoreInfo object
	CoreInfo* cInfo;

	//Iterate over unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get CoreInfo object
		cInfo = i->getData()->getData(*i);
		//Check if last k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the left side (i.e. the distance to the closest core k-mer on this unitig or the unitig's beginning) is not already too large
		if(!cInfo->sufBrdg && (cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && !lCrTooFar(i->len, cInfo->coreList, dlt)){
			//Do BFS on successive unitigs and check if we need to try a BFS on predecessors as well (which is the case only if there is a core k-mer on the current unitig or the BFS on successive unitigs was successful)
			if(!doSucBFS() && cInfo->coreList.empty()) continue;//TODO: Implement this function!
		}

		//Check if first k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the right side (i.e. the distance to the closest core kmer on this unitig or the unitig's end) is not already too large
		if()
	}
}