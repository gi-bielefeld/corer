#ifndef TRAVERSAL_HPP
#define TRAVERSAL_HPP

#define DEFAULT_DELTA 50

//A compare function to prioritize shortest paths
inline const bool prioShrtst(Path left, Path right){ return left.first > right.first; }

//This function returns the distance from the unitig border to the closest core k-mer depending on border and strand
inline const uint32_t getCoreDist(const neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false>& n, const bool& isSuc){
	return (isSuc ^ n->strand ? n->len - n->getData()->getData(*n)->coreList.back().second : n->getData()->getData(*n)->coreList.front().first + 1);
}

//This function returns the minium path length from a given list of paths. If list is empty UINT32_MAX is returned.
inline const uint32_t findMinPthLen(const list<Path>& pths){
	uint32_t min = UINT32_MAX;

	//Iterate over paths and update minimum
	for(list<Path>::const_iterator i = pths.begin(); i != pths.end(); ++i) min = min(min, i->first);

	return min;
}

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path>& queue, const uint32_t& dpth, list<Path>& res);

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths);

//This function performs a BFS of the given depths on all predecessors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doPredBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths);

#endif