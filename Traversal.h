#ifndef TRAVERSAL_HPP
#define TRAVERSAL_HPP

#include <queue>
#include <list>

#include "CoreInfo.h"
#include "TravTrack.h"

#define DEFAULT_DELTA 50

//A path through the graph is a list of unitigs
using Path = pair<uint32_t, list<UnitigColorMap<CoreInfo>>>;

//A priority queue for TravTrack elements
using TravTrackQueue = priority_queue<TravTrack, vector<TravTrack>, const bool (*)(const TravTrack&, const TravTrack&)>;

//A compare function to prioritize shortest paths
inline const bool prioShrtst(const Path& left, const Path& right){ return left.first > right.first; }

//A compare function to prioritize TravTracks belonging to short paths
inline const bool prioShrtst(const TravTrack& s, const TravTrack& t){ return s.cDist > t.cDist; }//TODO This function still needs to be tested!

//This function calculates an offset on a unitig depending on the given strand
inline const uint32_t calcOff(const uint32_t& refOff, const size_t& uLen, const bool& strand){ return (strand ? refOff : uLen - (refOff + 1)); }

//This function returns the distance from the unitig border to the closest core k-mer depending on border and strand
inline const uint32_t getCoreDist(const neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false>& n, const bool& isSuc){
	return (isSuc ^ n->strand ? n->len - n->getData()->getData(*n)->coreList.back().second : n->getData()->getData(*n)->coreList.front().first + 1);
}

//This function returns the minium path length from a given list of paths. If list is empty UINT32_MAX is returned.
inline const uint32_t findMinPthLen(const list<Path>& pths){
	uint32_t mini = UINT32_MAX;

	//Iterate over paths and update minimum
	for(list<Path>::const_iterator i = pths.begin(); i != pths.end(); ++i) mini = min(mini, i->first);

	return mini;
}

//This function iterates over all paths of the given path list and adds distance information to all unitigs involved representing the distance to the core k-mer at the end of each corresponding path
void addDists(const list<Path>& pthLst, const bool& isSucPth);

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res);

//This function extends the top priority path of the given priority queue on predecessive unitigs. It it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Qtherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expPredPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res);

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths);

//This function performs a BFS of the given depths on all predecessors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doPredBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths);

#endif