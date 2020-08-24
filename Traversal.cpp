#include "Traversal.h"
#include "CoreInfo.h"

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<pair<uint32_t, Path>>& queue, const uint32_t& dpth, list<Path>& res){
	//Iterator to iterate over successors
	ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it;
	//Variable to store an extended path
	pair<uint32_t, Path> extPth;	

	//Get iterator of last unitig in top priority path
	it = queue.top().second.back().getSuccessors();

	//Make sure there are successors to iterate over
	if(it.hasSuccessors()){
		//Iterate over successors
		for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> suc = it.begin(); suc != it.end(); ++suc){
			//Check if distance to next core k-mer (if known) is too large
			if((suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist) > dpth - queue.top().first){
				//Extending the path on this successor does not make sense
				continue;
			}
			//Check if there is a core k-mer on this successor and if it is close enough
			if(!suc->getData()->getData(*suc)->coreList.empty() && (suc->strand ? suc->getData()->getData(*suc)->coreList.front().first + 1 : suc->len - suc->getData()->getData(*suc)->coreList.back().second) <= dpth - queue.top().first){
				//Add path to results
				res.push_back(queue.top());
				//Add successor to path
				res.back.push_back(*suc);
				//Update path length
				res.back.first += (suc->strand ? suc->getData()->getData(*suc)->coreList.front().first + 1 : suc->len - suc->getData()->getData(*suc)->coreList.back().second);
			}
			//Check if adding all k-mers on successive unitig to path does not make it too long
			if(queue.top().first + suc->len < dpth){
				//Get path
				extPth = queue.top();
				//Add current successor to path
				extPth.second.push_back(*suc);
				//Update path length
				extPth.first += suc->len;
				//Insert path to queue
				queue.push(extPth);
			}
		}
	}

	//Remove top priority path from queue
	queue.pop();

	//Call function again if queue is not empty
	if(!queue.empty()) expSucPths(queue, dpth, res);
}

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//Compare function to prioritize shortest paths
	auto cmp = [](pair<uint32_t, Path> left, pair<uint32_t, Path> right) { return left.first > right.first; };
	//Priority queue to store explored paths
	priority_queue<pair<uint32_t, Path>, vector<pair<uint32_t, Path>>, decltype(cmp)> queue(cmp);

	//Add initial path to priority queue
	queue.push(pair(0, orig));
	//Explore paths
	expSucPths(queue, dpth, resPths);//TODO Implement this function!

	//Check if we could find valid paths
	if(!resPths.empty()){
		//Add core k-mer distance information for unitigs on result paths
		addDists(resPths);//TODO Implement this function!
		//Report success
		return true;
	} else{
		return false;
	}
}

//This function performs a BFS of the given depths on all predecessors of the given unitig. It returns true if a core k-mer could be reached.
const bool doPredBFS(){

}