#include "Traversal.h"
#include "CoreInfo.h"

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path>& queue, const uint32_t& dpth, list<Path>& res){
	//Iterator to iterate over successors
	ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it;
	//Variable to store an extended path
	Path extPth;	

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
			if(!suc->getData()->getData(*suc)->coreList.empty() && getCoreDist(suc, true) <= dpth - queue.top().first){
				//Add path to results
				res.push_back(queue.top());
				//Add successor to path
				res.back.second.push_back(*suc);
				//Update path length
				res.back.first += getCoreDist(suc, true);
				//Move on with next successor
				continue;
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

//This function extends the top priority path of the given priority queue on predecessive unitigs. It it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Qtherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expPredPths(priority_queue<Path>& queue, const uint32_t& dpth, list<Path>& res){
	//Iterator to iterate over predecessors
	BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it;
	//Variable to store an extended path
	Path extPth;

	//Get iterator of last unitig in top priority path
	it = queue.top().second.back().getPredecessors();

	//Make sure there are predecessors to iterate over
	if(it.hasPredecessors()){
		//Iterate over predecessors
		for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> pred = it.begin(); pred != it.end(); ++pred){
			//Check if distance to next core k-mer (if known) is too large
			if((pred->strand ? pred->getData()->getData(*pred)->predCoreDist : pred->getData()->getData(*pred)->sucCoreDist) > dpth - queue.top().first){
				//Extending the path on this predecessor does not make sense
				continue;
			}

			//Check if there is a core k-mer on this predecessor and if it is close enough
			if(!pred->getData()->getData(*pred)->coreList.empty() && getCoreDist(pred, false) <= dpth - queue.top().first){
				//Add path to results
				res.push_back(queue.top());
				//Add predecessor to path
				res.back.second.push_back(*pred);
				//Update path length
				res.back.first += getCoreDist(pred, false);<-
			}
		}
	}
}

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(Path, Path)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(pair(0, uniLst));
	//Explore paths
	expSucPths(queue, dpth, resPths);//TODO Implement this function!

	//Check if we could find valid paths
	if(!resPths.empty()){
		//Add core k-mer distance for unitigs on result paths
		addDists(resPths, true);//TODO Implement this function!
		//Report success
		return true;
	} else{
		return false;
	}
}

//This function performs a BFS of the given depths on all predecessors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doPredBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(Path, Path)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(pair(0, uniLst));
	//Explore paths
	expPredPths(queue, dpth, resPths);//TODO Implement this function!

	//Check if we could find valid paths
	if(!resPths.empty()){
		//Add core k-mer distance for unitigs on result paths
		addDists(resPths, false);
		//Report success
		return true;
	} else{
		return false;
	}
}