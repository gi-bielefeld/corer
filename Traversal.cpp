#include "Traversal.h"

//This function iterates over all paths of the given path list and adds distance information to all unitigs involved representing the distance to the core k-mer at the end of each corresponding path
void addDists(const list<Path>& pthLst, const bool& isSucPth){
	uint32_t lstDst;
	list<UnitigColorMap<CoreInfo>>::const_reverse_iterator j;

	//Iterate over all paths
	for(list<Path>::const_iterator i = pthLst.begin(); i != pthLst.end(); ++i){
		//Get reverse iterator for path
		j = i->second.rbegin();

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

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Get iterator of last unitig in top priority path
	ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getSuccessors();;
	//Variable to store an extended path
	Path extPth;

	//Iterate over successors
	for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> suc = it.begin(); suc != it.end(); ++suc){
		//Check if distance to next core k-mer (if known) is too large
		if((suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist) > dpth - queue.top().first){
			//Testing
			if(suc->strand){
				cout << "3 Option 1" << endl;
			} else{
				cout << "3 Option 2" << endl;
			}

			//Extending the path on this successor does not make sense
			continue;
		}

		//Testing
		if(suc->strand){
			cout << "4 Option 1" << endl;
		} else{
			cout << "4 Option 2" << endl;
		}
		bool wasCloseEnough = false;

		//Check if there is a core k-mer on this successor and if it is close enough
		if(!suc->getData()->getData(*suc)->coreList.empty() && getCoreDist(suc, true) <= dpth - queue.top().first){
			//Testing
			cout << "6 Option 1" << endl;
			wasCloseEnough = true;
			cout << "5 Option 1" << endl;

			//Add path to results
			res.push_back(queue.top());
			//Add successor to path
			res.back().second.push_back(*suc);
			//Update path length
			res.back().first += getCoreDist(suc, true);
			//Move on with next successor
			continue;
		}

		//Testing
		if(suc->getData()->getData(*suc)->coreList.empty()){
			cout << "5 Option 2" << endl;
		} else{
			cout << "5 Option 1" << endl;

			if(!wasCloseEnough) cout << "6 Option 2" << endl;
		}

		//Check if adding all k-mers of successive unitig to path does not make it too long
		if(queue.top().first + suc->len < dpth){
			//Testing
			cout << "7 Option 2" << endl;

			//Get path
			extPth = queue.top();
			//Add current successor to path
			extPth.second.push_back(*suc);
			//Update path length
			extPth.first += suc->len;
			//Insert path to queue
			queue.push(extPth);
		}
		
		//Testing
		cout << "7 Option 1" << endl;
	}

	//Remove top priority path from queue
	queue.pop();

	//Call function again if queue is not empty
	if(!queue.empty()){
		//Testing
		cout << "8 Option 2" << endl;

		expSucPths(queue, dpth, res);
	} else{
		//Testing
		cout << "8 Option 1" << endl;
	}
}

//This function extends the top priority path of the given priority queue on predecessive unitigs. It it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Qtherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expPredPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Get iterator of last unitig in top priority path
	BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getPredecessors();;
	//Variable to store an extended path
	Path extPth;

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
				res.back().second.push_back(*pred);
				//Update path length
				res.back().first += getCoreDist(pred, false);
				//Move on with next predecessor
				continue;
			}

			//Check if adding all k-mers of predecessive unitig to path does not make it too long
			if(queue.top().first + pred->len < dpth){
				//Get path
				extPth = queue.top();
				//Add current predecessor to path
				extPth.second.push_back(*pred);
				//Update path length
				extPth.first += pred->len;
				//Insert path to queue
				queue.push(extPth);
			}
		}
	}

	//Remove top priority path from queue
	queue.pop();

	//Call function again if queue is not empty
	if(!queue.empty()) expPredPths(queue, dpth, res);
}

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(Path(0, uniLst));
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
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(Path(0, uniLst));
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