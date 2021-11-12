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
			//Save distance for next unitig
			lstDst = j->len - j->getData()->getData(*j)->coreList.back().second;
			//Set distance
			j->getData()->getData(*j)->predCoreDist = lstDst;
		} else{
			//Save distance for next unitig
			lstDst = j->getData()->getData(*j)->coreList.front().first + 1;
			//Set distance
			j->getData()->getData(*j)->sucCoreDist = lstDst;
		}

		//Go to next unitig
		++j;

		//Iterate over all remaining unitigs in current path except for the last one
		for(uint32_t c = 2; c < i->second.size(); ++c, ++j){
			//Update saved distance for next unitig
			lstDst += j->len;

			//Check which distance we have to set
			if(isSucPth ^ j->strand){
				//Set distance
				j->getData()->getData(*j)->predCoreDist = min(j->getData()->getData(*j)->predCoreDist, lstDst);
			} else{
				//Set distance
				j->getData()->getData(*j)->sucCoreDist = min(j->getData()->getData(*j)->sucCoreDist, lstDst);
			}
		}
	}
}

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Distance to next core k-mer
	uint32_t coreDist;
	//Get iterator of last unitig in top priority path
	ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getSuccessors();
	//Variable to store an extended path
	Path extPth;

	//Iterate over successors
	for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> suc = it.begin(); suc != it.end(); ++suc){
		//Check if distance to next core k-mer (if known) is too large
		coreDist = (suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist);

		if(coreDist != UINT32_MAX && coreDist > dpth - queue.top().first){
			//Extending the path on this successor does not make sense
			continue;
		}

		//Check if there is a core k-mer on this successor and if it is close enough
		if(!suc->getData()->getData(*suc)->coreList.empty() && getCoreDist(suc, true) <= dpth - queue.top().first){
			//Add path to results
			res.push_back(queue.top());
			//Add successor to path
			res.back().second.push_back(*suc);
			//Update path length
			res.back().first += getCoreDist(suc, true);
			//Move on with next successor
			continue;
		}

		//Check if adding all k-mers of successive unitig to path does not make it too long
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

	//Remove top priority path from queue
	queue.pop();

	//Call function again if queue is not empty
	if(!queue.empty()) expSucPths(queue, dpth, res);
}

//This function extends the top priority path of the given priority queue on predecessive unitigs. It it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Qtherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expPredPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Distance to next core k-mer
	uint32_t coreDist;
	//Get iterator of last unitig in top priority path
	BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getPredecessors();
	//Variable to store an extended path
	Path extPth;

	//Iterate over predecessors
	for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> pred = it.begin(); pred != it.end(); ++pred){
		//Check if distance to next core k-mer (if known) is too large
		coreDist = (pred->strand ? pred->getData()->getData(*pred)->predCoreDist : pred->getData()->getData(*pred)->sucCoreDist);

		if(coreDist != UINT32_MAX && coreDist > dpth - queue.top().first){
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

	//Remove top priority path from queue
	queue.pop();

	//Call function again if queue is not empty
	if(!queue.empty()) expPredPths(queue, dpth, res);
}

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//Distance to next core k-mer
	uint32_t coreDist;
	//Some neighbor iterator
	neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> suc;
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Variable to store an extended path
	Path extPth;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(Path(0, uniLst));

	//Explore paths
	while(!queue.empty()){
		//Get iterator of last unitig in top priority path
		ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getSuccessors();

		//Iterate over successors
		for(suc = it.begin(); suc != it.end(); ++suc){
			//Check if distance to next core k-mer (if known) is too large
			coreDist = (suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist);

			if(coreDist != UINT32_MAX && coreDist > dpth - queue.top().first){
				//Extending the path on this successor does not make sense
				continue;
			}

			//Check if there is a core k-mer on this successor and if it is close enough
			if(!suc->getData()->getData(*suc)->coreList.empty() && getCoreDist(suc, true) <= dpth - queue.top().first){
				//Add path to results
				resPths.push_back(queue.top());
				//Add successor to path
				resPths.back().second.push_back(*suc);
				//Update path length
				resPths.back().first += getCoreDist(suc, true);
				//Move on with next successor
				continue;
			}

			//Check if adding all k-mers of a successive unitig to path does not make it too long
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

		//Remove top priority path from queue
		queue.pop();
	}

	//Check if we could find valid paths
	if(!resPths.empty()){
		//Add core k-mer distance for unitigs on result paths
		addDists(resPths, true);
		//Report success
		return true;
	} else{
		return false;
	}
}

//This function performs a BFS of the given depths on all predecessors of the given unitig. It returns true if a core k-mer could be reached by any path.
const bool doPredBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths){
	//Distance to next core k-mer
	uint32_t coreDist;
	//Some neighbor iterator
	neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> pred;
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Variable to store an extended path
	Path extPth;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(Path(0, uniLst));

	//Explore paths
	while(!queue.empty()){
		//Get iterator of last unitig in top priority path
		BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getPredecessors();

		//Iterate over predecessors
		for(pred = it.begin(); pred != it.end(); ++pred){
			//Check if distance to next core k-mer (if known) is too large
			coreDist = (pred->strand ? pred->getData()->getData(*pred)->predCoreDist : pred->getData()->getData(*pred)->sucCoreDist);

			if(coreDist != UINT32_MAX && coreDist > dpth - queue.top().first){	
				//Extending the path on this predecessor does not make sense
				continue;
			}

			//Check if there is a core k-mer on this predecessor and if it is close enough
			if(!pred->getData()->getData(*pred)->coreList.empty() && getCoreDist(pred, false) <= dpth - queue.top().first){	
				//Add path to results
				resPths.push_back(queue.top());
				//Add predecessor to path
				resPths.back().second.push_back(*pred);
				//Update path length
				resPths.back().first += getCoreDist(pred, false);
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

		//Remove top priority path from queue
		queue.pop();
	}

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

//This function traverses the given graph according to the order dictated by the given priority queue while adding distance information about the closest core k-mers to each
//processed node. 
void annotateDists(ColoredCDBG<CoreInfo>& cdbg, TravTrackQueue& queue, const uint32_t& dlt){
	//A flag stating if a distance on a unitig could be updated
	bool distUpdated;
	//A unitig mapping of the unitig we continue our traversal on
	UnitigColorMap<CoreInfo> uni;
	//A neighbor iterator for the traversal
	neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> nIt;

	//Process all unitigs in the queue
	while(!queue.empty()){
		//Find the unitig from which we have to continue
		uni = cdbg.find(queue.top().track);

		//Check if we are dealing with a successive or predecessive path
		if(queue.top().isSucTrav){
			ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> fIt = uni.getSuccessors();

			//Iterate over successors
			for(nIt = fIt.begin(); nIt != fIt.end(); ++nIt){
				//Yet, nothing has been updated
				distUpdated = false;

				//Ensure that this successor has not already been processed in this direction
				if(nIt->strand && nIt->getData()->getData(*nIt)->predCoreDist == UINT32_MAX){
					//Annotate this successor with distance to closest core k-mer
					nIt->getData()->getData(*nIt)->predCoreDist = queue.top().cDist;
					distUpdated = true;
				} else if(!nIt->strand && nIt->getData()->getData(*nIt)->sucCoreDist == UINT32_MAX){
					//Annotate this successor with distance to closest core k-mer
					nIt->getData()->getData(*nIt)->sucCoreDist = queue.top().cDist;
					distUpdated = true;
				}

				//If this successor contains a core k-mer, was already annotated, has no successors or would lead to a too long traversal, we do not need to continue with it
				if(distUpdated && nIt->getData()->getData(*nIt)->coreList.empty() && nIt->getSuccessors().hasSuccessors() && queue.top().cDist + nIt->len <= dlt){
					//Add this successor to the queue to continue the traversal from here
					queue.push(TravTrack(queue.top().cDist + nIt->len, Kmer(nIt->mappedSequenceToString().substr(0, cdbg.getK()).c_str()), true));
				}
			}
		} else{
			BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> bIt = uni.getPredecessors();
		
			//Iterate over predecessors
			for(nIt = bIt.begin(); nIt != bIt.end(); ++nIt){
				//Yet, nothing has been updated
				distUpdated = false;

				//Ensure that this predecessor has not already been processed in this direction
				if(nIt->strand && nIt->getData()->getData(*nIt)->sucCoreDist == UINT32_MAX){
					//Annotate this predecessor with distance to closest core k-mer
					nIt->getData()->getData(*nIt)->sucCoreDist =queue.top().cDist;
					distUpdated = true;
				} else if(!nIt->strand && nIt->getData()->getData(*nIt)->predCoreDist == UINT32_MAX){
					//Annotate this predecessor with distance to closest core k-mer
					nIt->getData()->getData(*nIt)->predCoreDist = queue.top().cDist;
					distUpdated = true;
				}

				//If this predecessor contains a core k-mer, was already annotated, has no predecessors or would lead to a too long traversal, we do not need to continue with it
				if(distUpdated && nIt->getData()->getData(*nIt)->coreList.empty() && nIt->getPredecessors().hasPredecessors() && queue.top().cDist + nIt->len <= dlt){
					//Add this predecessor to the queue to continue the traversal from here
					queue.push(TravTrack(queue.top().cDist + nIt->len, Kmer(nIt->mappedSequenceToString().substr(0, cdbg.getK()).c_str()), false));
				}
			}
		}

		//Remove processed TravTrack from queue
		queue.pop();
	}
}