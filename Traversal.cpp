#include "Traversal.h"

//Testing
extern bool report = false;
uint32_t callCnt = 0;

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
		for(uint32_t c = 2; c < i->second.size(); ++c, ++j){
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
		}
	}
}

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Get iterator of last unitig in top priority path
	ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getSuccessors();
	//Variable to store an extended path
	Path extPth;

	//Testing
	uint32_t c = 0;
	// if(!queue.top().second.back().mappedSequenceToString().compare("GCAAGCCCGACCTTCGCCGGTCCCTGGGCATCGTGCTGCAGGATGTGAATCTGTTCACCGGAACGGTCATGGACAACATTCGCTACGGCAAGCTCGACGCCTCCGACGAGGAGTGCATCGAAGCCGCCAAGCTCACGAATGCCGA")){
	if(report){
		report = true;
		cout << "expSucPths: Iterating over successors" << endl;
		cout << "expSucPths: callCnt: " << ++callCnt << endl;
		UnitigColorMap<CoreInfo> tmp = queue.top().second.back().getGraph()->find(Kmer("GCAAGCCCGACCTTCGCCGGTCCCTGGGCATCGTGCTGCAGGATGTGAATCTGTTCACCGGAACGGTCATGGACAACATTCGCTACGGCAAGCTCGACGCCTCCGACGAGGAGTGCATCGAAGCCGCCAAGCTCACGAATGCCGA"));
		cout << "tmp is " << (tmp.isEmpty ? "" : "not ") << "empty" << endl;
		cout << "tmp: " << tmp.mappedSequenceToString() << endl;
		it = queue.top().second.back().getSuccessors();
		cout << "expSucPths: Outgoing unitig is " << queue.top().second.back().mappedSequenceToString() << endl;
		cout << "expSucPths: The unitig does " << (queue.top().second.back().getSuccessors().hasSuccessors() ? "" : "not ") << "have successors" << endl;
		cout << "We do get here" << endl;
		it = queue.top().second.back().getSuccessors();
	}

	//Iterate over successors
	for(neighborIterator<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> suc = it.begin(); suc != it.end(); ++suc){
		//Testing
		if(report)
			cout << "expSucPths: Check distance to next core k-mer if possible" << endl;

		//Check if distance to next core k-mer (if known) is too large
		if((suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist) > dpth - queue.top().first){
			//Extending the path on this successor does not make sense
			continue;
		}

		//Testing
		if(report)
			cout << "expSucPths: Check if core k-mer is present on current node" << endl;

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

		//Testing
		if(report)
			cout << "expSucPths: Check if adding all k-mers of this node makes path too long" << endl;

		//Check if adding all k-mers of successive unitig to path does not make it too long
		if(queue.top().first + suc->len < dpth){
			//Get path
			extPth = queue.top();

			//Testing
			if(report){
				cout << "expSucPths: The top priority path is of length " << extPth.first << " and consists of " << extPth.second.size() << " nodes" << endl;
				cout << "expSucPths: Last nodes in the path are" << endl;
				c = 0;
				for(list<UnitigColorMap<CoreInfo>>::const_reverse_iterator r = extPth.second.rbegin(); r != extPth.second.rend(); ++r){
					cout << r->mappedSequenceToString() << endl;

					// if(++c > 3)
					// 	break;
				}
			}

			//Add current successor to path
			extPth.second.push_back(*suc);
			//Update path length
			extPth.first += suc->len;

			//Testing
			if(report){
				cout << "expSucPths: The updated path is of length " << extPth.first << " and consists of " << extPth.second.size() << " nodes" << endl;
				cout << "expSucPths: Last nodes in the path are" << endl;
				c = 0;
				for(list<UnitigColorMap<CoreInfo>>::const_reverse_iterator r = extPth.second.rbegin(); r != extPth.second.rend(); ++r){
					cout << r->mappedSequenceToString() << endl;

					if(++c > 3)
						break;
				}
			}

			//Insert path to queue
			queue.push(extPth);
		}
	}

	//Testing
	if(report){
		cout << "expSucPths: Iteration over successors done" << endl;
		cout << "Number of paths in queue before removal: " << queue.size() << endl;
	}


	//Remove top priority path from queue
	queue.pop();

	//Testing
	if(report){
		cout << "Number of paths in queue after removal: " << queue.size() << endl;
		cout << "expSucPths: After removal of the top priority path, the new top priority path is of length " << queue.top().first << " and consists of " << queue.top().second.size() << " nodes" << endl;
		cout << "expSucPths: Last nodes in the path are" << endl;
		c = 0;
		for(list<UnitigColorMap<CoreInfo>>::const_reverse_iterator r = queue.top().second.rbegin(); r != queue.top().second.rend(); ++r){
			cout << r->mappedSequenceToString() << endl;

			// if(++c > 3)
			// 	break;
		}	
	}

	//Call function again if queue is not empty
	if(!queue.empty()) expSucPths(queue, dpth, res);
}

//This function extends the top priority path of the given priority queue on predecessive unitigs. It it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Qtherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expPredPths(priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)>& queue, const uint32_t& dpth, list<Path>& res){
	//Get iterator of last unitig in top priority path
	BackwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getPredecessors();
	//Variable to store an extended path
	Path extPth;

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

	//Testing
	// cout << "doSucBFS: Start to explore paths" << endl;

	//Explore paths
	expSucPths(queue, dpth, resPths);

	//Testing
	// cout << "doSucBFS: Paths explored" << endl;

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
	//A list to store the first path
	list<UnitigColorMap<CoreInfo>> uniLst;
	//Priority queue to store explored paths
	priority_queue<Path, vector<Path>, const bool (*)(const Path&, const Path&)> queue(prioShrtst);

	//Add first unitig to list
	uniLst.push_back(orig);
	//Add initial path to priority queue
	queue.push(Path(0, uniLst));
	//Explore paths
	expPredPths(queue, dpth, resPths);

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