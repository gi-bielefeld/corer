#include "Traversal.h"

//This function iterates over all paths of the given path list and adds distance information to all unitigs involved representing the distance to the core k-mer at the end of each corresponding path
void addDists(const list<Path>& pthLst, const bool& isSucPth){
	uint32_t lstDst;
	list<UnitigColorMap<CoreInfo>>::const_reverse_iterator j;

	//Testing
	// bool appeared = false;
	cout << "1 Option " << (pthLst.empty() ? "1" : "2") << endl;

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

		//Testing
		cout << "2 Option " << (i->first > 2 ? "1" : "2") << endl;
		if(isSucPth){
			cout << "3 Option " << (!j->strand ? "1" : "2") << endl;
		} else{
			cout << "4 Option " << (!j->strand ? "1" : "2") << endl;
		}
		// if(!j->mappedSequenceToString().compare("GCTGTTTAACA") || !j->mappedSequenceToString().compare("TGTTAAACAGC")){
		// 	cerr << "addDists: Distance set for unitig " << j->mappedSequenceToString() << endl;
		// 	cerr << "addDists: The whole path is" << endl;
		// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator k = i->second.begin(); k != i->second.end(); ++k) cerr << k->mappedSequenceToString() << endl;
		// 	cerr << "addDists: Distance is " << j->getData()->getData(*j)->sucCoreDist << endl;
		// 	appeared = true;
		// }

		//Go to next unitig
		++j;

		//Iterate over all remaining unitigs in current path except for the last one
		for(uint32_t c = 2; c < i->second.size(); ++c, ++j){
			//Update saved distance for next unitig
			lstDst += j->len;

			//Testing
			if(isSucPth){
				cout << "5 Option " << (!j->strand ? "1" : "2") << endl;

				if(!j->strand){
					if(j->getData()->getData(*j)->predCoreDist != UINT32_MAX && j->getData()->getData(*j)->predCoreDist != lstDst) cout << "7 Option 1" << endl;
				} else{
					if(j->getData()->getData(*j)->sucCoreDist != UINT32_MAX && j->getData()->getData(*j)->sucCoreDist != lstDst) cout << "7 Option 2" << endl;
				}
			} else{
				cout << "6 Option " << (!j->strand ? "1" : "2") << endl;

				if(!j->strand){
					if(j->getData()->getData(*j)->sucCoreDist != UINT32_MAX && j->getData()->getData(*j)->sucCoreDist != lstDst) cout << "8 Option 2" << endl;
				} else{
					if(j->getData()->getData(*j)->predCoreDist != UINT32_MAX && j->getData()->getData(*j)->predCoreDist != lstDst) cout << "8 Option 1" << endl;
				}
			}

			//Check which distance we have to set
			if(isSucPth ^ j->strand){
				//Set distance
				j->getData()->getData(*j)->predCoreDist = min(j->getData()->getData(*j)->predCoreDist, lstDst);
			} else{
				//Set distance
				j->getData()->getData(*j)->sucCoreDist = min(j->getData()->getData(*j)->sucCoreDist, lstDst);
			}

			//Testing
			// if(!j->mappedSequenceToString().compare("GCTGTTTAACA") || !j->mappedSequenceToString().compare("TGTTAAACAGC")){
			// 	cerr << "addDists: Distance set for unitig " << j->mappedSequenceToString() << endl;
			// 	cerr << "addDists: The whole path is" << endl;
			// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator k = i->second.begin(); k != i->second.end(); ++k) cerr << k->mappedSequenceToString() << endl;
			// 	cerr << "addDists: Distance is " << j->getData()->getData(*j)->sucCoreDist << endl;
			// 	appeared = true;
			// }
		}
	}

	//Testing
	// if(appeared){
	// 	cerr << "addDists: We do get here" << endl;
	// 	exit(EXIT_SUCCESS);
	// }
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

	//Testing
	// string path[] = {"GATGCTGTTTA"};
	// uint32_t counter = 0;
	// if(!orig.mappedSequenceToString().compare("ATGCTGTTTAA")){
	// 	cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// 	// exit(EXIT_SUCCESS);
	// }
	// if(!orig.mappedSequenceToString().compare("GCTGTTTAACA")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("CTGTTTAACAGG")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("CTGTTTAACAC")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("TGTTTAACACACAA")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("TAACACACAAT")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("AACACACAATAA")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("CACACAATAAA")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;
	// if(!orig.mappedSequenceToString().compare("ACACAATAAAA")) cerr << "doSucBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;

	//Explore paths
	while(!queue.empty()){
		//Get iterator of last unitig in top priority path
		ForwardCDBG<DataAccessor<CoreInfo>, DataStorage<CoreInfo>, false> it = queue.top().second.back().getSuccessors();

		//Iterate over successors
		for(suc = it.begin(); suc != it.end(); ++suc){
			//Testing
			// if(!orig.mappedSequenceToString().compare("ATGCTGTTTAA")) cerr << "doSucBFS: We have found successor " << suc->mappedSequenceToString() << endl;

			//Check if distance to next core k-mer (if known) is too large
			coreDist = (suc->strand ? suc->getData()->getData(*suc)->sucCoreDist : suc->getData()->getData(*suc)->predCoreDist);

			if(coreDist != UINT32_MAX && coreDist > dpth - queue.top().first){
				//Testing
				// if(!orig.mappedSequenceToString().compare("ATGCTGTTTAA")){
				// 	cerr << "doSucBFS: Distance to next known core k-mer is too large" << endl;
				// 	cerr << "doSucBFS: Distance is " << suc->getData()->getData(*suc)->sucCoreDist << " and allowed is only " << dpth - queue.top().first << endl;
				// 	++counter;
				// }

				//Extending the path on this successor does not make sense
				continue;
			}

			//Check if there is a core k-mer on this successor and if it is close enough
			if(!suc->getData()->getData(*suc)->coreList.empty() && getCoreDist(suc, true) <= dpth - queue.top().first){
				//Testing
				// if(!orig.mappedSequenceToString().compare("ATGCTGTTTAA")){
				// 	cerr << "doSucBFS: There is a core k-mer on this unitig which is close enough" << endl;
				// 	++counter;
				// }

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

			//Testing
			// if(!orig.mappedSequenceToString().compare("ATGCTGTTTAA")){
			// 	cerr << "doSucBFS: We continue with this unitig" << endl;
			// 	++counter;
			// }
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

	//Testing
	// if(!orig.mappedSequenceToString().compare("GTTTAACAGGTACG")) cerr << "doPredBFS: BFS on unitig " << orig.mappedSequenceToString() << endl;

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