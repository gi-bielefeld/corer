#include "Bridging.h"
#include "Traversal.cpp"

//This function checks if the distance from the end of a unitig to either its closest core k-mer or its beginning is already delta or more
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt){
	if(!coreList.empty()) return ulen - coreList.back().second > dlt;

	return ulen >= dlt;
}

//This function checks if the distance from the beginning of a unitig to either its closest core k-mer or its end is already delta or more
const bool rCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt){
	if(!coreList.empty()) return coreList.front().first >= dlt;

	return ulen >= dlt;
}

//This function takes a list of paths obtained by a BFS and a flag indicating whether paths are successive, and marks all involved non-core k-mers as bridging
void markBrdg(const list<Path>& pths, const bool& sucPths){
	//Iterate over paths
	for(list<Path>::const_iterator p = pths.begin(); p != pths.end(); ++p){
		//The first unitig of a successor's path potenially has a bridging suffix and vice versa
		list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin();

		if(sucPths){
			u->getData()->getData(*u)->sufBrdg = true;
		} else{
			u->getData()->getData(*u)->preBrdg = true;
		}

		++u;

		//Iterate over unitigs in path
		for(; u != p->second.end(); ++u){
			//Check which bridging flag we have to set
			if(sucPths ^ u->strand){
				//Mark k-mers at unitig's end (reference strand orientation) as bridging
				u->getData()->getData(*u)->sufBrdg = true;
			} else{
				//Mark k-mers at unitig's beginning (reference strand orientation) as bridging
				u->getData()->getData(*u)->preBrdg = true;
			}
		}
	}
}

//This function detects and marks all bridging k-mers between core parts on different unitigs in the graph
void detectBrdg(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& dlt){
	//Length of the path from a core k-mer to the current unitigs beginning
	uint32_t exstPthLen;
	//Pointer to current unitig's CoreInfo object
	CoreInfo* cInfo;
	//List of paths leading to reachable core k-mers on successive unitigs
	list<Path> sucPaths;
	//List of paths leading to reachable core k-mers on predecessive unitigs
	list<Path> predPaths;

	//Iterate over unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get CoreInfo object
		cInfo = i->getData()->getData(*i);

		//Testing
		cout << "Current unitig: " << i->mappedSequenceToString() << endl;
		if(cInfo->coreList.empty()){
			cout << "1 Option 1" << endl;
		} else{
			cout << "1 Option 2" << endl;
		}
		if(cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1){
			cout << "2 Option 2" << endl;
		} else{
			cout << "2 Option 1" << endl;
		}
		// if(cInfo->coreList.empty() && cInfo->sufBrdg){
		// 	if(lCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "3 Option 1" << endl;
		// 	} else{
		// 		cout << "3 Option 2" << endl;
		// 	}
		// }
		if(cInfo->coreList.empty() && !cInfo->sufBrdg){
			if(lCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "6 Option 1" << endl;
			} else{
				cout << "6 Option 2" << endl;
			}
		}
		// if(!cInfo->coreList.empty() && cInfo->sufBrdg){
		// 	if(lCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "4 Option 1" << endl;
		// 	} else{
		// 		cout << "4 Option 2" << endl;
		// 	}
		// }
		if(!cInfo->coreList.empty() && !cInfo->sufBrdg){
			if(lCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "7 Option 1" << endl;
			} else{
				cout << "7 Option 2" << endl;
			}
		}
		// if(!cInfo->coreList.empty() && cInfo->coreList.back().second == i->len - 1 && cInfo->sufBrdg){
		// 	if(lCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "5 Option 1" << endl;
		// 	} else{
		// 		cout << "5 Option 2" << endl;
		// 	}
		// }
		if(!cInfo->coreList.empty() && cInfo->coreList.back().second == i->len - 1 && !cInfo->sufBrdg){
			if(lCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "8 Option 1" << endl;
			} else{
				cout << "8 Option 2" << endl;
			}
		}
		// if((cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && cInfo->sufBrdg){
		// 	if(lCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "6 Option 1" << endl;
		// 	} else{
		// 		cout << "6 Option 2" << endl;
		// 	}
		// }
		if((cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && !cInfo->sufBrdg){
			if(lCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "9 Option 1" << endl;
			} else{
				cout << "9 Option 2" << endl;
			}
		}
		if(cInfo->sufBrdg){
			cout << "5 Option 1" << endl;
		} else{
			cout << "5 Option 2" << endl;
		}

		//Check if last k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the left side (i.e. the distance to the closest core k-mer on this unitig or the unitig's beginning) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && !cInfo->sufBrdg && !lCrTooFar(i->len, cInfo->coreList, dlt)){
			//Clear path list
			sucPaths = list<Path>();

			//Do BFS on successive unitigs and check if we need to try a BFS on predecessors as well (which is the case only if there is a core k-mer on the current unitig or the BFS on successive unitigs was successful)
			if(!doSucBFS(*i, (dlt + 1) / 2, sucPaths) && cInfo->coreList.empty()){
				//Testing
				cout << "4 Option 2" << endl;

				continue;
			} else{
				//Testing
				if(sucPaths.empty()){
					cout << "4 Option 1" << endl;
				} else{
					if(cInfo->coreList.empty()){
						cout << "3 Option 2" << endl;
					} else{
						cout << "3 Option 1" << endl;
					}
				}
			}
		}

		//Testing
		if(cInfo->preBrdg){
			cout << "10 Option 1" << endl;
		} else{
			cout << "10 Option 2" << endl;
		}
		// if(cInfo->coreList.empty() && cInfo->preBrdg){
		// 	if(rCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "15 Option 1" << endl;
		// 	} else{
		// 		cout << "15 Option 2" << endl;
		// 	}
		// }
		// cout << "Core list is " << (cInfo->coreList.empty() ? "" : "not ") << "empty" << endl;
		// cout << "Prefix bridging flag is " << (cInfo->preBrdg ? "" : "not ") << "set" << endl;
		if(cInfo->coreList.empty() && !cInfo->preBrdg){
			if(rCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "11 Option 1" << endl;
			} else{
				cout << "11 Option 2" << endl;
			}
		}
		// if(!cInfo->coreList.empty() && cInfo->preBrdg){
		// 	if(rCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "17 Option 1" << endl;
		// 	} else{
		// 		cout << "17 Option 2" << endl;
		// 	}
		// }
		if(!cInfo->coreList.empty() && !cInfo->preBrdg){
			if(rCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "12 Option 1" << endl;
			} else{
				cout << "12 Option 2" << endl;
			}
		}
		// if(!cInfo->coreList.empty() && cInfo->coreList.front().first == 0 && cInfo->preBrdg){
		// 	if(rCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "19 Option 1" << endl;
		// 	} else{
		// 		cout << "19 Option 2" << endl;
		// 	}
		// }
		if(!cInfo->coreList.empty() && cInfo->coreList.front().first == 0 && !cInfo->preBrdg){
			if(rCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "13 Option 1" << endl;
			} else{
				cout << "13 Option 2" << endl;
			}
		}
		// if((cInfo->coreList.empty() || cInfo->coreList.front().first > 0) && cInfo->preBrdg){
		// 	if(rCrTooFar(i->len, cInfo->coreList, dlt)){
		// 		cout << "21 Option 1" << endl;
		// 	} else{
		// 		cout << "21 Option 2" << endl;
		// 	}
		// }
		if((cInfo->coreList.empty() || cInfo->coreList.front().first > 0) && !cInfo->preBrdg){
			if(rCrTooFar(i->len, cInfo->coreList, dlt)){
				cout << "14 Option 1" << endl;
			} else{
				cout << "14 Option 2" << endl;
			}
		}

		//Check if first k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the right side (i.e. the distance to the closest core kmer on this unitig or the unitig's end) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.front().first > 0) && !cInfo->preBrdg && !rCrTooFar(i->len, cInfo->coreList, dlt)){
			//Clear path list
			predPaths = list<Path>();

			//Check if there are no core k-mers on this unitig
			if(cInfo->coreList.empty()){
				//Find path to core k-mer on successive unitigs with minimum length
				exstPthLen = findMinPthLen(sucPaths);

				//Testing
				cout << "15 Option 2" << endl;
			} else{
				//A path has to reach the leftmost core k-mer on this unitig only
				exstPthLen = cInfo->coreList.front().first;

				//Testing
				cout << "15 Option 1" << endl;
			}

			//Testing
			// cout << "Delta for BFS on predecessors is " << min((dlt + 1) / 2, (uint32_t) (dlt - i->len - exstPthLen)) << endl;

			//Do BFS on predecessive unitigs and mark all bridging k-mers if necessary (+1, because otherwise we count in the k-mer we start at)
			if(doPredBFS(*i, min((dlt + 1) / 2, (uint32_t) (dlt - i->len - exstPthLen + 1)), predPaths) || !cInfo->coreList.empty()){
				//Testing
				if(cInfo->coreList.empty()){
					cout << "16 Option 2" << endl;
				}
				if(predPaths.empty()){
					cout << "16 Option 1" << endl;
					cout << "17 Option 1" << endl;
				}

				//Mark k-mers of successive result paths as bridging
				markBrdg(sucPaths, true);
				//Mark k-mers of predecessive result paths as bridging
				markBrdg(predPaths, false);
			} else{
				//Testing
				cout << "17 Option 2" << endl;
			}
		}
	}
}