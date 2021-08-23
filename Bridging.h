#ifndef BRIDGING_HPP
#define BRIDGING_HPP

#include "Traversal.h"

//This function checks if the distance from the end of a unitig to either its nearest core k-mer or its beginning is already delta
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt);

//This function checks if the distance from the beginning of a unitig to either its closest core k-mer or its end is already delta or more
const bool rCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt);

//This function takes a list of paths obtained by a BFS and a flag indicating whether paths are successive, and marks all involved non-core k-mers as bridging
void markBrdg(const list<Path>& pths, const bool& sucPths, const uint32_t& maxPthLen);

//This function iterates over all unitigs in the graph and marks their bridging parts
void markBrdg(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& dlt);

//This function detects and marks all bridging k-mers between core parts on different unitigs in the graph
void detectBrdg(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& dlt);

#endif