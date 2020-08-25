#ifndef BRIDGING_HPP
#define BRIDGING_HPP

//A path through the graph is a list of unitigs
using Path = list<pair<uint32_t, UnitigColorMap<CoreInfo>>>;

//This function checks if the distance from the end of a unitig to either its nearest core k-mer or its beginning is already delta
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt);

//This function checks if the distance from the beginning of a unitig to either its closest core k-mer or its end is already delta or more
const bool rCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt);

#endif