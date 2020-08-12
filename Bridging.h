#ifndef BRIDGING_HPP
#define BRIDGING_HPP

//This function checks if the distance from the end of a unitig to either its nearest core k-mer or its beginning is already delta
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt);

#endif