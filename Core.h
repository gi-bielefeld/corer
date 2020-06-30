#ifndef CORE_HPP
#define CORE_HPP

#include "CoreInfo.h"

#define DEFAULT_CORE_RATIO 0.9
#define MIN_QUORUM 1

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig
void markCore(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& qrm, const uint32_t& dlt);

//This function checks if the given unitig fulfills the given quorum and returns true in this case; false otherwise
//ATTENTION: This function only works correctly if the given unitig only consists of 1 k-mer!
const bool chkQrm(UnitigColorMap<CoreInfo> &u, const uint32_t& q);

#endif