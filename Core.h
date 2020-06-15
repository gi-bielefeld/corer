#ifndef CORE_HPP
#define CORE_HPP

#define DEFAULT_CORE_RATIO 0.9

//This function traverses the graph marking all core k-mers and all bridging k-mers connecting core k-mers within the same unitig
void markCore(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& qrm, const uint32_t& dlt)

#endif