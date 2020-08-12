#ifndef CORE_INFO_HPP
#define CORE_INFO_HPP

#include <list>

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

class CoreInfo : public CCDBG_Data_t<CoreInfo> {

	public:

		CoreInfo(): preBrdg(false), sufBrdg(false) {}
		
		//NOTE: Actually, methods join, serialize and sub should be overloaded here. Since the core info will only be used during run time this is not necessary here.

		//Flag indicating a bridging k-mer at the beginning of a unitig (reference strand orientation)
		bool preBrdg;
		//Flag indicating a bridging k-mer at the end of a unitig (reference strand orientation)
		bool sufBrdg;
		//Increasingly ordered list of core k-mer intervals on a unitig
		list<pair<uint32_t, uint32_t>> coreList;

};

#endif