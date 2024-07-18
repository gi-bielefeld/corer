#ifndef CORE_INFO_HPP
#define CORE_INFO_HPP

#include <list>

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

class CoreInfo : public CCDBG_Data_t<CoreInfo> {

	public:

		CoreInfo(): preBrdg(false), sufBrdg(false), sucCoreDist(UINT32_MAX), predCoreDist(UINT32_MAX) {}
		
		//NOTE: Actually, methods join, serialize and sub should be overloaded here. Since the core info will only be used during run time this is not necessary here.

		//Flag indicating a bridging k-mer at the beginning of a unitig (reference strand orientation)
		bool preBrdg;
		//Flag indicating a bridging k-mer at the end of a unitig (reference strand orientation)
		bool sufBrdg;
		//Distance to next core k-mer on this or successive unitigs measured from unitig's beginning (reference strand orientation; 0 if unknown)
		uint32_t sucCoreDist;
		//Distance to next core k-mer on this or predecessive unitigs measured from unitig's end (reference strand orientation; 0 if unknown)
		uint32_t predCoreDist;
		//Increasingly ordered list of core k-mer intervals on a unitig
		list<pair<uint32_t, uint32_t>> coreList;

		//This function updates a unitig's coreList by incorporating a new core k-mer interval starting with the k-mer at position 
		//start and ending with the k-mer at position end on the unitig sequence. Already existing core k-mer intervals are merged
		//also by considering internal bridging k-mers during this procedure
		void updateCoreList(uint32_t start, uint32_t end, const uint32_t& dlt);
};

#endif