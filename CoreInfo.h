#ifndef CORE_INFO_HPP
#define CORE_INFO_HPP

#include <list>

class CoreInfo : public CCDBG_Data_t<CoreInfo> {

	public:

		CoreInfo() {}
		
		//NOTE: Actually, methods join, serialize and sub should be overloaded here. Since the core info will only be used during run time this is not necessary here.

	private:

		//List of core k-mer intervals on a unitig
		list<pair<uint32_t, uint32_t>> coreList;

};

#endif