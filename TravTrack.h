#ifndef TRAV_TRACK_HPP
#define TRAV_TRACK_HPP

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

class TravTrack {

	public:

		TravTrack(): cDist(0) {}
		
		//Flag to indicate the traversal direction on the corresponding unitig
		bool isSucTrav;
		//Distance to the core k-mer at which this traversal was started
		uint32_t cDist;
		//K-mer to track the unitig (and strand) at which the traversal currently is
		Kmer track;

};

#endif