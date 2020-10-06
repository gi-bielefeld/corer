#ifndef BRIDGING_TEST_HPP
#define BRIDGING_TEST_HPP

#include <list>

#include "Traversal.h"

using namespace std;

class CrTooFarTest : public ::testing::Test {

	protected:

		CrTooFarTest(): d(42), len(42) {}

		//Maximum distance from unitig border
		uint32_t d;
		//Length of the hypothetical unitig
		size_t len;
		//List of core k-mer intervals
		list<pair<uint32_t, uint32_t>> ints;
};

class MarkBrdgTest : public ::testing::Test {

	protected:

		MarkBrdgTest(){}

		//A list of paths
		list<Path> l;
};

#endif