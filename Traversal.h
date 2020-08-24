#ifndef TRAVERSAL_HPP
#define TRAVERSAL_HPP

#define DEFAULT_DELTA 50

//This function extends the top priority path of the given priority queue on successive unitigs. If it reaches a core k-mer within an exceptable distance, the corresponding path is added to the result list. Otherwise, it is discarded (if the path exceeds the given limit) or is reinserted into the queue. The function calls itself recursively until the queue is empty. Initially, the priority queue must not be empty
void expSucPths(priority_queue<pair<uint32_t, Path>>& queue, const uint32_t& dpth, list<Path>& res);

//This function performs a BFS of the given depths on all successors of the given unitig. It returns true if a core k-mer could be reached.
const bool doSucBFS(const UnitigColorMap<CoreInfo> orig, const uint32_t dpth, list<Path>& resPths);

#endif