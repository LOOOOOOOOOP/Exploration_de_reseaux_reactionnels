#ifndef explore_synthetic_network_h
#define explore_synthetic_network_h

#include "../network.h"

bool stop_exploring(Network);

void find_compound_neighbour(Network&,System&);
void find_all_compound_neighbours(Network&,System&);

void find_class_neighbour(Network&,System&);
void find_all_class_neighbours(Network&,System&);

void find_hyperedge_from_pair(Network&,std::pair<System,System>&);
void find_splitting_hyperedge_from_system(Network&,System&);

#endif
