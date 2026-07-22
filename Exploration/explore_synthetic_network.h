#ifndef explore_synthetic_network_h
#define explore_synthetic_network_h

#include "../Structures/network.h"


void explore(Network&,const Network&);

void find_compound_neighbour(Network&,const Network&,System&);
void find_all_compound_neighbours(Network&,const Network&,System&);

void find_class_neighbour(Network&,const Network&,System&);
void find_all_class_neighbours(Network&,const Network&,System&);

void find_hyperedge_from_pair(Network&,const Network&,std::pair<System,System>&);
void find_splitting_hyperedge_from_system(Network&,const Network&,System&);

#endif
