#ifndef generate_synthetic_network_h
#define generate_synthetic_network_h

#include "../network.h"
#include "../system.h"
#include "parameters.h"
#include "../calculate_IDs.h"


std::multiset<Atom> generate_atoms(ProbabilityFunction);
size_t generate_number_of_electrons(std::multiset<Atom>,int);
std::string generate_system_ID(std::multiset<Atom>,size_t,int);

System generate_system(Parameters);

std::deque<System> generate_initial_systems(size_t);

Network generate_synthetic_network(Parameters);


#endif
