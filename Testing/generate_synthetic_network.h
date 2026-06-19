#ifndef generate_synthetic_network_h
#define generate_synthetic_network_h

#include "../network.h"
#include "../system.h"
#include "parameters.h"
#include "../calculate_IDs.h"


std::multiset<Atom> generate_atoms(Parameters);
int generate_charge(std::multiset<Atom>,Parameters);
size_t generate_number_of_electrons(std::multiset<Atom>,int);
std::string generate_system_ID(Parameters,std::multiset<Atom>,size_t,int);

System generate_system(Parameters);

std::deque<System> generate_initial_systems(size_t);

float generate_barrier_between_compound_neighbours(Parameters,System&,System&);
float generate_barrier_between_class_neighbours(Parameters,System&,System&);
float generate_hyperedge_barrier(Parameters,std::multiset<System>&,std::multiset<System>&);

void generate_compound_neighbour(Parameters,Network&,Class&,System&);
void generate_compound_neighbourhood(Parameters,Network&,System&);

void generate_class_neighbour(Parameters,Network&,System&);
void generate_class_neighbourhood(Parameters,Network&,System&);

std::vector<std::pair<System,System>>& choose_pairs_of_systems_for_hyperedges(Parameters,Network&);
std::vector<System>& choose_systems_for_splitting_hyperedges(Parameters,Network&);
void generate_hyperedge_from_pair(Parameters,Network&,std::pair<System,System>&);
void generate_splitting_hyperedge(Parameters,Network&,System&);

Network generate_synthetic_network(Parameters);


#endif
