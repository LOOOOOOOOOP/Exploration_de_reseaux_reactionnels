#ifndef generate_synthetic_network_h
#define generate_synthetic_network_h

#include "../network.h"
#include "../system.h"
#include "parameters.h"
#include "../calculate_IDs.h"


std::multiset<Atom> generate_atoms(Parameters);
int generate_charge(std::multiset<Atom>,Parameters);
size_t generate_number_of_electrons(std::multiset<Atom>,int);
void add_atoms_to_ID(std::string&,std::multiset<Atom>);
void add_compound_number_to_ID(Parameters,std::string&,std::multiset<Atom>);
void add_conformer_number_to_ID(Parameters,std::string&,std::multiset<Atom>);
std::string generate_system_ID(Parameters,std::multiset<Atom>,size_t,int);
std::string generate_system_ID(Parameters,std::multiset<Atom>,size_t,int,std::string);

System generate_system(Parameters);

std::deque<System> generate_initial_systems(Parameters);

float generate_barrier_between_compound_neighbours(Parameters,System&,System&);
float generate_barrier_between_class_neighbours(Parameters,System&,System&);
float generate_hyperedge_barrier(Parameters,std::multiset<System>&,std::multiset<System>&);

void generate_compound_neighbour(Parameters,Network&,Class&,System&);
void generate_compound_neighbourhood(Parameters,Network&,System&);

void generate_class_neighbour(Parameters,Network&,System&);
void generate_class_neighbourhood(Parameters,Network&,Class&,System&);

void choose_a_pair_of_systems_for_hyperedge(Parameters,Network&,std::multiset<System>&);
System& choose_system_for_splitting_hyperedge(Parameters,Network&,System&);
void generate_hyperedge_from_pair(Parameters,Network&,std::multiset<System>&);
void generate_splitting_hyperedge(Parameters,Network&,System&);

Network generate_synthetic_network(Parameters);
Network generate_synthetic_network(Parameters,std::deque<System>&);


#endif
