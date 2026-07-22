#ifndef generate_synthetic_network_h
#define generate_synthetic_network_h

#include "parameters.h"
#include "../Structures/calculate_IDs.h"
#include "../Structures/network.h"
#include "../Structures/system.h"


const std::multiset<Atom> generate_atoms(const Parameters&);
const int generate_charge(const std::multiset<Atom>&,const Parameters&);
const size_t generate_number_of_electrons(const std::multiset<Atom>&,const int&);
void add_atoms_to_ID(std::string&,const std::multiset<Atom>&);
void add_compound_number_to_ID(const Parameters&,std::string&,const std::multiset<Atom>&);
void add_conformer_number_to_ID(const Parameters&,std::string&,const std::multiset<Atom>&);
const std::string generate_system_ID(const Parameters&,const std::multiset<Atom>&,const size_t&,const int&);
const std::string generate_system_ID(const Parameters&,const std::multiset<Atom>&,const size_t&,const int&,const std::string&); // Génère l'ID avec un numéro de composé précis (pour créer un voisin intra-composé)
System generate_system(const Parameters&);

std::deque<System> generate_initial_systems(const Parameters&);

const float generate_barrier_between_compound_neighbours(const Parameters&,const System&,const System&);
const float generate_barrier_between_class_neighbours(const Parameters&,const System&,const System&);
const float generate_hyperedge_barrier(const Parameters&,const std::multiset<System>&,const std::multiset<System>&);

void generate_compound_neighbour(const Parameters&,Network&,Class&,const System&);
void generate_compound_neighbourhood(const Parameters&,Network&,const System&);

void generate_class_neighbour(const Parameters&,Network&,const System&);
void generate_class_neighbourhood(const Parameters&,Network&,Class&,const System&);

void choose_a_pair_of_systems_for_hyperedge(const Parameters&,const Network&,std::multiset<System>&);
const System& choose_system_for_splitting_hyperedge(const Parameters&,const Network&);
void generate_hyperedge_from_pair(const Parameters&,Network&,const std::multiset<System>&);
void generate_splitting_hyperedge(const Parameters&,Network&,const System&);

const Network generate_synthetic_network(const Parameters&);
const Network generate_synthetic_network(const Parameters&,std::deque<System>&);

#endif
