#ifndef parameters_h
#define parameters_h

#include "../system.h"

#include <cstddef>
#include <set>
#include <string>
#include <bits/stdc++.h>


class Parameters
{
    public:
    unsigned int seed;

    std::set<Atom> possible_atoms;
    size_t min_number_of_atoms_in_system;
    size_t max_number_of_atoms_in_system;

    size_t number_of_initial_systems;
    bool initial_systems_are_from_different_compounds;    // les systèmes initiaux appartiennent tous à des composés différents
    bool limit_number_of_conformers_to_5;

    size_t number_of_generation_rounds;
    float percentage_of_pairs_per_round;   // pourcentage de paires de systèmes sélectionnées pour générer un voisin extra-classe à chaque ronde
    float minimum_number_of_pairs_per_round; // nombre minimal de paires de systèmes sélectionnées pour générer un voisin extra-classe à chaque ronde
    float percentage_of_splittings_per_round;   // pourcentage des systèmes sélectionnés pour se séparer en plusieurs voisins extra-classe à chaque ronde
    float minimum_number_of_splittings_per_round; // nombre minimal de systèmes sélectionnés pour se séparer en plusieurs voisins extra-classe à chaque ronde


    Parameters();

    std::multiset<Atom> atoms_distribution_function();
    std::string InChI_connectivity_sublayer(std::multiset<Atom>);
    std::string InChI_hydrogen_sublayer(std::multiset<Atom>);
    int charge_distribution_function(std::multiset<Atom>);

    size_t size_of_compound(std::multiset<Atom>);
    size_t number_of_compound_neighbours_distribution(System);
    size_t size_of_class(std::multiset<Atom>);
    size_t number_of_class_neighbours_distribution(System);

    float generate_barrier_between_compound_neighbours(System,System);    // génère une valeur de barrière entre deux systèmes du même composé
    float generate_barrier_between_class_neighbours(System,System); // génère une valeur de barrière entre deux systèmes de la même classe (mais pas du même composé)

    float generate_hyperedge_barrier(std::multiset<System>,std::multiset<System>);
};


#endif
