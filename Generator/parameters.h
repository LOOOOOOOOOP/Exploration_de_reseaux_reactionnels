#ifndef parameters_h
#define parameters_h

#include "../Structures/system.h"

#include <bits/stdc++.h>
#include <cstddef>
#include <set>
#include <string>


class Parameters
{
    public:
    unsigned int seed;

    std::set<Atom> possible_atoms;                      // Liste des atomes qui peuvent se retrouver dans un système
    size_t min_number_of_atoms_in_initial_systems;
    size_t max_number_of_atoms_in_initial_systems;

    size_t number_of_initial_systems;
    bool initial_systems_are_from_different_compounds;
    bool limit_number_of_conformers_to_5;               // Les composés ont un maximum de 5 éléments

    size_t number_of_generation_rounds;
    float percentage_of_pairs_per_round;                // Pourcentage de paires de systèmes sélectionnées pour générer une hyperarête de combinaison à chaque ronde
    float minimum_number_of_pairs_per_round;            // Nombre minimal de paires de systèmes sélectionnées pour générer une hyperarête de combinaison à chaque ronde
    float percentage_of_splittings_per_round;           // Pourcentage des systèmes sélectionnés pour générer une hyperarête de fragmentation à chaque ronde
    float minimum_number_of_splittings_per_round;       // Nombre minimal de systèmes sélectionnés pour générer une hyperarête de fragmentation à chaque ronde

    Parameters();

    // Paramètres sous forme de fonction

    const std::multiset<Atom> generate_atoms() const;                                   // Génère un ensemble d'atomes pour la génération d'un système
    /* Pas utilisés
    const std::string& compound_ID_connectivity_sublayer(const std::multiset<Atom>&);   // Précise la connectivité entre les atomes, comme dans le InChI. Contribue à identifier le composé
    const std::string& compound_ID_hydrogen_sublayer(const std::multiset<Atom>&);       // Précise la connectivité des hydrogènes, comme dans le InChI. Contribue à identifier le composé
    */
    const int generate_charge(const std::multiset<Atom>&) const;                        // Génère une charge en fonction des atomes pour la génération d'un système.
                                                                                        // Le nombre d'électrons peut ensuite être déduit

    const size_t size_of_compound(const std::multiset<Atom>&) const;                    // Sert à attribuer un numéro de conformère lors de la génération d'un système
    const size_t number_of_compound_neighbours_generator(const System&) const;       // Retourne le nombre de voisins intra-composé d'un système lors de l'ajout de ses voisins
    const size_t size_of_class(const std::multiset<Atom>&) const;                       // Sert à attribuer un numéro de composé lors de la génération d'un système
    const size_t number_of_class_neighbours_generator(const System&) const;          // Retourne le nombre de voisins extra-composé d'un système lors de l'ajout de ses voisins

    const float generate_barrier_between_compound_neighbours(const System&,const System&) const;             // Génère une valeur de barrière entre deux voisins intra-composé
    const float generate_barrier_between_class_neighbours(const System&,const System&) const;                // Génère une valeur de barrière entre deux voisins extra-composé
    const float generate_hyperedge_barrier(const std::multiset<System>&,const std::multiset<System>&) const; // Génère une valeur de barrière entre deux voisins extra-classe
};

#endif
