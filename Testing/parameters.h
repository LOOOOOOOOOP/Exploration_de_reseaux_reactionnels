#ifndef parameters_h
#define parameters_h

#include "../system.h"

#include <cstddef>
#include <set>
#include <string>


class Atom
{
    public:
    std::string symbol;
    int atomic_number;
    float probability;

    Atom();
    Atom(std::string s,int n, float p): symbol(s), atomic_number(n), probability(p){}

    friend bool operator==(Atom a, Atom b)
    {
        return (a.symbol == b.symbol);
    }

    friend bool operator<(Atom a, Atom b)
    {
        return (a.symbol < b.symbol);
    }

    friend bool operator>(Atom a, Atom b)
    {
        return (a.symbol > b.symbol);
    }
};

/////////////////////////////////////////////////////////////////////

class Parameters
{
    public:
    unsigned int seed;

    std::set<Atom> possible_atoms;
    size_t min_number_of_atoms_in_system;
    size_t max_number_of_atoms_in_system;

    size_t number_of_initial_systems;
    bool initial_systems_are_from_different_compounds;    // les systèmes initiaux appartiennent tous à des composés différents

    size_t number_of_generation_rounds;
    float percentage_of_pairs_per_round;   // pourcentage des paires de systèmes sélectionnées pour générer un voisin extra-classe à chaque ronde
    float percentage_of_splittings_per_round;   // pourcentage des systèmes sélectionnés pour se séparer en plusieurs voisins extra-classe à chaque ronde
    size_t max_splittedness;    // nombre maximal de voisins extra-classe générés à partir d'un système


    Parameters();

    std::multiset<Atom> atoms_distribution_function();
    int charge_distribution_function();

    size_t conformer_degree_distribution(System);   // retourne le nombre de voisins intra-composés d'un système

    float compound_barrier_distribution();    // génère une valeur de barrière entre deux systèmes du même composé
};


#endif
