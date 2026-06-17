#include "parameters.h"

#include <random>
#include <iostream>

using namespace std;

Parameters::Parameters() :
    seed(14),

    possible_atoms({Atom("H",1,0.5),
                    Atom("C",6,0.3),
                    Atom("O",8,0.1),
                    Atom("N",7,0.05),
                    Atom("F",9,0.002),
                    Atom("Cl",17,0.02),
                    Atom("Br",35,0.005),
                    Atom("I",53,0.001),
                    Atom("P",15,0.001),
                    Atom("S",16,0.02),
                    Atom("B",5,0.001)}),
    min_number_of_atoms_in_system(1),
    max_number_of_atoms_in_system(15),

    number_of_initial_systems(4),
    initial_systems_are_from_different_compounds(true),

    number_of_generation_rounds(10),
    percentage_of_pairs_per_round(0.10),
    percentage_of_splittings_per_round(0.10),
    max_splittedness(3)
{}

/////////////////////////////////////////////////////////////////////

multiset<Atom> Parameters::atoms_distribution_function()
{
    multiset<Atom> atoms;
    mt19937 generator(seed);

    vector<float> probabilities;
    for (set<Atom>::iterator it = possible_atoms.begin(); it != possible_atoms.end(); it++)
    {
        probabilities.push_back(it->probability);
    }

    discrete_distribution<int> atoms_distribution(probabilities.begin(),probabilities.end());

    uniform_int_distribution<int> distribute(min_number_of_atoms_in_system,max_number_of_atoms_in_system);
    int number_of_atoms = distribute(generator);

    for (int i = 0; i < number_of_atoms; i++)
    {
        int index = atoms_distribution(generator);
        atoms.insert(*next(possible_atoms.begin(), index));
    }
    return atoms;
}

/////////////////////////////////////////////////////////////////////

int Parameters::charge_distribution_function()
{
    mt19937 generator(seed);
    discrete_distribution<> charge_distribution({0.0125,0.0375,0.075,0.75,0.075,0.0375,0.0125});
    int charge = charge_distribution(generator) - 3;
    return charge;
}

/////////////////////////////////////////////////////////////////////

size_t Parameters::conformer_degree_distribution(System S)
{
    size_t degree;
    size_t number_of_non_H = S.atoms.size() - S.atoms.count("H");
    size_t number_of_rotatable_bonds = number_of_non_H - 1;
    size_t average_degree = pow(3,number_of_rotatable_bonds);
    mt19937 generator(seed);
    uniform_int_distribution<int> distribute(average_degree - number_of_rotatable_bonds, average_degree + number_of_rotatable_bonds);
    degree = distribute(generator);
    return degree;
}

/////////////////////////////////////////////////////////////////////

float Parameters::compound_barrier_distribution()
{
    mt19937 generator(seed);
    uniform_real_distribution<float> distribute(-5,5);
    return distribute(generator);
}
