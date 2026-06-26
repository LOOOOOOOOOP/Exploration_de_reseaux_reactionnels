#include "parameters.h"

#include <random>
#include <iostream>

using namespace std;

Parameters::Parameters() :
    seed(14),

    possible_atoms({Atom("H",1,0.6),
                    Atom("C",6,0.25),
                    Atom("O",8,0.05),
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
    limit_number_of_conformers_to_5(true),

    number_of_generation_rounds(10),
    percentage_of_pairs_per_round(0.10),
    percentage_of_splittings_per_round(0.10),
    max_splittedness(3)
{}

/////////////////////////////////////////////////////////////////////

multiset<Atom> Parameters::atoms_distribution_function()
{
    multiset<Atom> atoms;
    static mt19937 generator(seed);

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

string Parameters::InChI_connectivity_sublayer(multiset<Atom> atoms)    // Pas utilisé
{
    static mt19937 generator(seed);
    string sublayer = "/c";

    Atom H("H",1,0.6);
    size_t number_of_non_H = atoms.size() - atoms.count(H);

    if (number_of_non_H < 2)
        return "";

    // Numéros des atomes à connecter (parfois avec répétition)
    multiset<int> numbers;
    for (size_t i = 1; i <= number_of_non_H; i++)
    {
        numbers.insert(i);
        uniform_int_distribution<int> distribute(1,100);
        int random_number = distribute(generator);
        if (random_number > 90 && random_number <= 97)
            numbers.insert(i);
        else if (random_number > 97)
        {
            numbers.insert(i);
            numbers.insert(i);
        }
    }

    // Ajout de la connectivité au sublayer
    while (numbers.empty() == false)
    {
        // On choisit un numéro au hasard
        uniform_int_distribution<int> distribute(0,numbers.size()-1);
        int number_index = distribute(generator);
        set<int>::iterator iterator_to_number = next(numbers.begin(),number_index);
        int number = *iterator_to_number;

        if (sublayer[sublayer.length()-1] == 'c' || sublayer[sublayer.length()-1] == ')')
            sublayer += to_string(number);
        else
        {
            // 25% de chances d'être entre parenthèses
            uniform_int_distribution<int> distribute(1,100);
            int random_number = distribute(generator);
            if (random_number < 75)
                sublayer = sublayer + "-" + to_string(number);
            else
                sublayer = sublayer + "(" + to_string(number) + ")";
        }
        numbers.erase(iterator_to_number);
    }

    return sublayer;
}

/////////////////////////////////////////////////////////////////////

string Parameters::InChI_hydrogen_sublayer(multiset<Atom> atoms)    // pas utilisé
{
    string sublayer = "/h";

    return sublayer;
}

/////////////////////////////////////////////////////////////////////

int Parameters::charge_distribution_function(multiset<Atom> atoms)
{
    static mt19937 generator(seed);
    discrete_distribution<> charge_distribution({0.0125,0.0375,0.075,0.75,0.075,0.0375,0.0125});

    int number_of_electrons = 0;
    for (multiset<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        number_of_electrons += it->atomic_number;
    }

    int charge;
    do
    {
        charge = charge_distribution(generator) - 3;
    }
    while (number_of_electrons - charge < 0);

    return charge;
}

/////////////////////////////////////////////////////////////////////

size_t Parameters::number_of_compound_neighbours_distribution(System S)
{
    int number_of_neighbours;
    int number_of_non_H = S.atoms.size() - S.atoms.count(Atom("H",1,0.6));

    if (S.atoms.size() < 3 || number_of_non_H < 2)
        return 0;

    if (limit_number_of_conformers_to_5 == true)
    {
        static mt19937 generator;
        uniform_int_distribution<int> distribute(0,5);
        return(distribute(generator));
    }

    int number_of_rotatable_bonds = max(number_of_non_H - 1,0);
    int average_number_of_neighbours = pow(3,number_of_rotatable_bonds);
    static mt19937 generator(seed);
    uniform_int_distribution<int> distribute(average_number_of_neighbours - number_of_rotatable_bonds, average_number_of_neighbours + number_of_rotatable_bonds);
    number_of_neighbours = distribute(generator);
    return number_of_neighbours;
}

/////////////////////////////////////////////////////////////////////

size_t Parameters::size_of_class(multiset<Atom> atoms)
{
    size_t number_of_non_H = atoms.size() - atoms.count(Atom("H",1,0.6));
    size_t number_of_stereocentres = number_of_non_H / 5;
    size_t size_of_class = pow(2,number_of_stereocentres);
    return size_of_class;
}

/////////////////////////////////////////////////////////////////////

size_t Parameters::number_of_class_neighbours_distribution(System S)
{
    multiset<Atom> atoms = S.atoms;
    size_t size_class = size_of_class(atoms);
    static mt19937 generator(seed);

    if (size_class == 1)
        return 0;

    size_t average_number_of_neighbours = size_class * 3 / 4;
    if (average_number_of_neighbours < 5)
    {
        uniform_int_distribution<int> distribute(0, min(average_number_of_neighbours + 5,size_class));
        return(distribute(generator));
    }
    else
    {
        uniform_int_distribution<int> distribute(average_number_of_neighbours - 5, min(average_number_of_neighbours + 5,size_class));
        return(distribute(generator));
    }
}

/////////////////////////////////////////////////////////////////////

float Parameters::generate_barrier_between_compound_neighbours(System R,System P)
{
    static mt19937 generator(seed);
    uniform_real_distribution<float> distribute(0,5);
    return distribute(generator);
}

/////////////////////////////////////////////////////////////////////

float Parameters::generate_barrier_between_class_neighbours(System R, System P)
{
    static mt19937 generator(seed);
    uniform_real_distribution<float> distribute(5,50);
    return distribute(generator);
}

/////////////////////////////////////////////////////////////////////

float Parameters::generate_hyperedge_barrier(multiset<System> R, multiset<System> P)
{
    static mt19937 generator(seed);
    uniform_real_distribution<float> distribute(50,250);
    return distribute(generator);
}

