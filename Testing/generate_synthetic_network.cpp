#include "generate_synthetic_network.h"

#include <cstdlib>
#include <deque>
#include <set>
#include <string>

using namespace std;


multiset<Atom> generate_atoms(ProbabilityFunction p)
{
    multiset<Atom> atoms;
    return atoms;
}

int generate_charge(multiset<Atom> atoms)
{
    int charge;
    return charge;
}

/////////////////////////////////////////////////////////////////////

size_t generate_number_of_electrons(multiset<Atom> atoms,int charge)
{
    size_t number_of_electrons;

    for (multiset<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        number_of_electrons += it->atomic_number;
    }
    number_of_electrons += charge;

    return number_of_electrons;
}

/////////////////////////////////////////////////////////////////////

string generate_system_ID(multiset<Atom> atoms,size_t n_electrons,int charge)
{
    string ID;  // atomes ordonnés selon la notation de Hill



    return ID;
}

/////////////////////////////////////////////////////////////////////

System generate_system(Parameters P)
{
    multiset<Atom> atoms = generate_atoms(P.atoms_distribution_function);
    multiset<string> atoms_names;
    for (multiset<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_names.insert(it->symbol);
    }

    int charge = generate_charge(atoms);
    size_t number_of_electrons = generate_number_of_electrons(atoms,charge);
    string system_ID = generate_system_ID(atoms,number_of_electrons,charge);

    System S(system_ID,atoms_names,number_of_electrons);
    return S;
}

/////////////////////////////////////////////////////////////////////

deque<System> generate_initial_systems(Parameters P)
{
    deque<System> initial_systems;

    for (size_t i = 0; i < P.number_of_initial_systems; i++)
    {
        System* pS;
        bool compound_already_exists = false;
        do
        {
            System S = generate_system(P);

            if (initial_systems.size() == 0)
                continue;
            else
            {
            for (size_t j = 0; j < initial_systems.size(); j++)
            {
                if (calculate_InChI(initial_systems[j]) == calculate_InChI(S))
                {
                    compound_already_exists = true;
                    break;
                }
            }
            }
            pS = &S;
        }
        while (compound_already_exists == true);

        initial_systems.push_back(*pS);
    }

    return initial_systems;
}

/////////////////////////////////////////////////////////////////////

Network generate_synthetic_network(Parameters P)
{
    deque<System> initial_systems = generate_initial_systems(P);
    Network N(initial_systems);

    for (size_t i = 0; i < P.number_of_generation_rounds; i++)
    {
        // Pour les systèmes non explorés
            // Génération des voisins intra-composé

            // Génération des voisins extra-composé

        // Pour un certain nombre de paires de systèmes, génération des voisins extra-classe

        // Pour un certain nombre de systèmes, génération des voisins extra-classe en se splittant

    }
    return N;
}
