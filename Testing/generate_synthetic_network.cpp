#include "generate_synthetic_network.h"

#include <cstdlib>
#include <deque>
#include <set>
#include <string>
#include <random>

using namespace std;


multiset<Atom> generate_atoms(Parameters param)
{
    multiset<Atom> atoms = param.atoms_distribution_function();
    return atoms;
}

/////////////////////////////////////////////////////////////////////

int generate_charge(multiset<Atom> atoms, Parameters param)
{
    int charge = param.charge_distribution_function(atoms);
    return charge;
}

/////////////////////////////////////////////////////////////////////

size_t generate_number_of_electrons(multiset<Atom> atoms,int charge)
{
    size_t number_of_electrons = 0;

    for (multiset<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        number_of_electrons += it->atomic_number;
    }
    number_of_electrons += -charge;

    return number_of_electrons;
}

/////////////////////////////////////////////////////////////////////

string generate_system_ID(Parameters param,multiset<Atom> atoms,size_t n_electrons,int charge)
{
    string ID = "InChI=1S/";

    // ajout des atomes ordonnés selon la notation de Hill
    Atom C("C",6,0.25);
    Atom H("H",1,0.6);

    pair<multiset<Atom>::iterator,multiset<Atom>::iterator> range_C = atoms.equal_range(C);
    if (range_C.first->symbol == C.symbol)    // S'il y a un carbone dans les atomes
    {
        // Ajout des C dans ID
        for (multiset<Atom>::iterator it = range_C.first; it != range_C.second; it++)
        {
            ID.append(it->symbol);
        }
        atoms.erase(C);

        pair<multiset<Atom>::iterator,multiset<Atom>::iterator> range_H = atoms.equal_range(H);
        if (range_H.first->symbol == H.symbol)  // S'il y a un hydrogène dans les atomes
        {
            // Ajout des H dans ID
            for (multiset<Atom>::iterator it = range_H.first; it != range_H.second; it++)
            {
                ID.append(it->symbol);
            }
            atoms.erase(H);
        }
    }

    for (multiset<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        ID.append(it->symbol);
    }
/*
    // ajout du sublayer de connectivité
    ID = ID + param.InChI_connectivity_sublayer(atoms);

    // ajout du sublayer des hydrogènes
    ID = ID + param.InChI_hydrogen_sublayer(atoms);
*/
    // spécification alternative du composé contenant le système
    static mt19937 generator(param.seed);
    size_t size_of_class = param.size_of_class(atoms);
    uniform_int_distribution<int> distribute(1,size_of_class);
    int compound_number = distribute(generator);
    ID = ID + "/compound#" + to_string(compound_number);

    // ajout du sublayer de charge
    if (charge > 0)
        ID = ID + "/q+" + to_string(charge);
    else if (charge < 0)
        ID = ID + "/q" + to_string(charge);

    return ID;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

System generate_system(Parameters param)
{
    multiset<Atom> atoms = generate_atoms(param);
    int charge = generate_charge(atoms,param);
    size_t number_of_electrons = generate_number_of_electrons(atoms,charge);
    string system_ID = generate_system_ID(param,atoms,number_of_electrons,charge);

    System S(system_ID,atoms,number_of_electrons,charge);
    return S;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

deque<System> generate_initial_systems(Parameters param)
{
    deque<System> initial_systems;

    for (size_t i = 0; i < param.number_of_initial_systems; i++)
    {
        if (param.initial_systems_are_from_different_compounds == true)
        {
            bool compound_already_exists;
            do
            {
                compound_already_exists = false;
                System S = generate_system(param);

                if (initial_systems.size() > 0)
                {
                    if (param.initial_systems_are_from_different_compounds == true)
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
                }
            if (compound_already_exists == false)
                initial_systems.push_back(S);
            }
            while (compound_already_exists == true);
        }
        else    // les systèmes initiaux peuvent être du même composé
        {
            System S = generate_system(param);
            initial_systems.push_back(S);
        }
    }

    return initial_systems;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

float generate_barrier_between_compound_neighbours(Parameters param, System& R, System& P)
{
    float barrier = param.generate_barrier_between_compound_neighbours(R,P);
    return barrier;
}

/////////////////////////////////////////////////////////////////////

float generate_barrier_between_class_neighbours(Parameters param, System& R, System& P)
{
    float barrier = param.generate_barrier_between_class_neighbours(R,P);
    return barrier;
}

/////////////////////////////////////////////////////////////////////

float generate_hyperedge_barrier(Parameters param, multiset<System>& R, multiset<System>& P)
{
    float barrier = param.generate_hyperedge_barrier(R,P);
    return barrier;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void generate_compound_neighbour(Parameters param, Network& N, Class& C, System& R)
{
    System P(calculate_InChI(R),R.atoms,R.n_electrons,R.charge);
    float RP_barrier = generate_barrier_between_compound_neighbours(param,R,P);
    float PR_barrier = generate_barrier_between_compound_neighbours(param,R,P);
    N.add_system_in_network_from_edge(C,R,P,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////

void generate_compound_neighbourhood(Parameters param,Network& N,System& S)
{
    N.compound_unexplored_systems.pop_back();

    string S_class_ID = calculate_class_ID(S);
    string S_InChI = calculate_InChI(S);
    Class& C = N.classes.find(S_class_ID)->second;
    size_t i = S.insertion_rank_in_class;
    Compound& D = C.class_compounds.find(S_InChI)->second;

    // Nombre de voisins intra-composé que S possède déjà
    size_t number_of_neighbours = 0;
    for (map<string,System>::iterator it = D.compound_systems.begin(); it != D.compound_systems.end(); it++)
    {
        size_t j = it->second.insertion_rank_in_class;
        if (C.edges[i][j] != -1)
        {
            number_of_neighbours++;
        }
    }

    size_t new_number_of_neighbours = param.number_of_compound_neighbours_distribution(S);
    if (new_number_of_neighbours > number_of_neighbours)
    {
        for (size_t k = 0; k < new_number_of_neighbours - number_of_neighbours; k++)
            generate_compound_neighbour(param,N,C,S);
    }

    S.compound_explored = true;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void generate_class_neighbour(Parameters param, Network& N, Class& C,System& R)
{
    multiset<Atom> atoms = R.atoms;
    size_t number_of_electrons = R.n_electrons;
    int charge = R.charge;
    string ID;
    do
    {
        ID = generate_system_ID(param,atoms,number_of_electrons,charge);
    }
    while (ID == calculate_InChI(R));
    System P(ID,atoms,number_of_electrons,charge);

    float RP_barrier = generate_barrier_between_class_neighbours(param,R,P);
    float PR_barrier = generate_barrier_between_class_neighbours(param,R,P);
    N.add_system_in_network_from_edge(C,R,P,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////

void generate_class_neighbourhood(Parameters param, Network& N, Class& C, System& S)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

vector<pair<System,System>>& choose_pairs_of_systems_for_hyperedges(Parameters param, Network& N)
{
    vector<pair<System,System>> R;
    return R;
}

/////////////////////////////////////////////////////////////////////

vector<System>& choose_systems_for_splitting_hyperedges(Parameters param, Network& N)
{
    vector<System> R;
    return R;
}

/////////////////////////////////////////////////////////////////////

void generate_hyperedge_from_pair(Parameters param, Network& N, pair<System,System>&)
{
    return;
}

/////////////////////////////////////////////////////////////////////

void generate_splitting_hyperedge(Parameters param, Network& N, System& S)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

Network generate_synthetic_network(Parameters param)
{
    deque<System> initial_systems = generate_initial_systems(param);
    Network N(initial_systems);
/*
    for (size_t i = 0; i < param.number_of_generation_rounds; i++)
    {
        // Pour les systèmes non classe-explorés
            // Génération des voisins extra-composé

        // Pour un certain nombre de paires de systèmes, génération des voisins extra-classe

        // Pour un certain nombre de systèmes, génération des voisins extra-classe en se splittant

        // Pour les systèmes non composé-explorés
            // Génération des voisins intra-composé
        while (N.compound_unexplored_systems.size() > 0)
        {
            generate_compound_neighbourhood(param,N,N.compound_unexplored_systems.back());
        }
    }
*/
    return N;
}

