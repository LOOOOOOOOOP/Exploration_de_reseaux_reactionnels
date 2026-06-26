#include "generate_synthetic_network.h"
#include "print_tests.h"

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

    multiset<Atom> atoms_copy = atoms;
    pair<multiset<Atom>::iterator,multiset<Atom>::iterator> range_C = atoms_copy.equal_range(C);
    if (range_C.first->symbol == C.symbol)    // S'il y a un carbone dans les atomes
    {
        // Ajout des C dans ID
        for (multiset<Atom>::iterator it = range_C.first; it != range_C.second; it++)
        {
            ID.append(it->symbol);
        }
        atoms_copy.erase(C);

        pair<multiset<Atom>::iterator,multiset<Atom>::iterator> range_H = atoms_copy.equal_range(H);
        if (range_H.first->symbol == H.symbol)  // S'il y a un hydrogène dans les atomes
        {
            // Ajout des H dans ID
            for (multiset<Atom>::iterator it = range_H.first; it != range_H.second; it++)
            {
                ID.append(it->symbol);
            }
            atoms_copy.erase(H);
        }
    }

    for (multiset<Atom>::iterator it = atoms_copy.begin(); it != atoms_copy.end(); it++)
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
    string S_class_ID = calculate_class_ID(S);
    string S_InChI = calculate_InChI(S);
    string S_system_ID = S.system_ID;
    Class& C = N.classes.find(S_class_ID)->second;
    size_t i = S.insertion_rank_in_class;
    Compound& D = C.class_compounds.find(S_InChI)->second;
    System& S_in_network = D.compound_systems.find(S_system_ID)->second;

    N.compound_unexplored_systems.pop_back();

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

    size_t new_number_of_neighbours = param.number_of_compound_neighbours_distribution(S_in_network);
    if (new_number_of_neighbours > number_of_neighbours)
    {
        size_t number_of_new_neighbours = new_number_of_neighbours - number_of_neighbours;

        // Certains des nouveaux voisins sont des systèmes qui existaient déjà dans le composé
        size_t number_of_new_already_existing_neighbours = min(number_of_new_neighbours,(D.compound_systems.size()-1) * 3 / 4);
        static mt19937 generator;
        uniform_int_distribution<int> distribute(0,D.compound_systems.size()-1);
        for (size_t k = 0; k < number_of_new_already_existing_neighbours; k++)
        {
            int neighbour_index = distribute(generator);
            System& P = next(D.compound_systems.begin(),neighbour_index)->second;
            N.add_system_in_network_from_edge(C,S_in_network,P,generate_barrier_between_compound_neighbours(param,S_in_network,P),generate_barrier_between_compound_neighbours(param,S_in_network,P));
        }
        // On génère le reste des voisins requis
        size_t number_of_new_non_existent_neighbours = number_of_new_neighbours - number_of_new_already_existing_neighbours;

        for (size_t k = 0; k < number_of_new_non_existent_neighbours; k++)
            generate_compound_neighbour(param,N,C,S_in_network);
    }

    S_in_network.compound_explored = true;
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

void generate_class_neighbourhood(Parameters param, Network& N, System& S)
{
    string S_class_ID = calculate_class_ID(S);
    string S_InChI = calculate_InChI(S);
    string S_system_ID = S.system_ID;
    size_t i = S.insertion_rank_in_class;
    Class& C = N.classes.find(S_class_ID)->second;
    Compound& D = C.class_compounds.find(S_InChI)->second;
    System& S_in_network = D.compound_systems.find(S_system_ID)->second;

    N.class_unexplored_systems.pop_back();

    // Nombre de voisins extra-composé que S possède déjà
    size_t number_of_neighbours = 0;
    size_t number_of_extra_compound_systems = 0;
    for (map<string,Compound>::iterator it = C.class_compounds.begin(); it != C.class_compounds.end(); it++)
    {
        if (it->second == D)
            continue;
        Compound& D2 = it->second;
        for (map<string,System>::iterator it2 = D2.compound_systems.begin(); it2 != D2.compound_systems.end(); it2++)
        {
            size_t j = it2->second.insertion_rank_in_class;
            if (C.edges[i][j] != -1)
            {
                number_of_neighbours++;
            }
            number_of_extra_compound_systems++;
        }
    }

    size_t new_number_of_neighbours = param.number_of_class_neighbours_distribution(S_in_network);
    if (new_number_of_neighbours > number_of_neighbours)
    {
        size_t number_of_new_neighbours = new_number_of_neighbours - number_of_neighbours;

        // Certains des nouveaux voisins sont des systèmes qui existaient déjà dans la classe
        size_t number_of_new_already_existing_neighbours = min(number_of_new_neighbours,number_of_extra_compound_systems * 3 / 4);
        static mt19937 generator;
        uniform_int_distribution<int> distribute_compound_index(0,C.class_compounds.size() - 1);
        for (size_t k = 0; k < number_of_new_already_existing_neighbours; k++)
        {
            int compound_index;
            do
                compound_index = distribute_compound_index(generator);
            while (next(C.class_compounds.begin(),compound_index)->second == D);

            Compound& D2 = next(C.class_compounds.begin(),compound_index)->second;

            uniform_int_distribution<int> distribute_neighbour_index(0,D2.compound_systems.size() - 1);
            int neighbour_index = distribute_neighbour_index(generator);

            System& P = next(D2.compound_systems.begin(),neighbour_index)->second;
            N.add_system_in_network_from_edge(C,S_in_network,P,generate_barrier_between_class_neighbours(param,S_in_network,P),generate_barrier_between_class_neighbours(param,S_in_network,P));
        }
        // On génère le reste des voisins requis
        size_t number_of_new_non_existent_neighbours = number_of_new_neighbours - number_of_new_already_existing_neighbours;

        for (size_t k = 0; k < number_of_new_non_existent_neighbours; k++)
            generate_class_neighbour(param,N,C,S_in_network);
    }

    S_in_network.class_explored = true;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void choose_a_pair_of_systems_for_hyperedge(Network& N,multiset<System>& R)
{
    static mt19937 generator;
    bool hyperedge_already_exists = false;

    // Premier réactif
    uniform_int_distribution<int> random_class_index(0,N.classes.size() - 1);
    int class_index_1 = random_class_index(generator);
    Class& C_1 = next(N.classes.begin(),class_index_1)->second;

    uniform_int_distribution<int> random_compound_index_1(0,C_1.class_compounds.size() - 1);
    int compound_index_1 = random_compound_index_1(generator);
    Compound& D_1 = next(C_1.class_compounds.begin(),compound_index_1)->second;

    uniform_int_distribution<int> random_system_index_1(0,D_1.compound_systems.size() - 1);
    int system_index_1 = random_system_index_1(generator);
    System& S_1 = next(D_1.compound_systems.begin(),system_index_1)->second;

    R.insert(S_1);

    do
    {
        // Deuxième réactif
        int class_index_2 = random_class_index(generator);
        Class& C_2 = next(N.classes.begin(),class_index_2)->second;

        uniform_int_distribution<int> random_compound_index_2(0,C_2.class_compounds.size() - 1);
        int compound_index_2 = random_compound_index_2(generator);
        Compound& D_2 = next(C_2.class_compounds.begin(),compound_index_2)->second;

        uniform_int_distribution<int> random_system_index_2(0,D_2.compound_systems.size() - 1);
        int system_index_2 = random_system_index_2(generator);
        System& S_2 = next(D_2.compound_systems.begin(),system_index_2)->second;

        R.insert(S_2);

        // Vérification qu'une hyperarête avec ces réactifs (ou produits) n'existe pas déjà
        for (set<Hyperedge>::iterator it = N.hyperedges.begin(); it != N.hyperedges.end(); it++)
        {
            if (it->reactants == R || it->products == R)
            {
                hyperedge_already_exists = true;
                R.erase(S_2);
                break;
            }
        }
    }
    while (hyperedge_already_exists == true);
}

/////////////////////////////////////////////////////////////////////

System& choose_system_for_splitting_hyperedge(Network& N)
{
    static mt19937 generator;
    bool hyperedge_already_exists = false;
    System* pR;

    do
    {
        uniform_int_distribution<int> random_class_index(0,N.classes.size() - 1);
        int class_index = random_class_index(generator);
        Class& C = next(N.classes.begin(),class_index)->second;

        uniform_int_distribution<int> random_compound_index(0,C.class_compounds.size() - 1);
        int compound_index = random_compound_index(generator);
        Compound& D = next(C.class_compounds.begin(),compound_index)->second;

        uniform_int_distribution<int> random_system_index(0,D.compound_systems.size() - 1);
        int system_index = random_system_index(generator);
        pR = &next(D.compound_systems.begin(),system_index)->second;

        // Vérification qu'une hyperarête avec ce réactif (ou produit) n'existe pas déjà
        for (set<Hyperedge>::iterator it = N.hyperedges.begin(); it != N.hyperedges.end(); it++)
        {
            multiset<System> Reactants = {*pR};
            if (it->reactants == Reactants || it->products == Reactants)
            {
                hyperedge_already_exists = true;
                break;
            }
        }
    }
    while (hyperedge_already_exists == true);

    System& R = *pR;
    return R;
}

/////////////////////////////////////////////////////////////////////

void generate_hyperedge_from_pair(Parameters param, Network& N, multiset<System>& R)
{
    // Génération du produit
    multiset<Atom> P_atoms = R.begin()->atoms;
    for (multiset<Atom>::iterator it = next(R.begin())->atoms.begin(); it != next(R.begin())->atoms.end(); it++)
        P_atoms.insert(*it);
    size_t P_n_electrons = R.begin()->n_electrons + next(R.begin())->n_electrons;
    int P_charge = R.begin()->charge + next(R.begin())->charge;
    string P_ID = generate_system_ID(param,P_atoms,P_n_electrons,P_charge);
    System P(P_ID,P_atoms,P_n_electrons,P_charge);
    multiset<System> Product = {P};

    // Génération de l'hyperarête
    float RP_barrier = generate_hyperedge_barrier(param,R,Product);
    float PR_barrier = generate_hyperedge_barrier(param,R,Product);
    N.add_system_in_network_from_hyperedge(R,Product,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////

void generate_splitting_hyperedge(Parameters param, Network& N, System& S)
{
    size_t total_number_of_atoms = S.atoms.size();
    multiset<Atom> available_atoms = S.atoms;
    static mt19937 generator;

    // Séparation des atomes
    uniform_int_distribution<int> distribute_number_of_atoms(1,total_number_of_atoms - 1);
    size_t number_of_atoms_P1 = distribute_number_of_atoms(generator);
    multiset<Atom> atoms_P1;
    for (size_t i = 0; i < number_of_atoms_P1; i++)
    {
        uniform_int_distribution<int> distribute_atoms_index(0,available_atoms.size() - 1);
        size_t index = distribute_atoms_index(generator);
        atoms_P1.insert(*next(available_atoms.begin(),index));
        available_atoms.erase(next(available_atoms.begin(),index));
    }
    multiset<Atom> atoms_P2 = available_atoms;

    // Séparation de la charge
    uniform_int_distribution<int> distribute_charge(min(S.charge,0),max(0,S.charge));
    int charge_P1 = distribute_charge(generator);
    int charge_P2 = S.charge - charge_P1;

    // Création des produits
    size_t n_electrons_P1 = generate_number_of_electrons(atoms_P1,charge_P1);
    string ID_P1 = generate_system_ID(param,atoms_P1,n_electrons_P1,charge_P1);
    System P1(ID_P1,atoms_P1,n_electrons_P1,charge_P1);

    size_t n_electrons_P2 = generate_number_of_electrons(atoms_P2,charge_P2);
    string ID_P2 = generate_system_ID(param,atoms_P2,n_electrons_P2,charge_P2);
    System P2(ID_P2,atoms_P2,n_electrons_P2,charge_P2);

    // Ajout de l'hyperarête
    multiset<System> R = {S};
    multiset<System> P = {P1,P2};
    float RP_barrier = generate_hyperedge_barrier(param,R,P);
    float PR_barrier = generate_hyperedge_barrier(param,R,P);
    N.add_system_in_network_from_hyperedge(R,P,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

Network generate_synthetic_network(Parameters param)
{
    deque<System> initial_systems = generate_initial_systems(param);
    Network N(initial_systems);

    for (size_t i = 0; i < param.number_of_generation_rounds; i++)
    {
        // Pour les systèmes non classe-explorés
            // Génération des voisins extra-composé
        while (N.class_unexplored_systems.size() > 0)
            generate_class_neighbourhood(param,N,N.class_unexplored_systems.back());

        // Pour un certain nombre de paires de systèmes, génération des voisins extra-classe
        size_t number_of_pairs = N.number_of_systems * param.percentage_of_pairs_per_round;
        for (size_t i = 0; i < number_of_pairs; i++)
        {
            multiset<System> R;
            choose_a_pair_of_systems_for_hyperedge(N,R);
            generate_hyperedge_from_pair(param,N,R);
        }
/*
        // Pour un certain nombre de systèmes, génération des voisins extra-classe en se splittant
        size_t number_of_splitting_systems = N.number_of_systems * param.percentage_of_splittings_per_round;
        for (size_t i = 0; i < number_of_splitting_systems; i++)
        {
            System R = choose_system_for_splitting_hyperedge(N);
            generate_splitting_hyperedge(param,N,R);
        }
*/  // infinite loop
        // Pour les systèmes non composé-explorés
            // Génération des voisins intra-composé
        while (N.compound_unexplored_systems.size() > 0)
            generate_compound_neighbourhood(param,N,N.compound_unexplored_systems.back());
    }

    return N;
}

