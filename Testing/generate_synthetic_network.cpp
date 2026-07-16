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

void add_atoms_to_ID(string& ID,multiset<Atom> atoms)
{
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
}

/////////////////////////////////////////////////////////////////////

void add_compound_number_to_ID(Parameters param,string& ID,multiset<Atom> atoms)
{
    static mt19937 generator(param.seed);
    size_t size_of_class = param.size_of_class(atoms);
    uniform_int_distribution<int> distribute(1,size_of_class);
    size_t compound_number = distribute(generator);
    ID = ID + "/compound#" + to_string(compound_number);
}

/////////////////////////////////////////////////////////////////////

void add_conformer_number_to_ID(Parameters param,string& ID,multiset<Atom> atoms)
{
    static mt19937 generator(param.seed);
    size_t size_of_compound = param.size_of_compound(atoms);
    uniform_int_distribution<int> distribute_conformer_number(1,size_of_compound);
    size_t conformer_number = distribute_conformer_number(generator);
    ID = ID + "/[" + to_string(conformer_number) + "]";
}

/////////////////////////////////////////////////////////////////////

string generate_system_ID(Parameters param,multiset<Atom> atoms,size_t n_electrons,int charge)
{
    string ID = "";
/*  // pas utilisé //
    ID += "InChI=1S/";
*/

    // ajout de la liste des atomes
    add_atoms_to_ID(ID,atoms);

/*  // pas utilisé //
    // ajout du sublayer de connectivité
    ID = ID + param.InChI_connectivity_sublayer(atoms);

    // ajout du sublayer des hydrogènes
    ID = ID + param.InChI_hydrogen_sublayer(atoms);
*/

    // ajout du sublayer de charge
    if (charge > 0)
        ID = ID + "/q+" + to_string(charge);
    else if (charge < 0)
        ID = ID + "/q" + to_string(charge);

    // spécification du composé contenant le système
    add_compound_number_to_ID(param,ID,atoms);

    // spécification du numéro du conformère dans le composé
    add_conformer_number_to_ID(param,ID,atoms);

    return ID;
}

/////////////////////////////////////////////////////////////////////

string generate_system_ID(Parameters param,multiset<Atom> atoms,size_t n_electrons,int charge,string compound_number)
{
    string ID = "";

/*  // pas utilisé //
    ID += "InChI=1S/";
*/

    // ajout de la liste des atomes
    add_atoms_to_ID(ID,atoms);

/*  // pas utilisé //
    // ajout du sublayer de connectivité
    ID = ID + param.InChI_connectivity_sublayer(atoms);

    // ajout du sublayer des hydrogènes
    ID = ID + param.InChI_hydrogen_sublayer(atoms);
*/

    // ajout du sublayer de charge
    if (charge > 0)
        ID = ID + "/q+" + to_string(charge);
    else if (charge < 0)
        ID = ID + "/q" + to_string(charge);

    // spécification du composé contenant le système
    ID += compound_number;

    // spécification du numéro du conformère dans le composé
    add_conformer_number_to_ID(param,ID,atoms);

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

                for (size_t j = 0; j < initial_systems.size(); j++)
                {
                    if (calculate_InChI(initial_systems[j]) == calculate_InChI(S))
                    {
                        compound_already_exists = true;
                        break;
                    }
                }
                if (compound_already_exists == false)
                    initial_systems.push_back(S);
            }
            while (compound_already_exists == true);
        }
        else    // les systèmes initiaux peuvent être du même composé
        {
            bool system_already_exists;
            do
            {
                system_already_exists = false;
                System S = generate_system(param);

                for (size_t j = 0; j < initial_systems.size(); j++)
                {
                    if (S == initial_systems[j])
                    {
                        system_already_exists = true;
                        break;
                    }
                }
                if (system_already_exists == false)
                    initial_systems.push_back(S);
            }
            while (system_already_exists == true);
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
    string ID = R.system_ID;
    ID.erase(ID.find("/["),string::npos);
    string compound_number = ID.substr(ID.find("/compound#"),string::npos);
    string P_ID = generate_system_ID(param,R.atoms,R.n_electrons,R.charge,compound_number);
    System P(P_ID,R.atoms,R.n_electrons,R.charge);

    if (P == R)
        return;

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
    Compound& D = C.class_compounds.find(S_InChI)->second;
    System& S_in_network = D.compound_systems.find(S_system_ID)->second;

    N.compound_unexplored_systems.pop_back();

    size_t number_of_neighbours = param.number_of_compound_neighbours_distribution(S_in_network);
    for (size_t i = 0; i < number_of_neighbours; i ++)
        generate_compound_neighbour(param,N,C,S_in_network);

    S_in_network.compound_explored = true;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void generate_class_neighbour(Parameters param, Network& N, Class& C,System& R)
{
    string ID = generate_system_ID(param,R.atoms,R.n_electrons,R.charge);
    System P(ID,R.atoms,R.n_electrons,R.charge);

    if (P == R)
        return;

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
    Class& C = N.classes.find(S_class_ID)->second;
    Compound& D = C.class_compounds.find(S_InChI)->second;
    System& S_in_network = D.compound_systems.find(S_system_ID)->second;

    N.class_unexplored_systems.pop_back();

    size_t number_of_neighbours = param.number_of_class_neighbours_distribution(S_in_network);
    for (size_t i = 0; i < number_of_neighbours; i ++)
        generate_class_neighbour(param,N,C,S_in_network);

    S_in_network.class_explored = true;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void choose_a_pair_of_systems_for_hyperedge(Parameters param,Network& N,multiset<System>& R)
{
    static mt19937 generator(param.seed);
    bool hyperedge_already_exists;

    uniform_int_distribution<int> random_class_index(0,N.classes.size() - 1);

    // Premier réactif
    int class_index_1 = random_class_index(generator);
    Class& C_1 = next(N.classes.begin(),class_index_1)->second;

    uniform_int_distribution<int> random_compound_index_1(0,C_1.class_compounds.size() - 1);
    int compound_index_1 = random_compound_index_1(generator);
    Compound& D_1 = next(C_1.class_compounds.begin(),compound_index_1)->second;

    uniform_int_distribution<int> random_system_index_1(0,D_1.compound_systems.size() - 1);
    int system_index_1 = random_system_index_1(generator);
    System& S_1 = next(D_1.compound_systems.begin(),system_index_1)->second;

    R.insert(S_1);

    // Deuxième réactif
    int class_index_2 = random_class_index(generator);
    Class& C_2 = next(N.classes.begin(),class_index_2)->second;

    uniform_int_distribution<int> random_compound_index_2(0,C_2.class_compounds.size() - 1);
    int compound_index_2 = random_compound_index_2(generator);
    Compound& D_2 = next(C_2.class_compounds.begin(),compound_index_2)->second;

    uniform_int_distribution<int> random_system_index_2(0,D_2.compound_systems.size() - 1);
    int system_index_2 = random_system_index_2(generator);
    System& S_2 = next(D_2.compound_systems.begin(),system_index_2)->second;

    multiset<System>::iterator S_2_iterator = R.insert(S_2);
}

/////////////////////////////////////////////////////////////////////

System& choose_system_for_splitting_hyperedge(Parameters param,Network& N)
{
    static mt19937 generator(param.seed);
    bool system_cant_be_split;
    System* pR;

    do
    {
        system_cant_be_split = false;

        uniform_int_distribution<int> random_class_index(0,N.classes.size() - 1);
        int class_index = random_class_index(generator);
        Class& C = next(N.classes.begin(),class_index)->second;

        uniform_int_distribution<int> random_compound_index(0,C.class_compounds.size() - 1);
        int compound_index = random_compound_index(generator);
        Compound& D = next(C.class_compounds.begin(),compound_index)->second;

        uniform_int_distribution<int> random_system_index(0,D.compound_systems.size() - 1);
        int system_index = random_system_index(generator);
        pR = &next(D.compound_systems.begin(),system_index)->second;

        // Vérification que le système peut être splitté
        if (pR->atoms.size() == 1)
        {
            system_cant_be_split = true;
        }
    }
    while (system_cant_be_split == true);

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
    static mt19937 generator(param.seed);

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
    Network N = generate_synthetic_network(param,initial_systems);
    return N;
}
/////////////////////////////////////////////////////////////////////

Network generate_synthetic_network(Parameters param, deque<System>& initial_systems)
{
    Network N(initial_systems);

    for (size_t i = 0; i < param.number_of_generation_rounds; i++)
    {
        // Pour les systèmes non classe-explorés, génération des voisins extra-composé
        while (N.class_unexplored_systems.size() > 0)
            generate_class_neighbourhood(param,N,N.class_unexplored_systems.back());

        // Pour un certain nombre de paires de systèmes, génération des voisins extra-classe en se combinant
        size_t number_of_pairs = max(N.number_of_systems * param.percentage_of_pairs_per_round,param.minimum_number_of_pairs_per_round);
        for (size_t j = 0; j < number_of_pairs; j++)
        {
            multiset<System> R;
            choose_a_pair_of_systems_for_hyperedge(param,N,R);
            generate_hyperedge_from_pair(param,N,R);
        }

        // Nombre de systèmes qui peuvent être splittés
        size_t number_of_splittable_systems = 0;
        for (map<string,Class>::iterator it = N.classes.begin(); it != N.classes.end(); it++)
        {
            Class& C = it->second;
            for (map<string,Compound>::iterator it2 = C.class_compounds.begin(); it2 != C.class_compounds.end(); it2++)
            {
                Compound& D = it2->second;
                for (map<string,System>::iterator it3 = D.compound_systems.begin(); it3 != D.compound_systems.end(); it3++)
                {
                    System& S = it3->second;
                    if (S.atoms.size() > 1)
                        number_of_splittable_systems++;
                }
            }
        }
        // Pour un certain nombre de systèmes, génération des voisins extra-classe en se splittant
        size_t number_of_systems_to_be_split = max(N.number_of_systems * param.percentage_of_splittings_per_round,param.minimum_number_of_splittings_per_round);
        for (size_t j = 0; j < min(number_of_splittable_systems,number_of_systems_to_be_split); j++)
        {
            System R = choose_system_for_splitting_hyperedge(param,N);
            generate_splitting_hyperedge(param,N,R);
        }

        // Pour les systèmes non composé-explorés, génération des voisins intra-composé
        while (N.compound_unexplored_systems.size() > 0)
            generate_compound_neighbourhood(param,N,N.compound_unexplored_systems.back());
    }
    return N;
}
