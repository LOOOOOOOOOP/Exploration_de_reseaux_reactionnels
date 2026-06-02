#include "network.h"
#include "Testing/print_tests.h"
#include "calculate_IDs.h"  // pour tester avec print_class_matrix

#include <iostream>
#include <queue>
//#include <tuple>


using namespace std;

int main()
{
    // Création de systèmes moléculaires
    multiset<string> atoms_S = {"H","H","C","C","C","H","H"};
    System S("H2CCCH2 - 0",atoms_S,16);

    multiset<string> atoms_T = {"C","H","C","H","C","H","H"};
    System T("CHCHCH2 - 0",atoms_T,16);

    multiset<string> atoms_U = {"C","H","C","H","C","H","H"};
    System U("CHCHCH2 - 1",atoms_U,16);

    multiset<string> atoms_V = {"H","H"};
    System V("H2 - 0",atoms_V,2);

    multiset<string> atoms_W = {"H","H","O"};
    System W("H2O - 0",atoms_W,8);

    multiset<string> atoms_X = {"H"};
    System X("H+ - 0",atoms_X,0);

    // Regroupement des systèmes dans un deque initial
    deque<System> initial_systems = {S,T,U,V,W,X};

    // Création du réseau initial
    Network reaction_network(initial_systems);

    // Récupération de T
    Class& C = reaction_network.classes.begin()->second;
    Compound& D = C.class_compounds.begin()->second;
    System& T_in_network = D.compound_systems.begin()->second;

    // Récupération de U
    System& U_in_network = D.compound_systems.find(U.system_ID)->second;

    // Ajout d'une arête entre T et U
    reaction_network.add_system_in_network_from_edge(C,T_in_network,U_in_network,12.74,5.8);

    // Ajout d'une hyperarête entre V et X
    multiset<System> reactants = {V};
    multiset<System> products = {X,X};
    reaction_network.add_system_in_network_from_hyperedge(reactants,products,24.0,2.34);

    // Découverte d'un conformère Y à partir de S
    System Y("H2CCCH2 - 1",atoms_S,16);
    reaction_network.add_system_in_network_from_edge(C,S,Y,0.3,0.2);

    // Découverte d'un isomère Alpha à partir de S
    System Alpha("CHCHCH2 - 2",atoms_U,16);
    reaction_network.add_system_in_network_from_edge(C,S,Alpha,12,6);

    // Découverte de produit Z et re-découverte de produit X à partir de W
    System Z("OH- - 0",{"O","H"},8);
    multiset<System> reactants2 = {W};
    multiset<System> products2 = {X,Z};
    reaction_network.add_system_in_network_from_hyperedge(reactants2,products2,4,9.5);

    // Découverte de Beta à partir de W et X
    multiset<string> atoms_Beta = {"H","H","H","O"};
    System Beta("H3O+ - 0",atoms_Beta,8);
    multiset<System> reactants3 = {W,X};
    multiset<System> products3 = {Beta};
    reaction_network.add_system_in_network_from_hyperedge(reactants3,products3,2.33,9.41);

    // Affichage
    print_network(reaction_network);
    print_class_matrix(C);
    print_hyperedges(reaction_network);

    cout << endl << "Hello reaction network!" << endl;
    return 0;
}
