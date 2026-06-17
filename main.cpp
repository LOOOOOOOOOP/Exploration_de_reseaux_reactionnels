#include "network.h"
#include "Testing/print_tests.h"
#include "Testing/parameters.h"
#include "Testing/generate_synthetic_network.h"
#include "calculate_IDs.h"  // pour tester avec print_class_matrix

#include <iostream>
#include <queue>
//#include <tuple>


using namespace std;

int main()
{
/////////////////// Tests de structure ///////////////////
/*
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
    Class& C_T = reaction_network.classes.begin()->second;
    Compound& D_T = C_T.class_compounds.begin()->second;
    System& T_in_network = D_T.compound_systems.begin()->second;

    // Ajout d'une arête entre T et U
    reaction_network.add_system_in_network_from_edge(C_T,T_in_network,U,12.74,5.8);

    // Récupération de V
    Class& C_V = reaction_network.classes.find(calculate_class_ID(V))->second;
    Compound& D_V = C_V.class_compounds.find(calculate_InChI(V))->second;
    System& V_in_network = D_V.compound_systems.find(V.system_ID)->second;

    // Ajout d'une hyperarête entre V et X,X
    multiset<System> reactants = {V_in_network};
    multiset<System> products = {X,X};
    reaction_network.add_system_in_network_from_hyperedge(reactants,products,24.0,2.34);

    // Récupération de S
    Class& C_S = reaction_network.classes.find(calculate_class_ID(S))->second;
    Compound& D_S = C_S.class_compounds.find(calculate_InChI(S))->second;
    System& S_in_network = D_S.compound_systems.find(S.system_ID)->second;

    // Découverte d'un conformère Y à partir de S
    System Y("H2CCCH2 - 1",atoms_S,16);
    reaction_network.add_system_in_network_from_edge(C_S,S_in_network,Y,0.3,0.2);

    // Redécouverte de T à partir de S
    reaction_network.add_system_in_network_from_edge(C_S,S_in_network,T,12,6);

    // Récupération de W
    Class& C_W = reaction_network.classes.find(calculate_class_ID(W))->second;
    Compound& D_W = C_W.class_compounds.find(calculate_InChI(W))->second;
    System& W_in_network = D_W.compound_systems.find(W.system_ID)->second;

    // Découverte de produit Z et re-découverte de produit X à partir de W
    System Z("OH- - 0",{"O","H"},8);
    multiset<System> reactants2 = {W_in_network};
    multiset<System> products2 = {X,Z};
    reaction_network.add_system_in_network_from_hyperedge(reactants2,products2,4,9.5);

    // Récupération de X
    Class& C_X = reaction_network.classes.find(calculate_class_ID(X))->second;
    Compound& D_X = C_X.class_compounds.find(calculate_InChI(X))->second;
    System& X_in_network = D_X.compound_systems.find(X.system_ID)->second;

    // Découverte de Alpha à partir de W et X
    multiset<string> atoms_Alpha = {"H","H","H","O"};
    System Alpha("H3O+ - 0",atoms_Alpha,8);
    multiset<System> reactants3 = {W_in_network,X_in_network};
    multiset<System> products3 = {Alpha};
    reaction_network.add_system_in_network_from_hyperedge(reactants3,products3,2.33,9.41);

    // Affichage
    print_network(reaction_network);
    print_class_matrix(C_T);
    print_hyperedges(reaction_network);
*/

//////////////// Tests de réseaux synthétiques ////////////////

    Parameters param;
    param.seed = time(0);
    param.initial_systems_are_from_different_compounds = false;
    //System S = generate_system(param);
    //print_system(S);

    Network reaction_network = generate_synthetic_network(param);
    print_network(reaction_network);


    cout << endl << "Hello reaction network!" << endl;
    return 0;
}
