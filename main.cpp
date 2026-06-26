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
    // Atomes
    Atom H("H",1,0.6);
    Atom C("C",6,0.25);
    Atom O("O",8,0.05);

    // Création de systèmes moléculaires
    multiset<Atom> atoms_S = {H,H,C,C,C,H,H};
    System S("H2CCCH2",atoms_S,16,0);

    multiset<Atom> atoms_T = {C,H,C,H,C,H,H};
    System T("CHCHCH2",atoms_T,16,0);

    multiset<Atom> atoms_U = {C,H,C,H,C,H,H};
    System U("CHCHCH2",atoms_U,16,0);

    multiset<Atom> atoms_V = {H,H};
    System V("H2",atoms_V,2,0);

    multiset<Atom> atoms_W = {H,H,O};
    System W("H2O",atoms_W,8,0);

    multiset<Atom> atoms_X = {H};
    System X("H+",atoms_X,0,1);

    // Regroupement des systèmes dans un deque initial
    deque<System> initial_systems = {S,T,U,V,W,X};

    // Création du réseau initial
    Network reaction_network(initial_systems);

    // Récupération de T
    Class& C_T = reaction_network.classes.begin()->second;
    Compound& D_T = C_T.class_compounds.begin()->second;
    System& T_in_network = D_T.compound_systems.begin()->second;

    // Récupération de U
    Class& C_U = reaction_network.classes.find(calculate_class_ID(U))->second;
    Compound& D_U = C_U.class_compounds.find(calculate_InChI(U))->second;
    System& U_in_network = next(D_U.compound_systems.begin(),1)->second;

    // Ajout d'une arête entre T et U
    reaction_network.add_system_in_network_from_edge(C_T,T_in_network,U_in_network,12.74,5.8);

    // Récupération de V
    Class& C_V = next(reaction_network.classes.begin(),2)->second;
    Compound& D_V = C_V.class_compounds.begin()->second;
    System& V_in_network = D_V.compound_systems.begin()->second;

    // Ajout d'une hyperarête entre V et Beta,Beta
    multiset<System> reactants = {V_in_network};
    System Beta("H",atoms_X,1,0);
    multiset<System> products = {Beta,Beta};
    reaction_network.add_system_in_network_from_hyperedge(reactants,products,24.0,2.34);

    // Récupération de S
    Class& C_S = reaction_network.classes.begin()->second;
    Compound& D_S = next(C_S.class_compounds.begin(),1)->second;
    System& S_in_network = D_S.compound_systems.begin()->second;

    // Découverte d'un conformère Y à partir de S
    System Y("H2CCCH2",atoms_S,16,0);
    reaction_network.add_system_in_network_from_edge(C_S,S_in_network,Y,0.3,0.2);

    // Redécouverte de T à partir de S
    reaction_network.add_system_in_network_from_edge(C_S,S_in_network,T_in_network,12,6);

    // Récupération de W
    Class& C_W = next(reaction_network.classes.begin(),4)->second;
    Compound& D_W = C_W.class_compounds.begin()->second;
    System& W_in_network = D_W.compound_systems.begin()->second;

    // Récupération de X
    Class& C_X = next(reaction_network.classes.begin(),1)->second;
    Compound& D_X = C_X.class_compounds.begin()->second;
    System& X_in_network = D_X.compound_systems.begin()->second;

    // Découverte de produit Z et re-découverte de produit X à partir de W
    System Z("OH-",{O,H},8,-1);
    multiset<System> reactants2 = {W_in_network};
    multiset<System> products2 = {X_in_network,Z};
    reaction_network.add_system_in_network_from_hyperedge(reactants2,products2,4,9.5);

    // Découverte de Alpha à partir de W et X
    multiset<Atom> atoms_Alpha = {H,H,H,O};
    System Alpha("H3O+",atoms_Alpha,8,1);
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
    param.initial_systems_are_from_different_compounds = true;
    param.max_number_of_atoms_in_system = 14;
    param.number_of_generation_rounds = 2;  // bug si plus que 2

    //System S = generate_system(param);
    //print_system(S);

    Network reaction_network = generate_synthetic_network(param);
    print_network(reaction_network,false);
    print_hyperedges(reaction_network);


    cout << endl << "Hello reaction network!" << endl;
    return 0;
}
