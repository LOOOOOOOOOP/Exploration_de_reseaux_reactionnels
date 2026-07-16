#include "network.h"
#include "Testing/print_tests.h"
#include "Testing/parameters.h"
#include "Testing/generate_synthetic_network.h"
#include "calculate_IDs.h"
#include "Testing/output_results.h"

#include <iostream>
#include <queue>

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
    multiset<Atom> atoms_CHCHCH2comp1conf1 = {H,H,C,C,C,H,H};
    System CHCHCH2comp1conf1("CHCHCH2/compound#1/[1]",atoms_CHCHCH2comp1conf1,16,0);

    multiset<Atom> atoms_CHCHCH2comp2conf1 = {C,H,C,H,C,H,H};
    System CHCHCH2comp2conf1("CHCHCH2/compound#2/[1]",atoms_CHCHCH2comp2conf1,16,0);

    multiset<Atom> atoms_CHCHCH2comp2conf2 = {C,H,C,H,C,H,H};
    System CHCHCH2comp2conf2("CHCHCH2/compound#2/[2]",atoms_CHCHCH2comp2conf2,16,0);

    multiset<Atom> atoms_H2 = {H,H};
    System H2("H2/compound#1/[1]",atoms_H2,2,0);

    multiset<Atom> atoms_H2O = {H,H,O};
    System H2O("H2O/compound#1/[1]",atoms_H2O,8,0);

    multiset<Atom> atoms_Hplus = {H};
    System Hplus("H/q+1/compound#1/[1]",atoms_Hplus,0,1);

    // Regroupement des systèmes dans un deque initial
    deque<System> initial_systems = {CHCHCH2comp1conf1,CHCHCH2comp2conf1,CHCHCH2comp2conf2,H2,H2O,Hplus};

    // Création du réseau initial
    Network reaction_network(initial_systems);

    // Récupération de CHCHCH2comp2conf1
    Class& C_CHCHCH2comp2conf1 = reaction_network.classes.begin()->second;
    Compound& D_CHCHCH2comp2conf1 = C_CHCHCH2comp2conf1.class_compounds.begin()->second;
    System& CHCHCH2comp2conf1_in_network = D_CHCHCH2comp2conf1.compound_systems.begin()->second;

    // Récupération de CHCHCH2comp2conf2
    Class& C_CHCHCH2comp2conf2 = reaction_network.classes.find(calculate_class_ID(CHCHCH2comp2conf2))->second;
    Compound& D_CHCHCH2comp2conf2 = C_CHCHCH2comp2conf2.class_compounds.find(calculate_InChI(CHCHCH2comp2conf2))->second;
    System& CHCHCH2comp2conf2_in_network = next(D_CHCHCH2comp2conf2.compound_systems.begin(),1)->second;

    // Ajout d'une arête entre CHCHCH2comp2conf1 et CHCHCH2comp2conf2
    reaction_network.add_system_in_network_from_edge(C_CHCHCH2comp2conf1,CHCHCH2comp2conf1_in_network,CHCHCH2comp2conf2_in_network,12.74,5.8);

    // Récupération de H2
    Class& C_H2 = next(reaction_network.classes.begin(),2)->second;
    Compound& D_H2 = C_H2.class_compounds.begin()->second;
    System& H2_in_network = D_H2.compound_systems.begin()->second;

    // Ajout d'une hyperarête entre H2 et H,H
    multiset<System> reactants = {H2_in_network};
    System H_system("H/compound#1/[1]",atoms_Hplus,1,0);
    multiset<System> products = {H_system,H_system};
    reaction_network.add_system_in_network_from_hyperedge(reactants,products,24.0,2.34);

    // Récupération de CHCHCH2comp1conf1
    Class& C_CHCHCH2comp1conf1 = reaction_network.classes.begin()->second;
    Compound& D_CHCHCH2comp1conf1 = next(C_CHCHCH2comp1conf1.class_compounds.begin(),1)->second;
    System& CHCHCH2comp1conf1_in_network = D_CHCHCH2comp1conf1.compound_systems.begin()->second;

    // Découverte d'un conformère CHCHCH2comp1conf2 à partir de CHCHCH2comp1conf1
    System CHCHCH2comp1conf2("CHCHCH2/compound#1/[2]",atoms_CHCHCH2comp1conf1,16,0);
    reaction_network.add_system_in_network_from_edge(C_CHCHCH2comp1conf1,CHCHCH2comp1conf1_in_network,CHCHCH2comp1conf2,0.3,0.2);

    // Redécouverte de CHCHCH2comp2conf1 à partir de CHCHCH2comp1conf1
    reaction_network.add_system_in_network_from_edge(C_CHCHCH2comp1conf1,CHCHCH2comp1conf1_in_network,CHCHCH2comp2conf1_in_network,12,6);

    // Récupération de H2O
    Class& C_H2O = next(reaction_network.classes.begin(),4)->second;
    Compound& D_H2O = C_H2O.class_compounds.begin()->second;
    System& H2O_in_network = D_H2O.compound_systems.begin()->second;

    // Récupération de Hplus
    Class& C_Hplus = next(reaction_network.classes.begin(),1)->second;
    Compound& D_Hplus = C_Hplus.class_compounds.begin()->second;
    System& Hplus_in_network = D_Hplus.compound_systems.begin()->second;

    // Découverte de produit OHminus et re-découverte de produit Hplus à partir de H2O
    System OHminus("OH/q-1/compound#1/[1]",{O,H},8,-1);
    multiset<System> reactants2 = {H2O_in_network};
    multiset<System> products2 = {Hplus_in_network,OHminus};
    reaction_network.add_system_in_network_from_hyperedge(reactants2,products2,4,9.5);

    // Découverte de H3Oplus à partir de H2O et Hplus
    multiset<Atom> atoms_H3Oplus = {H,H,H,O};
    System H3Oplus("H3O/q+1/compound#1/[1]",atoms_H3Oplus,8,1);
    multiset<System> reactants3 = {H2O_in_network,Hplus_in_network};
    multiset<System> products3 = {H3Oplus};
    reaction_network.add_system_in_network_from_hyperedge(reactants3,products3,2.33,9.41);

    // Affichage
    print_network(reaction_network,true);
    print_class_matrix(C_CHCHCH2comp2conf1);
    print_hyperedges(reaction_network);
*/

//////////////// Tests de réseaux synthétiques ////////////////

    Parameters param;
    Network reaction_network = generate_synthetic_network(param);
    print_network(reaction_network,false);
    print_hyperedges(reaction_network);
    output_all(param,reaction_network,"test");

    cout << endl << "Hello reaction network!" << endl;
    return 0;
}
