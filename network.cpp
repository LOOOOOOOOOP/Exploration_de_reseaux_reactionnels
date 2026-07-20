#include "assert.h"
#include "calculate_IDs.h"
#include "network.h"

#include <iostream>

using namespace std;


Network::Network() : number_of_systems(0){}

/////////////////////////////////////////////////////////////////////

Network::Network(deque<System>& initial_systems) : number_of_systems(0)
{
    while (!initial_systems.empty())
    {
        add_isolated_system_to_network(initial_systems.front());
        initial_systems.pop_front();
    }
}

/////////////////////////////////////////////////////////////////////

void Network::add_isolated_system_to_network(System& S)
{
    const string ID = calculate_class_ID(S);
    map<const string,Class>::iterator it = classes.find(ID);

    if (it == classes.end()) // Si la classe n'existe pas
        classes.insert(make_pair(ID,Class(number_of_systems,compound_unexplored_systems,class_unexplored_systems,ID,S)));
    else    // La classe existe
        it->second.add_isolated_system_in_class(number_of_systems,compound_unexplored_systems,class_unexplored_systems,S);
}

/////////////////////////////////////////////////////////////////////

void Network::add_system_in_network_from_edge(Class& C, const System& R, System& P, const float& RP_barrier, const float& PR_barrier)
{
    assert (calculate_class_ID(P) == C.class_ID);
    C.add_system_in_class_from_edge(number_of_systems,compound_unexplored_systems,class_unexplored_systems,R,P,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////

void Network::add_system_in_network_from_hyperedge(const multiset<System>& R, multiset<System>& P, const float& RP_barrier, const float& PR_barrier)
{
    // Vérification que l'hyperarête n'existe pas déjà
    bool combination_hyperedge = R.size() > 1;  // Vrai s'il s'agit d'une hyperarête de combinaison (et non de fragmentation)
    for (set<Hyperedge>::iterator it4 = hyperedges.begin(); it4 != hyperedges.end(); it4++)
    {
        if (combination_hyperedge == true && (it4->reactants == R || it4->products == R))
            return; // Un ensemble de systèmes ne peut s'unir que d'une seule façon
        if ((it4->reactants == R && it4->products == P) || (it4->reactants == P && it4->products == R))
            return; // La même hyperarête existe déjà
    }

    // Ajout des produits au réseau
    for (set<System>::iterator i = P.begin(); i != P.end(); i++)
    {
        System P_i = *i;
        add_isolated_system_to_network(P_i);
    }

    // Ajout de l'hyperarête
    hyperedges.insert(Hyperedge(R,P,make_pair(RP_barrier,PR_barrier)));
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Fonctions à implémenter (simulation numérique)

void Network::simulate(vector<float> initial_concentrations)
{
    update_systems_concentration_function(initial_concentrations);
    update_edges_flux();
}

/////////////////////////////////////////////////////////////////////

void Network::update_systems_concentration_function(vector<float> initial_concentrations)
{
    return;
}

/////////////////////////////////////////////////////////////////////

void Network::update_edges_flux()
{
    return;
}
