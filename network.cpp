#include "network.h"
#include "calculate_IDs.h"

#include "assert.h"
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
/////////////////////////////////////////////////////////////////////

void Network::add_isolated_system_to_network(System& S)
{
    const string ID = calculate_class_ID(S);
    map<const string,Class>::iterator it = classes.find(ID);

    if (it == classes.end()) // Si la classe n'existe pas
        classes.insert(make_pair(ID,Class(number_of_systems,compound_unexplored_systems,class_unexplored_systems,ID,S)));
    else
        it->second.add_isolated_system_in_class(number_of_systems,compound_unexplored_systems,class_unexplored_systems,S);
}

/////////////////////////////////////////////////////////////////////

void Network::add_system_in_network_from_edge(Class& C, System& R, System& P, float RP_barrier, float PR_barrier)
{
    assert (calculate_class_ID(P) == C.class_ID);

    C.add_system_in_class_from_edge(number_of_systems,compound_unexplored_systems,class_unexplored_systems,R,P,RP_barrier,PR_barrier);
}

/////////////////////////////////////////////////////////////////////

void Network::add_system_in_network_from_hyperedge(multiset<System>& R,multiset<System>& P,float RP_barrier,float PR_barrier)
{
    // Ajout des systèmes
    for (set<System>::iterator i = P.begin(); i != P.end(); i++) // Pour chaque produit
    {
        System Pi = *i;
        string class_ID_P = calculate_class_ID(Pi);
        string InChI_P = calculate_InChI(Pi);

        map<const string,Class>::iterator it = classes.find(class_ID_P);
        if (it == classes.end()) // Si la classe n'existe pas
        {
            Class C(number_of_systems,compound_unexplored_systems,class_unexplored_systems,class_ID_P,Pi);
            classes.insert(make_pair(class_ID_P,C));
            Compound D(InChI_P);
            //update_system_ID(Pi,D);
            D.compound_systems.insert(make_pair(Pi.system_ID,Pi));
            C.class_compounds.insert(make_pair(InChI_P,D));

            compound_unexplored_systems.push_back(Pi);
            class_unexplored_systems.push_back(Pi);

            P.erase(i);
            P.insert(Pi);
        }
        else // Si la classe existe déjà
        {
            Class& C = it->second;
            map<const string,Compound>::iterator it2 = C.class_compounds.find(InChI_P);
            if (it2 == C.class_compounds.end()) // Si le composé n'existe pas
            {
                (Pi).insertion_rank_in_class = C.number_of_systems_in_class;
                C.number_of_systems_in_class++;

                (Pi).insertion_rank_in_network = number_of_systems;
                number_of_systems++;

                Compound D(InChI_P);
                //update_system_ID(Pi,D);
                D.compound_systems.insert(make_pair(Pi.system_ID,Pi));
                C.class_compounds.insert(make_pair(InChI_P,D));

                compound_unexplored_systems.push_back(Pi);
                class_unexplored_systems.push_back(Pi);

                P.erase(i);
                P.insert(Pi);
            }
            else // Si le composé existe déjà
            {
                Compound &D = it2->second;

                /* À garder lorsque la fonction same_system sera fonctionnelle
                bool new_product_system = true;

                for (map<const string,System>::iterator it3 = D.compound_systems.begin(); it3 != D.compound_systems.end(); it3++)
                {
                    if (same_system(Pi,it3->second) == true)
                    {
                        new_product_system = false;

                        P.erase(i);
                        P.insert(it3->second);

                        break;
                    }
                }

                if (new_product_system == true) // Si le système n'existe pas déjà
                {
                    (Pi).insertion_rank_in_class = C.number_of_systems_in_class;
                    C.number_of_systems_in_class++;

                    (Pi).insertion_rank_in_network = number_of_systems;
                    number_of_systems++;

                    update_system_ID(Pi,D);

                    D.compound_systems.insert(make_pair(Pi.system_ID,Pi));

                    compound_unexplored_systems.push_back(Pi);
                    class_unexplored_systems.push_back(Pi);

                    P.erase(i);
                    P.insert(Pi);
                }
                */
                ////// À enlever lorsque same_system fonctionnera
                map<const string,System>::iterator it3 = D.compound_systems.find(Pi.system_ID);
                if (it3 == D.compound_systems.end()) // Si le système n'existe pas
                {
                    (Pi).insertion_rank_in_class = C.number_of_systems_in_class;
                    C.number_of_systems_in_class++;

                    (Pi).insertion_rank_in_network = number_of_systems;
                    number_of_systems++;

                    D.compound_systems.insert(make_pair(Pi.system_ID,Pi));

                    compound_unexplored_systems.push_back(Pi);
                    class_unexplored_systems.push_back(Pi);

                    P.erase(i);
                    P.insert(Pi);
                }
                else
                {
                    P.erase(i);
                    P.insert(it3->second);
                }
                //////
            }
        }
    }

    // Ajout de l'hyperarête
    hyperedges.insert(Hyperedge(R,P,make_pair(RP_barrier,PR_barrier)));
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

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
