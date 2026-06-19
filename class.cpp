# include "class.h"
# include "calculate_IDs.h"
using namespace std;
# include <iostream>

Class::Class() : n_electrons(0), number_of_systems_in_class(0){}

/////////////////////////////////////////////////////////////////////

Class::Class(size_t& number_of_systems,vector<System>& compound_unexplored_systems,vector<System>& class_unexplored_systems,const string ID,System& S) : class_ID(ID), atoms(S.atoms), n_electrons(S.n_electrons)
{
    number_of_systems_in_class = 0;
    add_isolated_system_in_class(number_of_systems,compound_unexplored_systems,class_unexplored_systems,S);
}

/////////////////////////////////////////////////////////////////////

void Class::add_system_in_class_from_edge(size_t& number_of_systems,vector<System>& compound_unexplored_systems,vector<System>& class_unexplored_systems,System& R, System& P, float RP_barrier, float PR_barrier)
{
    bool new_product_system = true; // Indique si le système P existait déjà dans le réseau
    System* already_existing_P; // si P existe déjà dans le réseau, pointeur vers P

// Ajout du système dans le bon composé
    const string ID = calculate_InChI(P);
    map<const string,Compound>::iterator it = class_compounds.find(ID);

    if (it == class_compounds.end()) // Si le composé n'existe pas
    {
        Compound C(ID);

        P.insertion_rank_in_class = number_of_systems_in_class;
        number_of_systems_in_class++;

        P.insertion_rank_in_network = number_of_systems;
        number_of_systems++;

        update_system_ID(P,C);

        C.compound_systems.insert(make_pair(P.system_ID,P));
        class_compounds.insert(make_pair(ID,C));

        compound_unexplored_systems.push_back(P);
        class_unexplored_systems.push_back(P);
    }
    else    // le composé existe
    {
        Compound& C = it->second;
        /* Version avec same_system()
        for (map<string,System>::iterator it = C.compound_systems.begin(); it != C.compound_systems.end(); it++)
        {
            if (same_system(P,it->second) == true)
            {
                new_product_system = false;
                already_existing_P = &(it->second);
                break;
            }
        }

        if (new_product_system == true) // Si le système n'existe pas déjà
        {
            P.insertion_rank_in_class = number_of_systems_in_class;
            number_of_systems_in_class++;

            P.insertion_rank_in_network = number_of_systems;
            number_of_systems++;

            update_system_ID(P,C);

            C.compound_systems.insert(make_pair(P.system_ID,P));

            compound_unexplored_systems.push_back(P);
            class_unexplored_systems.push_back(P);
        }

        */

        ////// Version sans same_system()
        map<const string,System>::iterator itt = C.compound_systems.find(P.system_ID);
        if (itt == C.compound_systems.end())
        {
            P.insertion_rank_in_class = number_of_systems_in_class;
            number_of_systems_in_class++;

            P.insertion_rank_in_network = number_of_systems;
            number_of_systems++;

            update_system_ID(P,C);

            C.compound_systems.insert(make_pair(P.system_ID,P));

            compound_unexplored_systems.push_back(P);
            class_unexplored_systems.push_back(P);
        }
        else
        {
            new_product_system = false;
            already_existing_P = &(itt->second);
        }
        //////
    }

// Ajout des arêtes dans la matrice d'adjacence
if (number_of_systems_in_class == 1)    // Si c'est le premier système dans la classe
{
    vector<float> v = {-1};
    edges.push_back(v);
    edges_flux.push_back(v);
}

else
{
    if (new_product_system == true) // Si un nouveau système a été ajouté
    {
        // Ajout de la dernière colonne
        for (size_t i = 0; i < edges.size(); i++)
        {
            if (i == R.insertion_rank_in_class)
            {
                edges[i].push_back(RP_barrier);
            }
            else
            {
                edges[i].push_back(-1);
            }
            edges_flux[i].push_back(-1);    ////// bad alloc
        }

        // Ajout de la dernière ligne
        vector<float> v;
        edges.push_back(v);
        edges_flux.push_back(v);
        for (size_t i = 0; i < edges.size(); i++)
        {
            if (i == R.insertion_rank_in_class)
                edges[edges.size()-1].push_back(PR_barrier);
            else
            {
                edges[edges.size()-1].push_back(-1);
            }
            edges_flux[edges.size()-1].push_back(-1);
        }
    }

    else // Si le système existait déjà, on ajoute l'arête
    {
        size_t i = R.insertion_rank_in_class;
        size_t j = already_existing_P->insertion_rank_in_class;
        edges[i][j] = RP_barrier;
        edges[j][i] = PR_barrier;
    }
}
}

/////////////////////////////////////////////////////////////////////

void Class::add_isolated_system_in_class(size_t& number_of_systems,vector<System>& compound_unexplored_systems,vector<System>& class_unexplored_systems,System& S)
{
    add_system_in_class_from_edge(number_of_systems,compound_unexplored_systems,class_unexplored_systems,S,S,-1,-1);
}
