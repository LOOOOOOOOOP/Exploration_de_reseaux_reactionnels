#include "explore_synthetic_network.h"
#include "../calculate_IDs.h"

using namespace std;


// À implémenter
void explore(Network& N, const Network& synth_N)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Non fonctionnel
void find_compound_neighbour(Network& N, const Network& synth_N, System& R)
{
    string system_ID = R.system_ID;
    string compound_ID = calculate_compound_ID(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    const Class& complete_class = synth_N.classes.find(class_ID)->second;
    const Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    for (map<string,System>::const_iterator it = complete_compound.compound_systems.begin(); it != complete_compound.compound_systems.end(); it++)
    {
        if (R_compound.compound_systems.find(it->second.system_ID) != R_compound.compound_systems.end())
        {
            if (R_class.edges[R.insertion_rank_in_class][it->second.insertion_rank_in_class] != -1)
                continue;   // Si le système a déjà été découvert et qu'il y a une arête entre les deux, on continue à chercher un nouveau voisin
        }

        // Le système n'a pas été découvert ou il y a une arête à ajouter
        System P = it->second;
        float barrier_RP = complete_class.edges[R.insertion_rank_in_class][P.insertion_rank_in_class];
        float barrier_PR = complete_class.edges[P.insertion_rank_in_class][R.insertion_rank_in_class];
        N.add_system_in_network_from_edge(R_class,R,P,barrier_RP,barrier_PR);
        break;
    }
}

/////////////////////////////////////////////////////////////////////

// Non fonctionnel
void find_all_compound_neighbours(Network& N, const Network& synth_N, System& R)
{
    string system_ID = R.system_ID;
    string compound_ID = calculate_compound_ID(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    const Class& complete_class = synth_N.classes.find(class_ID)->second;
    const Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    while (R_compound.compound_systems.size() < complete_compound.compound_systems.size())
        find_compound_neighbour(N,synth_N,R);

    N.compound_unexplored_systems.pop_back();
    R.compound_explored = true;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Non fonctionnel
void find_class_neighbour(Network& N, const Network& synth_N, System& R)
{
/*
    string system_ID = R.system_ID;
    string compound_ID = calculate_compound_ID(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    const Class& complete_class = synth_N.classes.find(class_ID)->second;
    const Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    bool stop = false;
    map<const string,Class>::const_iterator it = synth_N.classes.begin();
    while (it != synth_N.classes.end() && stop == false)
    {
        map<const string,Compound>::const_iterator it2 = it->second.class_compounds.begin();
        while (it2 != it->second.class_compounds.end() && stop == false)
        {
            if (it2->second != R_compound)
            {
                for (map<const string,System>::const_iterator it3 = complete_compound.compound_systems.begin(); it3 != complete_compound.compound_systems.end(); it3++)
                {
                    if ((R_compound.compound_systems.find(calculate_compound_ID(it3->second)) == R_compound.compound_systems.end()) || R_class.edges[R.insertion_rank_in_class][it3->second.insertion_rank_in_class] == -1))   // Si le voisin potentiel n'a pas déjà été découvert ou qu'il n'y a pas d'arête entre les deux
                    {
                        System& P = it3->second;
                        float barrier_RP = complete_class.edges[R.insertion_rank_in_class][P.insertion_rank_in_class];
                        float barrier_PR = complete_class.edges[P.insertion_rank_in_class][R.insertion_rank_in_class];
                        N.add_system_in_network_from_edge(R_class,R,P,barrier_RP,barrier_PR);
                        stop = true;
                        break;
                    }
                }
            }
            it2++;
        }
        it++;
    }
*/
}

/////////////////////////////////////////////////////////////////////

// À implémenter
void find_all_class_neighbours(Network& N, const Network& synth_N, System& R)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// À implémenter
void find_hyperedge_from_pair(Network& N, const Network& synth_N, std::pair<System,System>& R)
{
    return;
}

/////////////////////////////////////////////////////////////////////

// À implémenter
void find_splitting_hyperedge_from_system(Network& N, const Network& synth_N, System& S)
{
    return;
}
