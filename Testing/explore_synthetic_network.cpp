#include "explore_synthetic_network.h"
#include "../calculate_IDs.h"

using namespace std;
/*
void explore(Network& N,const Network synth_N)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void find_compound_neighbour(Network& N,const Network synth_N,System& R)
{
    string system_ID = R.system_ID;
    string compound_ID = calculate_InChI(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    Class& complete_class = synth_N.classes.find(class_ID)->second;
    Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    for (map<string,System>::iterator it = complete_compound.compound_systems.begin(); it != complete_compound.compound_systems.end(); it++)
    {
        if ((R_compound.compound_systems.find(calculate_InChI(it->second)) == R_compound.compound_systems.end()) || (R_class.edges[R.insertion_rank_in_class][it->second.insertion_rank_in_class] == -1))   // Si le voisin potentiel n'a pas déjà été découvert ou qu'il n'y a pas d'arête entre les deux
        {
            System& P = it->second;
            float barrier_RP = complete_class.edges[R.insertion_rank_in_class][P.insertion_rank_in_class];
            float barrier_PR = complete_class.edges[P.insertion_rank_in_class][R.insertion_rank_in_class];
            N.add_system_in_network_from_edge(R_class,R,P,barrier_RP,barrier_PR);
            break;
        }
    }
}

/////////////////////////////////////////////////////////////////////

void find_all_compound_neighbours(Network& N,const Network synth_N,System& R)
{
    string system_ID = R.system_ID;
    string compound_ID = calculate_InChI(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    Class& complete_class = synth_N.classes.find(class_ID)->second;
    Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    while (R_compound.compound_systems.size() < complete_compound.compound_systems.size())
        find_compound_neighbour(N,synth_N,R);
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void find_class_neighbour(Network& N,const Network synth_N,System& S)
{
    string system_ID = R.system_ID;
    string compound_ID = calculate_InChI(R);
    string class_ID = calculate_class_ID(R);

    Class& R_class = N.classes.find(class_ID)->second;
    Compound& R_compound = R_class.class_compounds.find(compound_ID)->second;
    Class& complete_class = synth_N.classes.find(class_ID)->second;
    Compound& complete_compound = complete_class.class_compounds.find(compound_ID)->second;

    bool stop = false;
    map<string,Class>::iterator it = synth_N.classes.begin();
    while (it != synth_N.classes.end() && stop == false)
    {
        map<string,Compound>::iterator it2 = it->second.class_compounds.begin();
        while (it2 != it->second.class_compounds.end() && stop == false)
        {
            if (it2->second != R_compound)
            {
                for (map<string,System>::iterator it3 = complete_compound.compound_systems.begin(); it3 != complete_compound.compound_systems.end(); it3++)
                {
                    if ((R_compound.compound_systems.find(calculate_InChI(it3->second)) == R_compound.compound_systems.end()) || R_class.edges[R.insertion_rank_in_class][it3->second.insertion_rank_in_class] == -1))   // Si le voisin potentiel n'a pas déjà été découvert ou qu'il n'y a pas d'arête entre les deux
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
}

/////////////////////////////////////////////////////////////////////

void find_all_class_neighbours(Network& N,const Network synth_N,System& S)
{
    return;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void find_hyperedge_from_pair(Network& N,const Network synth_N,std::pair<System,System>& R)
{
    return;
}

/////////////////////////////////////////////////////////////////////

void find_splitting_hyperedge_from_system(Network& N,const Network synth_N,System& S)
{
    return;
}
*/
