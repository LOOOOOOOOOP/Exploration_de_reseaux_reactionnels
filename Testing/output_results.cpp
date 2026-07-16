#include "output_results.h"
#include <ctime>

using namespace std;

void output_parameters(Parameters param, string filename)
{
    ofstream output_file;
    output_file.open(filename,std::ios_base::app);

    for (int i = 0; i < 19; i++)
        output_file << "-";
    output_file << " Parameters ";
    for (int i = 0; i < 19; i++)
        output_file << "-";
    output_file << endl;

    output_file << "Seed: " << param.seed << endl;
    output_file << "Possible atoms: ";
    for (set<Atom>::iterator it = param.possible_atoms.begin(); it != param.possible_atoms.end(); it++)
    {
        output_file << it->symbol << " ";
    }
    output_file << endl;
    output_file << "Min number of atoms in initial systems: " << param.min_number_of_atoms_in_initial_systems << endl;
    output_file << "Max number of atoms in initial systems: " << param.max_number_of_atoms_in_initial_systems << endl;
    output_file << "Number of initial systems: " << param.number_of_initial_systems << endl;
    if (param.initial_systems_are_from_different_compounds)
        output_file << "Initial systems are from different compounds: Yes" << endl;
    else
        output_file << "Initial systems are from different compounds: No" << endl;
    if(param.limit_number_of_conformers_to_5)
        output_file << "Limit number of conformers to 5: Yes" << endl;
    else
        output_file << "Limit number of conformers to 5: No" << endl;
    output_file << "Number of generation rounds: " << param.number_of_generation_rounds << endl;
    output_file << "Percentage of pairs per round: " << param.percentage_of_pairs_per_round << endl;
    output_file << "Minimum number of pairs per round: " << param.minimum_number_of_pairs_per_round << endl;
    output_file << "Percentage of splittings per round: " << param.percentage_of_splittings_per_round << endl;
    output_file << "Minimum number of splittings per round: " << param.minimum_number_of_splittings_per_round << endl;

    for (int i = 0; i < 50; i++)
        output_file << "-";
    output_file << endl << endl;
}

/////////////////////////////////////////////////////////////////////

void output_systems(Network N, string filename)
{
    ofstream output_file;
    output_file.open(filename,std::ios_base::app);

    for (int i = 0; i < 20; i++)
        output_file << "-";
    output_file << " Systems ";
    for (int i = 0; i < 21; i++)
        output_file << "-";
    output_file << endl;

    for (map<const string,Class>::iterator it = N.classes.begin(); it != N.classes.end(); it++)
    {
        output_file<< "Class: " << it->second.class_ID.substr(0,it->second.class_ID.find(" | ")) << endl;
        for (map<const string,Compound>::iterator it2 = it->second.class_compounds.begin(); it2 != it->second.class_compounds.end(); it2++)
        {
            size_t number_of_conformers = it2->second.compound_systems.size();
            output_file << "  Compound: " << it2->second.InChI;
            output_file << endl;
            for (map<const string,System>::iterator it3 = it2->second.compound_systems.begin(); it3 != it2->second.compound_systems.end(); it3++)
            {
                output_file << "    " << it3->second.system_ID << endl;
            }
        }
        for (int i = 0; i < 50; i++)
            output_file << "-";
        output_file << endl;
    }
    output_file << endl;
}

/////////////////////////////////////////////////////////////////////

void output_edges(Network N, string filename)
{
    ofstream output_file;
    output_file.open(filename,std::ios_base::app);

    for (int i = 0; i < 21; i++)
        output_file << "-";
    output_file << " Edges ";
    for (int i = 0; i < 22; i++)
        output_file << "-";
    output_file << endl;

    for (map<const string,Class>::iterator it = N.classes.begin(); it != N.classes.end(); it++)
    {
        Class& C = it->second;
        for (size_t i = 0; i < C.edges.size(); i++)
        {
            for (size_t j = 0; j < C.edges[0].size(); j++)
            {
                if (C.edges[i][j] != -1)
                {
                    // On trouve les deux systèmes
                    System* R;
                    System* P;
                    bool R_found = false;
                    bool P_found = false;
                    map<const string,Compound>::iterator it2 = C.class_compounds.begin();
                    while (R_found == false || P_found == false)
                    {
                        Compound& D = it2->second;
                        for (map<const string,System>::iterator it3 = D.compound_systems.begin(); it3 != D.compound_systems.end(); it3++)
                        {
                            System& S = it3->second;
                            if (S.insertion_rank_in_class == i)
                            {
                                R = &S;
                                R_found = true;
                            }
                            else if (S.insertion_rank_in_class == j)
                            {
                                P = &S;
                                P_found = true;
                            }
                            if (R_found == true && P_found == true)
                                break;
                        }
                        it2++;
                    }
                    output_file << R->system_ID << " | " << P->system_ID << " | " << C.edges[i][j] << endl;
                }
            }
        }
    }
    for (int i = 0; i < 50; i++)
        output_file << "-";
    output_file << endl << endl;
}

/////////////////////////////////////////////////////////////////////

void output_hyperedges(Network N, string filename)
{
    ofstream output_file;
    output_file.open(filename,std::ios_base::app);

    for (int i = 0; i < 19; i++)
        output_file << "-";
    output_file << " Hyperedges ";
    for (int i = 0; i < 19; i++)
        output_file << "-";
    output_file << endl;

    for (set<Hyperedge>::iterator it = N.hyperedges.begin(); it != N.hyperedges.end(); it++)
    {
        output_file << "Reactants: ";
        for (multiset<System>::iterator it2 = it->reactants.begin(); it2 != it->reactants.end(); it2++)
        {
            output_file << calculate_class_ID(*it2).substr(0,calculate_class_ID(*it2).find(" | ")) << " ";
        }
        output_file << "| Products: ";
        for (multiset<System>::iterator it2 = it->products.begin(); it2 != it->products.end(); it2++)
        {
            output_file << calculate_class_ID(*it2).substr(0,calculate_class_ID(*it2).find(" | ")) << " ";
        }
        output_file << "| Barriers: " << it->barrier.first << ", " << it->barrier.second << endl;
    }
}

/////////////////////////////////////////////////////////////////////

void output_network(Network N, string filename)
{
    output_systems(N,filename);
    output_edges(N,filename);
    output_hyperedges(N,filename);
}

/////////////////////////////////////////////////////////////////////

void output_all(Parameters param, Network N, string name_of_created_file)
{
    string filename = "./Output/" + name_of_created_file + ".txt";

    ofstream output_file;
    output_file.open(filename);

    time_t timestamp;
    time(&timestamp);
    output_file << ctime(&timestamp) << endl;

    output_file << N.number_of_systems << " systems discovered (";
    size_t number_of_compounds = 0;
    for (map<string,Class>::iterator it = N.classes.begin(); it != N.classes.end(); it++)
    {
        number_of_compounds += it->second.class_compounds.size();
    }
    output_file << number_of_compounds << " compounds)" << endl << endl;

    output_parameters(param,filename);
    output_network(N,filename);
}
