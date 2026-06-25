#include "print_tests.h"

#include <iostream>
using namespace std;

void print_system(System S)
{
    cout << S.system_ID << endl;
}

void print_network(Network N)
{
    cout << "Number of systems discovered: " << N.number_of_systems << endl << endl;

    for (map<const string,Class>::iterator it = N.classes.begin(); it != N.classes.end(); it++)
    {
        cout << "////////////////////////////////////////////////////////////////////////////////" << endl;
        cout << "Class: " << it->second.class_ID << endl;
        for (map<const string,Compound>::iterator it2 = it->second.class_compounds.begin(); it2 != it->second.class_compounds.end(); it2++)
        {
            cout << "  Compound: " << it2->second.InChI << endl;
            if (it2->second.compound_systems.size() > 5) // if compound is big
                cout << "    (" << it2->second.compound_systems.size() << " conformers)" << endl;
            else // list conformers
            {
                for (map<const string,System>::iterator it3 = it2->second.compound_systems.begin(); it3 != it2->second.compound_systems.end(); it3++)
                {
                    cout << "    " << it3->second.system_ID << " | rank in class: " << it3->second.insertion_rank_in_class << " | rank in network: " << it3->second.insertion_rank_in_network << endl;
                }
            }
        }
    }
    cout << endl << "////////////////////////////////////////////////////////////////////////////////" << endl;
}

/////////////////////////////////////////////////////////////////////

void print_class_matrix(Class C)
{
    cout << endl << "Adjacency matrix of class " << C.class_ID << ":" << endl;
    for (size_t i = 0; i < C.edges.size(); i++)
    {
        for (size_t j = 0; j < C.edges.size(); j++)
        {
            cout << " " << C.edges[i][j] << "  ";
        }
        cout << endl;
    }
}

/////////////////////////////////////////////////////////////////////

void print_hyperedges(Network N)
{
    cout << endl << "Hyperedges:" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    for (set<Hyperedge>::iterator it = N.hyperedges.begin(); it != N.hyperedges.end(); it++)
    {
        cout << "Reactants: ";
        for (multiset<System>::iterator it2 = it->reactants.begin(); it2 != it->reactants.end(); it2++)
        {
            cout << it2->system_ID << " ";
        }
        cout << "| Products: ";
        for (multiset<System>::iterator it2 = it->products.begin(); it2 != it->products.end(); it2++)
        {
            cout << it2->system_ID << " ";
        }
        cout << "| Barriers: " << it->barrier.first << ", " << it->barrier.second;

        cout << endl << "----------------------------------------------------------------------" << endl;
    }


}
