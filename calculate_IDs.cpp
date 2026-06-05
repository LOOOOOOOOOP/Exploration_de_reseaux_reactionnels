# include "calculate_IDs.h"

#include <string>
#include <set>
using namespace std;


void update_system_ID(System& S,const Compound& C)
{
    size_t compound_number = C.compound_systems.size();
    S.system_ID = S.system_ID + " - " + to_string(compound_number);
}

/////////////////////////////////////////////////////////////////////

const string calculate_InChI(const System S)
{
    // Calcul à partir de System_ID
    string ID = (S.system_ID.substr(0,S.system_ID.find(" - ")));
    return ID;
}

/////////////////////////////////////////////////////////////////////

const string calculate_class_ID(const System S)
{
    string class_ID = "";
    multiset<string> a = S.atoms;
    size_t n = S.n_electrons;

    for (multiset<string>::iterator i = a.begin(); i != a.end(); i++)
    {
        class_ID.append(*i);
    }
    class_ID.append(" - " + to_string(n));
    return class_ID;
}

/////////////////////////////////////////////////////////////////////

bool same_system(const System S, const System T)
{
    return false;   // placeholder
}
