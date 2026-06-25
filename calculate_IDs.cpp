# include "calculate_IDs.h"

#include <string>
#include <set>
using namespace std;


void update_system_ID(System& S,const Compound& C)
{
    size_t compound_number = C.compound_systems.size();
    S.system_ID = S.system_ID + " [" + to_string(compound_number) + "]";
}

/////////////////////////////////////////////////////////////////////

const string calculate_InChI(const System S)
{
    // Calcul à partir de System_ID
    int pos = S.system_ID.find(" [");
    if (pos == string::npos)
        return S.system_ID;
    else
        return S.system_ID.substr(0,pos);
}

/////////////////////////////////////////////////////////////////////

const string calculate_class_ID(const System S)
{
    string class_ID = "";
    multiset<Atom> a = S.atoms;
    size_t n = S.n_electrons;

    for (multiset<Atom>::iterator i = a.begin(); i != a.end(); i++)
    {
        class_ID.append(i->symbol);
    }

    if (S.charge > 0)
    {
        if (S.charge == 1)
            class_ID.append("+");
        else
            class_ID.append(" " + to_string(S.charge) + "+");
    }
    else if (S.charge < 0)
    {
        if (S.charge == -1)
            class_ID.append("-");
        else
            class_ID.append(" " + to_string(-S.charge) + "-");
    }

    class_ID.append(" | " + to_string(n) + " electrons");
    return class_ID;
}

/////////////////////////////////////////////////////////////////////

bool same_system(const System S, const System T)
{
    return false;   // placeholder
}
