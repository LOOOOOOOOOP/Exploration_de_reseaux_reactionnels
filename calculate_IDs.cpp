# include "calculate_IDs.h"

#include <set>
#include <string>

using namespace std;


const string calculate_compound_ID(const System& S)
{
    int pos = S.system_ID.find("/[");
    if (pos == string::npos)
        return S.system_ID;
    else
        return S.system_ID.substr(0,pos);
}

/////////////////////////////////////////////////////////////////////

const string calculate_class_ID(const System& S)
{
    string class_ID = S.system_ID;
    class_ID.erase(class_ID.find("/"),string::npos);

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

    return class_ID;
}
