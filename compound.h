#ifndef compound_h
#define compound_h

#include "system.h"

#include <map>

class Compound
{
    public:
    const std::string InChI;
    std::map<const std::string,System> compound_systems; // ordonnés selon leur system_ID

    Compound();
    Compound(const std::string);    // crée le compound (vide) à partir du InChI

/////////////////////////////////////////////////////////////////////

    friend bool operator==(const Compound& C1, const Compound& C2)
    {
        return (C1.InChI == C2.InChI);
    }

    friend bool operator<(const Compound& C1, const Compound& C2)
    {
        return (C1.InChI < C2.InChI);
    }

    friend bool operator>(const Compound& C1, const Compound& C2)
    {
        return (C1.InChI > C2.InChI);
    }
};

#endif
