#ifndef compound_h
#define compound_h

#include "system.h"

#include <map>


class Compound
{
    public:
    const std::string compound_ID;                          // Format "liste_des_atomes/q±x/compound#Y"
    std::map<const std::string,System> compound_systems;    // Ordonnés selon leur system_ID

    Compound();
    Compound(const std::string);    // Crée le composé (vide) à partir du compound_ID

//////////////////////////////
    // Ordonnés selon leur ID
    friend bool operator==(const Compound& C1, const Compound& C2)
    {
        return (C1.compound_ID == C2.compound_ID);
    }

    friend bool operator<(const Compound& C1, const Compound& C2)
    {
        return (C1.compound_ID < C2.compound_ID);
    }

    friend bool operator>(const Compound& C1, const Compound& C2)
    {
        return (C1.compound_ID > C2.compound_ID);
    }
};

#endif
