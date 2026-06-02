#ifndef class_h
#define class_h

#include "system.h"
#include "compound.h"

#include <string>
#include <set>
#include <vector>

class Class
{
    private:
    const std::multiset<std::string> atoms;
    const size_t n_electrons;
    public:
    size_t number_of_systems_in_class;
    std::map<const std::string,Compound> class_compounds; // ordonnés selon leur InChI
    std::vector<std::vector<float>> edges; // matrice d'adjacence avec les barrières. Si les deux systèmes ne sont pas adjacents, la barrière est de -1
                                            // premier vecteur: contient les vecteurs-lignes

    const std::string class_ID;

    Class();
    Class(size_t&,const std::string,System&);

    void add_system_in_class(size_t&,System&,System&,float,float); // Ajoute un système dans la même classe qu'un autre système, en créant l'arête. Update aussi le nombre de systèmes dans le réseau
    void add_isolated_system_in_class(size_t&,System&); // Ajoute un système qui n'est pas adjacent à un autre système dans la classe

    /////////////////////////////////////////////////////////////////////

    friend bool operator==(const Class& C1, const Class& C2)
    {
        return (C1.class_ID == C2.class_ID);
    }

    friend bool operator<(const Class& C1, const Class& C2)
    {
        return (C1.class_ID < C2.class_ID);
    }

    friend bool operator>(const Class& C1, const Class& C2)
    {
        return (C1.class_ID > C2.class_ID);
    }

};

#endif
