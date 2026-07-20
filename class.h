#ifndef class_h
#define class_h

#include "compound.h"
#include "system.h"

#include <set>
#include <string>
#include <vector>


class Class
{
    public:
    const std::multiset<Atom> atoms;
    const size_t n_electrons;

    size_t number_of_systems_in_class;
    std::map<const std::string,Compound> class_compounds;   // Ordonnés selon leur compound_ID
    std::vector<std::vector<float>> edges;  // Matrice d'adjacence avec les barrières. Si les deux systèmes ne sont pas adjacents, la barrière est de -1
                                            // Premier vecteur: contient les vecteurs-lignes
    std::vector<std::vector<float>> edges_flux;     // Flux des arêtes de la matrice d'adjacence. Initialisé à -1

    const std::string class_ID; // Format "liste_des_atomes x±"
                                // où   "liste_des_atomes" est la concaténation des symboles des atomes tels qu'il sont ordonnés dans le multiset
                                //      " x±" est   "" si la charge est 0
                                //                  "+" si la charge est +1
                                //                  "-" si la charge est -1
                                //                  "|x|+" si la charge est x>0
                                //                  "|x|-" si la charge est x<0

    Class();
    Class(size_t&,std::vector<System>&,std::vector<System>&,const std::string&,System&);   // number_of_systems (in network), compound_unexplored_systems, class_unexplored_systems, class_ID, système
                                                                                                // La classe est créée à partir du système

    void add_system_in_class_from_edge(size_t&,std::vector<System>&,std::vector<System>&,const System&,System&,const float&,const float&);  // number_of_systems (in network), compound_unexplored_systems, class_unexplored_systems, système à ajouter, barrière R vers P, barrière P vers R
                                                                                                                                    // Ajoute un système dans la même classe qu'un autre système, en créant l'arête
    void add_isolated_system_in_class(size_t&,std::vector<System>&,std::vector<System>&,System&);   // number_of_systems (n network), compound_unexplored_systems, class_unexplored_systems, système à ajouter
                                                                                                    // Ajoute un système qui n'est pas adjacent à un autre système dans la classe, sans arête

//////////////////////////////
    // Ordonnées selon leur ID
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
