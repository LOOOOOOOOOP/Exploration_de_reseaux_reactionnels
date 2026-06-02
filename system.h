#ifndef system_h
#define system_h

#include <string>
#include <set>


class System
{
    public:
    const std::string system_ID;    // InChI - [numéro du conformère]
    bool explored;

    const std::multiset<std::string> atoms;
    const size_t n_electrons;
    std::string concentration; // fonction de concentration obtenu lors de la dernière simulation

    size_t insertion_rank_in_class;  // rang d'insertion dans la classe, à partir de 0. Détermine l'indice dans la matrice d'adjacence
    size_t insertion_rank_in_network;  // rang d'insertion dans le réseau, à partir de 0.

    System();
    System(const std::string,const std::multiset<std::string>,const size_t);

/////////////////////////////////////////////////////////////////////

    friend bool operator==(const System& S1, const System& S2)
    {
        return (S1.system_ID == S2.system_ID);
    }

    friend bool operator<(const System& S1, const System& S2)
    {
        return (S1.system_ID < S2.system_ID);
    }

    friend bool operator>(const System& S1, const System& S2)
    {
        return (S1.system_ID > S2.system_ID);
    }
};

#endif
