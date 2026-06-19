#ifndef system_h
#define system_h

#include <string>
#include <set>
#include <vector>

class Atom
{
    public:
    std::string symbol;
    int atomic_number;
    float probability;

    Atom();
    Atom(std::string s,int n, float p): symbol(s), atomic_number(n), probability(p){}

    friend bool operator==(Atom a, Atom b)
    {
        return (a.symbol == b.symbol);
    }

    friend bool operator<(Atom a, Atom b)
    {
        return (a.symbol < b.symbol);
    }

    friend bool operator>(Atom a, Atom b)
    {
        return (a.symbol > b.symbol);
    }
};

/////////////////////////////////////////////////////////////////////

class System
{
    public:
    std::string system_ID;      // À l'instantiation, seulement le InChI
                                // Est ensuite modifié par update_system_ID pour devenir InChI - [numéro du conformère]
    bool compound_explored = false; // indique si les voisins intra-composé ont été ajoutés
    bool class_explored = false;    // indique si les voisins extra-composé ont été ajoutés

    const std::multiset<Atom> atoms;
    const size_t n_electrons;
    const int charge;
    std::vector<std::pair<float,float>> concentration;  // fonction de concentration obtenu lors de la dernière simulation.
                                                        // (i+1)ième paire: t_i, c(t_i)

    size_t insertion_rank_in_class;  // rang d'insertion dans la classe, à partir de 0. Détermine l'indice dans la matrice d'adjacence
    size_t insertion_rank_in_network;  // rang d'insertion dans le réseau, à partir de 0.

    System();
    System(const std::string,const std::multiset<Atom>,const size_t,const size_t);    // InChI, atomes, nombre d'électrons, charge

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
