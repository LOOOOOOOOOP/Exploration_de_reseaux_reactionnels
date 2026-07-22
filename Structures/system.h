#ifndef system_h
#define system_h

#include <set>
#include <string>
#include <vector>


class Atom
{
    public:
    const std::string symbol;
    const int atomic_number;
    const float probability;  // Probabilité que l'atome se retrouve dans un système moléculaire généré aléatoirement

    Atom();
    Atom(const std::string s,const int n, const float p): symbol(s), atomic_number(n), probability(p){}

//////////////////////////////
    // Ordonnés selon leur symbole atomique
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
    const std::string system_ID;    // Format "liste_des_atomes/q±x/compound#Y/[Z]"
                                    // où   liste_des_atomes est la concaténation des symboles des atomes selon la notation de Hill
                                    //      ±x est la charge du système (si x=0, il n'y a pas de "/q±x")
                                    //      Y est le numéro du composé (à partir de 1)
                                    //      Z est le numéro du conformère à l'intérieur du composé (à partir de 1)
    const std::multiset<Atom> atoms;
    const size_t n_electrons;
    const int charge;
    std::vector<std::pair<float,float>> concentration;  // Fonction de concentration obtenue lors de la dernière simulation numérique
                                                        // (i+1)ième paire: t_i, c(t_i)
                                                        // Pas de valeur par défaut

    size_t insertion_rank_in_class;     // Rang d'insertion dans la classe, à partir de 0. Détermine l'indice dans la matrice d'adjacence
    size_t insertion_rank_in_network;   // Rang d'insertion dans le réseau, à partir de 0.

    bool compound_explored = false; // Indique si les voisins intra-composé ont été ajoutés
    bool class_explored = false;    // Indique si les voisins extra-composé ont été ajoutés

    System();
    System(const std::string,const std::multiset<Atom>,const size_t,const size_t);    // system_ID, atoms, n_electrons, charge
    System(const System&);

//////////////////////////////

    // Ordonnés selon leur ID
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
