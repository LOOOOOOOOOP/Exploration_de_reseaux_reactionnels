#ifndef system_h
#define system_h

#include <string>
#include <set>
#include <vector>


class System
{
    public:
    std::string system_ID;      // À l'instantiation, system_ID = InChI.
                                // Lorsque le système est ajouté au réseau (et donc à un composé),
                                // le ID est modifié pour devenir InChI - [numéro du conformère]
    bool compound_explored = false; // indique si les voisins intra-composé ont été ajoutés
    bool class_explored = false;    // indique si les voisins extra-composé ont été ajoutés

    const std::multiset<std::string> atoms;
    const size_t n_electrons;
    std::vector<std::pair<float,float>> concentration;  // fonction de concentration obtenu lors de la dernière simulation.
                                                        // (i+1)ième paire: t_i, c(t_i)

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
