#ifndef network_h
#define network_h

#include "class.h"
#include "hyperedge.h"

#include <map>
#include <unordered_set>
#include <set>
#include <deque>

class Network
{
    public:
    std::map<const std::string,Class> classes;    // ordonnées selon leur class_ID
    std::set<Hyperedge> hyperedges;
    size_t number_of_systems;

    std::vector<System> compound_unexplored_systems;   // Systèmes dont les voisins intra-composé n'ont pas été ajoutés
    std::vector<System> class_unexplored_systems;  // Systèmes dont les voisins extra-composé n'otn pas été ajoutés

    Network();
    Network(std::deque<System>&);   // Les systèmes sont copiés du deque, puis ajoutés dans le réseau.

    void add_isolated_system_to_network(System&); // Ajoute un système sans ajouter d'arête ou d'hyperarête
    void add_system_in_network_from_edge(Class&,System&,System&,float,float);
    void add_system_in_network_from_hyperedge(std::multiset<System>&,std::multiset<System>&,float,float);

    void simulate(std::vector<float>);  // Paramètre: vecteur des concentrations des systèmes initiaux, dans l'ordre dans lequel ils ont été donnés
                                        // Update la fonction de concentration de chaque système et le flux de chaque arête et hyperarête
    void update_systems_concentration_function(std::vector<float>);
    void update_edges_flux();
};

#endif
