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

    Network();
    Network(std::deque<System>&);   // Les systèmes sont copiés du deque, puis ajoutés dans le réseau.

    void simulate(std::vector<float>);  // Paramètre: vecteur des concentrations des systèmes initiaux, dans l'ordre dans lequel ils ont été donnés
                                        // Update la fonction de concentration de chaque système
    bool take_decision();   // Retourne faux si la condition d'arrêt a été atteinte
    void generate_neighbours(std::set<System>);
    void add_isolated_system_to_network(System&); // Ajoute un système sans ajouter d'arête ou d'hyperarête
    void add_system_in_network_from_edge(Class&,System&,System&,float,float);
    void add_system_in_network_from_hyperedge(std::multiset<System>&,std::multiset<System>&,float,float);
};

#endif
