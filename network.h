#ifndef network_h
#define network_h

#include "class.h"
#include "hyperedge.h"

#include <deque>
#include <map>
#include <set>
#include <unordered_set>


class Network
{
    public:
    std::map<const std::string,Class> classes;    // Ordonnées selon leur class_ID
    std::set<Hyperedge> hyperedges;
    size_t number_of_systems;

    std::vector<System> compound_unexplored_systems;    // Systèmes dont les voisins intra-composé n'ont pas été ajoutés
    std::vector<System> class_unexplored_systems;       // Systèmes dont les voisins extra-composé n'ont pas été ajoutés

    Network();
    Network(std::deque<System>&);   // Les systèmes sont copiés du deque, puis ajoutés dans le réseau

    void add_isolated_system_to_network(System&);   // Ajoute un système sans ajouter d'arête ou d'hyperarête
    void add_system_in_network_from_edge(Class&,const System&,System&,const float&,const float&);   // Ajoute une arête entre un système connu et un système potentiellement inconnu
    void add_system_in_network_from_hyperedge(const std::multiset<System>&,std::multiset<System>&,const float&,const float&);   // Ajoute une hyperarête entre un ensemble de systèmes connu et un ensemble de systèmes potentiellement inconnus

//////////////////////////////
    // Fonctions à implémenter (simulation numérique)
    void simulate(std::vector<float>);  // Paramètre: vecteur des concentrations des systèmes initiaux, dans l'ordre dans lequel ils ont été donnés
                                        // Met à jour la fonction de concentration de chaque système et le flux de chaque arête et hyperarête
    void update_systems_concentration_function(std::vector<float>);
    void update_edges_flux();
};

#endif
