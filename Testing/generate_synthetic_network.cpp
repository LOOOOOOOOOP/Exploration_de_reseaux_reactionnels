#include "generate_synthetic_network.h"

#include <cstdlib>
#include <deque>

using namespace std;


System generate_system()
{

}

/////////////////////////////////////////////////////////////////////

std::deque<System> generate_initial_systems(size_t number_of_initial_systems)
{
    deque<System> initial_systems;
    for (size_t i = 0; i < number_of_initial_systems; i++)
    {
        System S = generate_system();
        initial_systems.push_back(S);
    }
    return initial_systems;
}

/////////////////////////////////////////////////////////////////////

Network generate_synthetic_network(size_t number_of_initial_systems,size_t network_depth)
{
    deque<System> initial_systems = generate_initial_systems(number_of_initial_systems);
    Network N(initial_systems);

    for (size_t i = 0; i < network_depth; i++)
    {
        // Pour les systèmes non explorés
            // Génération des voisins intra-composé

            // Génération des voisins extra-composé

        // Pour un certain nombre de paires de systèmes, génération des voisins extra-classe

        // Pour un certain nombre de systèmes, génération des voisins extra-classe en se splittant en 2 ou en 3

    }
    return N;
}
