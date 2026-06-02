#ifndef generate_synthetic_network_h
#define generate_synthetic_network_h

#include "../network.h"
#include "../system.h"
#include <deque>

System generate_system();
std::deque<System> generate_initial_systems(size_t);
Network generate_synthetic_network(size_t,size_t);  // Nombre de systèmes initiaux et profondeur du réseau synthétique


#endif
