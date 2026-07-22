#ifndef output_results_h
#define output_results_h

#include "parameters.h"
#include "../Structures/calculate_IDs.h"
#include "../Structures/network.h"


void output_parameters(const Parameters&,const std::string&);
void output_systems(Network&,const std::string&);
void output_edges(Network&,const std::string&);
void output_hyperedges(Network&,const std::string&);
void output_network(Network&,const std::string&);
void output_all(const Parameters&,Network&,const std::string&);    // crée un fichier dans ./Output

#endif
