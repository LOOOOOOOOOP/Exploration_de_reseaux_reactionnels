#ifndef output_results_h
#define output_results_h

#include "parameters.h"
#include "../network.h"
#include "../calculate_IDs.h"

void output_parameters(Parameters,std::string);
void output_systems(Network,std::string);
void output_edges(Network,std::string);
void output_hyperedges(Network,std::string);
void output_network(Network,std::string);
void output_all(Parameters,Network,std::string);

#endif
