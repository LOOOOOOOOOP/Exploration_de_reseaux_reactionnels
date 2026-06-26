#ifndef print_tests_h
#define print_tests_h

#include "../network.h"
#include "../class.h"

void print_system(System);
void print_network(Network,bool);   // booléen: détailler les petits composés ou pas
void print_class_matrix(Class);
void print_hyperedges(Network);

#endif
