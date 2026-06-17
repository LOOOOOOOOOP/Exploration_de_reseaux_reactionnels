#ifndef print_tests_h
#define print_tests_h

#include "../network.h"
#include "../class.h"

void print_system(System);
void print_network(Network); // liste les systèmes découverts avec leur ID et leur fonction de concentration
void print_class_matrix(Class);
void print_hyperedges(Network);

#endif
