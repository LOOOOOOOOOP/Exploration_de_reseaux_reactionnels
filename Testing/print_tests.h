#ifndef print_tests_h
#define print_tests_h

#include "../Structures/class.h"
#include "../Structures/network.h"


void print_system(const System&);
void print_network(const Network&,const bool&);   // booléen: vrai pour détailler le contenu des petits composés
void print_class_matrix(const Class&);
void print_hyperedges(const Network&);

#endif
