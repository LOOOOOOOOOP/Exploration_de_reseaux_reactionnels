#ifndef calculate_IDs_h
#define calculate_IDs_h

#include "system.h"
#include "compound.h"

void update_system_ID(System&,const Compound&); // Rajoute le numéro du conformère au system_ID
                                                // (à l'instantiation du système, son ID est seulement son InChI)

const std::string calculate_InChI(const System);    // Retrouve le InChI pour déterminer le compound auquel appartient le système

const std::string calculate_class_ID(const System);

bool same_system(const System, const System);   // Indique si les deux systèmes sont indentiques (même conformère)

#endif
