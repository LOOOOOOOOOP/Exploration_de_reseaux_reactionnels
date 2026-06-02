#ifndef calculate_IDs_h
#define calculate_IDs_h

#include "system.h"
#include "compound.h"

const std::string calculate_system_ID(const Compound);  // Calcul le system_ID à partir du Compound
                                                        // (avant l'instantiation du système)

const std::string calculate_InChI(const System);    // Retrouve le InChI pour déterminer le compound auquel appartient le système

const std::string calculate_class_ID(const System);

#endif
