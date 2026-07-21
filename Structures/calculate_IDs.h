#ifndef calculate_IDs_h
#define calculate_IDs_h

#include "compound.h"
#include "system.h"


const std::string calculate_compound_ID(const System&);  // Retrouve le compound_ID pour déterminer le compound auquel appartient le système
const std::string calculate_class_ID(const System&);     // Retrouve le class_ID pour déterminer la classe à laquelle appartient le système

#endif
