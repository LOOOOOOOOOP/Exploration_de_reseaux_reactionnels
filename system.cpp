#include "system.h"
using namespace std;

System::System() : compound_explored(false),class_explored(false),n_electrons(0),charge(0),insertion_rank_in_class(0),insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const string InChI,const multiset<Atom> a,const size_t n,const size_t c) :
    system_ID(InChI),compound_explored(false), class_explored(false), atoms(a), n_electrons(n), charge(c), insertion_rank_in_class(0), insertion_rank_in_network(0){}
