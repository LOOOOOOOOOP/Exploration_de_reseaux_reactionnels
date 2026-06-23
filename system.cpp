#include "system.h"
using namespace std;

System::System() : compound_explored(false),class_explored(false),n_electrons(0),charge(0),insertion_rank_in_class(0),insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const string InChI,const multiset<Atom> a,const size_t n,const size_t c) :
    system_ID(InChI),compound_explored(false), class_explored(false), atoms(a), n_electrons(n), charge(c), insertion_rank_in_class(0), insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const System& S) :
    system_ID(S.system_ID),
    compound_explored(S.compound_explored),
    class_explored(S.class_explored),
    atoms(S.atoms),
    n_electrons(S.n_electrons),
    charge(S.charge),
    concentration(S.concentration),
    insertion_rank_in_class(S.insertion_rank_in_class),
    insertion_rank_in_network(S.insertion_rank_in_network){}
