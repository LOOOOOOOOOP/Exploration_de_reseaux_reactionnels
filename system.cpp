#include "system.h"

using namespace std;


System::System() :
    n_electrons(0),
    charge(0),
    compound_explored(false),
    class_explored(false),
    insertion_rank_in_class(0),
    insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const string ID,const multiset<Atom> a,const size_t n,const size_t c) :
    system_ID(ID),
    atoms(a),
    n_electrons(n),
    charge(c),
    compound_explored(false),
    class_explored(false),
    insertion_rank_in_class(0),
    insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const System& S) :
    system_ID(S.system_ID),
    atoms(S.atoms),
    n_electrons(S.n_electrons),
    charge(S.charge),
    concentration(S.concentration),
    compound_explored(S.compound_explored),
    class_explored(S.class_explored),
    insertion_rank_in_class(S.insertion_rank_in_class),
    insertion_rank_in_network(S.insertion_rank_in_network){}
