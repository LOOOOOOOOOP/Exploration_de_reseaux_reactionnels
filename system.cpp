#include "system.h"
using namespace std;

System::System() : explored(false),n_electrons(0),insertion_rank_in_class(0),insertion_rank_in_network(0){}

/////////////////////////////////////////////////////////////////////

System::System(const string id,const multiset<string> a,const size_t n) :
    system_ID(id), explored(false), atoms(a), n_electrons(n), insertion_rank_in_class(0), insertion_rank_in_network(0){}
