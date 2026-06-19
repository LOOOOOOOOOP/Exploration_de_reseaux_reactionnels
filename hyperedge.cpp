# include "hyperedge.h"
using namespace std;

Hyperedge::Hyperedge(){}

Hyperedge::Hyperedge(multiset<System>& r, multiset<System>& p, pair<float, float> b) :
    reactants(r),
    products(p),
    barrier(b),
    flux(make_pair(-1,-1)){}
