# include "hyperedge.h"

using namespace std;


Hyperedge::Hyperedge(){}

Hyperedge::Hyperedge(const multiset<System>& r, const multiset<System>& p, const pair<float, float>& b) :
    reactants(r),
    products(p),
    barrier(b),
    flux(make_pair(-1,-1)){}
