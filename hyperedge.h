#ifndef hyperedge_h
#define hyperedge_h

#include "system.h"

#include <set>

class Hyperedge
{
    public:
    const std::multiset<System> reactants;
    const std::multiset<System> products;
    const std::pair<float,float> barrier;

    Hyperedge();
    Hyperedge(std::multiset<System>&, std::multiset<System>&, std::pair<float,float>);

/////////////////////////////////////////////////////////////////////

    friend bool operator==(Hyperedge H1,Hyperedge H2)
    {
        return (std::min(H1.barrier.first,H1.barrier.second) == std::min(H2.barrier.first,H2.barrier.second));
    }

    friend bool operator<(Hyperedge H1,Hyperedge H2)
    {
        return (std::min(H1.barrier.first,H1.barrier.second) < std::min(H2.barrier.first,H2.barrier.second));
    }

    friend bool operator>(Hyperedge H1,Hyperedge H2)
    {
        return (std::min(H1.barrier.first,H1.barrier.second) > std::min(H2.barrier.first,H2.barrier.second));
    }
};

#endif
