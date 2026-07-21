#ifndef hyperedge_h
#define hyperedge_h

#include "system.h"

#include <set>


class Hyperedge
{
    public:
    const std::multiset<System> reactants;
    const std::multiset<System> products;
    const std::pair<float,float> barrier;   // (barrière R vers P, barrière P vers R)
    std::pair<float,float> flux;            // Initialisé à (-1,-1)

    Hyperedge();
    Hyperedge(const std::multiset<System>&, const std::multiset<System>&, const std::pair<float,float>&);

//////////////////////////////
    // Ordonnées selon la valeur de leur plus petite barrière
    friend bool operator==(Hyperedge H1,Hyperedge H2)
    {
        return ((H1.reactants == H2.reactants) && (H1.products == H2.products));
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
