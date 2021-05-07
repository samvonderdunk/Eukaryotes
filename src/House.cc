#include "House.hh"

House::House() : Bead()
{
}

House::~House() {
}

House::House(const House &house) : Bead(house)
{
}

Bead* House::Clone() const
{
  return new House(*this);
}
