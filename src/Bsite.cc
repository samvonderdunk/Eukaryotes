#include "Bsite.hh"

Bsite::Bsite() : Bead()
{
  activity=0;
  for(int i=0; i<sequence_length; i++) sequence[i] = 0;
}

Bsite::Bsite(int act, bool seq[]) : Bead()
{
  activity=act;
  for(int i=0; i<sequence_length; i++) sequence[i] = seq[i];
}

Bsite::Bsite(const Bsite &bsite) : Bead(bsite)
{
  activity=bsite.activity;
  for(int i=0; i<sequence_length; i++) sequence[i] = bsite.sequence[i];
}

Bsite::~Bsite()
{
}

Bead* Bsite::Clone() const
{
  return new Bsite(*this);
}

void Bsite::RandomBsite()
{
  activity = (uniform()>0.5) ? -1 : 1;
  for (int i=0; i<sequence_length; k++)	sequence[k] = (uniform()>0.5) ? true : false;
}
