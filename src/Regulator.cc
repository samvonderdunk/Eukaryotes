#include "Regulator.hh"
#include "Header.hh"

Regulator::Regulator() : Bead()
{
  type = 0;
  threshold = 0;
  activity = 0;
  for(int i=0; i<sequence_length; i++) sequence[i] = 0;
  expression = 0;
	express = 0;
}

Regulator::Regulator(const Regulator &reg) : Bead(reg)
{
  type=reg.type;
  threshold=reg.threshold;
  activity=reg.activity;
  for(int i=0; i<sequence_length; i++) sequence[i] = reg.sequence[i];
  expression=reg.expression;
	express=reg.xpr;
}

Bead* Regulator::Clone() const
{
  return new Regulator(*this);
}

void Regulator::RandomRegulator()
{
  type = 0;	//This will be assigned somewhere else (in the genome)???
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
  for (int i=0; i<sequence_length; i++)	sequence[k] = (uniform()>0.5) ? true : false;
  expression = 0;
	express = 0;
}

Regulator::~Regulator()
{
}
