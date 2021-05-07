#include "Regulator.hh"
#include "Header.hh"

Regulator::Regulator() : Bead()
{
  type = 0;
  threshold = 0;
  activity = 0;
  for(int i=0; i<sequence_length; i++) sequence[i] = 0;
  expression = 0;
}

Regulator::Regulator(int t, int th, int act, bool bd_dom[], int expr) : Bead()
{
	type=t;
	threshold=th;
	activity=act;
	for(int i=0; i<binding_length; i++) binding_domain[i] = bd_dom[i];
	expression=expr;
}

Regulator::Regulator(const Regulator &reg) : Bead(reg)
{
  type=reg.type;
  threshold=reg.threshold;
  activity=reg.activity;
  for(int i=0; i<sequence_length; i++) sequence[i] = reg.sequence[i];
  expression=reg.expression;
}

Regulator::~Regulator()
{
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
}
