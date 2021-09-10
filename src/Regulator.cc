#include "Regulator.hh"

Regulator::Regulator() : Bead()
{
	int i;

	duplicate = false;
  type = 0;
  threshold = 0;
  activity = 0;
	for(i=0; i<signalp_length; i++)	signalp[i] = false;
  for(i=0; i<sequence_length; i++) sequence[i] = false;
  expression = 0;
	express = 0;
}

Regulator::Regulator(int typ, int thr, int act, bool sig[], bool seq[], int exp) : Bead()
{
	int i;

	duplicate = false;
	type = typ;
	threshold = thr;
	activity = act;
	for(i=0; i<signalp_length; i++)	signalp[i] = sig[i];
	for(i=0; i<sequence_length; i++) sequence[i] = seq[i];
  expression = exp;
	express = 0;
}


Regulator::Regulator(const Regulator &reg) : Bead(reg)
{
	int i;

	duplicate = reg.duplicate;
  type = reg.type;
  threshold = reg.threshold;
  activity = reg.activity;
	for(i=0; i<signalp_length; i++) signalp[i] = reg.signalp[i];
  for(i=0; i<sequence_length; i++) sequence[i] = reg.sequence[i];
  expression = reg.expression;
	express = reg.express;
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
	int i;

	duplicate = false;
	type = (int)(uniform()*nr_types);	//Useful if we do perfect_transport; otherwise will be redefined after InventRegulator anyway.
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
  for (i=0; i<sequence_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;
	for (i=0; i<signalp_length; i++)	signalp[i] = (uniform()>0.5) ? true : false;
  expression = 0;
	express = 0;
}
