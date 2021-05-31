#include "Regulator.hh"
#include "Header.hh"

Regulator::Regulator() : Bead()
{
	int i;

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

	type = typ;
	threshold = thr;
	acitivity = act;
	for(i=0; i<signalp_length; i++)	signalp[i] = sig[i];
	for(i=0; i<sequence_length; i++) sequence[i] = seq[i];
  expression = exp;
	express = 0;
}


Regulator::Regulator(const Regulator &reg) : Bead(reg)
{
	int i;

  type=reg.type;
  threshold=reg.threshold;
  activity=reg.activity;
  for(i=0; i<sequence_length; i++) sequence[i] = reg.sequence[i];
	for(i=0; i<signalp_length; i++) signalp[i] = reg.signalp[i];
  expression=reg.expression;
	express=reg.express;
}

Bead* Regulator::Clone() const
{
  return new Regulator(*this);
}

void Regulator::RandomRegulator()
{
	int i;

  type = 0;	//This will be assigned somewhere else (in the genome)???
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
  for (i=0; i<sequence_length; i++)	sequence[k] = (uniform()>0.5) ? true : false;
	for (i=0; i<signalp_length; i++)	signalp[k] = (uniform()>0.5) ? true : false;
  expression = 0;
	express = 0;
}

Regulator::~Regulator()
{
}
