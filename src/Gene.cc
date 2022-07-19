#include "Gene.hh"

Gene::Gene(int k) : Bead(k)		//This is the desirable constructor because it passes the kind (REGULATOR or EFFECTOR) on the to Bead constructor.
{
  type = 0;
  threshold = 0;
  expression = 0;
	express = 0;
}

Gene::Gene(int k, int typ, int thr, int exp) : Bead(k)
{
	type = typ;
	threshold = thr;
  expression = exp;
	express = 0;
}


Gene::Gene(const Gene &gene) : Bead(gene)
{
  type = gene.type;
  threshold = gene.threshold;
  expression = gene.expression;
	express = gene.express;
}

Gene::~Gene()
{
}

void Gene::Randomize(int organelle)
{
	type = 0;
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
  expression = 0;
	express = 0;
}

bool Gene::Mutate(int organelle)
{
	bool is_mutated = false;

	if ( MutateParameter(&threshold, mu[THRESHOLD][kind]) )			is_mutated = true;
	if ( MutateType(&type, mu[TYPE][kind]) )											is_mutated = true;

	return is_mutated;
}
