#include "Gene.hh"

// Gene::Gene() : Bead()		//Careful: with the default constructor we don't set the Bead kind.
// {
// 	int i;
//
//   type = 0;
//   threshold = 0;
// 	for(i=0; i<typeseq_length; i++)	typeseq[i] = false;
// 	for(i=0; i<signalp_length; i++)	signalp[i] = false;
//   expression = 0;
// 	express = 0;
// }

Gene::Gene(int k) : Bead(k)		//This is the desirable constructor because it passes the kind (REGULATOR or EFFECTOR) on the to Bead constructor.
{
	int i;

  type = 0;
  threshold = 0;
	for(i=0; i<signalp_length; i++)	signalp[i] = false;
  expression = 0;
	express = 0;
}

Gene::Gene(int k, int typ, int thr, bool sig[], int exp) : Bead(k)
{
	int i;

	type = typ;
	threshold = thr;
	for(i=0; i<signalp_length; i++)	signalp[i] = sig[i];
  expression = exp;
	express = 0;
}


Gene::Gene(const Gene &gene) : Bead(gene)
{
	int i;

  type = gene.type;
  threshold = gene.threshold;
	for(i=0; i<signalp_length; i++) signalp[i] = gene.signalp[i];
  expression = gene.expression;
	express = gene.express;
}

Gene::~Gene()
{
}

void Gene::Randomize()
{
	int i;

	type = 0;
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
	for (i=0; i<signalp_length; i++)	signalp[i] = (uniform()>0.5) ? true : false;	//Warning: genes with random localization are produced, even if mu[..][SIGNALP][..] is set to 0.
  expression = 0;
	express = 0;
}

bool Gene::Mutate(int organelle)
{
	bool is_mutated = false;

	if ( MutateParameter(&threshold, mu[organelle][THRESHOLD][kind]) )							is_mutated = true;
	if ( MutateBitstring(signalp, signalp_length, mu[organelle][SIGNALP][kind]) )		is_mutated = true;

	return is_mutated;
}
