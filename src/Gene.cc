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
	for(i=0; i<typeseq_length; i++)	typeseq[i] = false;
	for(i=0; i<signalp_length; i++)	signalp[i] = false;
  expression = 0;
	express = 0;
}

Gene::Gene(int k, int typ, int thr, bool tsq[], bool sig[], int exp) : Bead(k)
{
	int i;

	type = typ;
	threshold = thr;
	for(i=0; i<typeseq_length; i++)	typeseq[i] = tsq[i];
	for(i=0; i<signalp_length; i++)	signalp[i] = sig[i];
  expression = exp;
	express = 0;
}


Gene::Gene(const Gene &gene) : Bead(gene)
{
	int i;

  type = gene.type;
  threshold = gene.threshold;
	for(i=0; i<typeseq_length; i++)	typeseq[i] = gene.typeseq[i];
	for(i=0; i<signalp_length; i++) signalp[i] = gene.signalp[i];
  expression = gene.expression;
	express = gene.express;
}

Gene::~Gene()
{
}

Bead* Gene::Clone() const
{
  return new Gene(*this);
}

void Gene::Randomize()
{
	int i;

	cout << "Gene randomize" << endl;

	type = 0;
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
	for (i=0; i<typeseq_length; i++)	typeseq[i] = (uniform()>0.5) ? true : false;
	for (i=0; i<signalp_length; i++)	signalp[i] = (uniform()>0.5) ? true : false;
  expression = 0;
	express = 0;
}

bool Gene::Mutate(double mut_factor)
{
	bool is_mutated = false;

	if ( MutateParameter(&threshold, mu_threshold[kind]*mut_factor) )								is_mutated = true;
	if ( MutateBitstring(signalp, signalp_length, mu_signalp[kind]*mut_factor) )		is_mutated = true;

	return is_mutated;
}

void Gene::DefineTypeFromSeq()
{
	int i;
	bool found_type = false;

	for (i=1; i<6; i++)
	{
		if (BindingAffinity(typeseq, typeseq_defs[i-1], typeseq_length) <= typeseq_length - typeseq_hdist)
		{
			type = i;
			found_type = true;
			break;
		}
	}

	if (type == 0 || ((type > 0 && type < 6) && !found_type))
	{
		type = 6+(int)(uniform()*45);	//Type invention gets random type. Also for divergence from type 1-5, a new random type will be defined. In all other cases (type 1-5 to type 1-5, or type >6 to type >6), you don't have to do anything.
	}
}

string Gene::Show(bool terminal, bool type_only) const
{
	int i;
	string Content, color_prefix, color_suffix;
	std::stringstream ss;

	if (terminal)
	{
		color_prefix = "\033[94m";
		color_suffix = "\033[0m";
	}
	else
	{
		color_prefix = "";
		color_suffix = "";
	}

	ss << "(" << color_prefix << "R" << type << ":";
	if (!type_only) ss << threshold << ":";
	if (!type_only)
	{
		for(i=0; i<typeseq_length; i++)	ss << typeseq[i];
		ss << ":";
		for(i=0; i<signalp_length; i++)	ss << signalp[i];
		ss << ":";
	}
	ss << color_suffix << ")";

	Content = ss.str();
	ss.clear();

	return Content;
}
