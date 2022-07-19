#include "Regulator.hh"

Regulator::Regulator() : Gene(REGULATOR)
{
  activity = 0;
	sequence.reset();
}

Regulator::Regulator(int typ, int thr, int act, std::bitset<regulator_length>& seq, int exp) : Gene(REGULATOR, typ, thr, exp)
{
	activity = act;
	sequence = seq;
}


Regulator::Regulator(const Regulator &reg) : Gene(reg)
{
  activity = reg.activity;
	sequence = reg.sequence;
}

Regulator::~Regulator()
{
}

Bead* Regulator::Clone() const
{
  return new Regulator(*this);
}

void Regulator::Randomize(int organelle)
{
	int i;

	Gene::Randomize(organelle);
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
  for (i=0; i<regulator_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;

}

bool Regulator::Mutate(int organelle)
{
	bool is_mutated = false;

	Gene::Mutate(organelle);	//Mutate signalp and threshold
	if ( MutateParameter(&activity, mu[ACTIVITY][REGULATOR]) )		is_mutated = true;
	if ( MutateBitstring(sequence, mu[SEQUENCE][REGULATOR]) )		is_mutated = true;

	return is_mutated;
}

string Regulator::Show(bool terminal, bool type_only) const
{
	//There is quite some overlap between effector and regulatory Show(), so you'd think I just encode this overlap only in the Gene::Show() version. However, then I would have to slice in some elements here, because the regulatory activity preceeds the typeseq and signalp.
	string Content, color_prefix, color_suffix;
	std::stringstream ss;
	int i;

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
	ss << activity << ":";
	ss << ":";
	for (i=0; i<regulator_length; i++)
	{
		ss << sequence[i];
	}
	ss << color_suffix << ")";

	Content = ss.str();
	ss.clear();

	return Content;
}
