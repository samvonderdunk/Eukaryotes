#include "Regulator.hh"

Regulator::Regulator() : Gene(REGULATOR)
{
	int i;

  activity = 0;
  for(i=0; i<regulator_length; i++) sequence[i] = false;

}

Regulator::Regulator(int typ, int thr, int act, bool sig[], bool seq[], int exp) : Gene(REGULATOR, typ, thr, sig, exp)
{
	int i;

	activity = act;
	for(i=0; i<regulator_length; i++) sequence[i] = seq[i];

}


Regulator::Regulator(const Regulator &reg) : Gene(reg)
{
	int i;

  activity = reg.activity;
  for(i=0; i<regulator_length; i++) sequence[i] = reg.sequence[i];

}

Regulator::~Regulator()
{
}

Bead* Regulator::Clone() const
{
  return new Regulator(*this);
}

void Regulator::Randomize()
{
	int i;

	cout << "Regulator randomize" << endl;

	Gene::Randomize();
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
  for (i=0; i<regulator_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;

}

bool Regulator::Mutate(double mut_factor)
{
	bool is_mutated = false;

	Gene::Mutate(mut_factor);	//Mutate signalp and threshold
	if ( MutateParameter(&activity, mu_activity[REGULATOR]*mut_factor) )										is_mutated = true;
	if ( MutateBitstring(sequence, regulator_length, mu_sequence[REGULATOR]*mut_factor) )		is_mutated = true;

	return is_mutated;
}

string Regulator::Show(bool terminal, bool type_only) const
{
	//There is quite some overlap between effector and regulatory Show(), so you'd think I just encode this overlap only in the Gene::Show() version. However, then I would have to slice in some elements here, because the regulatory activity preceeds the typeseq and signalp.
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
	ss << activity << ":";
	if (!type_only)
	{
		for(i=0; i<signalp_length; i++)	ss << signalp[i];
		ss << ":";
	}
	for(i=0; i<regulator_length; i++)	ss << sequence[i];
	ss << color_suffix << ")";

	Content = ss.str();
	ss.clear();

	return Content;
}
