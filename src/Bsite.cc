#include "Bsite.hh"

Bsite::Bsite() : Bead(BSITE)
{
	int i;

  activity = 0;
  for(i=0; i<regulator_length; i++) sequence[i] = 0;
}

Bsite::Bsite(int act, bool seq[]) : Bead(BSITE)
{
	int i;

  activity = act;
  for(i=0; i<regulator_length; i++) sequence[i] = seq[i];
}

Bsite::Bsite(const Bsite &bsite) : Bead(bsite)
{
	int i;

  activity = bsite.activity;
  for(i=0; i<regulator_length; i++) sequence[i] = bsite.sequence[i];
}

Bsite::~Bsite()
{
}

Bead* Bsite::Clone() const
{
  return new Bsite(*this);
}

void Bsite::Randomize()
{
	int i;

  activity = (uniform()>0.5) ? -1 : 1;
  for (i=0; i<regulator_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;
}

bool Bsite::Mutate(double mut_factor)
{
	bool is_mutated = false;

	if ( MutateParameter(&activity, mu_activity[BSITE]*mut_factor) )										is_mutated = true;
	if ( MutateBitstring(sequence, regulator_length, mu_sequence[BSITE]*mut_factor) )		is_mutated = true;

	return is_mutated;
}

string Bsite::Show(bool terminal, bool type_only) const
{
	int i;
	string Content, color_prefix, color_suffix;
	std::stringstream ss;

	if(terminal){
		color_prefix = "\033[92m";
		color_suffix = "\033[0m";
	}
	else
	{
		color_prefix = "";
		color_suffix = "";
	}

	ss << "(" << color_prefix << activity << ":";
	for(i=0; i<regulator_length; i++)	ss << sequence[i];
	ss << color_suffix << ")";
	Content = ss.str();
	ss.clear();

	return Content;
}
