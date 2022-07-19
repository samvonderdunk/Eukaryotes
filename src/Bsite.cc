#include "Bsite.hh"

Bsite::Bsite() : Bead(BSITE)
{
  activity = 0;
	sequence.reset();
}

Bsite::Bsite(int act, std::bitset<regulator_length>& seq) : Bead(BSITE)
{
  activity = act;
	sequence = seq;
}

Bsite::Bsite(const Bsite &bsite) : Bead(bsite)
{
  activity = bsite.activity;
	sequence = bsite.sequence;
}

Bsite::~Bsite()
{
}

Bead* Bsite::Clone() const
{
  return new Bsite(*this);
}

void Bsite::Randomize(int organelle)
{
	int i;

  activity = (uniform()>0.5) ? -1 : 1;
  for (i=0; i<regulator_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;	//Is set
}

bool Bsite::Mutate(int organelle)
{
	bool is_mutated = false;

	if ( MutateParameter(&activity, mu[ACTIVITY][BSITE]) )		is_mutated = true;
	if ( MutateBitstring(sequence, mu[SEQUENCE][BSITE]) )		is_mutated = true;

	return is_mutated;
}

string Bsite::Show(bool terminal, bool type_only) const
{
	string Content, color_prefix, color_suffix;
	std::stringstream ss;
	int i;

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
	for (i=0; i<regulator_length; i++)
	{
		ss << sequence[i];
	}
	ss << color_suffix << ")";

	Content = ss.str();
	ss.clear();

	return Content;
}
