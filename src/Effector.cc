#include "Effector.hh"

Effector::Effector() : Gene(EFFECTOR)
{
	int i;
	for(i=0; i<effector_length; i++) sequence[i] = false;
}

Effector::Effector(int typ, int thr, bool sig[], bool seq[], int exp) : Gene(EFFECTOR, typ, thr, sig, exp)
{
	int i;
	for(i=0; i<effector_length; i++) sequence[i] = seq[i];
}

Effector::Effector(const Effector &eff) : Gene(eff)
{
	int i;

  for(i=0; i<effector_length; i++) sequence[i] = eff.sequence[i];
}

Effector::~Effector()
{
}

Bead* Effector::Clone() const
{
  return new Effector(*this);
}

bool Effector::Mutate(double mut_factor)
{
	bool is_mutated = false;

	Gene::Mutate(mut_factor);	//Mutate signalp and threshold
	if ( MutateBitstring(sequence, effector_length, mu_sequence[EFFECTOR]*mut_factor) )
	{
		is_mutated = true;
		DefineTypeFromSeq();		//Check if we have to change the type.
	}

	return is_mutated;
}

void Effector::DefineTypeFromSeq()
{
	int i;
	bool found_type = false;

	for (i=1; i<6; i++)
	{
		if (BindingAffinity(sequence, effector_types[i-1], effector_length) <= effector_length - eff_hdist)
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

string Effector::Show(bool terminal, bool type_only) const
{
	int i;
	string Content, color_prefix, color_suffix;
	std::stringstream ss;

	if (terminal)
	{
		color_prefix = "\033[93m";
		color_suffix = "\033[0m";
	}
	else
	{
		color_prefix = "";
		color_suffix = "";
	}

	ss << "(" << color_prefix << "E" << type << ":";
	if (!type_only) ss << threshold << ":";
	if (!type_only)
	{
		for(i=0; i<signalp_length; i++)	ss << signalp[i];
		ss << ":";
	}
	for(i=0; i<effector_length; i++)	ss << sequence[i];
	ss << color_suffix << ")";

	Content = ss.str();
	ss.clear();

	return Content;
}
