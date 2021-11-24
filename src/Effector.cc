#include "Effector.hh"

Effector::Effector() : Gene(EFFECTOR)
{
	int i;
	for(i=0; i<effector_length; i++) sequence[i] = false;
}

Effector::Effector(int typ, int thr, bool sig[], bool seq[], int exp) : Gene(REGULATOR, typ, thr, sig, exp)
{
	int i;
	for(i=0; i<effector_length; i++) sequence[i] = seq[i];
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
