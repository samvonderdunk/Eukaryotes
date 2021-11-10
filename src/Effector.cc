#include "Effector.hh"

Effector::Effector() : Gene(EFFECTOR)
{
}

Effector::Effector(int typ, int thr, bool tsq[], bool sig[], int exp) : Gene(REGULATOR, typ, thr, tsq, sig, exp)
{
}

bool Effector::Mutate(double mut_factor)
{
	bool is_mutated = false;

	Gene::Mutate(mut_factor);	//Mutate signalp and threshold
	if ( MutateBitstring(typeseq, typeseq_length, mu_typeseq[EFFECTOR]*mut_factor) )
	{
		is_mutated = true;
		DefineTypeFromSeq();		//Check if we have to change the type.
	}

	return is_mutated;
}
