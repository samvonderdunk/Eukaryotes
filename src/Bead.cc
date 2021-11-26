#include "Bead.hh"

Bead::Bead()
{
	kind=-1;
	duplicate=false;
}

Bead::Bead(int k)
{
	kind = k;
	duplicate = false;
}

Bead::~Bead()
{
}

Bead::Bead(const Bead &b)
{
	kind = b.kind;
	duplicate = b.duplicate;
}

bool Bead::MutateParameter(int* value, double m)
{
	int old_value = (*value);
	bool is_mutated = false;

	if (uniform() < m)
	{
		if (uniform() < 0.8)	*value = (uniform()>0.5) ? (*value)+1 : (*value)-1;
		else									*value = (int)(uniform()*(2*WeightRange+1) - WeightRange);
	}

	if ((*value) != old_value)		is_mutated = true;
	return is_mutated;
}

bool Bead::MutateBitstring(bool* bitstring, int bitstring_length, double m)
{
	int i;
	bool is_mutated = false;

	for(i=0; i<bitstring_length; i++)
	{
		if(uniform() < m)
		{
			if (bitstring[i] == false)			bitstring[i] = true;
			else if (bitstring[i] == true)	bitstring[i] = false;
			is_mutated = true;
		}
	}

	return is_mutated;
}
