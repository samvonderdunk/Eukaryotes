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

bool Bead::MutateType(int* value, double m)
{
	int old_value = (*value);
	bool is_mutated = false;

	if (uniform() < m)
	{
		*value = 1+(int)(uniform()*nr_gene_types);
	}

	if ((*value) != old_value)		is_mutated = true;
	return is_mutated;
}

bool Bead::MutateBitstring(std::bitset<regulator_length>& bitstring, double m)
{
	int i;
	bool is_mutated = false;

	for(i=0; i<regulator_length; i++)
	{
		if(uniform() < m)
		{
			bitstring.flip(i);
			is_mutated = true;
		}
	}

	return is_mutated;
}

bool Bead::MutateBitstring(std::bitset<effector_length>& bitstring, double m)
{
	int i;
	bool is_mutated = false;

	for(i=0; i<effector_length; i++)
	{
		if(uniform() < m)
		{
			bitstring.flip(i);
			is_mutated = true;
		}
	}

	return is_mutated;
}

bool Bead::MutateBitstring(std::bitset<signalp_length>& bitstring, double m)
{
	int i;
	bool is_mutated = false;

	for(i=0; i<signalp_length; i++)
	{
		if(uniform() < m)
		{
			bitstring.flip(i);
			is_mutated = true;
		}
	}

	return is_mutated;
}
