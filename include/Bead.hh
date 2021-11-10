#ifndef BeadHeader
#define BeadHeader

#include "Header.hh"

class Bead {
	public:
		int kind;	//Kind of bead: HOUSE, BSITE, REGULATOR, EFFECTOR.
		bool duplicate;

		Bead();				//Default/empty constructor.
		Bead(int k);	//Constructor with kind.
		explicit Bead(const Bead &b);
		virtual ~Bead();

		virtual Bead* Clone() const=0;

		virtual bool Mutate(double mut_factor)=0;	//No parameters to mutate for a plain bead, but function used by all derived classes.
		virtual void Randomize()=0;
		bool MutateParameter(int* value, double mu);
		bool MutateBitstring(bool* bitstring, int bitstring_length, double mu);
		int BindingAffinity(const bool* sequenceA, const bool* sequenceB, int seqlen = sequence_length) const;

		virtual string Show(bool terminal, bool type_only=false) const=0;

};

inline int Bead::BindingAffinity(const bool* sequenceA, const bool* sequenceB, int seqlen) const
{
	int affinity = 0;
	for (int i=0; i<seqlen; i++)
	{
		affinity += (int)(sequenceA[i]^sequenceB[i]);
	}
	return affinity;
}

#endif
