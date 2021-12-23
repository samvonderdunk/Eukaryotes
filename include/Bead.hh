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

		virtual bool Mutate(int organelle)=0;	//No parameters to mutate for a plain bead, but function used by all derived classes.
		virtual void Randomize()=0;
		bool MutateParameter(int* value, double m);

		//Try overloading the BindingAffinity and MutateBitstring functions with different length bitsets.
		bool MutateBitstring(std::bitset<regulator_length>& bitstring, double m);
		bool MutateBitstring(std::bitset<effector_length>& bitstring, double m);
		bool MutateBitstring(std::bitset<signalp_length>& bitstring, double m);

		int BindingAffinity(std::bitset<regulator_length> sequenceA, const std::bitset<regulator_length>& sequenceB) const;
		int BindingAffinity(std::bitset<effector_length> sequenceA, const std::bitset<effector_length>& sequenceB) const;
		int BindingAffinity(std::bitset<signalp_length> sequenceA, const std::bitset<signalp_length>& sequenceB) const;

		virtual string Show(bool terminal, bool type_only=false) const=0;

};

inline int Bead::BindingAffinity(std::bitset<regulator_length> sequenceA, const std::bitset<regulator_length>& sequenceB) const
{
	return (int)(sequenceA ^= sequenceB).count();
}

inline int Bead::BindingAffinity(std::bitset<effector_length> sequenceA, const std::bitset<effector_length>& sequenceB) const
{
	return (int)(sequenceA ^= sequenceB).count();
}

inline int Bead::BindingAffinity(std::bitset<signalp_length> sequenceA, const std::bitset<signalp_length>& sequenceB) const
{
	return (int)(sequenceA ^= sequenceB).count();
}

#endif
