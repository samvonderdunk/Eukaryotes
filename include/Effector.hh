#ifndef EffectorHeader
#define EffectorHeader

// #include "Bead.hh"
#include "Gene.hh"
#include "Header.hh"

// Effector genes, with some non-regulatory function.

class Effector : public Gene {
	public:

		bool sequence[effector_length];

		Effector();	//Still make a default constructor that tells the Gene (and Bead) constructors which kind of bead we are creating.
		Effector(int typ, int thr, bool sig[], bool seq[], int exp);
		explicit Effector(const Effector &eff);
		virtual ~Effector();

		virtual Bead* Clone() const;

		void DefineTypeFromSeq();
		bool Mutate(double mut_factor);
		string Show(bool terminal, bool type_only=false) const;

};

#endif
