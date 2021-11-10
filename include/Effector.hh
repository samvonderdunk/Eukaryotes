#ifndef EffectorHeader
#define EffectorHeader

// #include "Bead.hh"
#include "Gene.hh"
#include "Header.hh"

// Effector genes, with some non-regulatory function.

class Effector : public Gene {
	public:

		using Gene::Gene;	//Using the constructors and destructors of Gene class. All other functions of Gene class can be used by any Effector object anyway.
		Effector();	//Still make a default constructor that tells the Gene (and Bead) constructors which kind of bead we are creating.
		Effector(int typ, int thr, bool tsq[], bool sig[], int exp);


		bool Mutate(double mut_factor);

};

#endif
