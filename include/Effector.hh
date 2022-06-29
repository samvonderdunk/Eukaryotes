#ifndef EffectorHeader
#define EffectorHeader

// #include "Bead.hh"
#include "Gene.hh"
#include "Header.hh"

// Effector genes, with some non-regulatory function.

class Effector : public Gene {
	public:
		std::bitset<effector_length> sequence;

		Effector();	//Still make a default constructor that tells the Gene (and Bead) constructors which kind of bead we are creating.
		Effector(int typ, int thr, std::bitset<effector_length>& seq, int exp);
		explicit Effector(const Effector &eff);
		virtual ~Effector();

		virtual Bead* Clone() const;

		void DefineTypeFromSeq();
		bool Mutate(int organelle);
		string Show(bool terminal, bool type_only=false) const;

};

#endif
