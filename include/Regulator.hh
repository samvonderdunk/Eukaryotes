#ifndef RegulatorHeader
#define RegulatorHeader

// #include "Bead.hh"
#include "Gene.hh"
#include "Header.hh"

// Regulatory genes, i.e. transcript factors.

class Regulator : public Gene {
	public:
		int activity;
		bool sequence[regulator_length];

		Regulator();
		Regulator(int typ, int thr, int act, bool sig[], bool seq[], int exp);
		explicit Regulator(const Regulator &reg);
		virtual ~Regulator();

		virtual Bead* Clone() const;
		void Randomize();

		bool Mutate(int organelle);
		string Show(bool terminal, bool type_only=false) const;

};

#endif
