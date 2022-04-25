#ifndef RegulatorHeader
#define RegulatorHeader

// #include "Bead.hh"
#include "Gene.hh"
#include "Header.hh"

// Regulatory genes, i.e. transcript factors.

class Regulator : public Gene {
	public:
		int activity;
		std::bitset<regulator_length> sequence;

		Regulator();
		Regulator(int typ, int thr, int act, std::bitset<signalp_length>& sig, std::bitset<regulator_length>& seq, int exp);
		explicit Regulator(const Regulator &reg);
		virtual ~Regulator();

		virtual Bead* Clone() const;

		void Randomize(int organelle);
		bool Mutate(int organelle);
		string Show(bool terminal, bool type_only=false) const;

};

#endif
