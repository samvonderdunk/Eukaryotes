#ifndef RegulatorHeader
#define RegulatorHeader

#include "Bead.hh"

// Regulatory genes, i.e. transcription factors.

class Regulator : public Bead {
	public:
		int type;	//Used to determine phenotypic readout from expression.
		int threshold; //At what regulatory effect does the gene start to become expressed.
		int activity;
		bool sequence[sequence_length];
		int expression;

		Regulator();
		Regulator(int typ, int thr, int act, bool seq[], int exp);
		explicit Regulator(const Regulator &reg);
		~Regulator();

		virtual Bead* Regulator() const;
		void RandomRegulator();
};

#endif
