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
		int expression;	//Current expression state.
		int express;		//New expression state.

		Regulator();
		explicit Regulator(const Regulator &reg);
		virtual Bead* Regulator() const;
		void RandomRegulator();
		~Regulator();

};

#endif
