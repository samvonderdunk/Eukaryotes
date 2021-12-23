#ifndef GeneHeader
#define GeneHeader

#include "Bead.hh"
#include "Header.hh"

// Main gene class: genes are beads that can be expressed (threshold, express, expression), their products transported (signalp), and they usually will be given an identity (type, typeseq) even if this may not be superfluous in some cases.

class Gene : public Bead {
	public:
		int type;	//Used to determine phenotypic readout from expression.
		int threshold; //At what regulatory effect does the gene start to become expressed.
		// bool signalp[signalp_length];
		std::bitset<signalp_length> signalp;
		int expression;	//Current expression state.
		int express;		//New expression state.

		Gene(int k);	//Can we do without the default constructor?
		Gene(int k, int typ, int thr, std::bitset<signalp_length>& sig, int exp);
		explicit Gene(const Gene &gene);
		virtual ~Gene();

		virtual Bead* Clone() const=0;
		void Randomize();

		bool Mutate(int organelle);
		virtual string Show(bool terminal, bool type_only=false) const=0;
};

#endif
