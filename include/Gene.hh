#ifndef GeneHeader
#define GeneHeader

#include "Bead.hh"
#include "Header.hh"

// Main gene class: genes are beads that can be expressed (threshold, express, expression), their products transported (signalp), and they usually will be given an identity (type, typeseq) even if this may not be superfluous in some cases.

class Gene : public Bead {
	public:
		int type;	//Used to determine phenotypic readout from expression.
		int threshold; //At what regulatory effect does the gene start to become expressed.
		bool typeseq[typeseq_length];
		bool signalp[signalp_length];
		int expression;	//Current expression state.
		int express;		//New expression state.

		Gene(int k);	//Can we do without the default constructor?
		Gene(int k, int typ, int thr, bool tsq[], bool sig[], int exp);
		explicit Gene(const Gene &gene);
		virtual ~Gene();

		virtual Bead* Clone() const;
		void Randomize();

		bool Mutate(double mut_factor);
		void DefineTypeFromSeq();
		string Show(bool terminal, bool type_only=false) const;

};

#endif
