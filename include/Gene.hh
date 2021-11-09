#ifndef GeneHeader
#define GeneHeader

#include "Bead.hh"
#include "Header.hh"

// Geney genes, i.e. transcription factors.

class Gene : public Bead {
	public:
		int type;	//Used to determine phenotypic readout from expression.
		int threshold; //At what regulatory effect does the gene start to become expressed.
		int activity;
		bool typeseq[typeseq_length];
		bool sequence[sequence_length];
		bool signalp[signalp_length];
		int expression;	//Current expression state.
		int express;		//New expression state.

		Gene();
		Gene(int typ, int thr, int act, bool tsq[], bool sig[], bool seq[], int exp);
		explicit Gene(const Gene &gene);
		virtual ~Gene();

		virtual Bead* Clone() const;
		void RandomGene();

		string Show(bool terminal, bool type_only) const;

};

#endif
