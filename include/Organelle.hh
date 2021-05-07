#ifndef OrganelleHeader
#define OrganelleHeader

#include "Header.hh"
#include "Genome.hh"

class Organelle{
	public:
		int Stage;
		std::list<Bead*>* ExpressedGenes;
		Genome* G;

		typedef std::list<Bead*>::iterator i_bead;

		Organelle();
		~Organelle();

		void UpdateExpression();	//Mostly defers to Genome-level function.
		void UpdateState();
		int EvaluateState(int eval_state);

};

#endif
