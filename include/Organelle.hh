#ifndef OrganelleHeader
#define OrganelleHeader

#include "Header.hh"
#include "Genome.hh"

class Organelle{
	public:
		int Stage;
		bool privilige;
		std::list<Bead*>* ExpressedGenes;
		Genome* G;
		double fitness;

		typedef std::list<Bead*>::iterator i_bead;

		Organelle();
		~Organelle();

		void UpdateExpression();	//Mostly defers to Genome-level function.
		void UpdateState();
		int EvaluateState(int eval_state);

		void Mitosis(Organelle* parent, unsigned long long id_count);
		void Replicate(double resource);

};

#endif
