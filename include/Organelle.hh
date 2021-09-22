#ifndef OrganelleHeader
#define OrganelleHeader

#include "Header.hh"
#include "Genome.hh"

class Organelle{
	public:
		int Stage;
		bool privilige;
		std::list<Bead*>* ExpressedGenes;
		int nr_native_expressed;
		Genome* G;
		double fitness;
		double nutrient_claim;

		//Fossil info.
		bool alive;
		bool mutant;
		int time_of_appearance;
		unsigned long long fossil_id;
		Organelle* Ancestor;

		typedef std::list<Bead*>::iterator i_bead;

		Organelle();
		~Organelle();

		void UpdateState();
		int EvaluateState(int eval_state, int* readout);

		void Mitosis(Organelle* parent, unsigned long long id_count);
		void Replicate(double resource);
		void Abort();

		void InitialiseOrganelle(string genome, string expression);
		void CloneOrganelle(Organelle* ImageO, unsigned long long id_count);

		string Show(bool backup);
		string Output(bool backup);

};

#endif
