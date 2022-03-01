#ifndef OrganelleHeader
#define OrganelleHeader

#include "Header.hh"
#include "Genome.hh"

class Organelle{
	public:
		int Stage;
		bool privilige;
		Genome* G;
		int nr_houses;
		double fitness;
		double nutrient_claim;

		//Fossil info.
		bool alive;
		bool mutant;
		bool lifetime_mutant;	//Explicitly set whenever mutant is false (i.e. no mutations during the organelle's own division, but later during its life due to transfer).
		bool exp_gene_transfer;	//Checks if an expressed gene was transferred through cut-and-paste, so that ExpressedGenes should be reset.
		int time_of_appearance;
		unsigned long long fossil_id;
		Organelle* Ancestor;

		typedef std::list<Bead*>::iterator i_bead;	//Not used anywhere currently.
		typedef std::list<Regulator*>::iterator i_reg;

		Organelle();
		~Organelle();

		void UpdateState();
		int EvaluateState(int eval_state, int* readout);

		double CalculateFitness(int target_nr, double real_nr);

		void Mitosis(Organelle* parent, unsigned long long id_count);
		void Replicate(double resource);
		void Abort();

		void InitialiseOrganelle(string genome, string expression, string definition);
		void CloneOrganelle(Organelle* ImageO, unsigned long long id_count);

		string Show(bool backup);
		string Output(bool backup);

};

#endif
