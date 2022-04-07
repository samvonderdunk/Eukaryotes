#ifndef EukaryoteHeader
#define EukaryoteHeader

#include "Header.hh"
#include "Organelle.hh"
#include "Cell.hh"

class Eukaryote : public Cell{
	public:
		Organelle* Host;
		std::vector<Organelle*>* Symbionts;
		int nr_symbionts;

		Eukaryote();
		virtual ~Eukaryote();

		void CalculateCellFitness();
		bool BasalDeath();
		void DeathOfCell();

		bool LostSymbionts(bool output);
		void DeathOfSymbiont(int s);
		void DeathOfHost();

		bool FailedDivision();
		Cell* Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count);
		void CloneSymbiont(int s, Eukaryote* SourceE, Eukaryote* TargetE, Fossils* FP, unsigned long long* pid_count);
		void SymbiontDivisions(Fossils* FP, unsigned long long* pid_count);

		void UpdateOrganelles();
		void Replication(nuts n);

		void GeneTransport();	//Including leakage & all transport at this point.

		void DNATransferToHost();
		void DNATransfertoSymbiont(Organelle* Symbiont);

		void InitialiseCell(int input_nr);
		void CloneCell(Cell* ImageC, unsigned long long* pid_count);

		string Show(bool include_organelles, bool include_genomes);
};

#endif
