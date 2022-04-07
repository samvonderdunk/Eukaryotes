#ifndef ProkaryoteHeader
#define ProkaryoteHeader

#include "Header.hh"
#include "Organelle.hh"
#include "Cell.hh"

class Prokaryote : public Cell {
	public:
		Organelle* Vesicle;

		Prokaryote();
		virtual ~Prokaryote();

		void CalculateCellFitness();
		bool BasalDeath();
		void DeathOfCell();

		bool FailedDivision();
		Cell* Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count);
		void UpdateOrganelles();
		void Replication(nuts n);

		void InitialiseCell(int input_nr);
		void CloneCell(Cell* ImageC, unsigned long long* pid_count);

		string Show(bool include_organelles, bool include_genomes);
};

#endif
