#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Organelle.hh"
#include "Fossils.hh"

class Cell{
	public:
		Organelle* Vesicle;
		int barcode;	//Track different strains or individuals.

		typedef std::list<Bead*>::iterator i_bead;
		typedef std::list<Regulator*>::iterator i_reg;
		typedef std::vector<Organelle*>::iterator i_org;

		Cell();	//Default constructor should not be used (instead call with kind, as below).
		Cell(int k);
		virtual ~Cell();

		void CalculateCellFitness();
		bool BasalDeath();
		void DeathOfCell();
		bool FailedDivision();
		Cell* Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count);
		void UpdateOrganelles();
		void Replication(double nuts);

		i_bead TransferGene(i_bead it, Organelle* Source, Organelle* Target, bool include_distal, bool cut_and_paste);
		void TransferBead(i_bead it, Organelle* Target);

		void InitialiseCell(int input_nr);
		void CloneCell(Cell* ImageC, unsigned long long* pid_count);

		string Show();
		string Show(bool include_organelles, bool include_genomes);
};

#endif

//NOTE: here we probably need to define gene types, because they can go from one organelle to another.
