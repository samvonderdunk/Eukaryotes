#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Organelle.hh"
#include "Fossils.hh"

class Cell{
	public:
		int kind;	//Prokaryote (0) or Eukaryote (1)
		int barcode;

		typedef std::tuple<double,double,double> nuts;
		typedef std::list<Bead*>::iterator i_bead;
		typedef std::list<Regulator*>::iterator i_reg;
		typedef std::vector<Organelle*>::iterator i_org;

		Cell();	//Default constructor should not be used (instead call with kind, as below).
		Cell(int k);
		virtual ~Cell();

		virtual void CalculateCellFitness()=0;
		virtual bool BasalDeath()=0;
		virtual void DeathOfCell()=0;
		virtual bool FailedDivision()=0;
		virtual Cell* Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count)=0;
		virtual void UpdateOrganelles()=0;
		virtual void Replication(nuts n)=0;

		i_bead TransferGene(i_bead it, Organelle* Source, Organelle* Target, bool include_distal, bool cut_and_paste);
		void TransferBead(i_bead it, Organelle* Target);

		virtual void InitialiseCell(int input_nr)=0;
		virtual void CloneCell(Cell* ImageC, unsigned long long* pid_count)=0;

		string Show();
		virtual string Show(bool include_organelles, bool include_genomes)=0;
};

#endif

//NOTE: here we probably need to define gene types, because they can go from one organelle to another.
