#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Organelle.hh"

class Cell{
	public:
		Organelle* Host;
		std::vector<Organelle*>* Symbionts;
		int nr_symbionts;	//Perhaps I don't need to keep track of this, and could just determine Symbionts on the go always?
		//Perhaps will also get its own list of molecules present, and perhaps a hash table for sequence matching...
		int barcode;

		typedef std::list<Bead*>::iterator i_bead;
		typedef std::vector<Organelle*>::iterator i_org;

		Cell();
		~Cell();

		void UpdateOrganelles();
		void RegulatorTransport();	//Including leakage & all transport at this point.

		void DNATransferToHost();
		void DNATransfertoSymbiont(Organelle* Symbiont);
		void TransferGene(i_bead it, Organelle* Source, Organelle* Target);
		void TransferBead(i_bead it, Organelle* Target);

		void InitialiseCell(int input_nr);
		void CloneCell(Cell* ImageC, unsigned long long* pid_count);

		bool CheckCellDeath(bool output);
		void DeathOfSymbiont(int s);
		void DeathOfHost();
		void DeathOfCell();
		void SingleCellOutput(bool death_event);

};

#endif

//NOTE: here we probably need to define gene types, because they can go from one organelle to another.
