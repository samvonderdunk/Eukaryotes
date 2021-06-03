#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Organelle.hh"

class Cell{
	public:
		Organelle* Host;
		std::list<Organelle*>* SymbiontList;
		int nr_symbionts;	//Perhaps I don't need to keep track of this, and could just determine SymbiontList on the go always?
		//Perhaps will also get its own list of molecules present, and perhaps a hash table for sequence matching...

		typedef std::list<Bead*>::iterator i_bead;
		typedef std::list<Organelle*>::iterator i_symbiont;

		Cell();
		~Cell();

		void UpdateOrganelles();
		void RegulatorTransport();	//Including leakage.
		bool ActiveTransport(i_bead it, list<Bead*>* SourceCompartment, list<Bead*>* TargetCompartment);

		void DNATransferToHost();
		void DNATransfertoSymbiont(Organelle* Symbiont);
		void TransferGene(i_bead it, Organelle* Source, Organelle* Target);
		void TransferBead(i_bead it, Organelle* Target);

		void InitialiseCell();
		void CloneCell(Cell* ImageC, unsigned long long id_count);

};

#endif

//NOTE: here we probably need to define gene types, because they can go from one organelle to another.
