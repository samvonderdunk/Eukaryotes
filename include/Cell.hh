#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Organelle.hh"

class Cell{
	public:
		Organelle* Host;
		Organelle* Symbionts[HS];
		int nr_symbionts;
		//Perhaps will also get its own list of molecules present, and perhaps a hash table for sequence matching...

		typedef std::list<Bead*>::iterator i_bead;

		Cell();
		~Cell();

		void UpdateOrganelles();
		void RegulatorTransport();	//Including leakage.
		bool ActiveTransport(i_bead it, list<Bead*>* SourceCompartment, list<Bead*>* TargetCompartment);

};

#endif

//NOTE: here we probably need to define gene types, because they can go from one organelle to another.
