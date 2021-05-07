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

		Cell();
		~Cell();

};

#endif
