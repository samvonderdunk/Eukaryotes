#ifndef OrganelleHeader
#define OrganelleHeader

#include "Header.hh"
#include "Genome.hh"

class Organelle{
	public:
		int ExpressedGenes;
		int Stage;
		Genome* G;

		Organelle();
		~Organelle();

};

#endif
