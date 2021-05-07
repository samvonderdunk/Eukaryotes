#include "Genome.hh"

Genome::Genome() {
	BeadList=NULL;
	g_length=0;
	gnr_genes=0;
	gnr_bsites=0;
	gnr_transporters=0;
	gnr_houses=0;
	fork_position=0;
	terminus_position=0;
	is_mutated=false;
}

Genome::~Genome() {
	iter i;

	if(BeadList!=NULL) {
		i=BeadList->begin();
		while(i!=BeadList->end()) {
			delete (*i);
			i++;
		}

		i=BeadList->erase(BeadList->begin(),BeadList->end());
		delete BeadList;
		BeadList=NULL;
	}
}
