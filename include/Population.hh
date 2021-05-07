#ifndef PopulationHeader
#define PopulationHeader

#include "Header.hh"
#include "Cell.hh"
#include "FossilRecord.hh"
// #include <cstdio>
// #include <stdlib.h>

class Population
{
	public:
		Cell* Space[NR][NC];
		// FossilRecord* HostFossils;
		// FossilRecord* SymbiontFossils;

		unsigned long long cell_count;
		unsigned long long organelle_count;

		Population();
		~Population();
};

#endif
