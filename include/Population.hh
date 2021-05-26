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

		unsigned long long id_count;

		typedef std::pair<int,int> coords;

		Population();
		~Population();

		void UpdatePopulation();
		void DeathOfSymbiont(int i, int j, int s);
		void DeathOfHost(int i, int j);
		coords PickNeighbour(int i, int j);
		double CollectNutrients(int i, int j);

		void InitialisePopulation();
};

#endif
