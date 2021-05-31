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
		typedef std::list<Organelle*>::iterator i_symbiont;

		Population();
		~Population();

		void UpdatePopulation();
		void DeathOfSymbiont(int i, int j, i_symbiont iS);
		void DeathOfHost(int i, int j);
		coords PickNeighbour(int i, int j);
		double CollectNutrients(int i, int j);

		void InitialisePopulation();

		void ShowGeneralProgress();
};

#endif
