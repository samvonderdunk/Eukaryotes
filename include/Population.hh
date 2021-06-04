#ifndef PopulationHeader
#define PopulationHeader

#include "Header.hh"
#include "Cell.hh"
#include "Fossils.hh"
// #include <cstdio>
// #include <stdlib.h>

class Population
{
	public:
		Cell* Space[NR][NC];
		Fossils* FossilSpace;

		unsigned long long id_count;

		typedef std::pair<int,int> coords;
		typedef std::vector<Organelle*>::iterator i_org;
		typedef std::list<Organelle*>::iterator i_fos;
		typedef std::list<unsigned long long>::iterator i_ull;

		Population();
		~Population();

		void UpdatePopulation();

		bool CheckCellDeath(int i, int j);
		void DeathOfSymbiont(int i, int j, int s);
		void DeathOfHost(int i, int j);
		void DeathOfCell(int i, int j);
		coords PickNeighbour(int i, int j);
		double CollectNutrients(int i, int j);

		void InitialisePopulation();

		void PruneFossilRecord();
		void OutputBackup();
		void ShowGeneralProgress();
};

#endif
