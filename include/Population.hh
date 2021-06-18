#ifndef PopulationHeader
#define PopulationHeader

#include "Header.hh"
#include "Cell.hh"
#include "Fossils.hh"

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

		void FollowSingleCell();
		void ResetSingleCell(Cell** CP, Cell** CP_reset);

		void UpdatePopulation();

		coords PickNeighbour(int i, int j);
		double CollectNutrients(int i, int j);

		void InitialisePopulation();
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();
		Organelle* FindInFossilRecord(unsigned long long AncID);

		void PruneFossilRecord();
		void OutputGrid(bool backup);
		void ShowGeneralProgress();
};

#endif
