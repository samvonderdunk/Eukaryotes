#ifndef PopulationHeader
#define PopulationHeader

#include "Header.hh"
#include "Cell.hh"
#include "Fossils.hh"

class Population
{
	public:
		Cell* Space[NR][NC];
		double NutrientSpace[NR][NC];
		Fossils* FossilSpace;
		std::vector<unsigned long long> Lineage;

		unsigned long long id_count;

		typedef std::pair<int,int> coords;
		typedef std::vector<Organelle*>::iterator i_org;
		typedef std::list<Organelle*>::iterator i_fos;
		typedef std::list<unsigned long long>::iterator i_ull;
		typedef std::vector<unsigned long long>::iterator i_lin;

		Population();
		~Population();

		void FollowSingleCell();
		void ResetSingleCell(Cell** CP, Cell** CP_reset);

		void UpdatePopulation();

		bool CheckLineage(int i, int j);
		void LogLineage(int i, int j);
		coords PickNeighbour(int i, int j);
		void CollectNutrientsFromSite(int i, int j);

		void InitialisePopulation();
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();
		Organelle* FindInFossilRecord(unsigned long long AncID);
		void ReadLineageFile();

		void PruneFossilRecord(int m, int n);
		void OutputGrid(bool backup);
		void OutputLineage(int i, int j);
		void ShowGeneralProgress();
};

#endif
