#ifndef PopulationHeader
#define PopulationHeader

#include "Header.hh"
#include "Cell.hh"
#include "Fossils.hh"

class Population
{
	public:
		Cell* Space[NR_max][NC_max];
		double NutrientSpace[NR_max][NC_max];
		Fossils* FossilSpace;
		std::vector<unsigned long long> Lineage;

		unsigned long long id_count;
		int nr_strains;

		typedef std::pair<int,int> coords;
		typedef std::tuple<double,double,double> nuts;
		typedef std::vector<Organelle*>::iterator i_org;
		typedef std::list<Organelle*>::iterator i_fos;
		typedef std::list<unsigned long long>::iterator i_ull;
		typedef std::vector<unsigned long long>::iterator i_lin;

		Population();
		~Population();

		void UpdatePopulation();

		void CloneSymbiont(int source_i, int source_j, int s, Cell* NewCell);
		void WellMix();
		bool CheckLineage(int i, int j);
		void LogLineage(int i, int j);
		coords PickNeighbour(int i, int j);
		void CollectNutrientsFromSite(int i, int j);
		nuts HandleNutrientClaims(int i, int j);

		void InitialisePopulation();
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();
		Organelle* FindInFossilRecord(unsigned long long AncID);

		void PruneFossilRecord();
		void OutputGrid(bool backup);
		void OutputLineage(int i, int j);
		void ShowGeneralProgress();
};

#endif
