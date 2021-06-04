#include "Population.hh"

Population::Population()
{
	int i,j;
	id_count = 0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		Space[i][j]=NULL;
	}
	FossilSpace = new Fossils();
}

Population::~Population()	//Skips deconstructor of Cell, because we need to check whether organelles are also in the FossilSpace.
{
	int i,j,s;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if((Space[i][j]) != NULL)
		{
			if (Space[i][j]->Host->mutant)	FossilSpace->EraseFossil(Space[i][j]->Host->fossil_id);
			delete Space[i][j]->Host;
			Space[i][j]->Host = NULL;
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if (Space[i][j]->Symbionts->at(s)->mutant)	FossilSpace->EraseFossil(Space[i][j]->Symbionts->at(s)->fossil_id);
				delete Space[i][j]->Symbionts->at(s);
				Space[i][j]->Symbionts->at(s) = NULL;
			}
		}
	}

	delete FossilSpace;
	FossilSpace = NULL;
}

void Population::UpdatePopulation()
{
	int update_order[NR*NC];
	int u, i, j, s;
	double nutrients;
	Organelle* SymbiontCopy;

	// if(Time%TimeSaveGrid==0)	PrintFieldToFile();
	if(Time%TimePruneFossils==0 && Time!=0)	PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)	FossilSpace->ExhibitFossils();
	// if(Time%TimeSaveBackup==0 && Time!=0)	OutputBackup();
	if(Time%TimeTerminalOutput==0)	ShowGeneralProgress();

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;
		j = update_order[u]%NC;

		if (Space[i][j] != NULL)
		{
			random_shuffle(Space[i][j]->Symbionts->begin(), Space[i][j]->Symbionts->end());

			/* Basal death */
			if (uniform() < death_rate_host)
			{
				DeathOfCell(i, j);
				continue;
			}
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if (uniform() < death_rate_symbiont)
				{
					DeathOfSymbiont(i, j, s);
					Space[i][j]->nr_symbionts--;
				}
			}
			if (CheckCellDeath(i, j))	continue;
			/* ~Basal death */

			/* Failed division */
			if (Space[i][j]->Host->Stage == 5)
			{
				DeathOfCell(i, j);
				continue;
			}
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if ( Space[i][j]->Symbionts->at(s)->Stage == 5)
				{
					DeathOfSymbiont(i, j, s);
					Space[i][j]->nr_symbionts--;
				}
			}
			if (CheckCellDeath(i, j))	continue;
			/* ~Failed division */

			/* Division of host */
			if (Space[i][j]->Host->Stage == 4)
			{
				coords neigh = PickNeighbour(i, j);
				if (Space[neigh.first][neigh.second] != NULL)		DeathOfCell(neigh.first, neigh.second);
				Space[neigh.first][neigh.second] = new Cell();
				Space[neigh.first][neigh.second]->Host = new Organelle();
				id_count++;
				Space[neigh.first][neigh.second]->Host->Mitosis(Space[i][j]->Host, id_count);

				for (s=0; s<Space[i][j]->nr_symbionts; s++)
				{
					if (uniform() < 0.5)	//Transfer the symbiont to the new cell (just a matter of using the right pointers).
					{
						Space[neigh.first][neigh.second]->Symbionts->push_back(Space[i][j]->Symbionts->at(s));
						Space[i][j]->Symbionts->erase(Space[i][j]->Symbionts->begin()+s);
						Space[neigh.first][neigh.second]->nr_symbionts++;
						Space[i][j]->nr_symbionts--;
					}
				}
				if (CheckCellDeath(i, j))	continue;
				if (!CheckCellDeath(neigh.first, neigh.second))	//Don't continue, since we have killed a cell at a different site.
				{
					Space[neigh.first][neigh.second]->DNATransferToHost();	//Instead, use the cell death check to know whether we need to do DNA transfer.
					if (Space[neigh.first][neigh.second]->Host->mutant)
					{
						FossilSpace->BuryFossil(Space[neigh.first][neigh.second]->Host);	//Bury fossil only after potential transfer events (also count as mutations).
					}
				}
			}
			/* ~Division of host */

			/* Division of symbionts */
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if (Space[i][j]->Symbionts->at(s)->Stage == 4)
				{
					SymbiontCopy = new Organelle();
					id_count++;
					SymbiontCopy->Mitosis(Space[i][j]->Symbionts->at(s), id_count);
					Space[i][j]->Symbionts->push_back(SymbiontCopy);
					Space[i][j]->DNATransfertoSymbiont(SymbiontCopy);
					Space[i][j]->nr_symbionts++;
					if (Space[i][j]->Symbionts->at(Space[i][j]->nr_symbionts-1)->mutant)
					{
						FossilSpace->BuryFossil(Space[i][j]->Symbionts->at(Space[i][j]->nr_symbionts-1));
					}
				}
			}
			/* ~Division of symbionts */

			Space[i][j]->UpdateOrganelles();	//Expression dynamics within cell.

			/* Replication */
			nutrients = CollectNutrients(i, j);
			if (Space[i][j]->Host->Stage == 2   &&   Space[i][j]->Host->privilige)
			{
				Space[i][j]->Host->Replicate(nutrients);
			}
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if ( Space[i][j]->Symbionts->at(s)->Stage == 2   &&   Space[i][j]->Symbionts->at(s)->privilige)
				{
					Space[i][j]->Symbionts->at(s)->Replicate(nutrients);
				}
			}
			/* ~Replication */

		}
	}
}

bool Population::CheckCellDeath(int i, int j)
{
	if (Space[i][j]->nr_symbionts == 0)
	{
		DeathOfCell(i, j);
		return true;				//The cell died.
	}
	else
	{
		return false;				//The cell lives.
	}
}

void Population::DeathOfSymbiont(int i, int j, int s)
{
	if (!Space[i][j]->Symbionts->at(s)->mutant)
	{
		delete Space[i][j]->Symbionts->at(s);
	}
	else
	{
		Space[i][j]->Symbionts->at(s)->alive = false;
	}

	Space[i][j]->Symbionts->erase(Space[i][j]->Symbionts->begin()+s);	//We erase the pointer from the Symbionts vector. Similar to setting Space[i][j]->Host to NULL for host death (see below).
}

void Population::DeathOfHost(int i, int j)
{
	if (!Space[i][j]->Host->mutant)
	{
		delete Space[i][j]->Host;
	}
	else
	{
		Space[i][j]->Host->alive = false;
	}

	Space[i][j]->Host = NULL;
}

void Population::DeathOfCell(int i, int j)	//Called when host died or when last symbiont died, responsible for cleaning of remaining organelles.
{
	int s;

	if (Space[i][j]->Host != NULL)	DeathOfHost(i, j);
	for (s=0; s<Space[i][j]->nr_symbionts; s++)
	{
		DeathOfSymbiont(i, j, s);
		Space[i][j]->nr_symbionts--;
	}
	Space[i][j] = NULL;
}

Population::coords Population::PickNeighbour(int i, int j)
{
	int nrow = i, ncol = j, ni, nj, random_neighbour;

	while (nrow == i && ncol == j)	//Try again if you pick yourself.
	{
		random_neighbour = (int)(uniform()*9);
		ni = random_neighbour/3;
		nj = random_neighbour%3;

		//Wrap grid boundaries
		nrow = i+ni-1;
		if(nrow < 0)	nrow += NR;
		else if(nrow >= NR)	nrow -= NR;

		ncol = j+nj-1;
		if(ncol < 0)	ncol += NC;
		else if(ncol >= NC)	ncol -= NC;
	}
	return std::make_pair(nrow, ncol);
}

double Population::CollectNutrients(int i, int j)
{
	//Collects nutrients per organelle at site (i,j).
	//
	//	n_ij = n_ext / x_i (1 - sum(x)/k)
	//
	//	I don't know whether I should cast things to double only during calculations or just always have them as doubles.
	//
	int ii, jj, nrow, ncol;
	int organelle_density = 0;

	//First obtain organelle density in 3x3 neighbourhood.
	for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
	{
		if (ii == i && jj == j)	continue;

		if (ii < 0)	nrow = ii + NR;	//-1 becomes -1+100=99, i.e. the last index of a row with 100 sites (0-99).
		else if (ii >= NR)	nrow = ii - NR;
		else	nrow = ii;
		if (jj < 0)	ncol = jj + NC;
		else if (jj >= NC)	ncol = jj - NC;
		else	ncol = jj;

		if (Space[nrow][ncol] != NULL)
		{
			organelle_density += Space[nrow][ncol]->nr_symbionts + 1;		//Host always counts for one.
		}
	}

	// return (nutrient_abundance / (double)Space[i][j]->nr_symbionts)  *  (1. - (double)organelle_density / max_organelle_density);
	return (nutrient_abundance - (double)organelle_density) / (double)Space[i][j]->nr_symbionts;
}


void Population::InitialisePopulation()
{
	Cell* CellZero;
	Cell* CellOne;
	int i, j, s;

	//First create one Cell.
	CellZero = new Cell();
	CellZero->InitialiseCell();

	//Now fill the field with this cell.
	for(i=0; i<NR; i++) for(j=0; j<NC; j++){
		// if(row<20 && col<20)	//Initialise square
		if (uniform() < 0.1)	//Initialise lower density
		{
			CellOne = new Cell();
			CellOne->CloneCell(CellZero, &id_count);
			Space[i][j] = CellOne;
			FossilSpace->BuryFossil(Space[i][j]->Host);
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				FossilSpace->BuryFossil(Space[i][j]->Symbionts->at(s));
			}
		}
	}

	delete CellZero;	//I cannot delete PP_Copy, because each is actually turned into one of grid spaces. I can however delete this single bit of memory.
	CellZero = NULL;
}

/* ######################################################################## */
/*				FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS						*/
/* ######################################################################## */

void Population::PruneFossilRecord()
{
	std::list<unsigned long long> AllFossilIDs;
	int s, fossil_record_size;
	unsigned long long fossilID;
	i_fos iF;
	i_ull findit;
	Organelle* lastCA;

	//Trace all living individuals back to the beginning, storing all individuals along their lineage.
	for(int i=0; i<NR; i++)	for(int j=0; j<NC; j++)
	{
		if(Space[i][j] != NULL)
		{
			lastCA = Space[i][j]->Host->Ancestor;
			while(lastCA != NULL)
			{
				AllFossilIDs.push_back(lastCA->fossil_id);
				lastCA = lastCA->Ancestor;
			}

			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				lastCA = Space[i][j]->Symbionts->at(s)->Ancestor;
				while(lastCA != NULL)
				{
					AllFossilIDs.push_back(lastCA->fossil_id);
					lastCA = lastCA->Ancestor;
				}
			}

		}
	}
	// Delete duplicates (e.g. Agent 9 in example asci will be located 4 times. Agent 11 two times, etc.)
	AllFossilIDs.sort();
	AllFossilIDs.unique();

	// Delete all in FossilList that are not in AllFossilIDs (unless they are still living):
	fossil_record_size = (*FossilSpace).FossilRecord.size();
	iF = (*FossilSpace).FossilRecord.begin();
	while(iF != (*FossilSpace).FossilRecord.end())
	{
		// Search if stored agent was also found by tracing back:
		fossilID = (*iF)->fossil_id;
		findit = std::find(AllFossilIDs.begin(),AllFossilIDs.end(),fossilID);
		// If not, delete the fossil unless it is still alive or is still saved in the graveyard. If a prokaryote dies, the graveyard-flag remains for one ShowGeneralProgress() cycle at most, so that the fossil can be deleted at the next pruning step. If ShowGeneralProgress() precedes PruneFossilRecord(), this is issue is even avoided, because flags are already removed off dead prokaryotes.
		if(findit==AllFossilIDs.end() && !(*iF)->alive)
		{
			delete *iF;
			iF = FossilSpace->FossilRecord.erase(iF);
		}
		else
		{
			++iF;
		}
	}

	AllFossilIDs.clear();
	cout << "ID count (" << id_count << ")\tFossil record (pruned from " << fossil_record_size << " to " << (*FossilSpace).FossilRecord.size() << ")" << endl;
}

/* ######################################################################## */
/*				OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT						*/
/* ######################################################################## */

void Population::ShowGeneralProgress()
{
	int i, j, host_count=0, symbiont_count=0;
	int stage_counts[6] = {0, 0, 0, 0, 0, 0};
	i_org iS;

	for (i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		if ( Space[i][j] != NULL )
		{
			host_count++;
			symbiont_count += Space[i][j]->nr_symbionts;
			stage_counts[Space[i][j]->Host->Stage]++;
			iS = Space[i][j]->Symbionts->begin();
			while (iS != Space[i][j]->Symbionts->end())
			{
				stage_counts[(*iS)->Stage]++;
				iS++;
			}
		}
	}

	cout << "Time " << Time;
	cout << "\tHosts " << host_count;
	cout << "\tSymbionts " << symbiont_count;
	for (i=0; i<6; i++)	cout << "\tS(" << i << ") " << stage_counts[i];
	cout << endl;

}
