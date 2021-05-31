#include "Population.hh"

Population::Population()
{
	int i,j;
	id_count=0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		Space[i][j]=NULL;
	}
}

Population::~Population()
{
	int i,j;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if((Space[i][j])!=NULL)
		{
			delete (Space[i][j]);
			Space[i][j]=NULL;
		}
	}
}

void Population::UpdatePopulation()
{
	//Possible speedups by sorting the symbionts in between certain dynamics...
	//Currently, update order of cells is random, but per cell it is fixed (first host, then all symbionts in order); this should be changed in the future.

	int update_order[NR*NC];
	int u, i, j, s, x;
	double nutrients;
	i_symbiont iS;
	Organelle* SymbiontCopy;

	// if(Time%TimeSaveGrid==0)	PrintFieldToFile();
	// if(Time%TimePruneFossils==0 && Time!=0)	PruneFossilRecord();
	// if(Time%TimeOutputFossils==0 && Time!=0)	Fossils->ExhibitFossils();
	// if(Time%TimeSaveBackup==0 && Time!=0)	OutputBackup();
	if(Time%TimeTerminalOutput==0)	ShowGeneralProgress();

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;	//Row index.
		j = update_order[u]%NC;	//Column index.

		if (Space[i][j] != NULL)	//There is life at this site; so there is a host -- we do host dynamics here, symbiont dynamics are done mostly in Cell.cc??
		{

			//Here mix the symbionts if any, and them put them all in the first few spaces. So that we get fair order of updating but do not waste time on empty spots.

			//Basal death events for hosts and symbionts.
			if (uniform() < death_rate_host)
			{
				DeathOfHost(i, j);
				continue;
			}

			iS = Space[i][j]->SymbiontList->begin();
			while(iS != Space[i][j]->SymbiontList->end())
			{
				if (uniform() < death_rate_symbiont)
				{
					DeathOfSymbiont(i, j, iS);
				}
				iS++;
			}
			if (Space[i][j]->Host == NULL)	continue;	//After potential symbiont death, the host might have died.

			//Failed division events for hosts and endosymbionts.
			if (Space[i][j]->Host.Stage == 5)
			{
				DeathOfHost(i, j);
				continue;
			}

			iS = iS=Space[i][j]->SymbiontList->begin();
			while(iS != Space[i][j]->SymbiontList->end())
			{
				if ( (*iS).Stage == 5)
				{
					DeathOfSymbiont(i, j, iS);
				}
				iS++;
			}
			if (Space[i][j]->Host == NULL)	continue;	//After potential symbiont death, the host might have died.

			//Successfull division events for hosts and endosymbionts.
			if (Space[i][j]->Host.Stage == 4)
			{
				coords neigh = PickNeighbour(i, j);
				if (Space[neigh.first][neigh.second] != NULL)		DeathOfHost(neigh.first, neigh.second);
				id_count++;
				Space[neigh.first][neigh.second] = new Cell();
				Space[neigh.first][neigh.second]->Host = new Organelle();
				Space[neigh.first][neigh.second]->Host->Mitosis(Space[i][j]->Host, id_count);
				Space[neigh.first][neigh.second]->DNATransferToHost();

				iS = Space[i][j]->SymbiontList->begin();
				while (iS != Space[i][j]->SymbiontList->end())
				{
					if (uniform() < 0.5)	//Transfer the symbiont to the new cell (just a matter of using the right pointers).
					{
						SymbiontCopy = new Organelle();
						SymbiontCopy->CloneOrganelle(*iS);	//Does this work?
						Space[neigh.first][neigh.second]->SymbiontList->push_back(SymbiontCopy);
						Space[neigh.first][neigh.second]->nr_symbionts++;
						DeathOfSymbiont(i, j, iS);
					}
					else
					{
						iS++;	//See above, check that the erase does what I expected.
					}
				}
			}

			iS = Space[i][j]->SymbiontList->begin();
			while (iS != Space[i][j]->SymbiontList->end())
			{
				if ((*iS).Stage == 4)
				{
					id_count++;
					SymbiontCopy = new Organelle();
					SymbiontCopy->Mitosis((*iS), id_count);
					SymbiontList->push_back(SymbiontCopy);
					Space[i][j]->DNATransfertoSymbiont(SymbiontCopy);
					Space[i][j]->nr_symbionts++;
				}
				iS++;
			}

			//Expression dynamics within cell.
			Space[i][j]->UpdateOrganelles();

			//Do replication of host and symbiont genomes.
			nutrients = CollectNutrients(i, j);
			if (Space[i][j]->Host.Stage == 2   &&   Space[i][j]->Host.privilige == true)
			{
				Space[i][j]->Host->Replicate(nutrients);
			}

			iS = Space[i][j]->SymbiontList->begin();
			while (iS != Space[i][j]->SymbiontList->end())
			{
				if ( (*iS).Stage == 2   &&   (*iS).privilige == true)
				{
					(*iS)->Replicate(nutrients);
				}
				iS++;
			}

		}
	}
}

void Population::DeathOfSymbiont(int i, int j, i_symbiont iS)
{
	delete (*iS);
	iS = Space[i][j]->SymbiontList->erase(iS);
	Space[i][j]->nr_symbionts--;
	if (Space[i][j]->nr_symbionts == 0)
	{
		DeathOfHost(i, j);
	}
}

void Population::DeathOfHost(int i, int j)
{
	delete Space[i][j];	//All internal objects are deleted by deconstructors.
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

	return (nutrient_abundance / (double)Space[i][j]->nr_symbionts)  *  (1. - (double)organelle_density / max_organelle_density)
}

void Population::InitialisePopulation()
{
	Cell* CellZero;
	Cell* CellOne;
	int i, j;

	//First create one Cell.
	CellZero = new Cell();
	CellZero->InitialiseCell(id_count);

	//Now fill the field with this cell.
	for(i=0; i<NR; i++) for(j=0; j<NC; j++){
		// if(row<20 && col<20)	//Initialise square
		if (uniform() < 0.1)	//Initialise lower density
		{
			id_count_++;	//Make sure the first individual gets p_id_count_ of 1.
			CellOne = new Cell();
			CellOne->CloneCell(CellZero, id_count);
			CellOne->Ancestor = NULL;
			Space[i][j] = CellOne;
			//Fossils something...
		}
	}

	delete CellZero;	//I cannot delete PP_Copy, because each is actually turned into one of grid spaces. I can however delete this single bit of memory.
	CellZero = NULL;
}


void Population::ShowGeneralProgress()
{
	int i, j, host_count=0, symbiont_count=0;
	int stage_counts[6] = {0, 0, 0, 0, 0, 0};
	i_symbiont iS;

	for (i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		if ( Space[i][j] != NULL )
		{
			host_count++;
			symbiont_count += Space[i][j]->nr_symbionts;
			stage_counts[Space[i][j]->Host.Stage]++;
			iS = Space[i][j]->SymbiontList->begin();
			while (iS != Space[i][j]->SymbiontList->end())
			{
				stage_counts[(*iS).Stage]++;
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
