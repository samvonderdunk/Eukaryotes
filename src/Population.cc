#include "Population.hh"

Population::Population()
{
	int i,j;
	cell_count=0;
	organelle_count=0;

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

	int update_order[NR*NC];
	int u, i, j, s, x;
	double nutrients;

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;	//Row index.
		j = update_order[u]%NC;	//Column index.

		if (Space[i][j] != NULL)	//There is life at this site; so there is a host -- we do host dynamics here, symbiont dynamics are done mostly in Cell.cc??
		{

			//Basal death events for hosts and symbionts.
			if (uniform() < death_rate_host)
			{
				DeathOfHost(i, j);
				continue;
			}

			for (s=0;s<HS;s++)
			{
				if (Space[i][j]->Symbionts[s] != NULL)
				{
					if (uniform() < death_rate_symbiont)
					{
						DeathOfSymbiont(i, j, s);
					}
				}
			}
			if (Space[i][j]->Host == NULL)	continue;	//After potential symbiont death, the host might have died.

			//Failed division events for hosts and endosymbionts.
			if (Space[i][j]->Host.Stage == 5)
			{
				DeathOfHost(i, j);
				continue;
			}

			for (s=0;s<HS;s++)
			{
				if (Space[i][j]->Symbionts[s] != NULL)
				{
					if (Space[i][j]->Symbionts[s].Stage == 5)
					{
						DeathOfSymbiont(i, j, s);
					}
				}
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
				for (s=0;s<HS;s++)
				{
					if (Space[i][j]->Symbionts[s] != NULL && uniform() < 0.5)
					{
						Space[neigh.first][neigh.second]->Symbionts[s] = new Organelle();
						Space[neigh.first][neigh.second]->Symbionts[s]->Mitosis(Space[i][j]->Symbionts[s], id_count)
					}
				}
			}

			for (s=0;s<HS;s++)
			{
				if (Space[i][j]->Symbionts[s] != NULL)
				{
					if (Space[i][j]->Symbionts[s].Stage == 4)
					{
						while (x == s)		x = (int)(uniform()*HS);	//Don't pick yourself!
						if (Space[i][j]->Symbionts[x] != NULL)		DeathOfSymbiont(i, j, x);
						id_count++;
						Space[i][j]->Symbionts[x] = new Organelle();
						Space[i][j]->Symbionts[x]->Mitosis(Space[i][j]->Symbionts[s], id_count);
					}
				}
			}

			//Expression dynamics within cell.
			Space[i][j]->UpdateOrganelles();

			//Do replication of host and symbiont genomes.
			nutrients = CollectNutrients(i, j);
			if (Space[i][j]->Host.Stage == 2   &&   Space[i][j]->Host.privilige == true)
			{
				Space[i][j]->Host->Replicate(nutrients);
			}

			for (s=0;s<HS;s++)
			{
				if (Space[i][j]->Symbionts[s] != NULL)
				{
					if (Space[i][j]->Symbionts[s].Stage == 2   &&   Space[i][j]->Symbionts[s].privilige == true)
					{
						Space[i][j]->Symbionts[s]->Replicate(nutrients);
					}
				}
			}

		}
	}
}

void Population::DeathOfSymbiont(int i, int j, int s)
{
	delete (Space[i][j]->Symbionts[s]);
	Space[i][j]->Symbionts[s] = NULL;
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
