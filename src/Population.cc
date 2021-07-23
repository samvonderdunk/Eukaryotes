#include "Population.hh"

Population::Population()
{
	int i,j;
	id_count = 0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		Space[i][j]=NULL;
		NutrientSpace[i][j]=0.;
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
			if (Space[i][j]->Host->mutant || trace_lineage || log_lineage)	FossilSpace->EraseFossil(Space[i][j]->Host->fossil_id);
			delete Space[i][j]->Host;
			Space[i][j]->Host = NULL;
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if (Space[i][j]->Symbionts->at(s)->mutant || trace_lineage || log_lineage)	FossilSpace->EraseFossil(Space[i][j]->Symbionts->at(s)->fossil_id);
				delete Space[i][j]->Symbionts->at(s);
				Space[i][j]->Symbionts->at(s) = NULL;
			}
		}
	}

	delete FossilSpace;
	FossilSpace = NULL;
}

void Population::FollowSingleCell()
{
	Cell* C_init, * C, * C_copy;
	Organelle* S_copy;
	int s, nutrients;

	C_init = new Cell();
	C_init->InitialiseCell();
	C = new Cell();
	C->CloneCell(C_init, &id_count);

	for (Time=TimeZero; Time<SimTime+1; Time++)
	{
		C->SingleCellOutput(false);	//Print all relevant data.

		//Failed division.
		if (C->Host->Stage == 5)
		{
			if (follow_with_fixed_symbionts)
			{
				//Just reset the host.
				delete C->Host;
				C->Host = NULL;
				C->Host = new Organelle();
				C->Host->CloneOrganelle(C_init->Host, id_count);
			}
			else
			{
				//Actually kill the host, and therewith the entire cell.
				C->SingleCellOutput(true);
				ResetSingleCell(&C, &C_init);
				continue;
			}
		}
		for (s=C->nr_symbionts-1; s>=0; s--)
		{
			if ( C->Symbionts->at(s)->Stage == 5)
			{
				if (follow_with_fixed_symbionts)
				{
					//Just reset the symbiont.
					delete C->Symbionts->at(s);
					C->Symbionts->at(s) = NULL;
					C->Symbionts->at(s) = new Organelle();
					C->Symbionts->at(s)->CloneOrganelle(C_init->Symbionts->at(s), id_count);
				}
				else
				{
					//Actually kill the symbiont.
					C->DeathOfSymbiont(s);
					C->nr_symbionts--;
				}
			}
		}
		if (C->CheckCellDeath(true))
		{
			ResetSingleCell(&C, &C_init);
			continue;
		}

		//Division of host.
		if (C->Host->Stage == 4   &&   (uniform() < C->Host->fitness))
		{
			C_copy = new Cell();
			C_copy->Host = new Organelle();
			id_count++;
			C_copy->Host->Mitosis(C->Host, id_count);

			if (!follow_with_fixed_symbionts)
			{
				for (s=C->nr_symbionts-1; s>=0; s--)
				{
					if (uniform() < 0.5)	//Transfer the symbiont to the new cell (just a matter of using the right pointers).
					{
						C_copy->Symbionts->push_back(C->Symbionts->at(s));
						C->Symbionts->erase(C->Symbionts->begin()+s);
						C_copy->nr_symbionts++;
						C->nr_symbionts--;
					}
				}
			}
			delete C_copy;	//We only track one of the offspring every time.
			C_copy = NULL;

			if (C->CheckCellDeath(true))
			{
				ResetSingleCell(&C, &C_init);
				continue;
			}
		}

		//Division of symbionts.
		for (s=0; s<C->nr_symbionts; s++)
		{
			if (C->Symbionts->at(s)->Stage == 4   &&   (uniform() < C->Symbionts->at(s)->fitness))
			{
				S_copy = new Organelle();
				id_count++;
				S_copy->Mitosis(C->Symbionts->at(s), id_count);
				if (follow_with_fixed_symbionts)
				{
					delete S_copy;
					S_copy = NULL;
				}
				else
				{
					S_copy->Ancestor = C->Symbionts->at(s);
					C->Symbionts->push_back(S_copy);
					C->nr_symbionts++;
				}
			}
		}

		C->UpdateOrganelles();	//Expression dynamics within cell.

		//Replication.
		nutrients = nutrient_abundance/(double)(C->nr_symbionts+1);
		if (C->Host->Stage == 2   &&   C->Host->privilige)
		{
			C->Host->Replicate(nutrients);
		}
		for (s=0; s<C->nr_symbionts; s++)
		{
			if ( C->Symbionts->at(s)->Stage == 2   &&   C->Symbionts->at(s)->privilige)
			{
				C->Symbionts->at(s)->Replicate(nutrients);
			}
		}

	}
}

void Population::ResetSingleCell(Cell** CP, Cell** CP_reset)
{
	delete *CP;
	*CP = NULL;

	*CP = new Cell();
	(*CP)->CloneCell(*CP_reset, &id_count);
}



void Population::UpdatePopulation()
{
	int update_order[NR*NC];
	int u, i, j, s;
	Organelle* SymbiontCopy;

	for(i=0; i<NR; i++)	for(j=0; j<NC; j++)	NutrientSpace[i][j] = 0.;
	for(i=0; i<NR; i++)	for(j=0; j<NC; j++)	CollectNutrientsFromSite(i,j);

	if(Time%TimeTerminalOutput==0)															ShowGeneralProgress();
	if(Time%TimeSaveGrid==0)																		OutputGrid(false);
	if(Time%TimePruneFossils==0 && Time!=0)											PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)										FossilSpace->ExhibitFossils();
	if(Time%TimeSaveBackup==0 && Time!=0)												OutputGrid(true);

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC], uniform_shuffle);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;
		j = update_order[u]%NC;

		if (Space[i][j] != NULL)
		{
			random_shuffle(Space[i][j]->Symbionts->begin(), Space[i][j]->Symbionts->end(), uniform_shuffle);

			if (log_lineage)	LogLineage(i,j);

			/* Basal death */
			if (uniform() < death_rate_host)
			{
				Space[i][j]->DeathOfCell();
				Space[i][j] = NULL;
				continue;
			}
			for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
			{
				if (uniform() < death_rate_symbiont)
				{
					Space[i][j]->DeathOfSymbiont(s);
					Space[i][j]->nr_symbionts--;
				}
			}
			if (Space[i][j]->CheckCellDeath(false))
			{
				Space[i][j] = NULL;
				continue;
			}
			/* ~Basal death */

			/* Failed division */
			if (Space[i][j]->Host->Stage == 5)
			{
				Space[i][j]->DeathOfCell();
				Space[i][j] = NULL;
				continue;
			}
			for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
			{
				if ( Space[i][j]->Symbionts->at(s)->Stage == 5)
				{
					Space[i][j]->DeathOfSymbiont(s);
					Space[i][j]->nr_symbionts--;
				}
			}
			if (Space[i][j]->CheckCellDeath(false))
			{
				Space[i][j] = NULL;
				continue;
			}
			/* ~Failed division */

			/* Division of host */
			if (Space[i][j]->Host->Stage == 4   &&   (uniform() < Space[i][j]->Host->fitness))
			{
				coords neigh = PickNeighbour(i, j);
				if (Space[neigh.first][neigh.second] != NULL)
				{
					Space[neigh.first][neigh.second]->DeathOfCell();
					Space[neigh.first][neigh.second] = NULL;
				}
				Space[neigh.first][neigh.second] = new Cell();
				Space[neigh.first][neigh.second]->Host = new Organelle();
				id_count++;
				Space[neigh.first][neigh.second]->Host->Mitosis(Space[i][j]->Host, id_count);

				for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
				{
					if (uniform() < 0.5)	//Transfer the symbiont to the new cell (just a matter of using the right pointers).
					{
						Space[neigh.first][neigh.second]->Symbionts->push_back(Space[i][j]->Symbionts->at(s));
						Space[i][j]->Symbionts->erase(Space[i][j]->Symbionts->begin()+s);
						Space[neigh.first][neigh.second]->nr_symbionts++;
						Space[i][j]->nr_symbionts--;
					}
				}
				if (!Space[neigh.first][neigh.second]->CheckCellDeath(false))	//Don't continue, since we have killed the child which lies at a different site.
				{
					Space[neigh.first][neigh.second]->DNATransferToHost();	//Instead, use the cell death check to know whether we need to do DNA transfer.
					if (Space[neigh.first][neigh.second]->Host->mutant || trace_lineage || log_lineage)
					{
						FossilSpace->BuryFossil(Space[neigh.first][neigh.second]->Host);	//Bury fossil only after potential transfer events (also count as mutations).
					}

				}
				else	Space[neigh.first][neigh.second] = NULL;
				if (Space[i][j]->CheckCellDeath(false))
				{
					Space[i][j] = NULL;
					continue;
				}
			}
			/* ~Division of host */

			/* Division of symbionts */
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if (Space[i][j]->Symbionts->at(s)->Stage == 4   &&   (uniform() < Space[i][j]->Symbionts->at(s)->fitness))
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
			if (Space[i][j]->Host->Stage == 2   &&   Space[i][j]->Host->privilige)
			{
				Space[i][j]->Host->Replicate(NutrientSpace[i][j]);
			}
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				if ( Space[i][j]->Symbionts->at(s)->Stage == 2   &&   Space[i][j]->Symbionts->at(s)->privilige)
				{
					Space[i][j]->Symbionts->at(s)->Replicate(NutrientSpace[i][j]);
				}
			}
			/* ~Replication */

		}
	}
}

bool Population::CheckLineage(int i, int j)
{
	i_lin check_id;
	check_id = lower_bound(Lineage.begin(), Lineage.end(), Space[i][j]->Host->fossil_id);
	if (*check_id == Space[i][j]->Host->fossil_id)	return true;
	else																						return false;
}

void Population::LogLineage(int i, int j)
{
	if (Lineage.size() > (size_t)0)
	{
		if (CheckLineage(i,j))
		{
			OutputLineage(i,j);
		}
	}
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

void Population::CollectNutrientsFromSite(int i, int j)
{
	int ii, jj, nrow, ncol, cell_density=0, organelle_density=0;
	double nutrient_share;

	//First obtain organelle density in 3x3 neighbourhood.
	for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
	{
		if (nutrient_competition == 1 && (ii == i && jj == j))	continue;		//Only skipped in classic nutrient function.

		if (ii < 0)					nrow = ii + NR;
		else if (ii >= NR)	nrow = ii - NR;
		else								nrow = ii;
		if (jj < 0)					ncol = jj + NC;
		else if (jj >= NC)	ncol = jj - NC;
		else								ncol = jj;

		if (Space[nrow][ncol] != NULL)
		{
			cell_density++;
			organelle_density += Space[nrow][ncol]->nr_symbionts + 1;		//Host always counts for one.
		}
	}

	if (nutrient_competition == 1)	//Finish classic nutrient function here.
	{
		if (Space[i][j] == NULL)	NutrientSpace[i][j] += nutrient_abundance - (double)organelle_density;
		else											NutrientSpace[i][j] += (nutrient_abundance - (double)organelle_density) / (double)(Space[i][j]->nr_symbionts+1);
	}

	else
	{
		if (organelle_density == 0)					nutrient_share = nutrient_abundance;
		else if (nutrient_competition == 2)	nutrient_share = nutrient_abundance / (double)organelle_density;	//Used in first smooth nutrient function.
		else																nutrient_share = nutrient_abundance / (double)cell_density;	//Used in second smooth nutrient function.

		for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
		{
			if (ii < 0)					nrow = ii + NR;
			else if (ii >= NR)	nrow = ii - NR;
			else								nrow = ii;
			if (jj < 0)					ncol = jj + NC;
			else if (jj >= NC)	ncol = jj - NC;
			else								ncol = jj;

			if (nutrient_competition == 2)		NutrientSpace[nrow][ncol] += nutrient_share;	//Used in first nutrient function.
			else	//Used in second nutrient function.
			{
				if (Space[nrow][ncol] == NULL)	NutrientSpace[nrow][ncol] += nutrient_share;
				else														NutrientSpace[nrow][ncol] += nutrient_share / (double)(Space[nrow][ncol]->nr_symbionts+1);
			}
		}
	}
}

/* ######################################################################## */
/*			INITIALISATION	INITIALISATION	INITIALISATION	INITIALISATION			*/
/* ######################################################################## */

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
		if (uniform() < 0.6)	//Initialise lower density
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

void Population::ContinuePopulationFromBackup()
{
	int i, j, s;

	ReadBackupFile();
	if(anctrace_reboot != "")	ReadAncestorFile();
	else	//Reset the fossil_ids if we are not using an existing fossil record.
	{
		id_count = 0;
		for (i=0; i<NR; i++)	for(j=0; j<NC; j++)
		{
			if (Space[i][j]!=NULL)
			{
				id_count++;
				Space[i][j]->Host->fossil_id = id_count;	//Other things such as time_of_appearance and Ancestor are set to zero by the EmptyProkaryote function.
				for (s=0; Space[i][j]->nr_symbionts; s++)
				{
					id_count++;
					Space[i][j]->Symbionts->at(s)->fossil_id = id_count;
				}
			}
		}
	}

	//Only do these things if we're starting at an irregular timepoint (where we wouldn't already store info inside UpdatePop.).
	if (TimeZero%TimeTerminalOutput != 0)										ShowGeneralProgress();
	if (TimeZero%TimeSaveGrid != 0)													OutputGrid(false);
	if (TimeZero%TimePruneFossils != 0 && !trace_lineage)		PruneFossilRecord();
	if (TimeZero%TimeOutputFossils != 0 && !trace_lineage)	FossilSpace->ExhibitFossils();
	if (TimeZero%TimeSaveBackup != 0)												OutputGrid(true);
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_reboot.c_str());
	string line, data;
	char* data_element;
	string::iterator sit;
	Genome::i_bead it;
	int r, c, s, begin_data, end_data, success, stage, pfork, pterm, nr=NR, nc=NC, init_seed, read_header = 0, count_lines = 0;	//For now, set nr and nc for backup-file to NR and NC parameters (i.e. if backup-file does not contain header; for old backups).
	unsigned long long org_id, anc_id, sdraws, iull;
	char temp_is_mutant[20], temp_priv[20];
	double fit;
	size_t pos;
	Organelle* O, * OF, * SaveO;
	Cell* C;
	i_fos iF;

	if (!infile.is_open())
	{
		printf("Backup-file could not be opened.\n");
		exit(1);
	}

	printf("Reading backup from file: %s\n", backup_reboot.c_str());
	while(getline(infile,line))
	{

		//Read "Header" of backup file.
		if (line=="### Header ###")
		{
			read_header = 1;
			continue;
		}
		else if (line=="#### Main ####")
		{
			read_header = -1;
			continue;
		}
		else if (read_header==1)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "NR:%d\tNC:%d", &nr, &nc);
			if(success != 2)
			{
				cerr << "Could not read NR and NC from backup-file.\n" << endl;
				exit(1);
			}
			else	cout << "nr=" << nr << ", nc=" << nc << endl;
			read_header = 2;
			continue;
		}
		else if (read_header==2)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "Initial seed:%d\tSeed draws:%llu", &init_seed, &sdraws);
			if (success != 2)
			{
				cerr << "Could not read seed status from backup-file.\n" << endl;
				exit(1);
			}
			else
			{
				dsfmt_init_gen_rand(&dsfmt, init_seed);	//Used to seed uniform().
				srand(init_seed);	//Used to seed random_shuffle(...).
				cout << "Simulating " << sdraws << " random draws, initial seed=" << init_seed << endl;
				for(iull=0; iull<sdraws; iull++){
					uniform();
				}
			}
			read_header = 3;	//Not doing anything with this yet.
			continue;
		}

		//Start reading "Main" data from backup file.
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);	//Print some progress.


		pos = line.find("NULL");
		if (pos != string::npos)	//Empty site.
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "%d %d NULL\n", &r, &c);
			if (success != 2)
			{
				cerr << "Could not read row and column from empty site. Backup file potentially corrupt.\n" << endl;
			}
			Space[r][c] = NULL;
		}
		else											//Non-empty site.
		{
			//Make the new organelle, we'll decide at the end (?) whether it is a host and requires a cell to be made or whether it is a symbiont.
			O = new Organelle();

			//Read BeadList. Some variables are here initiated, but will be overwritten by information from the backup file below.
			begin_data = line.find_first_of("(");
			end_data = line.find_last_of(")");
			data = line.substr(begin_data, end_data-begin_data+1);
			O->G->ReadGenome(data);

			begin_data = line.find_first_of("{");
			end_data = line.find_first_of("}");
			data = line.substr(begin_data, end_data-begin_data+1);	//Brackets are stripped off in ReadExpression function.
			O->G->ReadExpression(data);

			//Read organelle data.
			begin_data = line.find_last_of("[");
			end_data = line.find_last_of("]");
			data = line.substr(begin_data+1, end_data-begin_data-1);
			data_element = strtok((char*)data.c_str(),"\t");
			while(data_element != NULL)
			{
				success = sscanf(data_element, "%llu %llu %s %lf %d %s %d %d", &org_id, &anc_id, temp_is_mutant, &fit, &stage, temp_priv, &pfork, &pterm);
				if(success != 8)
				{
					cerr << "Could not find sufficient information for this prokaryote. Backup file potentially corrupt.\n" << endl;
					exit(1);
				}

				data_element = strtok(NULL, "\t");
				O->Stage = stage;
				O->fitness = fit;
				O->G->fork_position = pfork;
				O->G->terminus_position = pterm;
				O->fossil_id = org_id;
				O->mutant = (strcmp( temp_is_mutant, "Y") == 0);
				O->privilige = (strcmp( temp_priv, "Y") == 0);
				O->alive = true;
				if (org_id > id_count) id_count = org_id;

				SaveO = O;	//SaveO takes the pointer value of O.
				O = FindInFossilRecord(SaveO->fossil_id);	//Check the current organelle was already put in the fossil record; O is reassigned as pointer.
				if (O == NULL)	//It is not yet in the fossil record. Act like nothing happened.
				{
					O = SaveO;	//Continue with the SaveO pointer.
					SaveO = NULL;
				}
				else	//It is in the fossil record.
				{
					O->CloneOrganelle(SaveO, SaveO->fossil_id);	//Continue with the previously made pointer that was already in the fossil record, but give it all the data we've read in so far.
					delete SaveO;
					SaveO = NULL;
					break;	//Although it is a mutant, we have already buried its fossil, so no need to take the otherwise standard action for mutant (below).
				}

				if (O->mutant)
				{
					//If this live organelle is a mutant, we will store the organelle itself in the fossil record.
					FossilSpace->BuryFossil(O);
				}
				else
				{
					//If this live organelle is NOT a mutant, we will check if its direct ancestor is already known (in which case we can simply point to that); otherwise, we construct the ancestor only with its ID, completing its construction in ReadAncestorFile.
					O->Ancestor = FindInFossilRecord(anc_id);
					if (O->Ancestor == NULL)	//i.e. Ancestor not currently present in fossil record.
					{
						OF = new Organelle();	//Make the ancestor with the right name...
						OF->fossil_id = anc_id;
						FossilSpace->BuryFossil(OF);
						O->Ancestor = OF;	//...and point to it, so that our current organelle (O) now has a defined ancestor.
					}
				}
			}

			data_element = (char*)line.c_str();
			success = sscanf(data_element, "%d %d %d\t", &r, &c, &s);
			if (success != 3)
			{
				cerr << "Could not read r, c, s from non-empty site. Backup file potentially corrupt.\n" << endl;
			}
			if (s == -1)
			{
				C = new Cell();
				C->Host = O;
				Space[r][c] = C;
			}
			else
			{
				Space[r][c]->Symbionts->push_back(O);
				Space[r][c]->nr_symbionts++;
			}

		}
		count_lines++;
	}

	if ((r != nr-1) || (c != nc-1))
	{
		cerr << "Backup file incomplete or too large\nRows: " << r << " of " << nr << "\nCols: " << c << " of " << nc << endl;
		exit(1);
	}
}

void Population::ReadAncestorFile()
{
	ifstream infile(anctrace_reboot.c_str());
	string line, data;
	int begin_data, end_data, TimeOA, count_alive = 0, count_fossils = 0, count_lines = 0;
	unsigned long long ID, AncID;
	i_fos iF;
	Organelle* O;

	if (!infile.is_open())
	{
		printf("Ancestor-file could not be opened.\n");
		exit(1);
	}

	printf("Reading ancestors from file: %s\n", anctrace_reboot.c_str());
	while(getline(infile,line))
	{
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);
		end_data = line.find("\t");
		data = line.substr(0,end_data);
		stringstream(data) >> ID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> AncID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> TimeOA;

		begin_data = end_data;
		end_data = line.size();
		data = line.substr(begin_data+1, end_data-begin_data);

		// 3 OPTIONS:
		// 1) This line corresponds to an organelle that is still alive, we only need to know the time_of_appearance and where to have the Ancestor point to.
		// 2) This line corresponds to an organelle that is not alive, but which was created only with a fossil_id because it is the ancestor of an extant (non-mutant) organelle; and because we wanted to store the Ancestor of that extant, non-mutant organelle, we pre-made the ancestor, but we can now give it a genome and few other variables.
		// 3) This line corresponds to a dead ancestor (mutant) that is not the direct ancestor of any of the living organelles, and so we will create it here (similar to case 2, but now it's not yet created).


		O = FindInFossilRecord(ID);
		if (O != NULL)	//This line was already in the fossil record.
		{
			if (O->alive)	count_alive++;	//Live mutants; most info read from backup.
			else	//Dead mutants, must have been put in the fossil record because it is the direct ancestor of a live and non-mutant individual
			{
				count_fossils++;
				O->mutant = true;
				O->G->ReadGenome(data);
				O->G->terminus_position = O->G->g_length;	//Relevant for printing.
			}
			O->time_of_appearance = TimeOA;

			if (AncID == 0)	O->Ancestor = NULL;
			else
			{
				O->Ancestor = FindInFossilRecord(AncID);
				if (O->Ancestor == NULL)
				{
					cerr << "Error: requested ancestor not available in fossil record.\n" << endl;
					exit(1);
				}
			}
		}
		else	//We did not break out of the loop, so we have apparently not encountered this ID among the current list of fossils.
		{
			count_fossils++;
			O = new Organelle();
			O->fossil_id = ID;
			O->mutant = true;
			O->time_of_appearance = TimeOA;
			if(AncID == 0)	O->Ancestor = NULL;
			else
			{
				O->Ancestor = FindInFossilRecord(AncID);
				if (O->Ancestor == NULL)
				{
					cerr << "Error: requested ancestor not available in fossil record.\n" << endl;
					exit(1);
				}
			}

			//Set up its ghost genome.	It only has beads.
			O->G->ReadGenome(data);
			O->G->terminus_position = O->G->g_length;	//Relevant for printing.
			FossilSpace->BuryFossil(O);
		}
		count_lines++;
	}
	assert (count_lines == count_alive+count_fossils);
	assert (FossilSpace->FossilRecord.size() == (size_t)count_lines);
}

Organelle* Population::FindInFossilRecord(unsigned long long AncID)
{
	i_fos iF;

	iF = FossilSpace->FossilRecord.begin();
	while(iF != FossilSpace->FossilRecord.end())
	{
		if ((*iF)->fossil_id == AncID)
		{
			return *iF;
		}
		iF++;
	}

	return NULL;
}

void Population::ReadLineageFile()
{
	ifstream infile(lineage_record.c_str());
	string line, data;
	int end_data, count_lines = 0;
	unsigned long long ID;

	if (!infile.is_open())
	{
		printf("Lineage-file could not be opened.\n");
		exit(1);
	}

	printf("Reading lineage from file: %s\n", lineage_record.c_str());
	while(getline(infile,line))
	{
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);
		end_data = line.find("\t");
		data = line.substr(0,end_data);

		stringstream(data) >> ID;
		Lineage.push_back(ID);
		count_lines++;
	}
	assert (Lineage.size() == (size_t)count_lines);
	sort(Lineage.begin(), Lineage.end());
}

/* ######################################################################## */
/*				FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS						*/
/* ######################################################################## */

void Population::PruneFossilRecord()
{
	std::list<unsigned long long> AllFossilIDs;
	int i, j, s, fossil_record_size;
	unsigned long long fossilID;
	i_fos iF;
	i_ull findit;
	Organelle* lastCA;
	i_lin iL;

	if (trace_lineage && Time == SimTime)
	{
		//Only trace hosts that are in the lineage record.

		iL = Lineage.begin();
		while (iL != Lineage.end() && *iL <= id_count)
		{
			AllFossilIDs.push_back(*iL);

			//Find location of this individual in the field.
			lastCA = NULL;
			for(i=0; i<NR; i++)	for(j=0; j<NC; j++)
			{
				if (Space[i][j] != NULL)
				{
					if (Space[i][j]->Host->fossil_id == *iL)
					{
						lastCA = Space[i][j]->Host->Ancestor;
						break;
					}
				}
			}

			//Store entire lineage of this individual.
			while(lastCA != NULL)
			{
				AllFossilIDs.push_back(lastCA->fossil_id);
				lastCA = lastCA->Ancestor;
			}
			iL++;
		}
	}

	else
	{
		//Trace all living individuals back to the beginning, storing all individuals along their lineage.
		for(i=0; i<NR; i++)	for(j=0; j<NC; j++)
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
	}
	// Delete duplicates.
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
		if(findit==AllFossilIDs.end() && (!(*iF)->alive || (trace_lineage && Time == SimTime)) )	//Fossil not needed/found in relevant branches AND either the organelle is already dead OR we're tracing a lineage (in which case we want to clean our tree more quickly, live individuals will not suddenly produce offspring that are in the lineage record).
		{
			if (!(*iF)->alive)	delete *iF;	//Only delete things that are not alive anymore.
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
/*				BACKUP	BACKUP	BACKUP	BACKUP	BACKUP	BACKUP	BACKUP						*/
/* ######################################################################## */

void Population::OutputGrid(bool backup)
{
	FILE* f;
	char OutputFile[800];
	int i, j, s;

	if (backup)
	{
		sprintf(OutputFile, "%s/backups/backup%08d.txt", folder.c_str(), Time);
		f=fopen(OutputFile, "w");
		if (f == NULL)	printf("Failed to open file for writing the backup.\n");

		//Print NR and NC to file.
		fprintf(f, "### Header ###\nNR:%d\tNC:%d\nInitial seed:%d\tSeed draws:%llu\n#### Main ####\n", NR, NC, initial_seed, seed_draws);
	}
	else
	{
		sprintf(OutputFile, "%s/snapsamples/field%08d.txt", folder.c_str(), Time);
		f=fopen(OutputFile, "w");
		if (f == NULL)	printf("Failed to open file for writing the snapshot.\n");
	}

	for (i=0; i<NR; i++) for(j=0; j<NC; j++) {
		if(Space[i][j] == NULL)
		{
		 	fprintf(f, "%d %d NULL %f\n", i, j, NutrientSpace[i][j]);
		}
		else	//Print internal state and genome of prokaryote to file.
		{
			fprintf(f, "%d %d -1 %f\t%s\n", i, j, NutrientSpace[i][j], Space[i][j]->Host->Output(backup).c_str());
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				fprintf(f, "%d %d %d %f\t%s\n", i, j, s, NutrientSpace[i][j], Space[i][j]->Symbionts->at(s)->Output(backup).c_str());
			}
		}
	}

	fclose(f);
}

void Population::OutputLineage(int i, int j)
{
	FILE* f;
	char OutputFile[800];
	int s;

	sprintf(OutputFile, "%s/lineage.out", folder.c_str());
	f=fopen(OutputFile, "a");
	if (f == NULL)	printf("Failed to open file for writing the backup.\n");

	fprintf(f, "%d %d %d -1 %f\t%s\n", Time, i, j, NutrientSpace[i][j], Space[i][j]->Host->Output(true).c_str());
	for (s=0; s<Space[i][j]->nr_symbionts; s++)
	{
		fprintf(f, "%d %d %d %d %f\t%s\n", Time, i, j, s, NutrientSpace[i][j], Space[i][j]->Symbionts->at(s)->Output(true).c_str());
	}

	fclose(f);
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

	if (host_count==0)
	{
		cout << "And since there is no more life, we will stop the experiment here.\n" << endl;
		exit(1);
	}

}
