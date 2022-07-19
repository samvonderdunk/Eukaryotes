#include "Population.hh"

Population::Population()
{
	int i,j;
	id_count = 0;
	nr_strains = 0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		Space[i][j]=NULL;
		NutrientSpace[i][j]=0.;
	}
	FossilSpace = new Fossils();
}

Population::~Population()	//Skips deconstructor of Cell, because we need to check whether organelles are also in the FossilSpace.
{
	int i,j;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if ((Space[i][j]) != NULL)
		{
			if (Space[i][j]->Vesicle->mutant) FossilSpace->EraseFossil(Space[i][j]->Vesicle->fossil_id);
			delete Space[i][j];
			Space[i][j] = NULL;
		}
	}

	delete FossilSpace;
	FossilSpace = NULL;
}

void Population::UpdatePopulation()
{
	int update_order[NR*NC];
	int u, i, j;
	Cell* C;

	if (well_mixing)	WellMix();	//Well-mixing (before determining nutrient levels).

	for(i=0; i<NR; i++)	for(j=0; j<NC; j++)	NutrientSpace[i][j] = 0.;
	for(i=0; i<NR; i++)	for(j=0; j<NC; j++)	CollectNutrientsFromSite(i,j);

	if(Time%TimeTerminalOutput==0)																						ShowGeneralProgress();
	if(Time%TimeSaveGrid==0)																									OutputGrid(false);
	if(Time%TimePruneFossils==0 && Time!=0)																		PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)																	FossilSpace->ExhibitFossils();
	if(Time%TimeSaveBackup==0 && Time!=0)
	{
		if (cell_fitness)
		{
			for(i=0; i<NR; i++) for(j=0; j<NC; j++)	if (Space[i][j]!=NULL)	Space[i][j]->CalculateCellFitness();
		}
		OutputGrid(true);
	}

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC], uniform_shuffle);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;
		j = update_order[u]%NC;

		if (Space[i][j] != NULL)
		{
			/* Basal death */
			if (Space[i][j]->BasalDeath())
			{
				delete Space[i][j];		//Note: collect these things in a function?
				Space[i][j] = NULL;
				continue;
			}

			/* Failed division */
			if (Space[i][j]->FailedDivision())
			{
				delete Space[i][j];
				Space[i][j] = NULL;
				continue;
			}

			/* Successful host division */
			coords neigh = PickNeighbour(i, j);
			C = Space[i][j]->Division(&Space[neigh.first][neigh.second], FossilSpace, &id_count);
			if (C != NULL)	//Division returns the newborn eukaryote, and we are to kill the current cell (which cannot be done from inside the cell class).
			{
				delete Space[i][j];
				if (!empty_division_killing) //If we don't want to kill a neighbour, swap parent and daugther cell. Don't kill a neibhour; we're done here.
				{
					Space[i][j] = C;
					C = NULL;
				}
				else	//We do want to kill a neighbour, so we're not going to swap; newcell will be put in the neighbour spot, parent simply dies.
				{
					Space[i][j] = NULL;
					if (Space[neigh.first][neigh.second] != NULL)
					{
						Space[neigh.first][neigh.second]->DeathOfCell();
						delete Space[neigh.first][neigh.second];
						Space[neigh.first][neigh.second] = NULL;
					}
					Space[neigh.first][neigh.second] = C;
					C = NULL;
					continue;	//We end up with an empty parent cell, so we cannot do any regulation updating, etc.
				}
			}

			/* Update expression */
			Space[i][j]->UpdateOrganelles();	//Expression dynamics within cell.

			/* Replication */
			Space[i][j]->Replication(NutrientSpace[i][j]);

		}
	}
}

void Population::WellMix()
{
	int Space1D[NR*NC];
	int u, i, j;
	Cell* SpaceMirror[NR][NC];

	for(u=0; u<NR*NC; u++)	Space1D[u]=u;
	random_shuffle(&Space1D[0], &Space1D[NR*NC], uniform_shuffle);
	for(u=0; u<NR*NC; u++)	SpaceMirror[u/NC][u%NC] = Space[Space1D[u]/NC][Space1D[u]%NC];
	for(i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		Space[i][j] = SpaceMirror[i][j];
		SpaceMirror[i][j] = NULL;	//Should not be really necessary.
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
	double nut_condition, nutrient_share, starting_nuts, claim_density=0., orgs_at_site;

	nut_condition = nutrient_condition[j/(NC/nr_sectors)];

	//First obtain organelle and cell density in 3x3 neighbourhood.
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
			organelle_density++;
			claim_density += Space[nrow][ncol]->Vesicle->nutrient_claim;
		}
	}

	if (Space[i][j] == NULL)										orgs_at_site = 0.;
	else																				orgs_at_site = 1.;

	//If we don't have wrapped columns, border pixels get fewer nutrients to start with.
	starting_nuts = nut_condition;

	if (nutrient_competition == 0)	//Constant nutrient level.
	{
		NutrientSpace[i][j] += nut_condition;
	}

	else if (nutrient_competition == 1)	//Finish classic nutrient function here.
	{
		if (Space[i][j] == NULL)	NutrientSpace[i][j] += starting_nuts - (double)organelle_density;
		else											NutrientSpace[i][j] += (starting_nuts - (double)organelle_density) / orgs_at_site;
	}

	else
	{
		if (organelle_density == 0)					nutrient_share = starting_nuts;
		else if (nutrient_competition == 2)	nutrient_share = starting_nuts / (double)organelle_density;	//Used in first smooth nutrient function.
		else if (nutrient_competition == 6)	nutrient_share = starting_nuts / claim_density;
		else																nutrient_share = starting_nuts / (double)cell_density;	//Used in second smooth nutrient function and in third nutshare_evolve protocol.

		//Here we actually add the nutrient share to each site.
		for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
		{
			if (ii < 0)					nrow = ii + NR;
			else if (ii >= NR)	nrow = ii - NR;
			else								nrow = ii;
			if (jj < 0)					ncol = jj + NC;
			else if (jj >= NC)	ncol = jj - NC;
			else								ncol = jj;

			if (nutrient_competition == 2 || nutrient_competition == 4 || nutrient_competition == 6)		NutrientSpace[nrow][ncol] += nutrient_share;	//Used in first smooth nutrient function. But also when we evolve nutrient_sharing, the whole site gets the unshared total nutrients. These will be further depleted upon replication.
			else	//Used in second nutrient function.
			{
				if (Space[nrow][ncol] == NULL)	NutrientSpace[nrow][ncol] += nutrient_share;
				else
				{
					NutrientSpace[nrow][ncol] += nutrient_share;
				}
			}
		}
	}
}


/* ######################################################################## */
/*			INITIALISATION	INITIALISATION	INITIALISATION	INITIALISATION			*/
/* ######################################################################## */

void Population::InitialisePopulation()
{
	Cell* InitCells[max_input_files] = {NULL};
	int k, i, j, blocks, r, c;
	double b, b_dbl, b_int;
	Cell* C;

	//Each input file (genome & expression) will be made into a cell.
	for (k=0; k<max_input_files; k++)
	{
		if (genome_files[k] == "")
		{
			nr_strains = k;
			break;
		}

		InitCells[k] = new Cell();
		InitCells[k]->InitialiseCell(k);
	}

	if (strain_competition == 3)	//Determine how many squares will divide the grid (hopefully 2x2, 3x3 or 4x4 and no more).
	{
		if (NR != NC)
		{
			printf("Warning: Might be difficult to initialise blocks of strains on a non-square field.\n");
		}

		b = sqrt((double)nr_strains);
		b_dbl = modf(b, &b_int);	//Split double into its integer and fractional parts.
		blocks = (int) b_int;
		if (b_dbl > 0.)	blocks++;
	}

	for (i=0; i<NR; i++) for(j=0; j<NC; j++)
	{

		//Determine which strain we will put here.
		if (strain_competition == 1)
		{
			k = j / (NC/nr_strains);
		}
		else if (strain_competition == 2)
		{
			k = (int)(uniform()*nr_strains);
		}
		else if (strain_competition == 3)
		{
			r = i / (NR/blocks);
			c = j / (NC/blocks);

			if (r < blocks && c < blocks)	//Otherwise leave empty, outside square blocks.
			{
				k = r*blocks + c;
			}
		}

		if (k >= nr_strains)	continue;	//Leave empty, if NC is not divisible by nr_strains (i.e. when you don't get an integer).

		//Take InitCell[k] and copy to Space[i][j].
		//Still we might not want to fill the entire grid, so roll another die.

		if (uniform() < 0.6)
		{
			C = new Cell();
			C->CloneCell(InitCells[k], &id_count);
			FossilSpace->BuryFossil(C->Vesicle);
			Space[i][j] = C;
		}
	}

	for (k=0; k<max_input_files; k++)
	{
		if (InitCells[k] != NULL)	delete InitCells[k];
		InitCells[k] = NULL;
	}
}

void Population::ContinuePopulationFromBackup()
{
	int i, j;

	ReadBackupFile();
	if(anctrace_file != "")	ReadAncestorFile();
	else	//Reset the fossil_ids if we are not using an existing fossil record.
	{
		id_count = 0;
		for (i=0; i<NR; i++)	for(j=0; j<NC; j++)
		{
			if (Space[i][j]!=NULL)
			{
				id_count++;
				Space[i][j]->Vesicle->fossil_id = id_count;
			}
		}
	}

	//Only do these things if we're starting at an irregular timepoint (where we wouldn't already store info inside UpdatePop.).
	if (TimeZero%TimeTerminalOutput != 0)										ShowGeneralProgress();
	if (TimeZero%TimeSaveGrid != 0)													OutputGrid(false);
	if (TimeZero%TimePruneFossils != 0 && mutations_on)			PruneFossilRecord();
	if (TimeZero%TimeOutputFossils != 0 && mutations_on)		FossilSpace->ExhibitFossils();
	if (TimeZero%TimeSaveBackup != 0)												OutputGrid(true);
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_file.c_str());
	string line, data;
	char* data_element, *number;
	string::iterator sit;
	Genome::i_bead it;
	int r, c, s, begin_data, end_data, success, stage, pfork, pterm, nr, nc, count_lines = 0, nn = 0, idx_primary, idx_secondary, it_cntr;
	unsigned long long org_id, anc_id;
	char temp_is_mutant[20], temp_priv[20];
	double fit, nutcl;
	size_t pos;
	Organelle* O, * OF, * SaveO;
	Cell* C;
	i_fos iF;

	if (!infile.is_open())
	{
		printf("Backup-file could not be opened.\n");
		exit(1);
	}

	printf("Reading backup from file: %s\n", backup_file.c_str());
	while(getline(infile,line))
	{

		if (line.substr(0,1)=="#")	continue;

		else if (line.substr(0,2)=="NR")
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "NR:%d\tNC:%d", &nr, &nc);
			if(success != 2)
			{
				cerr << "Could not read NR and NC from backup-file.\n" << endl;
				exit(1);
			}
			else if(nr != NR || nc != NC)
			{
				cerr << "NR and NC of backup did not match parameter settings.\n" << endl;
				exit(1);
			}
		}

		else if(line.substr(0,9)=="dsfmt.idx")
		{
			data_element = (char*)line.c_str();

			success = sscanf(data_element, "dsfmt.idx:%d", &dsfmt.idx);
			if (success != 1)
			{
				cerr << "Could not read dsfmt.idx from backup-file.\n" << endl;
				exit(1);
			}
		}

		else if(line.substr(0,12)=="dsfmt.status")
		{
			line = line.substr(13);

			number = strtok((char*)line.c_str(),",");
			while (number != NULL)
			{
				idx_primary = nn/2;
				idx_secondary = nn%2;

				success = sscanf(number, "%" PRIu64 , &dsfmt.status[idx_primary].u[idx_secondary]);//&dsfmt_element);
				if (success != 1)
				{
					cerr << "Could not read dsfmt element dsfmt.status[" << idx_primary << "].u[" << idx_secondary << "] from backup-file.\n" << endl;
					exit(1);
				}

				number = strtok(NULL, ",");
				nn ++;
			}
		}

		else			//Start reading "Main" data from backup file.
		{
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
				end_data = line.find("\t", begin_data+1);
				data = line.substr(begin_data, end_data-begin_data);
				O->G->ReadDefinition(data);

				begin_data = end_data+1;
				end_data = line.find_last_of(")");
				data = line.substr(begin_data, end_data-begin_data+1);
				O->G->ReadGenome(data);

				begin_data = line.find_first_of("{");
				end_data = line.find_first_of("}");
				data = line.substr(begin_data, end_data-begin_data+1);	//Brackets are stripped off in ReadExpression function.
				O->G->ReadExpression(data);
				O->G->NativeExpression();	//Only native expression has to be set; once we get to UpdateOrganelles, we will set the nr_native_expressed, do gene transport, and update the cell-cycle stage.

				//Read organelle data.
				begin_data = line.find_last_of("[");
				end_data = line.find_last_of("]");
				data = line.substr(begin_data+1, end_data-begin_data-1);
				data_element = strtok((char*)data.c_str(),"\t");
				while(data_element != NULL)
				{
					success = sscanf(data_element, "%llu %llu %s %lf %lf %d %s %d %d", &org_id, &anc_id, temp_is_mutant, &fit, &nutcl, &stage, temp_priv, &pfork, &pterm);
					if(success != 9)
					{
						cerr << "Could not find sufficient information for this prokaryote. Backup file potentially corrupt.\n" << endl;
						exit(1);
					}

					data_element = strtok(NULL, "\t");
					O->Stage = stage;
					O->fitness = fit;
					O->nutrient_claim = nutcl;
					O->G->fork_position = pfork;
					O->G->terminus_position = pterm;
					O->fossil_id = org_id;
					O->mutant = (strcmp( temp_is_mutant, "Y") == 0);
					O->privilige = (strcmp( temp_priv, "Y") == 0);
					O->alive = true;
					if (org_id > id_count) id_count = org_id;

					//nr_houses is not stored in the backup, so we calculate it. Note that this might be a burden for workflows that use backups (although runtime is probably still the bigger burden).
					it = O->G->BeadList->begin();
					it_cntr = 0;
					while (it != O->G->BeadList->end())
					{
						if ((*it)->kind == HOUSE)	O->nr_houses++;
						it++;
						it_cntr++;
						if (it_cntr == O->G->terminus_position)	break;
					}

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
				C = new Cell();
				C->Vesicle = O;
				Space[r][c] = C;
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
	ifstream infile(anctrace_file.c_str());
	string line, data, data2;
	int begin_data, end_data, TimeOA, count_alive = 0, count_fossils = 0, count_lines = 0;
	unsigned long long ID, AncID;
	double NutCL;
	i_fos iF;
	Organelle* O;

	if (!infile.is_open())
	{
		printf("Ancestor-file could not be opened.\n");
		exit(1);
	}

	printf("Reading ancestors from file: %s\n", anctrace_file.c_str());
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
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> NutCL;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data+1, end_data-begin_data-1);

		begin_data = end_data;
		end_data = line.size();
		data2 = line.substr(begin_data+1, end_data-begin_data);

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
				O->G->ReadDefinition(data);
				O->G->ReadGenome(data2);
				O->G->terminus_position = O->G->g_length;	//Relevant for printing.
			}
			O->time_of_appearance = TimeOA;	//Only this was not in the backup.

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
			O->nutrient_claim = NutCL;
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

			//Set up its ghost genome.	It only has beads and type definitions.
			O->G->ReadDefinition(data);
			O->G->ReadGenome(data2);
			O->G->terminus_position = O->G->g_length;	//Relevant for printing.
			FossilSpace->BuryFossil(O);
		}
		count_lines++;
	}

	FossilSpace->SortFossils();

	assert (count_lines == count_alive+count_fossils);
	FossilSpace->ExhibitFossils();
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

/* ######################################################################## */
/*				FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS						*/
/* ######################################################################## */

void Population::PruneFossilRecord()
{
	std::list<unsigned long long> AllFossilIDs;
	int i, j, fossil_record_size;
	unsigned long long fossilID;
	i_fos iF;
	i_ull findit;
	Organelle* lastCA;
	i_lin iL;

	//Trace all living individuals back to the beginning, storing all individuals along their lineage.
	for(i=0; i<NR; i++)	for(j=0; j<NC; j++)
	{
		if(Space[i][j] != NULL)
		{
			lastCA = Space[i][j]->Vesicle->Ancestor;
			while(lastCA != NULL)
			{
				AllFossilIDs.push_back(lastCA->fossil_id);
				lastCA = lastCA->Ancestor;
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
		if(findit==AllFossilIDs.end() && (!(*iF)->alive))	//Fossil not needed/found in relevant branches AND either the organelle is already dead OR we're tracing a lineage (in which case we want to clean our tree more quickly, live individuals will not suddenly produce offspring that are in the lineage record).
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
	int i, j;

	if (backup)
	{
		sprintf(OutputFile, "%s/backups/backup%08d.txt", folder.c_str(), Time);
		f=fopen(OutputFile, "w");
		if (f == NULL)	printf("Failed to open file for writing the backup.\n");

		fprintf(f, "### Header ###\nNR:%d\tNC:%d\n", NR, NC);
		fprintf(f, "dsfmt.idx:%d\ndsfmt.status:", dsfmt.idx);
		for (i=0; i<DSFMT_N; i++)
		{
			fprintf(f, "%" PRIu64 ",%" PRIu64 ",", dsfmt.status[i].u[0], dsfmt.status[i].u[1]);
		}
		fprintf(f, "%" PRIu64 ",%" PRIu64 "\n### Main ###\n", dsfmt.status[DSFMT_N].u[0], dsfmt.status[DSFMT_N].u[1]);
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
			fprintf(f, "%d %d -2 %f\t%s\n", i, j, NutrientSpace[i][j], Space[i][j]->Vesicle->Output(backup).c_str());
		}
	}

	fclose(f);
}

/* ######################################################################## */
/*				OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT	OUTPUT						*/
/* ######################################################################## */

void Population::ShowGeneralProgress()
{
	int i, j;
	int cell_count[nr_sectors] = {0};
	int stage_counts[3][6] = {0};
	i_org iS;

	for (i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		if ( Space[i][j] != NULL )
		{
			cell_count[j/(NC/nr_sectors)]++;	// sector = j/(NC/nr_sectors)
			stage_counts[0][Space[i][j]->Vesicle->Stage]++;
		}
	}

	cout << "Time " << Time << "\tCells ";
	for (i=0; i<nr_sectors; i++)
	{
		cout << cell_count[i];
		if (i != nr_sectors-1)	cout << "|";
	}
	cout << endl;

	if (cell_count[0]+cell_count[1]+cell_count[2]+cell_count[3]+cell_count[4]+cell_count[5] == 0)	//Easier to use symbiont_count as we did not split this over sectors.
	{
		cout << "And since there is no more life, we will stop the experiment here.\n" << endl;
		exit(1);
	}

}
