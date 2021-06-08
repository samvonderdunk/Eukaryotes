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

	if(Time%TimeTerminalOutput==0)						ShowGeneralProgress();
	// if(Time%TimeSaveGrid==0)	PrintFieldToFile();
	if(Time%TimePruneFossils==0 && Time!=0)		PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)	FossilSpace->ExhibitFossils();
	if(Time%TimeSaveBackup==0 && Time!=0)			OutputBackup();

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC], uniform_shuffle);	//Is also set by initial_seed through srand(initial_seed), see World.cc

	for(u=0; u<NR*NC; u++)
	{
		i = update_order[u]/NC;
		j = update_order[u]%NC;

		if (Space[i][j] != NULL)
		{
			random_shuffle(Space[i][j]->Symbionts->begin(), Space[i][j]->Symbionts->end(), uniform_shuffle);

			/* Basal death */
			if (uniform() < death_rate_host)
			{
				DeathOfCell(i, j);
				continue;
			}
			for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
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
			for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
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
				if (!CheckCellDeath(neigh.first, neigh.second))	//Don't continue, since we have killed the child which lies at a different site.
				{
					Space[neigh.first][neigh.second]->DNATransferToHost();	//Instead, use the cell death check to know whether we need to do DNA transfer.
					if (Space[neigh.first][neigh.second]->Host->mutant)
					{
						FossilSpace->BuryFossil(Space[neigh.first][neigh.second]->Host);	//Bury fossil only after potential transfer events (also count as mutations).
					}
				}
				if (CheckCellDeath(i, j))	continue;
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
	for (s=Space[i][j]->nr_symbionts-1; s>=0; s--)
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
				id_count ++;
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
	if (TimeZero%TimePruneFossils != 0)		PruneFossilRecord();
	if (TimeZero%TimeOutputFossils != 0)	FossilSpace->ExhibitFossils();
	if (TimeZero%TimeSaveBackup != 0)			OutputBackup();
	if (TimeZero%TimeTerminalOutput != 0)	ShowGeneralProgress();
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_reboot.c_str());
	string line, data;
	char* data_element;
	string::iterator sit;
	Genome::i_bead it;
	int i, r, c, s, begin_data, end_data, success, stage, pfork, pterm, temp_is_mutant, temp_priv, nr=NR, nc=NC, init_seed, read_header = 0, count_lines = 0;	//For now, set nr and nc for backup-file to NR and NC parameters (i.e. if backup-file does not contain header; for old backups).
	unsigned long long org_id, anc_id, sdraws;
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
				cout << "Simulating " << sdraws << " random draws, initial seed=" << init_seed << endl;
				for(i=0; (unsigned long long)i<sdraws; i++){
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
				success = sscanf(data_element, "%d %lf %d %d %llu %llu %d %d", &stage, &fit, &pfork, &pterm, &org_id, &anc_id, &temp_is_mutant, &temp_priv);
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
				O->mutant = (temp_is_mutant==1);
				O->privilige = (temp_priv==1);
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
			lastCA = Space[i][j]->Host->Ancestor; //HERE?
			while(lastCA != NULL)
			{
				AllFossilIDs.push_back(lastCA->fossil_id);
				lastCA = lastCA->Ancestor;
			}

			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				lastCA = Space[i][j]->Symbionts->at(s)->Ancestor; //HERE?
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
/*				BACKUP	BACKUP	BACKUP	BACKUP	BACKUP	BACKUP	BACKUP						*/
/* ######################################################################## */

void Population::OutputBackup()
{
	FILE* f;
	char OutputFile[800];
	int i, j, s;

	sprintf(OutputFile, "%s/backups/backup%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL)	printf("Failed to open file for writing the backup.\n");

	//Print NR and NC to file.
	fprintf(f, "### Header ###\nNR:%d\tNC:%d\nInitial seed:%d\tSeed draws:%llu\n#### Main ####\n", NR, NC, initial_seed, seed_draws);

	for (i=0; i<NR; i++) for(j=0; j<NC; j++) {
		if(Space[i][j]==NULL){
		 	fprintf(f, "%d %d NULL\n", i, j);
		}
		else	//Print internal state and genome of prokaryote to file.
		{
			fprintf(f, "%d %d %d\t%s\n", i, j, -1, Space[i][j]->Host->OutputBackup().c_str());
			for (s=0; s<Space[i][j]->nr_symbionts; s++)
			{
				fprintf(f, "%d %d %d\t%s\n", i, j, s, Space[i][j]->Symbionts->at(s)->OutputBackup().c_str());
			}
		}
	}

	fclose(f);
	return;
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
