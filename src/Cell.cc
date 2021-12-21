#include "Cell.hh"

Cell::Cell()
{
	nr_symbionts = 0;
	barcode = -1;	//Not initialised here.
	Host = NULL;
	Host = new Organelle();
	Symbionts = NULL;
	Symbionts = new vector<Organelle*>();
}

Cell::~Cell()
{
	int s;

	delete Host;
	Host = NULL;

	for (s=0; s<nr_symbionts; s++)	delete Symbionts->at(s);

	Symbionts->erase(Symbionts->begin(),Symbionts->end());
	delete Symbionts;
	Symbionts = NULL;
}

void Cell::UpdateOrganelles()
{
	i_bead it;
	int s;

	//First determine gene expression (by making the ExpressedGenes list), then use this to determine movement of expressed genes, then update the organelle state, and lastly update gene expression on the actual genes. This updated gene expression will change the organelle state the next timestep during UpdateOrganelles.

	//Set up the expression lists.
	Host->G->NativeExpression(Host->ExpressedGenes);
	Host->nr_native_expressed = (int)Host->ExpressedGenes->size();
	for (s=0; s<nr_symbionts; s++)
	{
		Symbionts->at(s)->G->NativeExpression(Symbionts->at(s)->ExpressedGenes);
		Symbionts->at(s)->nr_native_expressed = (int) Symbionts->at(s)->ExpressedGenes->size();
	}


	GeneTransport();	//Move around expression products.

	if (!(host_growth == 1 && Host->Stage == 4))	//If waiting is free (option 1), host stays in M 'forever'.
	{
		Host->UpdateState();	//Update organelle state. If foreign products should not interfere with organelle state, put that in the UpdateState() function.
	}
	Host->G->UpdateGeneExpression(Host->ExpressedGenes);	//Update gene expression states.
	it = Host->ExpressedGenes->erase(Host->ExpressedGenes->begin(), Host->ExpressedGenes->end());	//Erase ExpressedGenes to avoid conflicts with pointers during cell dynamics.
	Host->nr_native_expressed = 0;

	for (s=0; s<nr_symbionts; s++)
	{
		Symbionts->at(s)->UpdateState();
		Symbionts->at(s)->G->UpdateGeneExpression(Symbionts->at(s)->ExpressedGenes);
		it = Symbionts->at(s)->ExpressedGenes->erase(Symbionts->at(s)->ExpressedGenes->begin(), Symbionts->at(s)->ExpressedGenes->end());
		Symbionts->at(s)->nr_native_expressed = 0;
	}
}

void Cell::GeneTransport()
{
	i_bead it, it2;
	int s, it_cntr;
	Gene* gene;

	for (s=0; s<nr_symbionts; s++)
	{
		//Movement of expressed regulators from symbionts to host.
		it = Symbionts->at(s)->ExpressedGenes->begin();
		while (it != Symbionts->at(s)->ExpressedGenes->end())
		{
			gene = dynamic_cast<Gene*>(*it);
			if (perfect_transport && gene->signalp[0] == false)
			{
				it2 = it;
				it--;
				Host->ExpressedGenes->splice(Host->ExpressedGenes->end(), *Symbionts->at(s)->ExpressedGenes, it2);
				Symbionts->at(s)->nr_native_expressed--;
			}
			else if (uniform() < leakage_to_host)
			{
				Host->ExpressedGenes->push_back(*it);
			}
			it++;
		}

		//Movement of expressed regulators from host to symbionts.
		it = Host->ExpressedGenes->begin();
		it_cntr = 0;
		while (it_cntr < Host->nr_native_expressed)
		{
			gene = dynamic_cast<Gene*>(*it);
			if (perfect_transport && gene->signalp[0] == true)
			{
				if (s == nr_symbionts-1)	//Only erase the expressed gene from the host if we get to the last symbiont (already transported to all other symbionts).
				{
					it2 = it;
					it--;
					it_cntr--;	//Important, bug in the first attempt of M. commu V.
					Symbionts->at(s)->ExpressedGenes->splice(Symbionts->at(s)->ExpressedGenes->end(), *Host->ExpressedGenes, it2);
					Host->nr_native_expressed--;	//I think this will be fine (although it defines the while-loop), because the iterator takes a step back.
				}
				else
				{
					Symbionts->at(s)->ExpressedGenes->push_back(*it);
				}
			}
			else if (uniform() < leakage_to_symbiont)
			{
				Symbionts->at(s)->ExpressedGenes->push_back(*it);
			}
			it++;
			it_cntr++;
		}
	}
}


void Cell::DNATransferToHost()
{
	i_bead it, insertsite;
	int s, k;
	double uu;

	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s) != NULL)
		{
			it = Symbionts->at(s)->G->BeadList->begin();
			while (it != Symbionts->at(s)->G->BeadList->end())
			{
				if ( (*it)->kind == HOUSE || (*it)->kind == BSITE )	//Let's simply always do copy-paste.
				{
					if ( uniform() < muT[SYMBIONT][(*it)->kind] )	TransferBead(it, Host);
					it++;
					break;
				}
				else if ( (*it)->kind == REGULATOR || (*it)->kind == EFFECTOR )
				{
					uu = uniform();
					if (Symbionts->at(s)->Stage >= 2)
					{
						if (uu < muT[SYMBIONT][(*it)->kind])					it = TransferGene(it, Symbionts->at(s), Host, false, false);
						else if (uu < 2*muT[SYMBIONT][(*it)->kind])		it = TransferGene(it, Symbionts->at(s), Host, true, false);
						else	it++;
					}
					else
					{
						if (uu < muT[SYMBIONT][(*it)->kind])					it = TransferGene(it, Symbionts->at(s), Host, false, true);
						else if (uu < 2*muT[SYMBIONT][(*it)->kind])		it = TransferGene(it, Symbionts->at(s), Host, true, true);
						else	it++;
					}
				}
			}

			for (k=0; k<4; k++)
			{
				assert(Symbionts->at(s)->G->gnr[k] == Symbionts->at(s)->G->CountBeads(k));
				assert(Host->G->gnr[k] == Host->G->CountBeads(k));
			}
		}
	}
}

void Cell::DNATransfertoSymbiont(Organelle* Symbiont)
{
	int k;
	double uu;

	i_bead it = Host->G->BeadList->begin();
	while( it != Host->G->BeadList->end() )
	{
		if ( (*it)->kind == HOUSE || (*it)->kind == BSITE )
		{
			if (uniform() < muT[HOST][(*it)->kind])	TransferBead(it, Symbiont);
			it++;
		}
		else if ( (*it)->kind == REGULATOR || (*it)->kind == EFFECTOR )
		{
			uu = uniform();
			if (Host->Stage >= 2)
			{
				if (uu < muT[HOST][(*it)->kind]) 					it = TransferGene(it, Host, Symbiont, false, false);
				else if (uu < 2*muT[HOST][(*it)->kind])		it = TransferGene(it, Host, Symbiont, true, false);
				else	it++;
			}
			else
			{
				if (uu < muT[HOST][(*it)->kind])					it = TransferGene(it, Host, Symbiont, false, true);
				else if (uu < 2*muT[HOST][(*it)->kind])		it = TransferGene(it, Host, Symbiont, true, true);
				else	it++;
			}
		}
	}

	for (k=0; k<4; k++)
	{
		assert(Host->G->gnr[k] == Host->G->CountBeads(k));
		assert(Symbiont->G->gnr[k] == Symbiont->G->CountBeads(k));
	}
}


Cell::i_bead Cell::TransferGene(i_bead it, Organelle* Source, Organelle* Target, bool include_distal, bool cut_and_paste)
{
	int copy_length;
	i_bead insertsite, first, last, ii;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	last = it;		//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=Source->G->FindFirstBsiteInFrontOfGene(it, include_distal);	//First tfbs in front of gene.
	copy_length=distance(first, last);

	Source->G->CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	if (include_distal)
	{
		ii = BeadListTemp.begin();	//Don't transfer the household genes, because these will result in immediate fitness penalties.
		while (ii != BeadListTemp.end())
		{
			if ((*ii)->kind == HOUSE)
			{
				delete (*ii);
				ii = BeadListTemp.erase(ii);
				copy_length--;
			}
			else
			{
				ii++;
			}
		}
	}

	if (cut_and_paste)
	{
		ii = first;
		while (ii != last)
		{
			if ((*ii)->kind != HOUSE)	//Houses are never transferred, also not with include_distal.
			{
				Source->G->g_length--;
				Source->G->terminus_position--;
				if ((*ii)->kind == BSITE)						Source->G->gnr[BSITE]--;
				else if ((*ii)->kind == REGULATOR)	Source->G->gnr[REGULATOR]--;
				else if ((*ii)->kind == EFFECTOR)		Source->G->gnr[EFFECTOR]--;
				delete (*ii);
				ii = Source->G->BeadList->erase(ii);
			}
			else
			{
				ii++;
			}
		}
		if (!Source->mutant)	//Cut-and-paste results in mutations during the lifetime of an organelle; these organelles (if not yet denoted as mutants) are flagged and added to the fossil record inside Population.cc.
		{
			Source->lifetime_mutant = true;
		}
	}

	//Splice the temporary chromosome into the full genome.
	insertsite=Target->G->FindRandomGenePosition(true,true);			//Find position to insert gene (may be the end of the genome).
	insertsite=Target->G->FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	Target->G->BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	insertsite--;	//Go to the just inserted gene.
	if ((*insertsite)->kind == REGULATOR)	Target->G->PotentialTypeChange(insertsite);	//Regulator will be redefined based on the new genomic context, using the new organelle's RegTypeList.

	//Increment the number of beads and the number of genes.
	Target->G->g_length+=copy_length;
	if ((*insertsite)->kind == REGULATOR)			Target->G->gnr[REGULATOR]++;
	else if((*insertsite)->kind == EFFECTOR)	Target->G->gnr[EFFECTOR]++;
	Target->G->gnr[BSITE]+=copy_length-1;	//The length of the whole transferred piece except for the gene (i.e. you will always transfer 1 gene with x bsites and nothing else).
	Target->G->is_mutated = true;
	Target->mutant = true;
	Target->G->terminus_position = Target->G->g_length;
	return last;
}

void Cell::TransferBead(i_bead it, Organelle* Target)
{
	Bead* bead = (*it)->Clone();
	i_bead insertsite = Target->G->FindRandomPosition(true);
	Target->G->BeadList->insert(insertsite, bead);
	Target->G->g_length++;
	Target->G->terminus_position = Target->G->g_length;
	switch ( bead->kind )
	{
		case BSITE:
			Target->G->gnr[BSITE]++;
			break;
		case HOUSE:
			Target->G->gnr[HOUSE]++;
			break;
	}
	Target->G->is_mutated = true;
	Target->mutant = true;
}

void Cell::InitialiseCell(int input_nr)
{
	int i=0;
	ifstream in_genome(genome_files[input_nr].c_str());
	ifstream in_expression(expression_files[input_nr].c_str());
	ifstream in_definition(definition_files[input_nr].c_str());
	string genome;
	string expression;
	string definition;
	Organelle* Symbiont;

	barcode = input_nr;

	if ( !in_genome.is_open() )
	{
		printf("Genome file %s could not be opened.\n", genome_files[input_nr].c_str());
		exit(1);
	}

	if ( !in_expression.is_open() )
	{
		printf("Expression file %s could not be opened.\n", expression_files[input_nr].c_str());
		exit(1);
	}

	if ( !in_definition.is_open() )
	{
		printf("Expression file %s could not be opened.\n", definition_files[input_nr].c_str());
		exit(1);
	}

	while ( !in_genome.eof() & !in_expression.eof() & !in_definition.eof() )
	{
		//We start reading below, because when we reach the end of the file, there will still be one round through the while-loop.
		if (i==1)
		{
			Host->InitialiseOrganelle(genome, expression, definition);
		}

		else if (i>1)
		{
			Symbiont = new Organelle();
			Symbiont->InitialiseOrganelle(genome, expression, definition);
			Symbiont->G->organelle = SYMBIONT;
			Symbionts->push_back(Symbiont);
			nr_symbionts++;
		}

		in_genome >> genome;
		in_expression >> expression;
		in_definition >> definition;


		i++;
	}
}


void Cell::CloneCell(Cell* ImageC, unsigned long long* pid_count)
{
	int s;
	Organelle* Symbiont;

	nr_symbionts = ImageC->nr_symbionts;
	barcode = ImageC->barcode;

	(*pid_count)++;
	Host->CloneOrganelle(ImageC->Host, *pid_count);

	for (s=0; s<ImageC->nr_symbionts; s++)
	{
		Symbiont = new Organelle();
		(*pid_count)++;
		Symbiont->CloneOrganelle(ImageC->Symbionts->at(s), *pid_count);
		Symbionts->push_back(Symbiont);
	}
}

bool Cell::CheckCellDeath(bool output)
{
	if (nr_symbionts == 0)
	{
		if (output)	SingleCellOutput(true);	//Used for FollowSingleCell.
		DeathOfCell();
		return true;				//The cell died.
	}
	else
	{
		return false;				//The cell lives.
	}
}

void Cell::DeathOfSymbiont(int s)
{
	if (Symbionts->at(s)->mutant || trace_lineage || log_lineage)
	{
		Symbionts->at(s)->alive = false;
	}
	else
	{
		delete Symbionts->at(s);
	}

	Symbionts->erase(Symbionts->begin()+s);	//We erase the pointer from the Symbionts vector. Similar to setting C->Host to NULL for host death (see below).
}

void Cell::DeathOfHost()
{
	if (Host->mutant || trace_lineage || log_lineage)
	{
		Host->alive = false;
	}
	else
	{
		delete Host;
	}

	Host = NULL;
}

void Cell::DeathOfCell()	//Called when host died or when last symbiont died, responsible for cleaning of remaining organelles.
{
	int s;

	if (Host != NULL)	DeathOfHost();
	for (s=nr_symbionts-1; s>=0; s--)
	{
		DeathOfSymbiont(s);
		nr_symbionts--;
	}
}

void Cell::CalculateCellFitness()
{
	int s;
	double fitness, nr_houses = 0.;

	nr_houses += Host->nr_houses;
	for (s=0; s<nr_symbionts; s++)
	{
		nr_houses += (double)Symbionts->at(s)->nr_houses / nr_symbionts;
	}

	fitness = Host->CalculateFitness(nr_household_genes, nr_houses);

	Host->fitness = fitness;
	for (s=0; s<nr_symbionts; s++)
	{
		Symbionts->at(s)->fitness = fitness;
	}
}

void Cell::SingleCellOutput(bool death_event)
{
	int s;
	unsigned long long anc_id;


	cout << "T " << Time << endl;

	//Host
	if (Host->Ancestor == NULL)	anc_id = 0;
	else												anc_id = Host->Ancestor->fossil_id;
	cout << "O -1\tI " << Host->fossil_id << "\tA " << anc_id << "\tS " << Host->Stage << "\tP " << Host->privilige << "\tF " << Host->G->fork_position << "\tE " << Host->G->ShowExpression(NULL, false) << endl;

	//Symbionts
	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s)->Ancestor == NULL)		anc_id = 0;
		else																			anc_id = Symbionts->at(s)->Ancestor->fossil_id;
		cout << "O " << s << "\tI " << Symbionts->at(s)->fossil_id << "\tA " << anc_id << "\tS " << Symbionts->at(s)->Stage << "\tP " << Symbionts->at(s)->privilige << "\tF " << Symbionts->at(s)->G->fork_position << "\tE " << Symbionts->at(s)->G->ShowExpression(NULL, false) << endl;
	}

}
