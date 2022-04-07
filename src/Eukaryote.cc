#include "Eukaryote.hh"

Eukaryote::Eukaryote() : Cell(EUKARYOTE)
{
	nr_symbionts = 0;
	Host = NULL;
	Host = new Organelle();
	Symbionts = NULL;
	Symbionts = new vector<Organelle*>();
}

Eukaryote::~Eukaryote()
{
	int s;

	if (Host != NULL)
	{
		delete Host;
		Host = NULL;
	}

	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s) != NULL)	delete Symbionts->at(s);
	}

	Symbionts->clear();
	delete Symbionts;
	Symbionts = NULL;
}

void Eukaryote::CalculateCellFitness()
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


bool Eukaryote::BasalDeath()
{
	int s;
	if (uniform() < death_rate_host)
	{
		DeathOfCell();
		return true;	//Host dies, cell dies.
	}
	if (!moran_symbionts)	//Otherwise, if no moran process, then symbionts cannot die.
	{
		for (s=nr_symbionts-1; s>=0; s--)
		{
			if (uniform() < death_rate_symbiont)
			{
				DeathOfSymbiont(s);
				nr_symbionts--;
			}
		}
		if (LostSymbionts(false))
		{
			DeathOfCell();
			return true;	//Last symbiont dies, cell dies.
		}
	}
	return false;	//Cell did not die.
}

void Eukaryote::DeathOfCell()	//Called when host died or when last symbiont died, responsible for cleaning of remaining organelles.
{
	int s;

	if (Host != NULL)	DeathOfHost();
	for (s=nr_symbionts-1; s>=0; s--)
	{
		if (Symbionts->at(s) != NULL)	DeathOfSymbiont(s);
		nr_symbionts--;
	}
}

bool Eukaryote::LostSymbionts(bool output)
{
	if (nr_symbionts == 0)	return true;				//The cell died.
	else										return false;				//The cell lives.
}

void Eukaryote::DeathOfSymbiont(int s)
{
	if (Symbionts->at(s)->mutant || trace_lineage || log_lineage)			Symbionts->at(s)->alive = false;
	else																															delete Symbionts->at(s);

	Symbionts->erase(Symbionts->begin()+s);	//We erase the pointer from the Symbionts vector. Similar to setting C->Host to NULL for host death (see below).
}

void Eukaryote::DeathOfHost()
{
	if (Host->mutant || trace_lineage || log_lineage)									Host->alive = false;
	else																															delete Host;

	Host = NULL;
}

bool Eukaryote::FailedDivision()
{
	int s;
	if (Host->Stage == 5)
	{
		DeathOfCell();
		return true;
	}
	for (s=nr_symbionts-1; s>=0; s--)
	{
		if ( Symbionts->at(s)->Stage == 5)
		{
			if (moran_symbionts)	Symbionts->at(s)->Abort();	//No death but just abortion.
			else
			{
				DeathOfSymbiont(s);
				nr_symbionts--;
			}
		}
	}

	if (LostSymbionts(false))
	{
		DeathOfCell();
		return true;
	}
	else	return false;
}



Cell* Eukaryote::Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count)
{
	Eukaryote* NewEuk;
	int s, pick_s;

	if (cell_fitness)	CalculateCellFitness();

	//Divide host.
	if (Host->Stage == 4 && (uniform() < Host->fitness) && (host_growth == 0 || (host_growth == 1 && *NewSite == NULL)))
	{
		NewEuk = new Eukaryote();
		NewEuk->barcode = barcode;
		(*pid_count)++;
		NewEuk->Host->Mitosis(Host, *pid_count);

		//Distribute symbionts.
		for (s=nr_symbionts-1; s>=0; s--)
		{
			if (moran_symbionts)	CloneSymbiont(s, this, NewEuk, FP, pid_count);
			else
			{
				if (uniform() < 0.5)	//Transfer the symbiont to the new cell (just a matter of using the right pointers).
				{
					NewEuk->Symbionts->push_back(Symbionts->at(s));
					Symbionts->erase(Symbionts->begin()+s);
					NewEuk->nr_symbionts++;
					nr_symbionts--;
				}
			}
		}

		//Safe distribution of symbionts.
		if (safe_symbiont_distribution)	//If we do not allow division to remove all symbionts from one of the daughters.
		{
			if (NewEuk->nr_symbionts == 0)
			{
				pick_s = (int) (uniform() * (double)nr_symbionts);
				CloneSymbiont(pick_s, this, NewEuk, FP, pid_count);
			}
			else if (nr_symbionts == 0)
			{
				pick_s = (int) (uniform() * (double)NewEuk->nr_symbionts);
				CloneSymbiont(pick_s, NewEuk, this, FP, pid_count);
			}
		}

		if (!NewEuk->LostSymbionts(false))	//NewEuk did not die (it received at least 1 symbiont.)
		{
			if (mutations_on)
			{
				NewEuk->DNATransferToHost();
			}
			if (NewEuk->Host->mutant || trace_lineage || log_lineage)
			{
				FP->BuryFossil(NewEuk->Host);	//Bury fossil only after potential transfer events (also count as mutations).
				for (s=0; s<NewEuk->nr_symbionts; s++)	//Check whether cut-and-paste led to new mutants among symbionts.
				{
					if (NewEuk->Symbionts->at(s)->lifetime_mutant && !NewEuk->Symbionts->at(s)->mutant)
					{
						FP->BuryFossil(NewEuk->Symbionts->at(s));
						NewEuk->Symbionts->at(s)->G->is_mutated = true;
						NewEuk->Symbionts->at(s)->mutant = true;
					}
				}
			}
		}
		else
		{
			delete NewEuk;
			NewEuk = NULL;
		}

		if (!LostSymbionts(false))	//Parent did not die (received at least 1 symbiont). If empty_division_killing == false, and NewEuk did die, don't overgrow NewSite after all. If empty_division_killing == true, we always want to overgrow the neighbour.
		{
			if (empty_division_killing || NewEuk != NULL)	//Either we always want to kill the neighbour, or we note that we have a viable child (and a viable parent, so put the child in the spot of the neighbour).
			{
				if (*NewSite != NULL)
				{
					(*NewSite)->DeathOfCell();
					delete *NewSite;
					*NewSite = NULL;
				}
				*NewSite = NewEuk;
				NewEuk = NULL;
			}
		}
		else
		{
			DeathOfCell();	//Make sure potential fossils of host and symbionts are stored properly.
			return NewEuk;	//Swap in Population.cc
		}
		return NULL;	//Nothing else needs to happen in Population.cc
	}
	return NULL;
}

void Eukaryote::CloneSymbiont(int s, Eukaryote* SourceE, Eukaryote* TargetE, Fossils* FP, unsigned long long* pid_count)
{
	Organelle* SymbiontCopy;

	SymbiontCopy = new Organelle();
	(*pid_count)++;
	SymbiontCopy->CloneOrganelle(SourceE->Symbionts->at(s), *pid_count);
	if (SymbiontCopy->mutant || trace_lineage || log_lineage)		FP->BuryFossil(SymbiontCopy);
	TargetE->Symbionts->push_back(SymbiontCopy);
	TargetE->nr_symbionts++;
}

void Eukaryote::SymbiontDivisions(Fossils* FP, unsigned long long* pid_count)
{
	Organelle* SymbiontCopy;
	int s, pick_s;

	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s)->Stage == 4   &&   (uniform() < Symbionts->at(s)->fitness))
		{
			SymbiontCopy = new Organelle();
			(*pid_count)++;
			SymbiontCopy->Mitosis(Symbionts->at(s), *pid_count);
			if (mutations_on)	DNATransfertoSymbiont(SymbiontCopy);
			if (moran_symbionts || (symbiont_overgrowth>0 && uniform()<(double)(nr_symbionts-1)/(symbiont_overgrowth-1)) )	//If there is overgrowth, the chance of killing a colleague is proportional to the number of colleagues divided by the total room available.
			{
				do	pick_s = (int) (uniform() * (double)nr_symbionts);
				while(pick_s == s);
				DeathOfSymbiont(pick_s);
			}
			else
			{
				nr_symbionts++;
			}
			Symbionts->push_back(SymbiontCopy);
			if (Symbionts->at(nr_symbionts-1)->mutant)
			{
				FP->BuryFossil(Symbionts->at(nr_symbionts-1));
			}
		}
	}
	if (Host->lifetime_mutant && !Host->mutant)
	{
		FP->BuryFossil(Host);
		Host->G->is_mutated = true;
		Host->mutant = true;
	}
}

void Eukaryote::UpdateOrganelles()
{
	i_bead it;
	int s;

	//Set number of native expressed genes. This is used in GeneTransport & UpdateState, so right after here. This avoids having to record native expressed genes anywhere after UpdateGeneExpression until here.
	Host->G->nr_native_expressed = (int)Host->G->ExpressedGenes->size();
	for (s=0; s<nr_symbionts; s++)
	{
		Symbionts->at(s)->G->nr_native_expressed = (int) Symbionts->at(s)->G->ExpressedGenes->size();
	}

	GeneTransport();	//Move around expression products.

	if (!(host_growth == 1 && Host->Stage == 4))	//If waiting is free (option 1), host stays in M 'forever'.
	{
		Host->UpdateState();	//Update organelle state. If foreign products should not interfere with organelle state, put that in the UpdateState() function.
	}
	Host->G->UpdateGeneExpression();	//Update gene expression states.
	for (s=0; s<nr_symbionts; s++)
	{
		Symbionts->at(s)->UpdateState();
		Symbionts->at(s)->G->UpdateGeneExpression();
	}
}

void Eukaryote::Replication(nuts n)
{
	int s;
	if (Host->Stage == 2 && Host->privilige)
	{
		Host->Replicate(Host->nutrient_claim * std::get<0>(n));
	}
	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s)->Stage == 2 && Symbionts->at(s)->privilige)
		{
			Symbionts->at(s)->Replicate(Symbionts->at(s)->nutrient_claim * std::get<1>(n));
		}
	}
}

void Eukaryote::GeneTransport()
{
	i_reg ir, ir2;
	int s, it_cntr;

	/*
	signal peptides:

	10 -> HOST
	01 -> SYMBIONT
	11 -> dual
	00 -> no transport (remains in compartment where produced).
	*/

	for (s=0; s<nr_symbionts; s++)
	{
		//Movement of expressed regulators from symbionts to host.
		ir = Symbionts->at(s)->G->ExpressedGenes->begin();
		while (ir != Symbionts->at(s)->G->ExpressedGenes->end())
		{
			if (perfect_transport && (*ir)->signalp.test(0) && !(*ir)->signalp.test(1))
			{	//Protein translocated to host.
				ir2 = ir;
				ir--;
				Host->G->ExpressedGenes->splice(Host->G->ExpressedGenes->end(), *Symbionts->at(s)->G->ExpressedGenes, ir2);
				Symbionts->at(s)->G->nr_native_expressed--;
			}
			else if ((perfect_transport && (*ir)->signalp.test(0) && (*ir)->signalp.test(1)) || uniform() < leakage_to_host)
			{	//Protein also transported to host.
				Host->G->ExpressedGenes->push_back(*ir);
			}
			ir++;
		}

		//Movement of expressed regulators from host to symbionts.
		ir = Host->G->ExpressedGenes->begin();
		it_cntr = 0;
		while (it_cntr < Host->G->nr_native_expressed)
		{
			if (perfect_transport && !(*ir)->signalp.test(0) && (*ir)->signalp.test(1))
			{	//Protein translocated to symbiont.
				if (s == nr_symbionts-1)	//Only erase the expressed gene from the host if we get to the last symbiont (already transported to all other symbionts).
				{
					ir2 = ir;
					ir--;
					it_cntr--;	//Important, bug in the first attempt of M. commu V.
					Symbionts->at(s)->G->ExpressedGenes->splice(Symbionts->at(s)->G->ExpressedGenes->end(), *Host->G->ExpressedGenes, ir2);
					Host->G->nr_native_expressed--;	//I think this will be fine (although it defines the while-loop), because the iterator takes a step back.
				}
				else
				{
					Symbionts->at(s)->G->ExpressedGenes->push_back(*ir);
				}
			}
			else if ((perfect_transport && (*ir)->signalp.test(0) && (*ir)->signalp.test(1)) || uniform() < leakage_to_symbiont)
			{
				Symbionts->at(s)->G->ExpressedGenes->push_back(*ir);
			}
			ir++;
			it_cntr++;
		}
	}
}


void Eukaryote::DNATransferToHost()
{
	i_bead it, insertsite;
	int s;
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
						if (uu < 0.5*muT[SYMBIONT][(*it)->kind])			it = TransferGene(it, Symbionts->at(s), Host, false, false);
						else if (uu < muT[SYMBIONT][(*it)->kind])			it = TransferGene(it, Symbionts->at(s), Host, true, false);
						else	it++;
					}
					else
					{
						if (uu < 0.5*muT[SYMBIONT][(*it)->kind])			it = TransferGene(it, Symbionts->at(s), Host, false, true);
						else if (uu < muT[SYMBIONT][(*it)->kind])			it = TransferGene(it, Symbionts->at(s), Host, true, true);
						else	it++;
					}
				}
			}

			Symbionts->at(s)->G->CheckBeadCounts();
			Host->G->CheckBeadCounts();
		}
	}
}

void Eukaryote::DNATransfertoSymbiont(Organelle* Symbiont)
{
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
				if (uu < 0.5*muT[HOST][(*it)->kind]) 			it = TransferGene(it, Host, Symbiont, false, false);
				else if (uu < muT[HOST][(*it)->kind])			it = TransferGene(it, Host, Symbiont, true, false);
				else	it++;
			}
			else
			{
				if (uu < 0.5*muT[HOST][(*it)->kind])			it = TransferGene(it, Host, Symbiont, false, true);
				else if (uu < muT[HOST][(*it)->kind])			it = TransferGene(it, Host, Symbiont, true, true);
				else	it++;
			}
		}
	}

	Host->G->CheckBeadCounts();
	Symbiont->G->CheckBeadCounts();
}

void Eukaryote::InitialiseCell(int input_nr)
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

void Eukaryote::CloneCell(Cell* ImageC, unsigned long long* pid_count)
{
	Eukaryote* ImageE = dynamic_cast<Eukaryote*>(ImageC);
	int s;
	Organelle* Symbiont;

	nr_symbionts = ImageE->nr_symbionts;
	barcode = ImageE->barcode;

	(*pid_count)++;
	Host->CloneOrganelle(ImageE->Host, *pid_count);

	for (s=0; s<ImageE->nr_symbionts; s++)
	{
		Symbiont = new Organelle();
		(*pid_count)++;
		Symbiont->CloneOrganelle(ImageE->Symbionts->at(s), *pid_count);
		Symbionts->push_back(Symbiont);
	}
}

string Eukaryote::Show(bool include_organelles, bool include_genomes)
{
	int s;
	string Content = Cell::Show();

	if (include_organelles)
	{
		Content += Host->Output(include_genomes);
		Content += "\n";
		for (s=0; s<nr_symbionts; s++)
		{
			Content += Symbionts->at(s)->Output(include_genomes);
			Content += "\n";
		}
	}

	return Content;
}
