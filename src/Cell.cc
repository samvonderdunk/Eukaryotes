#include "Cell.hh"

Cell::Cell()
{
	nr_symbionts = 0;
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


	RegulatorTransport();	//Move around expression products.

	Host->UpdateState();	//Update organelle state. If foreign products should not interfere with organelle state, put that in the UpdateState() function.
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

void Cell::RegulatorTransport()
{
	i_bead it;
	int s;

	for (s=0; s<nr_symbionts; s++)
	{
		//Movement of expressed regulators from symbionts to host.
		it = Symbionts->at(s)->ExpressedGenes->begin();
		while (it != Symbionts->at(s)->ExpressedGenes->end())
		{
			if (ActiveTransport(it, Symbionts->at(s)->ExpressedGenes, Host->ExpressedGenes) || uniform() < leakage_to_host)
			{
				Host->ExpressedGenes->push_back(*it);
			}
			it++;
		}

		//Movement of expressed regulators from host to symbionts.
		it = Host->ExpressedGenes->begin();
		while (distance(Host->ExpressedGenes->begin(), it) < Host->nr_native_expressed)
		{
			if (ActiveTransport(it, Host->ExpressedGenes, Symbionts->at(s)->ExpressedGenes) || uniform() < leakage_to_symbiont)
			{
				Symbionts->at(s)->ExpressedGenes->push_back(*it);
			}
			it++;
		}
	}
}

bool Cell::ActiveTransport(i_bead it, list<Bead*>* SourceCompartment, list<Bead*>* TargetCompartment)
{
	//For future use: based on signal peptides of proteins, more protein movement can occur.
	//Maybe passing BeadLists is a bit cumbersome as well...
	return false;	//For now, there will be no active transport.
}


void Cell::DNATransferToHost()
{
	i_bead it, insertsite;
	int s;

	for (s=0; s<nr_symbionts; s++)
	{
		if (Symbionts->at(s) != NULL)
		{
			it = Symbionts->at(s)->G->BeadList->begin();
			while (it != Symbionts->at(s)->G->BeadList->end())
			{
				switch ( Symbionts->at(s)->G->WhatBead(*it) )
				{
					case REGULATOR:
						if (uniform() < regulator_transfer_mu_StoH)	TransferGene(it, Symbionts->at(s), Host);
						break;
					case BSITE:
						if (uniform() < bsite_transfer_mu_StoH) TransferBead(it, Host);
						break;
					case HOUSE:
						if (uniform() < house_transfer_mu_StoH)	TransferBead(it, Host);
						break;
				}
				it++;
			}
		}
	}
}

void Cell::DNATransfertoSymbiont(Organelle* Symbiont)
{
	i_bead it = Host->G->BeadList->begin();
	while( it != Host->G->BeadList->end() )
	{
		switch ( Host->G->WhatBead(*it) )
		{
			case REGULATOR:
				if (uniform() < regulator_transfer_mu_HtoS) TransferGene(it, Host, Symbiont);
				break;
			case BSITE:
				if (uniform() < bsite_transfer_mu_HtoS) TransferBead(it, Symbiont);
				break;
			case HOUSE:
				if (uniform() < house_transfer_mu_HtoS) TransferBead(it, Symbiont);
				break;
		}
		it++;
	}
}


void Cell::TransferGene(i_bead it, Organelle* Source, Organelle* Target)
{
	int copy_length;
	i_bead insertsite, first, last;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	last = it;		//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=Source->G->FindFirstBsiteInFrontOfGene(it);	//First tfbs in front of gene.
	copy_length=distance(first, last);
	Source->G->CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	//Splice the temporary chromosome into the full genome.
	insertsite=Target->G->FindRandomGenePosition(true,true);			//Find position to insert gene (may be the end of the genome).
	insertsite=Target->G->FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	Target->G->BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	//Increment the number of beads and the number of genes.
	Target->G->g_length+=copy_length;
	Target->G->gnr_regulators++;
	Target->G->gnr_bsites+=copy_length-1;	//The length of the whole transferred piece except for the gene (i.e. you will always transfer 1 gene with x bsites and nothing else).
	Target->G->is_mutated = true;
	Target->mutant = true;
	Target->G->terminus_position = Target->G->g_length;
}



void Cell::TransferBead(i_bead it, Organelle* Target)
{
	Bead* bead = (*it)->Clone();
	i_bead insertsite = Target->G->FindRandomPosition(true);
	Target->G->BeadList->insert(insertsite, bead);
	Target->G->g_length++;
	Target->G->terminus_position = Target->G->g_length;
	switch ( Target->G->WhatBead(bead) )
	{
		case BSITE:
			Target->G->gnr_bsites++;
			break;
		case HOUSE:
			Target->G->gnr_houses++;
			break;
	}
	Target->G->is_mutated = true;
	Target->mutant = true;
}

void Cell::InitialiseCell()
{
	int i=0;
	ifstream in_genome(genome_initialisation.c_str());
	ifstream in_expression(expression_initialisation.c_str());
	string genome;
	string expression;
	Organelle* Symbiont;

	if ( !in_genome.is_open() )
	{
		printf("Genome file %s could not be opened.\n", genome_initialisation.c_str());
		exit(1);
	}

	if ( !in_expression.is_open() )
	{
		printf("Expression file %s could not be opened.\n", expression_initialisation.c_str());
		exit(1);
	}

	while ( !in_genome.eof() & !in_expression.eof() )
	{
		in_genome >> genome;
		in_expression >> expression;

		if (i==0)
		{
			Host->InitialiseOrganelle(genome, expression);
		}

		else
		{
			Symbiont = new Organelle();
			Symbiont->InitialiseOrganelle(genome, expression);
			Symbionts->push_back(Symbiont);
			nr_symbionts++;
		}

		i++;
	}

	//Print first individuals...
}


void Cell::CloneCell(Cell* ImageC, unsigned long long* pid_count)
{
	int s;
	Organelle* Symbiont;

	nr_symbionts = ImageC->nr_symbionts;

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
