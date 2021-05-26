#include "Cell.hh"

Cell::Cell()
{
	Host=NULL;
	for(int i=0;i<HS;i++)
	{
		Symbionts[i]=NULL;
	}
}

Cell::~Cell()
{
	delete Host;
	Host=NULL;
	for(int i=0;i<HS;i++)
	{
		if (Symbionts[i]!=NULL)
		{
			delete (Symbionts[i]);
			Symbionts[i]=NULL;
		}
	}
}

void Cell::UpdateOrganelles()
{
	int i;

	Host->UpdateExpression();
	for (i=0;i<HS;i++)	Symbionts[i]->UpdateExpression();

	RegulatorTransport();

	Host->UpdateState();
	for (i=0;i<HS;i++)	Symbionts[i]->UpdateState();

	//RegulatorTransport();	//You can do it here if leaking proteins should not affect cell-cycle stage (or before UpdateExpression(), which would be identical).
}

void Cell::RegulatorTransport()
{
	i_bead it;
	int i, nr_native_host_genes;

	nr_native_host_genes = (int) Host->ExpressedGenes->size();
	for (i=0;i<HS;i++)
	{
		//Movement of expressed regulators from symbionts to host.
		it = Symbionts[i]->ExpressedGenes->begin();
		while (it != Symbionts[i]->ExpressedGenes->end())
		{
			if (ActiveTransport(it, Symbionts[i]->ExpressedGenes, Host->ExpressedGenes) || uniform() < leakage_to_host)
			{
				Host->ExpressedGenes->push_back(*it);
			}
			it++;
		}

		//Movement of expressed regulators from host to symbionts.
		it = Host->ExpressedGenes->begin();
		while (distance(Host->ExpressedGenes->begin(), it) < nr_native_host_genes)
		{
			if (ActiveTransport(it, Host->ExpressedGenes, Symbionts[i]->ExpressedGenes) || uniform() < leakage_to_symbiont)
			{
				Symbionts[i]->ExpressedGenes->push_back(*it);
			}
			it++;
		}
	}
}

bool Cell::ActiveTransport(i_bead it, list<Bead*>* SourceCompartment, list<Bead*>* TargetCompartment)
{
	//For future use: if transporters are expressed, they can here lead to very likely transport of genes.
	//Maybe passing BeadLists is a bit cumbersome as well...
	return false;	//For now, there will be no active transport.
}


void Cell::DNATransferToHost()
{
	int i;
	i_bead it, insertsite;

	for (i=0;i<HS;i++)
	{
		if (Symbionts[i] != NULL)
		{
			it = Symbionts[i]->G->BeadList->begin();
			while (it != Symbionts[i]->G->BeadList->end())
			{
				switch Symbionts[i]->G->WhatBead(*it)
				{
					case 'R':
						if (uniform() < regulator_transfer_mu_StoH)	TransferGene(it, Host);
					case 'B':
						if (uniform() < bsite_transfer_mu_StoH) TransferBead(it, Host)
					case 'H':
						if (uniform() < house_transfer_mu_StoH)	TransferBead(it, Host);
				}
				it++;
			}
		}
	}
}

void Cell::DNATransfertoSymbiont(int s)
{
	i_bead it = Host->G->BeadList->begin();
	while( it != Host->G->BeadList->end() )
	{
		switch Host->G->WhatBead(*it)
		{
			case 'R':
				if (uniform() < regulator_transfer_mu_HtoS) TransferGene(it, Symbionts[s]);
			case 'B':
				if (uniform() < bsite_transfer_mu_HtoS) TransferGene(it, Symbionts[s]);
			case 'H':
				if (uniform() < house_transfer_mu_HtoS) TransferGene(it, Symbionts[s]);
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
}



void Cell::TransferBead(i_bead it, Organelle* Target)
{
	Bead* bead = (*it)->Clone();
	i_bead insertsite = Target->G->FindRandomPosition(true);
	Target->G->BeadList->insert(insertsite, bead);
	Target->G->g_length++;
	switch Target->G->WhatBead(bead)
	{
		case 'B':	Target->G->gnr_bsites++;
		case 'H':	Target->G->gnr_houses++;
	}
	Target->G->is_mutated = true;
}

void Cell::InitialiseCell(unsigned long long id_count)
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
			Host = new Organelle();
			Host->InitialiseOrganelle(genome, expression);
		}

		else
		{
			Symbiont = new Organelle();
			Symbiont->InitialiseOrganelle(genome, expression);
			Symbionts->push_back(Symbiont);
		}

		id_count++;
		i++;
	}

	//Print first individuals...
}


void CellOne::CloneCell(Cell* ImageC, unsigned long long id_count)
{
	int i;
	Host = new Organelle();
	Host->CloneOrganelle(ImageC->Host);
	id_count++;
	for (i=0;i<HS;i++)
	{
		Symbionts[i] = new Organelle();
		Symbionts[i]->CloneOrganelle(ImageC->Symbionts[i]);
		id_count++;
	}
}
