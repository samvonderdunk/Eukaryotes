#include "Cell.hh"

Cell::Cell()
{
	Vesicle = NULL;
	Vesicle = new Organelle();
}

Cell::Cell(int k)
{
	Vesicle = NULL;
	Vesicle = new Organelle();
	barcode = k;
}

Cell::~Cell()
{
	delete Vesicle;
	Vesicle = NULL;
}

void Cell::CalculateCellFitness()
{
	Vesicle->fitness = Vesicle->CalculateFitness(nr_household_genes, Vesicle->nr_houses);
}

bool Cell::BasalDeath()
{
	if (uniform() < death_rate)
	{
		DeathOfCell();
		return true;	//Cell dies.
	}
	else	return false;	//Cell lives.
}

void Cell::DeathOfCell()
{
	if (Vesicle != NULL)
	{
		if (Vesicle->mutant)		Vesicle->alive = false;
		else										delete Vesicle;

		Vesicle = NULL;
	}
}

bool Cell::FailedDivision()
{
	if (Vesicle->Stage == 5)
	{
		DeathOfCell();
		return true;
	}
	else	return false;
}

Cell* Cell::Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count)
{
	//Divide vesicle.
	Cell* NewCell;

	CalculateCellFitness();

	if (Vesicle->Stage == 4 && (uniform() < Vesicle->fitness) && (host_growth == 0 || (host_growth == 1 && *NewSite == NULL)))
	{
		NewCell = new Cell();
		NewCell->barcode = barcode;
		(*pid_count)++;
		NewCell->Vesicle->Mitosis(Vesicle, *pid_count);
		if (NewCell->Vesicle->mutant)		FP->BuryFossil(NewCell->Vesicle);

		if (*NewSite != NULL)
		{
			(*NewSite)->DeathOfCell();
			delete *NewSite;
		}
		*NewSite = NewCell;
	}

	return NULL;	//Unlike eukaryotic division, nothing can go wrong here; the newborn is always viable.
}

void Cell::UpdateOrganelles()
{
	//Set number of native expressed genes. This is used in GeneTransport & UpdateState, so right after here. This avoids having to record native expressed genes anywhere after UpdateGeneExpression until here.
	Vesicle->G->nr_native_expressed = (int)Vesicle->G->ExpressedGenes->size();

	if (!(host_growth == 1 && Vesicle->Stage == 4))	//Use host_growth also for prokaryotes!
	{
		Vesicle->UpdateState();	//Update organelle state. If foreign products should not interfere with organelle state, put that in the UpdateState() function.
	}
	Vesicle->G->UpdateGeneExpression();	//Update gene expression states.
}

void Cell::Replication(double nuts)
{
	if (Vesicle->Stage == 2 && Vesicle->privilige)
	{
		Vesicle->Replicate(Vesicle->nutrient_claim * nuts);
	}
}

Cell::i_bead Cell::TransferGene(i_bead it, Organelle* Source, Organelle* Target, bool include_distal, bool cut_and_paste)
{
	int copy_length;
	i_bead insertsite, first, last, ii;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	//If an expressed gene gets tranferred via cut-and-paste, we will have to flag the organelle, so that we can reset its gene expression (inside Population.cc).
	if (cut_and_paste)	//This has to happen before "it" is transferred and not reachable on source anymore...
	{
		Gene* gene = dynamic_cast<Gene*>(*it);
		if (gene->expression > 0)	Source->exp_gene_transfer = true;	//An expressed gene has been removed from the source genome.
	}

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
				Source->G->gnr[(*ii)->kind]--;
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
			Vesicle->InitialiseOrganelle(genome, expression, definition);
		}

		else if (i>1)
		{
			printf("Too many organelles inside prokaryote; check input files.\n");
			exit(1);
		}

		in_genome >> genome;
		in_expression >> expression;
		in_definition >> definition;

		i++;
	}
}

void Cell::CloneCell(Cell* ImageC, unsigned long long* pid_count)
{
	Cell* ImageP = dynamic_cast<Cell*>(ImageC);

	barcode = ImageP->barcode;

	(*pid_count)++;
	Vesicle->CloneOrganelle(ImageP->Vesicle, *pid_count);
}

string Cell::Show()
{
	string Content="";
	std::stringstream ss;

	ss << "<" << barcode << ">" << endl;

	Content += ss.str();
	ss.clear();

	return Content;
}

string Cell::Show(bool include_organelles, bool include_genomes)
{
	string Content = Show();

	if (include_organelles)
	{
		Content += Vesicle->Output(include_genomes);
		Content += "\n";
	}

	return Content;
}
