#include "Prokaryote.hh"

Prokaryote::Prokaryote() : Cell(PROKARYOTE)
{
	Vesicle = NULL;
	Vesicle = new Organelle();
}

Prokaryote::~Prokaryote()
{
	if (Vesicle != NULL)
	{
		delete Vesicle;
		Vesicle = NULL;
	}
}

void Prokaryote::CalculateCellFitness()
{
	Vesicle->fitness = Vesicle->CalculateFitness(nr_household_genes, Vesicle->nr_houses);
}


bool Prokaryote::BasalDeath()
{
	if (uniform() < death_rate_prok)
	{
		DeathOfCell();
		return true;	//Cell dies.
	}
	else	return false;	//Cell lives.
}

void Prokaryote::DeathOfCell()	//For Prokaryote, it basically merges DeathOfCell() and DeathOfHost(), because it is all so simple.
{
	if (Vesicle != NULL)
	{
		if (Vesicle->mutant || trace_lineage || log_lineage)		Vesicle->alive = false;
		else																										delete Vesicle;

		Vesicle = NULL;
	}
}

bool Prokaryote::FailedDivision()
{
	if (Vesicle->Stage == 5)
	{
		DeathOfCell();
		return true;
	}
	else	return false;
}

Cell* Prokaryote::Division(Cell** NewSite, Fossils* FP, unsigned long long* pid_count)
{
	//Divide vesicle.
	Prokaryote* NewProk;

	if (Vesicle->Stage == 4 && (uniform() < Vesicle->fitness) && (host_growth == 0 || (host_growth == 1 && *NewSite == NULL)))
	{
		NewProk = new Prokaryote();
		NewProk->barcode = barcode;
		(*pid_count)++;
		NewProk->Vesicle->Mitosis(Vesicle, *pid_count);
		if (NewProk->Vesicle->mutant || trace_lineage || log_lineage)		FP->BuryFossil(NewProk->Vesicle);

		if (*NewSite != NULL)
		{
			(*NewSite)->DeathOfCell();
			delete *NewSite;
		}
		*NewSite = NewProk;
	}

	return NULL;	//Unlike eukaryotic division, nothing can go wrong here; the newborn is always viable.
}

void Prokaryote::UpdateOrganelles()
{
	//Set number of native expressed genes. This is used in GeneTransport & UpdateState, so right after here. This avoids having to record native expressed genes anywhere after UpdateGeneExpression until here.
	Vesicle->G->nr_native_expressed = (int)Vesicle->G->ExpressedGenes->size();

	if (!(host_growth == 1 && Vesicle->Stage == 4))	//Use host_growth also for prokaryotes!
	{
		Vesicle->UpdateState();	//Update organelle state. If foreign products should not interfere with organelle state, put that in the UpdateState() function.
	}
	Vesicle->G->UpdateGeneExpression();	//Update gene expression states.
}

void Prokaryote::Replication(nuts n)
{
	if (Vesicle->Stage == 2 && Vesicle->privilige)
	{
		Vesicle->Replicate(Vesicle->nutrient_claim * std::get<2>(n));
	}
}

void Prokaryote::InitialiseCell(int input_nr)
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

void Prokaryote::CloneCell(Cell* ImageC, unsigned long long* pid_count)
{
	Prokaryote* ImageP = dynamic_cast<Prokaryote*>(ImageC);

	barcode = ImageP->barcode;

	(*pid_count)++;
	Vesicle->CloneOrganelle(ImageP->Vesicle, *pid_count);
}

string Prokaryote::Show(bool include_organelles, bool include_genomes)
{
	string Content = Cell::Show();

	if (include_organelles)
	{
		Content += Vesicle->Output(include_genomes);
		Content += "\n";
	}

	return Content;
}
