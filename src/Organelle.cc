#include "Organelle.hh"

Organelle::Organelle()
{
	Stage=0;
	privilige=false;
	nr_houses=0;
	fitness=1.;
	nutrient_claim=1.;
	G = NULL;
	G = new Genome();

	alive = false;
	mutant = false;
	lifetime_mutant = false;
	exp_gene_transfer = false;
	time_of_appearance = 0;
	fossil_id = 0;
	Ancestor = NULL;
}

Organelle::~Organelle()
{
	delete G;
	G = NULL;
}

void Organelle::UpdateState()
{
	i_reg ir;
	int i, it_cntr, eval_state = Stage;
	std::bitset<5> readout;	//States of the five cell-cycle regulators. Overloaded with the states of effector genes if they are functional.

	if (Stage == 2   &&   G->fork_position != G->terminus_position)	eval_state--;	//You cannot proceed to G2 without finishing replication.
	privilige = true;

	//Fill in the readout.
	ir = G->ExpressedGenes->begin();
	it_cntr = 0;
	while (ir != G->ExpressedGenes->end())
	{
		if (it_cntr < G->nr_native_expressed || mu[TYPE][REGULATOR]!=0.0)	//Native genes from the organelle itself; just look at the type. OR, gene types are mutating independently so we also directly consider gene types 1-5.
		{
			if ((*ir)->type < 6)			readout.set((*ir)->type-1);	//Exp===1.
		}
		else	//Foreign genes.
		{
			for (i=0; i<5; i++)
			{
				if (   (*ir)->BindingAffinity((*ir)->sequence, G->RegTypeList[i]->sequence) + abs((*ir)->activity - G->RegTypeList[i]->activity) <= seq_hdist   )
				{
					readout.set(i);	//Exp===1
				}
				break;	//We assume that only 1 gene type will be matched. This means that when two type definitions come to close by mutation, only the lower gene type will be scored, probably preventing types from coming to close. NEW: actually check out PotentialTypeChange() to see that the type defs can never come within the regtype_hdist.
			}
		}
		ir++;
		it_cntr++;
	}

	//Compare readout (expression states) with cell-cycle states.
	if (EvaluateState(StageTargets[eval_state], readout) == 5) Stage = eval_state + 1;	//Success!
	else if (EvaluateState(StageTargets[3], readout) == 5)	Stage = 5;	//Marked for death during attempted mitosis.
	else if (Stage == 2   &&   G->fork_position != G->terminus_position)	privilige = false;
}


double Organelle::CalculateFitness(int target_nr, double real_nr)
{
	if (minimum_houses)
	{
		return min(1., 1. - (target_nr - real_nr) / (float)10); //Hard-coded version where we only require a minimum number of household genes.
	}
	else
	{
		return 1. - abs(target_nr - real_nr) / (float)10;	//There are so many mutation rates to check that it probably is faster to just calculate the fitness always.
	}
}


void Organelle::Mitosis(Organelle* parent, unsigned long long id_count)
{
	//Main genome-level algorithm.
	G->organelle = parent->G->organelle;
	G->SplitGenome(parent->G);

	//Just recreate the Expression list.
	G->ExpressedGenes->clear();
	G->NativeExpression();
	parent->G->ExpressedGenes->clear();
	parent->G->NativeExpression();

	parent->Stage = 0;
	parent->privilige = false;

	nr_houses = G->gnr[HOUSE];	//Determined at birth; gnr[HOUSE] will change during replication.
	if (!cell_fitness)		fitness = CalculateFitness(nr_household_genes, G->gnr[HOUSE]);

	time_of_appearance = Time;
	fossil_id = id_count;
	alive = true;

	nutrient_claim = parent->nutrient_claim;
	if (nutshare_evolve && uniform() < nutrient_claim_mu)
	{
		nutrient_claim += nutrient_claim_mu_delta*(uniform()-0.5);
	}

	if (G->is_mutated)	mutant = true;

	if (parent->mutant)		Ancestor = parent;
	else									Ancestor = parent->Ancestor;
}

void Organelle::Replicate(double resource)
{
	if (G->fork_position < G->terminus_position)
	{
		G->ReplicateStep(resource);
	}
}

void Organelle::Abort()
{
	G->AbortChildGenome();
	G->ExpressedGenes->clear();	//Also renew expression state.
	G->NativeExpression();

	Stage = 0;
	privilige = false;
}

void Organelle::InitialiseOrganelle(string genome, string expression, string definition)
{
	mutant = true;
	alive = true;
	time_of_appearance = 0;	//We don't do anything with id_count, because we will clone the initial organelle.
	Stage = init_stage;

	G->ReadGenome(genome);
	G->ReadExpression(expression);
	G->ReadDefinition(definition);

	G->NativeExpression();	//Set up gene expression list.

	nr_houses = G->gnr[HOUSE];
	if (!cell_fitness)	fitness = CalculateFitness(nr_household_genes, G->gnr[HOUSE]);
	nutrient_claim = init_nutrient_claim;

	cout << G->Show(NULL, true, false) << endl;	//This should just do the Show() function.
}

void Organelle::CloneOrganelle(Organelle* ImageO, unsigned long long id_count)
{
	fossil_id = id_count;
	Stage = ImageO->Stage;
	privilige = ImageO->privilige;
	nr_houses = ImageO->nr_houses;
	fitness = ImageO->fitness;
	nutrient_claim = ImageO->nutrient_claim;
	mutant = ImageO->mutant;
	lifetime_mutant = ImageO->lifetime_mutant;
	alive = ImageO->alive;
	time_of_appearance = ImageO->time_of_appearance;

	Ancestor = ImageO->Ancestor;

	G->CloneGenome(ImageO->G);
	G->NativeExpression();
}

string Organelle::Show(bool backup)
{
	string Content="";
	unsigned long long AncestorID;
	string is_mutant, has_privilige;
	std::stringstream ss;

	has_privilige = (privilige)?"Y":"N";
	if (backup)
	{
		if (Ancestor==NULL)	AncestorID = 0;
		else	AncestorID = Ancestor->fossil_id;

		is_mutant = (mutant)?"Y":"N";

		ss << "[" << fossil_id << " " << AncestorID << " " << is_mutant << " " << fitness << " " << nutrient_claim << " " << Stage << " " << has_privilige << " " << G->fork_position << " " << G->terminus_position << "]";
	}
	else
	{
		ss << fossil_id << "\t" << fitness << "\t" << nutrient_claim << "\t" << Stage << "\t" << has_privilige << "\t" << G->terminus_position << "\t" << G->g_length;
	}

	Content += ss.str();
	ss.clear();

	return Content;
}

string Organelle::Output(bool backup)
{
	string Content="";

	Content += Show(backup);

	if (backup)
	{
		Content += "\t";
		Content += G->ShowExpression(NULL, false);
		Content += "\t";
		Content += G->ShowDefinition(false);
		Content += "\t";
		Content += G->Show(NULL, false, false);
	}

	return Content;
}
