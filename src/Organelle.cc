#include "Organelle.hh"

Organelle::Organelle()
{
	Stage=0;
	privilige=false;
	fitness=1.;
	nutrient_claim=1.;
	G = NULL;
	G = new Genome();
	ExpressedGenes = NULL;
	ExpressedGenes = new list<Bead*>();
	nr_native_expressed = 0;

	alive = false;
	mutant = false;
	lifetime_mutant = false;
	time_of_appearance = 0;
	fossil_id = 0;
	Ancestor = NULL;
}

Organelle::~Organelle()
{
	i_bead it;

	//We only delete the list of pointers, they point to the things in the genome G, so we will remove the actual genes below (or they might actually be genes from different genomes in which case we also don't want to delete them).
	if (ExpressedGenes != NULL)	it = ExpressedGenes->erase(ExpressedGenes->begin(), ExpressedGenes->end());
	delete ExpressedGenes;
	ExpressedGenes = NULL;

	delete G;
	G = NULL;
}

void Organelle::UpdateState()
{
	int readout[5] = {0, 0, 0, 0, 0};	//States of the five cell-cycle regulators.
	i_bead it, it2;
	Gene* gene, *reg2;
	int i;

	int eval_state = Stage, it_cntr;
	if (Stage == 2   &&   G->fork_position != G->terminus_position)	eval_state--;	//You cannot proceed to G2 without finishing replication.
	privilige = true;

	//Fill in the readout.
	it = ExpressedGenes->begin();
	it_cntr = 0;
	while (it != ExpressedGenes->end())
	{
		if (G->WhatBead(*it)==REGULATOR)
		{
			gene = dynamic_cast<Gene*>(*it);
			if (it_cntr < nr_native_expressed || independent_regtypes)	//Native genes from the organelle itself; just look at the type.
			//For now, under independent gene types, we assume no complication in the definition of regulatory type between compartments; something assigned to type 1 in the host, will also be counted as type 1 in the symbiont.
			{
				if (gene->type < 6)
				{
					readout[gene->type-1] += gene->expression;	//In case expression can ever be something beside 0 and 1 (otherwise you could just do ++).
				}
			}
			else	//Leaked/transported genes from other organelles; see if they match one of the 5 key regulator types in the current organelle; if not, not interesting for the state of the cell (i.e. do nothing). If you don't want to involve foreign expression in organelle state, exclude this entire part.
			{
				if (regtypes_by_regulation)
				{
					for (i=0; i<5; i++)
					{
						if (   G->BindingAffinity(gene->sequence, G->RegTypeList[i]->sequence) + abs(gene->activity - G->RegTypeList[i]->activity) <= regtype_hdist   )
						{
							readout[i] += gene->expression;
						}
						break;	//We assume that only 1 gene type will be matched. This means that when two type definitions come to close by mutation, only the lower gene type will be scored, probably preventing types from coming to close. NEW: actually check out PotentialTypeChange() to see that the type defs can never come within the regtype_hdist.
					}
				}
				else	//Primitive carry-over of Prokaryotes: check if you resemble one of the other genes currently in the genome that you end up in. This is primitive, because you could never fully replace a gene type, only function like already existing gene types.
				{
					it2 = G->BeadList->begin();
					while (it2 != G->BeadList->end())
					{
						if (G->WhatBead(*it2)==REGULATOR)
						{
							reg2 = dynamic_cast<Gene*>(*it2);
							if ( G->BindingAffinity(gene->sequence, reg2->sequence) == 0   &&   gene->activity == reg2->activity )
							{
								if (reg2->type < 6)
								{
									readout[reg2->type-1] += gene->expression;	//Add the expression of this type to the native type that it resembles.
								}
								break;	//We have found a matching gene, and possibly also read its expression (if type 1-5), so we don't have to search any further in the genome.
							}
						}
						it2++;
					}
				}
			}
		}
		it++;
		it_cntr++;
	}

	//Compare readout (expression states) with cell-cycle states.
	if (EvaluateState(eval_state, readout) == 5) Stage = eval_state + 1;	//Success!
	else if (EvaluateState(3, readout) == 5)	Stage = 5;	//Marked for death during attempted mitosis.
	else if (Stage == 2   &&   G->fork_position != G->terminus_position)	privilige = false;
}

int Organelle::EvaluateState(int eval_state, int* readout)
{
	int i, match=0;

	for (i=0;i<5;i++)
	{
		if ((readout[i]>0) == StageTargets[eval_state][i])	match++;
	}

	return match;
}

void Organelle::Mitosis(Organelle* parent, unsigned long long id_count)
{
	G->is_symbiont = parent->G->is_symbiont;
	G->SplitGenome(parent->G);

	parent->Stage = 0;
	parent->privilige = false;

	if (house_duplication_mu > 0.0 || house_deletion_mu > 0.0)	fitness = 1. - abs(nr_household_genes - G->gnr_houses) / (float)10;

	time_of_appearance = Time;
	fossil_id = id_count;
	alive = true;

	nutrient_claim = parent->nutrient_claim;
	if (nutshare_evolve && uniform() < nutrient_claim_mu)
	{
		nutrient_claim += nutrient_claim_mu_delta*(uniform()-0.5);
	}

	if (G->is_mutated)	mutant = true;

	if (parent->mutant || trace_lineage || log_lineage)		Ancestor = parent;
	else																									Ancestor = parent->Ancestor;
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

	fitness = 1. - abs(nr_household_genes - G->gnr_houses) / (float)10;
	nutrient_claim = init_nutrient_claim;

	cout << G->Show(NULL, true, false) << endl;	//This should just do the Show() function.
}

void Organelle::CloneOrganelle(Organelle* ImageO, unsigned long long id_count)
{
	fossil_id = id_count;
	Stage = ImageO->Stage;
	privilige = ImageO->privilige;
	fitness = ImageO->fitness;
	nutrient_claim = ImageO->nutrient_claim;
	mutant = ImageO->mutant;
	lifetime_mutant = ImageO->lifetime_mutant;
	alive = ImageO->alive;
	time_of_appearance = ImageO->time_of_appearance;

	Ancestor = ImageO->Ancestor;

	G->CloneGenome(ImageO->G);
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
