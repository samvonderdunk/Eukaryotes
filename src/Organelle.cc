#include "Organelle.hh"

Organelle::Organelle()
{
	Stage=0;
	privilige=false;
	fitness=1.;
	mutant = false;
	G = NULL;
	G = new Genome();
	ExpressedGenes = NULL;
	ExpressedGenes = new list<Bead*>();
	nr_native_expressed = 0;
}

Organelle::~Organelle()
{
	i_bead it;

	//We only delete the list of pointers, they point to the things in the genome G, so we will remove the actual genes below (or they might actually be genes from different genomes in which case we also don't want to delete them).
	it = ExpressedGenes->erase(ExpressedGenes->begin(), ExpressedGenes->end());
	delete ExpressedGenes;
	ExpressedGenes = NULL;

	delete G;
	G = NULL;
}

void Organelle::UpdateState()
{
	int readout[5] = {0, 0, 0, 0, 0};	//States of the five cell-cycle regulators.
	i_bead it, it2;
	Regulator* reg, *reg2;

	int eval_state = Stage;
	if (Stage == 2   &&   G->fork_position != G->terminus_position)	eval_state--;	//You cannot proceed to G2 without finishing replication.
	privilige = true;

	//Fill in the readout.
	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		if (G->WhatBead(*it)==REGULATOR)
		{
			reg = dynamic_cast<Regulator*>(*it);
			if (distance(ExpressedGenes->begin(),it) < nr_native_expressed)	//Native genes from the organelle itself; just look at the type.
			{
				if (reg->type < 6)
				{
					readout[reg->type-1] += reg->expression;	//In case expression can ever be something beside 0 and 1 (otherwise you could just do ++).
				}
			}
			else	//Leaked/transported genes from other organelles; see if they match one of the 5 key regulator types in the current organelle; if not, not interesting for the state of the cell (i.e. do nothing). If you don't want to involve foreign expression in organelle state, exclude this entire part.
			{
				it2 = G->BeadList->begin();
				while (it2 != G->BeadList->end())
				{
					if (G->WhatBead(*it2)==REGULATOR)
					{
						reg2 = dynamic_cast<Regulator*>(*it2);
						if ( G->BindingAffinity(reg->sequence, reg2->sequence) == 0   &&   reg->activity == reg2->activity)
						{
							readout[reg2->type-1] += reg->expression;	//Add the expression of this type to the native type that it resembles.
						}
						break;	//WARNING: check where this puts you, should be at the next ExpressedGene.
					}
					it2++;
				}
			}
		}
		it++;
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
		if ((readout[i]>0) == StageTargets[eval_state][i])
			match++;
	}

	return match;
}

void Organelle::Mitosis(Organelle* parent, unsigned long long id_count)
{
	G->SplitGenome(parent->G);

	parent->Stage = 0;
	privilige = false;

	if (house_duplication_mu > 0.0 || house_deletion_mu > 0.0)	//Only if they can actually be gained or lost.
	{
		fitness = 1. - abs(nr_household_genes - G->gnr_houses) / (float)10;
	}
}

void Organelle::Replicate(double resource)
{
	if (G->fork_position < G->terminus_position)
	{
		G->ReplicateStep(resource);
	}
}

void Organelle::InitialiseOrganelle(string genome, string expression)
{
	mutant = true;
	G->ReadGenome(genome);
	G->ReadExpression(expression);
	cout << G->Show(NULL, true, false) << endl;
	// G->SetExpression(ExpressedGenes, false);
}

void Organelle::CloneOrganelle(Organelle* ImageO)
{
	Stage = ImageO->Stage;
	privilige = ImageO->privilige;
	fitness = ImageO->fitness;
	mutant = ImageO->mutant;

	//First copy the genome, then the expression list can be updated (because we want pointers to the new genome!).
	G->CloneGenome(ImageO->G);
	// G->SetExpression(ExpressedGenes, false);
}

string Organelle::ShowExpression()
{
	string ExpressionContent="[";
	i_bead it;
	Regulator* reg;

	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		std::stringstream ss;
		reg = dynamic_cast<Regulator*>(*it);
		if (distance(ExpressedGenes->begin(),it) != 0)	ss << ",";
		ss << reg->type << ":" << reg->expression;
		ExpressionContent += ss.str();
		ss.clear();
		it++;
	}
	ExpressionContent += "]";
	return ExpressionContent;
}
