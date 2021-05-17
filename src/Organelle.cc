#include "Organelle.hh"

Organelle::Organelle()
{
	Stage=0;
	privilige=false;
	ExpressedGenes=NULL;
	G=NULL;
}

Organelle::~Organelle()
{
	i_bead it;

	if(ExpressedGenes!=NULL)
	{
		it=ExpressedGenes->begin();
		while(it!=ExpressedGenes->end())
		{
			delete (*it);
			it++;
		}

		it=ExpressedGenes->erase(ExpressedGenes->begin(), ExpressedGenes->end());
		delete ExpressedGenes;
		ExpressedGenes=NULL;
	}

	if (G!=NULL)
	{
		delete (G);
		G=NULL;
	}
}

void Organelle::UpdateExpression()
{
	G->UpdateExpression(ExpressedGenes);	//Starts after regulator transport and ends with only native organelle expression. Updates gene "express" variables using ExpressedGenes.
	G->EraseExpression(ExpressedGenes);
	G->SetExpression(ExpressedGenes, true);	//Transfer new expression states to ExpressedGenes (iterators of expressed genes).
}

void Organelle::UpdateState()
{
	int* readout[5] = {0, 0, 0, 0, 0};	//States of the five cell-cycle regulators.
	i_bead it, it2;
	Regulator* reg, reg2;

	int eval_state = Stage;
	if (Stage == 2   &&   G->fork_position != G->terminus_position)	eval_state--;	//You cannot proceed to G2 without finishing replication.
	privilige = true;

	//Fill in the readout.
	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		if (G->WhatBead(*it)=='R')
		{
			reg = dynamic_cast<Regulator*>(*it);
			if (distance(ExpressedGenes,it) < G->gnr_regulators)	//Native genes from the organelle itself; just look at the type.
			{
				if (reg->type < 5)
				{
					readout[reg->type-1] += reg->expression;	//In case expression can ever be something beside 0 and 1 (otherwise you could just do ++).
				}
			}
			else	//Leaked/transported genes from other organelles; see if they match one of the 5 key regulator types in the current organelle; if not, not interesting for the state of the cell (i.e. do nothing).
			{
				it2 = G->BeadList->begin();
				while (it2 != G->BeadList->end())
				{
					if (G->WhatBead(*it2)=='R')
					{
						reg2 = dynamic_cast<Regulator*>(*it2);
						if ( G->BindingAffinity(reg, reg2) == 0   &&   reg->activity == reg2->activity)
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
		if (readout[i] == StageTargets[eval_state][i])
			match++;
	}

	return match;
}

void Organelle::Mitosis(Organelle* parent, unsigned long long id_count)
{
	G->SplitGenome(parent->G);

	parent->Stage = 0;
	privilige = false;

	//Reset expression in parent and daughter by reading gene states.
	parent->G->EraseExpression(parent->ExpressedGenes);
	parent->G->SetExpression(parent->ExpressedGenes, false);
	G->EraseExpression(ExpressedGenes);
	G->SetExpression(ExpressedGenes, false);

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
