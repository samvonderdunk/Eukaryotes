#include "Organelle.hh"

Organelle::Organelle()
{
	stage=0;
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
	G->UpdateExpression(ExpressedGenes);	//Starts after regulator transport and ends with only native organelle expression.
}

void Organelle::UpdateState()
{
	int* readout[5] = {0, 0, 0, 0, 0};	//States of the five cell-cycle regulators.
	i_bead it;
	int eval_state = Stage;
	if (Stage == 2   &&   G->fork_position != G->terminus_position)	eval_state--;	//You cannot proceed to G2 without finishing replication.

	//Fill in the readout.
	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		if (G->WhatBead(*it)=='R')
		{
			Regulator* reg = dynamic_cast<Regulator*>(*it);
			if (reg->type < 5)
			{
				readout[reg->type-1] += reg->expression;	//In case expression can ever be something beside 0 and 1 (otherwise you could just do ++).
			}
		}
		it++;
	}

	//Compare readout (expression states) with cell-cycle states.
	if (EvaluateState(eval_state, readout) == 5) Stage = eval_state + 1;	//Success!
	else if (EvaluateState(3, readout) == 5)	Stage = 5;	//Marked for death during attempted mitosis.
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
