#include "Genome.hh"

Genome::Genome()
{
	BeadList=NULL;
	g_length=0;
	gnr_genes=0;
	gnr_bsites=0;
	gnr_transporters=0;
	gnr_houses=0;
	fork_position=0;
	terminus_position=0;
	is_mutated=false;
}

Genome::~Genome()
{
	i_bead it;

	if(BeadList!=NULL) {
		it=BeadList->begin();
		while(it!=BeadList->end()) {
			delete (*it);
			it++;
		}

		it=BeadList->erase(BeadList->begin(),BeadList->end());
		delete BeadList;
		BeadList=NULL;
	}
}

void Genome::UpdateExpression(list<Bead*>* ExpressedGenes)
{
	i_bead it;
	char wb;

	//Determine regulatory dynamics and put result into the express variable of each gene (Regulator or Transporter).
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		wb = WhatBead(*it)
		if (wb=='R' || wb=='T')
		{
			if(wb=='R')					Regulator* gene = dynamic_cast<Regulator*>(*it);			//Maybe recasting is not allowed, in which case I should separate the regulators and transporters over if statements.
			else if(wb=='T')		Transporter* gene = dynamic_cast<Transporter*>(*it);

			cum_effects -= gene->threshold;
			gene->express = max(min((int)cum_effects+1,1),0);	//For the +1, see explanation of the gene threshold in Regulator.hh
			cum_effects = 0;
		}
		else if (wb=='B')
		{
			i_gene = RegulatorCompetition(it, ExpressedGenes);
			if (i_gene != BeadList->end())
			{
				Regulator* reg = dynamic_cast<Regulator*>(*i_gene);
				Bsite* bs = dynamic_cast<Bsite*>(*it);
				cum_effects += bs->activity * reg->activity;
			}
		}
		it++;
		if (distance(BeadList->begin(), it) == terminus_position)	cum_effects = 0.;	//WARNING: check later how I define terminus_position; should be the first double bead (i.e. the one replicated first, lying directly after the original genome).
	}

	//Transfer express value to expression state of each gene, and store pointers to expressed genes in ExpressedGenes (object from the Organelle class).
	EraseExpressedGenes(ExpressedGenes);
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		wb = WhatBead(*it);
		if (wb=='R' || wb=='T')
		{
			if(wb=='R')					Regulator* gene = dynamic_cast<Regulator*>(*it);			//See similar potential issue in UpdateExpression().
			else if(wb=='T')		Transporter* gene = dynamic_cast<Transporter*>(*it);

			gene->expression = gene->express;
			if (gene->expression > 0)
			{
				ExpressedGenes->push_back(gene);
			}
		}
	}
}

void Genome::EraseExpressedGenes(list<Bead*>* ExpressedGenes)
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
		// delete ExpressedGenes;	//I don't think I want to delete the object itself, because we'll be rewriting it immediately.
		// ExpressedGenes=NULL;
	}
}

Genome::i_bead Genome::RegulatorCompetition(i_bead i_bsite, list<Bead*>* ExpressedGenes)
{
	//Optimise this, many possibilities:
	//	--> inline functions.
	//  --> reduce casting operations.
	//	--> only calculate binding affinity once.
	//	--> obtain binding affinity from some stored variable or look-up in hash-table?
	//	--> sort the affinities before rolling the die to shorten the second part.


	Bsite* bsite = dynamic_cast<Bsite*>(*i_bsite);
	i_bead it;
	double affinity, p_bind, z_partition = 1.;

	//Calculate partition function.
	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		if (WhatBead(*it)=='R')
		{
			Regulator* reg = dynamic_cast<Regulator*>(*it);
			affinity = (double)BindingAffinity(i_bsite, it);
			z_partition += reg->expression * k_zero * exp(affinity * epsilon);
		}
		it++;
	}

	//Pick one of the probabilities.
	double die_roll = uniform();
	die_roll -= 1 / z_partition;
	if (die_roll <= 0.)	return BeadList->end();

	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		if (WhatBead(*it)=='R')
		{
			Regulator* reg = dynamic_cast<Regulator*>(*it);
			affinity = (double)BindingAffinity(i_bsite, it);
			p_bind = ( reg->expression * k_zero * exp(affinity * epsilon) ) / z_partition;
			die_roll -= p_bind;
			if (die_roll <= 0.)
			{
				return it;
			}
		}
		it++;
	}

	printf("Error: partition function out of bounds.\n");
	exit(1);
}
