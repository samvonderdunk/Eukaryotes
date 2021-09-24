#include "Genome.hh"

Genome::Genome()
{
	BeadList=NULL;
	BeadList = new list<Bead*>();
	g_length=0;
	gnr_regulators=0;
	gnr_bsites=0;
	gnr_houses=0;
	fork_position=0;
	terminus_position=0;
	is_mutated=false;
}

Genome::~Genome()
{
	i_bead it;

	it=BeadList->begin();
	while(it!=BeadList->end())
	{
		delete (*it);
		it++;
	}

	it=BeadList->erase(BeadList->begin(),BeadList->end());
	delete BeadList;
	BeadList=NULL;
}



void Genome::UpdateGeneExpression(list<Bead*>* ExpressedGenes)
{
	i_bead it, i_gene;
	int wb, it_cntr;
	double cum_effects = 0.;
	Regulator* reg;
	Bsite* bs;

	//Determine regulatory dynamics and put result into the express variable of each gene (Regulator or otherwise...).
	it = BeadList->begin();
	it_cntr = 0;
	while (it != BeadList->end())
	{
		wb = WhatBead(*it);
		if (wb==REGULATOR)
		{
			reg = dynamic_cast<Regulator*>(*it);

			cum_effects -= (double)reg->threshold;
			reg->express = max(min((int)cum_effects+1,1),0);	//For the +1, see explanation of the gene threshold in Regulator.hh
			cum_effects = 0.;
		}
		else if (wb==BSITE)
		{
			i_gene = RegulatorCompetition(it, ExpressedGenes);
			if (i_gene != BeadList->end())
			{
				reg = dynamic_cast<Regulator*>(*i_gene);
				bs = dynamic_cast<Bsite*>(*it);
				cum_effects += bs->activity * reg->activity;
			}
		}
		it++;
		it_cntr++;
		if (it_cntr == terminus_position)	cum_effects = 0.;
	}

	//Realise the just calculated regulatory dynamics.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if (WhatBead(*it)==REGULATOR)
		{
			reg = dynamic_cast<Regulator*>(*it);
			reg->expression = reg->express;
		}
		it++;
	}
}



void Genome::NativeExpression(list<Bead*>* ExpressedGenes)
{
	i_bead it;

	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if (WhatBead(*it)==REGULATOR)
		{
			//See similar potential issue in UpdateExpression() and in BindingAffinity().
			Regulator* reg = dynamic_cast<Regulator*>(*it);
			if(reg->expression > 0)	ExpressedGenes->push_back(reg);	//Native genes are always stored in ExpressedGenes.
		}
		it++;
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
	Regulator* reg;
	i_bead it;
	double affinity, p_bind, z_partition = 1.;

	//Calculate partition function.
	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		reg = dynamic_cast<Regulator*>(*it);
		affinity = (double)BindingAffinity(bsite->sequence, reg->sequence);
		z_partition += reg->expression * k_zero * exp(affinity * epsilon);
		it++;
	}

	//Pick one of the probabilities.
	double die_roll = uniform();
	die_roll -= 1 / z_partition;
	if (die_roll <= 0.)	return BeadList->end();

	it = ExpressedGenes->begin();
	while (it != ExpressedGenes->end())
	{
		reg = dynamic_cast<Regulator*>(*it);
		affinity = (double)BindingAffinity(bsite->sequence, reg->sequence);
		p_bind = ( reg->expression * k_zero * exp(affinity * epsilon) ) / z_partition;
		die_roll -= p_bind;
		if (die_roll <= 0.)
		{
			return it;
		}
		it++;
	}

	printf("Error: partition function out of bounds.\n");
	exit(1);
}



void Genome::ReplicateStep(double resource)
{
	i_bead it, start, end;
	Bead* bead;
	int gene_length = 0, repl_remaining_steps, wb;
	double res_int, res_fract, fract_repl_remaining;

	if (relative_replication > 0)	//Modify resource to represent relative replication. Still, the fractional part of the resulting (normalised) resource is seen as a probability.
	{
		fract_repl_remaining = resource / relative_replication;
		resource = fract_repl_remaining * terminus_position;
	}

	res_fract = modf(resource, &res_int);	//Split double into its integer and fractional parts.
	repl_remaining_steps = (int) res_int;
	if (uniform() < res_fract)	repl_remaining_steps++;	//Regard fractional replication step as probability.

	if (repl_remaining_steps <= 0)	return;	//Nothing to do if no beads are allowed to replicate.

	start = BeadList->begin();
	advance(start, fork_position);	//Now it points to the the first bead to be replicated in this replication step.
	end = start;	//End starts at start.

	//This loop sets "end".
	while ((distance(BeadList->begin(),end) < terminus_position) && repl_remaining_steps > 0)	//The maximal position of end is defined by pos_anti_ori - 1 (pointing to the last bead of the parental genome). pos_anti_ori holds the number of genes in the parental genome, so if your distance to the first bead is pos_anti_ori, you are actually one past the last bead of the parental genome. The first replication step, this will point to NULL, but in consecutive steps it will point to a child bead; we want to point to an end point that does not change, hence the last bead of the parental genome.
	{
		if(WhatBead(*end)==REGULATOR)	repl_remaining_steps--;	//Genes for sure count for the repl_step_size.
		else if(!gene_replication)	repl_remaining_steps--;	//TFBSs count if replicate_entire_genes is set to false.
		gene_length++;
		end++;
	}
	end--;

	bool last_round = false;
	it = start;	//Loop over a number of beads (how many we will replicate in one step).
	if(start == end)	last_round = true;	//We go straight into the last round.
	while (it != BeadList->end())
	{
		//Copy bead.
		wb = WhatBead(*it);
		bead=(*it)->Clone();
		(*BeadList).push_back(bead);
		g_length++;
		switch (wb)
		{
			case REGULATOR:
				gnr_regulators++;
				break;
			case BSITE:
				gnr_bsites++;
				break;
			case HOUSE:
				gnr_houses++;
				break;
		}

		it++;
		if (last_round)	break;
		if (it == end)	last_round = true;	//We have apparently hit the last bead of the parental genome, so time for one final replication step.
	}


	if (g_length > 500)	//We only check after the full replication step, not after each replicated bead.
	{
		printf("Warning: genome sizes reached extravagant size (%d) during replication.\nExiting just to be safe...\n", g_length);
		cout << Show(NULL, true, false) << endl;
		exit(1);
	}

	fork_position += gene_length;	//Move the fork to the right.
	if (fork_position >= terminus_position)	//If this is the final replication step.
	{
		fork_position = terminus_position;
		assert(g_length == 2*terminus_position);
	}

}



void Genome::SplitGenome(Genome* parentG)	//Used to split a genome upon division
{
	//Find the fork with i_split.
	i_bead i_split = parentG->BeadList->begin();
	advance(i_split, parentG->terminus_position);	//terminus_position points to the end of the parental genome, which is now the first bead of the child genome.

	BeadList->splice(BeadList->begin(), *parentG->BeadList, i_split, parentG->BeadList->end());

	DevelopChildrenGenomes(parentG);
}



void Genome::AbortChildGenome()
{
	i_bead it = BeadList->begin();
	advance(it, terminus_position);

	while (it != BeadList->end())
	{
		switch ( WhatBead(*it) )
		{
			case REGULATOR:
				gnr_regulators--;
				break;
			case BSITE:
				gnr_bsites--;
				break;
			case HOUSE:
				gnr_houses--;
				break;
		}
		g_length--;
		delete(*it);
		it++;
	}

	it = BeadList->begin();
	advance(it, terminus_position);
	it = BeadList->erase(it, BeadList->end());

	fork_position = 0;
}



void Genome::DevelopChildrenGenomes(Genome* parentG)	//Function gets iterators of parental genome, but copies it to child and then acts on variables of child genome.
{
	i_bead it;
	vector<bool>* MutationList;
	int wb, del_length, dup_length, index, g_length_before_mut;
	int* pdup_length, * pdel_length;

	g_length = BeadList->size();
	g_length_before_mut = g_length;

	//Clean up the variables of the parental genome left behind.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		switch ( WhatBead(*it) )
		{
			case REGULATOR:
				parentG->gnr_regulators--;
				break;
			case BSITE:
				parentG->gnr_bsites--;
				break;
			case HOUSE:
				parentG->gnr_houses--;
				break;
		}
		parentG->g_length--;
		it++;
	}

	parentG->fork_position = 0;
	//Copy gene numbers, mutations happen next.
	gnr_regulators = parentG->gnr_regulators;
	gnr_bsites = parentG->gnr_bsites;
	gnr_houses = parentG->gnr_houses;

	if (mutations_on)	//START mutations.
	{
		del_length = 0;
		dup_length = 0;
		pdup_length = &dup_length;
		pdel_length = &del_length;
		index = 0;

		MutationList = new vector<bool>(g_length, false);
		//Mutate beads and set gnr_genes.
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			switch ( WhatBead(*it) )
			{
				case REGULATOR:
					MutationList->at(index) = true;
			 		it = MutateRegulator(it, pdel_length);
					break;
				case BSITE:
					MutationList->at(index) = true;
					it = MutateBsite(it, pdel_length);
					break;
				case HOUSE:
					MutationList->at(index) = true;
					it = MutateHouse(it, pdel_length);
					break;
			}
			index++;
		}

		//Check that all beads had the chance to mutate.
		vector<bool>::iterator mit = MutationList->begin();
		while (mit != MutationList->end())
		{
			assert((*mit));
			mit++;
		}
		delete MutationList;
		MutationList = NULL;

		//Look for duplicated genes and tfbs's, which we will actually duplicate here.
		//g_length is updated inside functions, gnr_regulators etc. are updated here on the outside; except gnr_bsites inside DuplicateGene.
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			if ((*it)->duplicate)
			{
				switch ( WhatBead(*it) )
				{
					case REGULATOR:
						it=DuplicateGene(it, pdup_length);;
						gnr_regulators++;
						break;
					case BSITE:
						it=DuplicateBsite(it);
						gnr_bsites++;
						(*pdup_length)++;
						break;
					case HOUSE:
						it=DuplicateHouse(it);
						gnr_houses++;
						(*pdup_length)++;
						break;
				}
			}
			else	it++;
		}

		//Check that no more genes are tagged for duplication.
		it = BeadList->begin();
		while(it != BeadList->end())
		{
			assert (!(*it)->duplicate);
			it++;
		}

		//Innovations.
		//Note that g_length is updated in DuplicateGene() and DeleteGene() so that any position along the genome is allowed for the novel bead. Here all counters are also updated within the functions.
		if(uniform() < regulator_innovation_mu)
		{
			it = InventRegulator();
			if (independent_regtypes)		DetermineRegType(it);
			else												PotentialTypeChange(it);
			(*pdup_length)++;
		}
		if(uniform() < bsite_innovation_mu)
		{
			InventBsite();
			(*pdup_length)++;
		}
		if(uniform() < house_innovation_mu)
		{
			InventHouse();
			(*pdup_length)++;
		}

		//Shuffling mutations.
		it = BeadList->begin();
		while(it != BeadList->end())
		{
			wb = WhatBead(*it);
			if(wb==REGULATOR && uniform() < regulator_shuffle_mu)	it=ShuffleGene(it);
			else if(wb==BSITE && uniform() < bsite_shuffle_mu)	it=ShuffleBead(it);
			else if(wb==HOUSE && uniform() < house_shuffle_mu) it=ShuffleBead(it);
			else	it++;
		}

	}	//END of mutations.

	if(mutations_on)
	{
		assert(g_length == g_length_before_mut + (*pdup_length) - (*pdel_length));
		assert(g_length == gnr_regulators + gnr_bsites + gnr_houses);
		assert((size_t)g_length == BeadList->size());
	}

	terminus_position = g_length;

}



void Genome::DetermineRegType(i_bead it)
{
	int i;
	bool type_change = false;
	Regulator* reg;

	reg = dynamic_cast<Regulator*>(*it);

	for (i=1; i<=5; i++)
	{
		if (BindingAffinity(reg->typeseq, regtype[i-1], typeseq_length) <= typeseq_length - 2)
		{
			reg->type = i;
			type_change = true;
			break;
		}
	}

	if (reg->type == 0 || ((reg->type >= 1 && reg->type <= 5) && !type_change))
	{
		reg->type = 6+(int)(uniform()*45);	//Type invention gets random type. Also for divergence from type 1-5, a new random type will be defined. In all other cases (type 1-5 to type 1-5, or type >6 to type >6), you don't have to do anything.
	}
}



void Genome::PotentialTypeChange(i_bead it)
{
	i_bead it2;
	Regulator* reg, *reg2;
	list<int> UsedTypes;
	bool found_matching_type = false;
	int type_abundance;

	reg = dynamic_cast<Regulator*>(*it);
	type_abundance = CountTypeAbundance(reg->type);

	it2 = BeadList->begin();	//The other genes in the genome
	while(it2 != BeadList->end())
	{
		if( WhatBead(*it2)==REGULATOR && it2!=it )	//We don't convert genes to themselves.
		{
			reg2 = dynamic_cast<Regulator*>(*it2);
			UsedTypes.push_back(reg2->type);

			if ( BindingAffinity(reg->sequence, reg2->sequence) == 0   &&   reg->activity == reg2->activity )
			{
				reg->type = reg2->type;		//Convert to existing gene type.
				found_matching_type = true;
				return;		//We have found a match, converted the gene; time to try mutation of the next bead.
			}
		}
		it2++;
	}
	if(found_matching_type == false && (type_abundance > 1 || reg->type==0))	//We haven't been able to convert it to an existing type, so let's define it as a new type. type_abundance should always be more than 1 because there is always an original copy on the parental section of the genome (I think that statement is no longer true, because we mutate after genomes have been split in two).
	{
		list<int>::iterator ix;
		int x=0;

		while (ix != UsedTypes.end())		//WARNING: check that this is not true at the start. Otherwise initialise ix as something else then end().
		{
			x++;
			ix = find(UsedTypes.begin(), UsedTypes.end(), x);
		}
		reg->type = x;
		return;
	}
}

int Genome::CountTypeAbundance(int type)
{
	i_bead it;
	Regulator* reg;
	int count_type = 0;

	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if (WhatBead(*it)==REGULATOR)
		{
			reg = dynamic_cast<Regulator*>(*it);
			if (reg->type == type)	count_type++;
		}
		it++;
	}

	return count_type;
}

/*
###########################################################################
###########################################################################
							 |\ /\  |  | ----  /\  ---- o / \  |\ |
							|  V  \ l_J   |   /- \  |   | L_J  | \|
###########################################################################
###########################################################################
*/

Genome::i_bead Genome::MutateRegulator(i_bead it, int* pdel_length)
{
	Regulator* reg;
	reg = dynamic_cast<Regulator*>(*it);
	bool potential_type_change = false;
	int i;

	double uu = uniform();
	if(uu < regulator_duplication_mu)
	{
		(*it)->duplicate = true;	//Mark for duplication during divison.
		is_mutated = true;
		it++;
	}

	else if(uu < regulator_deletion_mu+regulator_duplication_mu)
	{
		it = DeleteGene(it, pdel_length);
		gnr_regulators--;
		is_mutated = true;
	}

	else
	{
		if (uniform() < regulator_threshold_mu)	//Parameters mutate independently.
		{
			reg->threshold = ChangeParameter(reg->threshold);
			is_mutated = true;
		}

		if (uniform() < regulator_activity_mu)
		{
			reg->activity = ChangeParameter(reg->activity);
			potential_type_change = true;
			is_mutated = true;
		}

		if (independent_regtypes)
		{
			for(i=0; i<typeseq_length; i++)
			{
				if (uniform() < regulator_typeseq_mu)
				{
					if (reg->typeseq[i] == false)			reg->typeseq[i] = true;
					else if(reg->typeseq[i] == true)	reg->typeseq[i] = false;
					DetermineRegType(it);
					is_mutated = true;
				}
			}
		}

		for(i=0; i<sequence_length; i++)
		{
			if(uniform() < regulator_sequence_mu)
			{
				if (reg->sequence[i] == false)			reg->sequence[i] = true;
				else if (reg->sequence[i] == true)	reg->sequence[i] = false;
				potential_type_change = true;
				is_mutated = true;
			}
		}

		if (potential_type_change && !independent_regtypes){
			PotentialTypeChange(it);	//Check type change, activity or sequence has been mutated.
		}
		it++;
	}

	return it;
}

Genome::i_bead Genome::MutateBsite(i_bead it, int* pdel_length)
{
	Bsite* bsite;
	bsite = dynamic_cast<Bsite*>(*it);
	int i;

	double uu = uniform();
	if(uu < bsite_duplication_mu)
	{
		(*it)->duplicate = true;
		is_mutated = true;
		it++;
	}

	else if(uu < bsite_duplication_mu+bsite_deletion_mu)
	{
		it = DeleteBead(it);
		gnr_bsites--;
		(*pdel_length)++;
		is_mutated = true;
	}

	else
	{
		for (i=0; i<sequence_length; i++)
		{
			if(uniform() < bsite_sequence_mu)
			{
				if (bsite->sequence[i] == false) bsite->sequence[i] = true;
				else if (bsite->sequence[i] == true) bsite->sequence[i] = false;
				is_mutated = true;
			}
		}
		if(uniform() < bsite_activity_mu)
		{
			bsite->activity = ChangeParameter(bsite->activity);
			is_mutated = true;
		}
		it++;
	}
	return it;
}

Genome::i_bead Genome::MutateHouse(i_bead it, int* pdel_length)
{
	double uu = uniform();

	if(uu < house_duplication_mu)
	{
		(*it)->duplicate = true;
		is_mutated = true;
		it++;
	}
	else if(uu < house_duplication_mu+house_deletion_mu)
	{
		it = DeleteBead(it);
		gnr_houses--;
		(*pdel_length)++;
		is_mutated = true;
	}
	else	it++;

	return it;
}

int Genome::ChangeParameter(int value)
{
	int new_value = value;

	if (uniform()>0.8)
	{
		new_value = (uniform()>0.5) ? new_value+1 : new_value-1;
	}
	else
	{
		while (new_value == value)
		{
			new_value = (int)(uniform()*(2*WeightRange+1) - WeightRange);
		}
	}

	return new_value;
}

Genome::i_bead Genome::DeleteGene(i_bead it, int* pdel_length)
{
	i_bead first, last, it2;
	int del_length;

	last=it;	//gene position
	last++;	//one further than the gene position
	first=FindFirstBsiteInFrontOfGene(it);	//first tfbs in front of gene

	//Decrement the number of beads and the number of genes.
	del_length = distance(first, last);
	g_length -= del_length;
	gnr_bsites -= del_length-1;
	(*pdel_length) += del_length;

	it2=first;
	while( it2 != last )
	{
		delete *it2;
		it2++;
	}
	it=(*BeadList).erase(first, last);
	return it;
}

Genome::i_bead Genome::DeleteBead(i_bead it)
{
	delete (*it);
	it=(*BeadList).erase(it);

	g_length--;
	return it;
}

Genome::i_bead Genome::DuplicateGene(i_bead it, int* pdup_length)
{
	int copy_length;
	i_bead insertsite, first, last;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	(*it)->duplicate = false;	//Remove the duplication flag.

	//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last = it;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first = FindFirstBsiteInFrontOfGene(it);	//First tfbs in front of gene.
	copy_length = distance(first, last);
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	//Splice the temporary chromosome into the full genome.
	insertsite = FindRandomGenePosition(true,true);			//Find position to insert gene (may be the end of the genome).
	insertsite = FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	//Increment the number of beads and the number of genes.
	g_length += copy_length;
	gnr_bsites += copy_length-1;
	(*pdup_length) += copy_length;

	it = last;	//Make sure ii points to one position further than just duplicated gene.
	return it;
}

Genome::i_bead Genome::DuplicateBsite(i_bead it)
{
	i_bead tt, upstream;
	Bsite* bsite = dynamic_cast<Bsite*>(*it);

	bsite->duplicate = false;
	Bsite* bsite_new = new Bsite(*bsite);

	tt = FindRandomPosition(true);
	tt = (*BeadList).insert(tt, bsite_new);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	g_length++;
	it++;
	return it;
}

Genome::i_bead Genome::DuplicateHouse(i_bead it)
{
	i_bead tt, upstream;

	(*it)->duplicate = false;
	House* house_new = new House();

	tt = FindRandomPosition(true);	//Including the position beyond the last current bead.
	tt = (*BeadList).insert(tt, house_new);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	g_length++;
	it++;
	return it;
}

Genome::i_bead Genome::InventRegulator()
{
	Regulator* reg;
	i_bead insertsite;

	reg = new Regulator();
	reg->RandomRegulator();

	insertsite=FindRandomGenePosition(true,true);
	insertsite=FindFirstBsiteInFrontOfGene(insertsite);
	insertsite=BeadList->insert(insertsite, reg);

	g_length++;
	gnr_regulators++;
	is_mutated = true;
	return insertsite;
}

void Genome::InventBsite()
{
	Bsite* bsite;
	i_bead insertsite;

	bsite = new Bsite();
	bsite->RandomBsite();

	insertsite = FindRandomPosition(true);
	BeadList->insert(insertsite, bsite);

	g_length++;
	gnr_bsites++;
	is_mutated = true;
}

void Genome::InventHouse()
{
	House* house;
	i_bead insertsite;

	house = new House();
	insertsite = FindRandomPosition(true);
	BeadList->insert(insertsite, house);

	g_length++;
	gnr_houses++;
	is_mutated = true;
}

Genome::i_bead Genome::ShuffleGene(i_bead it)
{
	i_bead insertsite, first, last;

	last = it;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=FindFirstBsiteInFrontOfGene(it);	//First tfbs in front of gene.

	//Splice the temporary chromosome into the full genome.
	insertsite=FindRandomGenePosition(true,true);			//Find position to insert gene (may be end of genome).
	insertsite=FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.

	if ( (distance(BeadList->begin(), insertsite) >= distance(BeadList->begin(), first))   &&   (distance(BeadList->begin(), insertsite) < distance(BeadList->begin(), last)) )	//We are attempting to splice to same location, so don't do it, does not count as mutation. Although distance is a time-costly operation in UpdateGeneExpression, we probably do not bother with mutations.
	{
		it++;
		return it;
	}
	else
	{
		BeadList->splice(insertsite, *BeadList, first, last);
		is_mutated = true;
		return last;
	}
}

Genome::i_bead Genome::ShuffleBead(i_bead it)
{
	i_bead tt, nit;

	is_mutated = true;
	nit = it;
	nit++;
	tt = FindRandomPosition(true);
	BeadList->splice(tt, *BeadList, it);	//When we don't splice a whole range, the above problem for ShuffleGene should not be relevant.

	return nit;
}

Genome::i_bead Genome::FindFirstBsiteInFrontOfGene(i_bead it) const
{
	ri_bead rit(it);	//BEHIND in this case means an element to the left in the list (i.e if ii points to the 6th element, rii(ii) will make rii point to the 5th element). The function .base() as used below will make rii.base() point to the 6th element again.

	ri_bead rit2 = (*BeadList).rend();//search should be bounded
	while(rit != rit2)//begin not yet reached
	{
		if(WhatBead(*rit)==BSITE)	rit++;
		else	rit2 = rit;
	}
	return rit2.base();
}

Genome::i_bead Genome::FindRandomGenePosition(bool include_houses, bool include_end) const
{
	std::list<i_bead> pos;
	std::list<i_bead>::iterator ipos;
	i_bead it, it2;
	int randpos, add_houses=0, end=0;
	int wb;

	if(include_houses)	add_houses=gnr_houses;
	if(include_end)	end=1;

	if(gnr_regulators+add_houses==0)	return (*BeadList).end();
	else
	{
		it=(*BeadList).begin();
		while(it != (*BeadList).end())
		{
			wb = WhatBead(*it);
			if(wb==REGULATOR || (wb==HOUSE && include_houses))	pos.push_back(it);
			it++;
		}
		pos.push_back(it);
		randpos=(int)(uniform()*(gnr_regulators+add_houses+end));	//If you found the first gene, randpos will be 0 (no need to advance to another gene); if you find the last gene, randpos will be gnr_genes-1 which is enough to reach the gnr_genes'th gene.
		it2 = (*boost::next(pos.begin(),randpos));	//Possibly try advance(ii, randpos) instead.
		return it2;
	}
}

Genome::i_bead Genome::FindRandomPosition(bool include_end) const
{
	int randpos, end=0;
	i_bead it;

	if (include_end)	end=1;	//Insertion can also happen beyond the last bead; deletion only of the current beads.

	randpos = (int)(uniform()* (g_length+end));
	it = BeadList->begin();
	advance(it, randpos);

	return it;
}

void Genome::CopyPartOfGenomeToTemplate(i_bead begin, i_bead end, list<Bead*>* template_beadlist)
{
	i_bead it;
	Bead* bead;

	it=begin;
	while(it!=end)
	{
		bead=(*it)->Clone();
		(*template_beadlist).push_back(bead);
		it++;
	}
}


void Genome::CloneGenome(const Genome* ImageG)
{
	// BeadList = new list<Bead*>();
	CopyPartOfGenome(ImageG->BeadList->begin(),ImageG->BeadList->end());

	g_length = ImageG->g_length;
	gnr_regulators = ImageG->gnr_regulators;
	gnr_bsites = ImageG->gnr_bsites;
	gnr_houses = ImageG->gnr_houses;
	fork_position = ImageG->fork_position;
	terminus_position = ImageG->terminus_position;
	is_mutated = ImageG->is_mutated;
}


void Genome::CopyPartOfGenome(i_bead begin, i_bead end)
{
	//Note that no genome and bead counters are set in this function (see e.g. CloneGenome).
	i_bead it;
	Bead* bead;
	it=begin;
	while(it!=end)
	{
		bead=(*it)->Clone();
		(*BeadList).push_back(bead);
		it++;
	}
}

void Genome::ReadGenome(string genome)
{
	//Blueprint for this function is ReadBeadsFromString() in Prokaryotes. See that function for detailed comments.

	char* bead, *buffer;
	int success, type, threshold, activity;
	bool typeseq[typeseq_length], signalp[signalp_length], sequence[sequence_length];
	Regulator* reg;
	House* house;
	Bsite* bsite;

	bead = strtok((char*)genome.c_str(),".");
	while (bead != NULL)
	{
		if(bead[1] == 'R')
		{
			buffer = new char[typeseq_length+signalp_length+sequence_length+2];
			success = sscanf(bead, "(R%d:%d:%d:%s)", &type, &threshold, &activity, buffer);
			if(success != 4) cerr << "Could not find sufficient information for this regulatory gene. Genome file potentially corrupt. \n" << endl;

			ReadBuffer(buffer, typeseq, 'X', ':');
			ReadBuffer(buffer, signalp, ':', ':', 1, 2);
			ReadBuffer(buffer, sequence, ':', ')', 2);

			reg = new Regulator(type, threshold, activity, typeseq, signalp, sequence, 0);
			(*BeadList).push_back(reg);
			gnr_regulators++;
			g_length++;
			delete [] buffer;
			buffer = NULL;
		}
		else if(bead[1] == 'H')
		{
			house = new House();
			(*BeadList).push_back(house);
			gnr_houses++;
			g_length++;
		}
		else
		{
			buffer = new char[sequence_length+2];
			success = sscanf(bead, "(%d:%s)", &activity, buffer);
			if(success != 2) cerr << "Could not find sufficient information for this binding site. Genome file potentially corrupt. \n" << endl;

			ReadBuffer(buffer, sequence, 'X', ')');

			bsite = new Bsite(activity, sequence);
			(*BeadList).push_back(bsite);
			gnr_bsites++;
			g_length++;
			delete [] buffer;
			buffer = NULL;
		}

		bead = strtok(NULL, ".");
	}

	fork_position = 0;
	terminus_position = g_length;
}

void Genome::ReadBuffer(string buffer, bool* array, char start_sign, char stop_sign, int ith_start_sign, int ith_stop_sign)
{
	int q = 0;
	int start_reading = -1;
	int nstart = 0, nstop = 0;

	if (start_sign == 'X')					start_reading = 0;
	while(nstop < ith_stop_sign)
	{
		if (start_reading != -1)			array[q-start_reading] = (buffer[q]=='1');
		if (buffer[q] == start_sign)
		{
			nstart++;
			if (nstart == ith_start_sign){
				start_reading = q+1;
			}
		}
		q++;
		if (buffer[q] == stop_sign)		nstop++;
	}
}

void Genome::ReadExpression(string expression)
{
	string data;
	int i;
	i_bead it;
	Regulator* reg;

	data = expression.substr(1,expression.size()-2);	//Strip accolades or square brackets.

	it = BeadList->begin();
	for(i=0; i<gnr_regulators; i++)
	{
		while (it != BeadList->end() && WhatBead(*it) != REGULATOR) it++;

		if ((size_t) i > data.length())
		{
			cerr << "Could not find sufficient expression states. Expression file potentially corrupt.\n" << endl;
			exit(1);
		}

		reg = dynamic_cast<Regulator*>(*it);
		reg->expression = (data[i]=='1');
		it++;
	}
}

string Genome::Show(list<Bead*>* chromosome, bool terminal, bool only_parent)
{
	string GenomeContent="", expressed_prefix, reg_color_prefix, reg_color_suffix, bsite_color_prefix, bsite_color_suffix, house_color_prefix, house_color_suffix, prefix;
	i_bead it, end;
	Regulator* reg;
	Bsite* bsite;
	int i;


	if(terminal){
		// expressed_prefix = "\033[43m";
		reg_color_prefix = "\033[94m";
		reg_color_suffix = "\033[0m";
		bsite_color_prefix = "\033[92m";
		bsite_color_suffix = "\033[0m";
		house_color_prefix = "\033[95m";
		house_color_suffix = "\033[0m";
	}
	else
	{
		reg_color_prefix = "";
		reg_color_suffix = "";
		bsite_color_prefix = "";
		bsite_color_suffix = "";
		house_color_prefix = "";
		house_color_suffix = "";
	}

	if(chromosome == NULL) chromosome = this->BeadList;
	it = chromosome->begin();
	if (only_parent)
	{
		end = chromosome->begin();
		advance(end, terminus_position);
	}
	else	end = chromosome->end();

	while(it != end)
	{
		if(it != chromosome->begin()) GenomeContent += ".";
		GenomeContent += "(";

		std::stringstream ss;
		switch (WhatBead(*it))
		{
			case REGULATOR:
				reg=dynamic_cast<Regulator*>(*it);
				// if (reg->expression > 0)	prefix = expressed_prefix;
				// else									prefix = reg_color_prefix;
				ss << reg_color_prefix << "R" << reg->type << ":" << reg->threshold << ":" << reg->activity << ":";
				for(i=0; i<typeseq_length; i++)	ss << reg->typeseq[i];
				ss << ":";
				for(i=0; i<signalp_length; i++)	ss << reg->signalp[i];
				ss << ":";
				for(i=0; i<sequence_length; i++)	ss << reg->sequence[i];
				ss << reg_color_suffix;
				GenomeContent += ss.str();
				ss.clear();
				break;
			case BSITE:
				bsite=dynamic_cast<Bsite*>(*it);
				ss << bsite_color_prefix << bsite->activity << ":";
				for(i=0; i<sequence_length; i++)	ss << bsite->sequence[i];
				ss << bsite_color_suffix;
				GenomeContent += ss.str();
				ss.clear();
				break;
			case HOUSE:
				ss << house_color_prefix << "H" << house_color_suffix;
				GenomeContent += ss.str();
				ss.clear();
				break;
		}

		GenomeContent += ")";
		it++;
	}
	return GenomeContent;
}

string Genome::ShowExpression(list<Bead*>* chromosome, bool only_parent)
{
	i_bead it, end;
	Regulator* reg;
	string ExpressionContent="{";

	if(chromosome == NULL) chromosome = this->BeadList;
	it = chromosome->begin();
	if (only_parent)
	{
		end = chromosome->begin();
		advance(end, terminus_position);
	}
	else	end = chromosome->end();

	while (it != end)
	{
		if (WhatBead(*it)==REGULATOR)
		{
			std::stringstream ss;
			reg = dynamic_cast<Regulator*>(*it);
			ss << reg->expression;
			ExpressionContent += ss.str();
			ss.clear();
		}
		it++;
	}

	ExpressionContent += "}";
	return ExpressionContent;
}
