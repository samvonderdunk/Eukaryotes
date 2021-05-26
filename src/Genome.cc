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
}



void Genome::EraseExpression(list<Bead*>* ExpressedGenes)
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



void Genome::SetExpression(list<Bead*>* ExpressedGenes, bool Updating)
{
	i_bead it;

	it = BeadList->begin();
	while (it != BeadList->end())
	{
		wb = WhatBead(*it);
		if (wb=='R' || wb=='T')
		{
			//See similar potential issue in UpdateExpression() and in BindingAffinity().
			if(wb=='R')					Regulator* gene = dynamic_cast<Regulator*>(*it);
			else if(wb=='T')		Transporter* gene = dynamic_cast<Transporter*>(*it);

			if (Updating)								gene->expression = gene->express;		//If we have just updated gene states, expression should be read from "express". Otherwise we might just want to update ExpressedGenes (as in Mitosis).
			if (gene->expression > 0)		ExpressedGenes->push_back(gene);
		}
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
			affinity = (double)BindingAffinity(bsite, reg);
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
			affinity = (double)BindingAffinity(bsite, reg);
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



void Genome::ReplicateStep(double resource)
{
	i_bead it, start, end;
	Bead* bead;
	int gene_length = 0;
	int repl_remaining_steps;
	double res_int, res_fract;
	double fract_repl_remaining;
	char wb;

	if (relative_replication)	//Modify resource to represent relative replication. Still, the fractional part of the resulting (normalised) resource is seen as a probability.
	{
		fract_repl_remaining = resource / rel_repl_full;
		resource = fract_repl_remaining * terminus_position;
	}

	res_fract = modf(resource, &res_int);	//Split double into its integer and fractional parts.
	repl_remaining_steps = (int) res_int;
	if (uniform() < res_fract)	repl_remaining_steps++;	//Regard fractional replication step as probability.

	if (repl_remaining_steps==0)	return;	//Nothing to do if no beads are allowed to replicate.

	start = BeadList->begin();
	advance(start, fork_position);	//Now it points to the the first bead to be replicated in this replication step.
	end = start;	//End starts at start.

	//This loop sets "end".
	while ((distance(BeadList->begin(),end) < terminus_position) && repl_remaining_steps > 0)	//The maximal position of end is defined by pos_anti_ori - 1 (pointing to the last bead of the parental genome). pos_anti_ori holds the number of genes in the parental genome, so if your distance to the first bead is pos_anti_ori, you are actually one past the last bead of the parental genome. The first replication step, this will point to NULL, but in consecutive steps it will point to a child bead; we want to point to an end point that does not change, hence the last bead of the parental genome.
	{
		wb = WhatBead(*end);
		if(wb=='R' || wb=='T')	repl_remaining_steps--;	//Genes for sure count for the repl_step_size.
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
		switch wb
		{
			case 'R':	gnr_regulators++;
			case 'T': gnr_transporters++;
			case 'B': gnr_bsites++;
			case 'H': gnr_houses++;
		}

		it++;
		if (last_round)	break;
		if (it == end)	last_round = true;	//We have apparently hit the last bead of the parental genome, so time for one final replication step.
	}


	if (g_length > 500)	//We only check after the full replication step, not after each replicated bead.
	{
		printf("Warning: genome sizes reached extravagant size (%d) during replication.\nExiting just to be safe...\n", g_length);
		cout << PrintContent(NULL, true, false) << endl;
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

	BeadList=new list<Bead*>();
	BeadList->splice(BeadList->begin(), *parentG->BeadList, i_split, parentG->BeadList->end());

	DevelopChildrenGenomes(parentG);
}



void Genome::DevelopChildrenGenomes(Genome* parentG)	//Function gets iterators of parental genome, but copies it to child and then acts on variables of child genome.
{
	i_bead it;
	vector<bool>* MutationList;
	int del_length, dup_length, index;
	int* pdup_length, * pdel_length;
	char wb;

	g_length = BeadList->size();

	//Clean up the variables of the parental genome left behind.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		switch ( WhatBead(*it) )
		{
			case 'R': parentG->gnr_regulators--;
			case 'T':	parentG->gnr_transporters--;
			case 'B': parentG->gnr_bsites--;
			case 'H': parentG->gnr_houses--;
		}
		parentG->g_length--;
		it++;
	}

	parentG->fork_position = 0;
	//Copy gene numbers, mutations happen next.
	gnr_regulators = parentG->gnr_regulators;
	gnr_transporters = parentG->gnr_transporters;
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
				case 'R':
					MutationList->at(index) = true;
			 		it = MutateRegulator(it, pdel_length);
				case 'T':
					MutationList->at(index) = true;
					it = MutateTransporter(it, pdel_length);
				case 'B':
					MutationList->at(index) = true;
					it = MutateBsite(it, pdel_length);
				case 'H':
					MutationList->at(index) = true;
					it = MutateHouse(it, pdel_length);
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
		//g_length is updated inside functions, gnr_regulators etc. are updated here on the outside.
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			if (duplicate)
			{
				switch ( WhatBead(*it) )
				{
					case 'R':
						it=DuplicateGene(it, pdup_length);;
						gnr_regulators++;
					case 'T':	//WARNING: Check that this now works for both 'R' and 'T'.
						it=DuplicateGene(it, pdup_length);
						gnr_transporters++;
					case 'B':
						it=DuplicateBsite(it);
						gnr_bsites++;
						(*pdup_length)++;
					case 'H':
						it=DuplicateHouse(it);
						gnr_houses++;
						(*pdup_length)++;
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
			PotentialTypeChange(it);
			(*pdup_length)++;
		}
		if(uniform() < transporter_innovation_mu)
		{
			it = InventTransporter();
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
			if(wb=='R' && uniform() < regulator_shuffle_mu)	it=ShuffleGene(it);
			else if(wb=='T' && uniform() < transporter_shuffle_mu)	it=ShuffleGene(it);
			else if(wb=='B' && uniform() < bsite_shuffle_mu)	it=ShuffleBead(it);
			else if(wb=='H' && uniform() < house_shuffle_mu) it=ShuffleBead(it);
			else	it++;
		}

	}	//END of mutations.
}



void Genome::PotentialTypeChange(i_bead it)
{
	i_bead it2;
	Regulator* reg2, reg=dynamic_cast<Regulator*>(*it);
	bool genes_are_the_same;
	list<int> UsedTypes;
	bool found_matching_type = false;
	int type_abundance = CountTypeAbundance(reg->type);

	it2 = BeadList->begin();	//The other genes in the genome
	while(it2 != BeadList->end())
	{
		if( WhatBead(*it)=='R' && it2!=it)	//We don't convert genes to themselves.
		{
			reg2 = dynamic_cast<Regulator*>(*it2);
			UsedTypes.push_back(reg2->type);

			if ( BindingAffinity(reg, reg2) == 0   &&   reg->activity == reg2->activity )
			{
				reg->type = reg2->type;		//Convert to existing gene type.
				found_matching_type = true;
				return;		//We have found a match, converted the gene; time to try mutation of the next bead.
			}
		}
		it2++;
	}
	if(found_matching_type == false && (type_abundance > 1 || reg->type==0))	//We haven't been able to convert it to an existing type, so let's define it as a new type. type_abundance should always be more than 1 because there is always an original copy on the parental section of the genome.
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
		if(uniform() < regulator_threshold_mu)	//Parameters mutate independently.
		{
			reg->threshold = ChangeGeneParameter(reg->threshold);
			is_mutated = true;
		}

		if (uniform() < regulator_activity_mu)
		{
			reg->activity = ChangeGeneParameter(reg->activity);
			potential_type_change = true;
			is_mutated = true;
		}

		for(i=0; i<sequence_length; i++)
		{
			if(uniform() < regulator_sequence_mu)
			{
				if (reg->sequence[i] == false) reg->sequence[i] = true;
				else if (reg->sequence[i] == true) reg->sequence[i] = false;
				potential_type_change = true;
				is_mutated = true;
			}
		}

		if (potential_type_change){
			PotentialTypeChange(it);	//Check type change, activity or sequence has been mutated.
		}
		it++;
	}

	return it;
}

Genome::i_bead Genome::MutateTransporter(i_bead it, int* pdel_length)
{
	//Future use, so for now I don't care about the types. Maybe later this will be useful, similar to regulators.
	Transporter* tp;
	tp = dynamic_cast<Transporter*>(*it);
	int i;

	double uu = uniform();
	if(uu < transporter_duplication_mu)
	{
		(*it)->duplicate = true;	//Mark for duplication during divison.
		is_mutated = true;
		it++;
	}

	else if(uu < transporter_deletion_mu+transporter_duplication_mu)
	{
		it = DeleteGene(it, pdel_length);
		gnr_transporters--;
		is_mutated = true;
	}

	else
	{
		if(uniform() < transporter_threshold_mu)	//Parameters mutate independently.
		{
			tp->threshold = ChangeGeneParameter(tp->threshold);
			is_mutated = true;
		}

		for(i=0; i<sequence_length; i++)
		{
			if(uniform() < transporter_sequence_mu)
			{
				if (tp->sequence[i] == false) tp->sequence[i] = true;
				else if (tp->sequence[i] == true) tp->sequence[i] = false;
				is_mutated = true;
			}
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
		it = DeleteBsite(it);
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

Genome::iter Genome::MutateHouse(i_bead it, int* pdel_length)
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
		it = DeleteHouse(it);
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
		new_value = (uniform()>0.5) ? new_value++ : new_value--;
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
	first=FindFirstBsiteInFrontOfGene(it);	//First tfbs in front of gene.
	copy_length=distance(first, last);
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	//Splice the temporary chromosome into the full genome.
	insertsite=FindRandomGenePosition(true,true);			//Find position to insert gene (may be the end of the genome).
	insertsite=FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	//Increment the number of beads and the number of genes.
	g_length+=copy_length;
	(*pdup_length)+=copy_length;

	it=last;	//Make sure ii points to one position further than just duplicated gene.
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

Genome::i_bead Genome::InventTransporter()
{
	Transporter* tp;
	i_bead insertsite;

	tp = new Transporter();
	tp->RandomTransporter();

	insertsite=FindRandomGenePosition(true,true);
	insertsite=FindFirstBsiteInFrontOfGene(insertsite);
	insertsite=BeadList->insert(insertsite, tp);

	g_length++;
	gnr_transporters++;
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
	house->RandomHouse();

	insertsite = FindRandomPosition(true);
	BeadList->insert(insertsite, house);

	g_length++;
	gnr_houses++;
	is_mutated = true;
}

Genome::i_bead Genome::ShuffleGene(i_bead it)
{
	i_bead insertsite, first, last, it2;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	is_mutated = true;
	//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last = it;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=FindFirstBsiteInFrontOfGene(it);	//First tfbs in front of gene.
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.
	last--;	//This is important, because if we move the virtual copy of the gene directly downstream of its original (i.e. it does not really move), then if last is still pointing to the bead that was originally adjacent to the gene, both the virtual copy and the original will be removed.

	//Splice the temporary chromosome into the full genome.
	insertsite=FindRandomGenePosition(true,true);			//Find position to insert gene (may be end of genome).
	insertsite=FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	last++;	//Now make last point to the next bead after the original gene again (see above), so that we will only remove the original gene with its preceding binding sites.
	//Remove the bead from its original position, taken from GeneDeletion().
	it2=first;
	while( it2 != last )
	{
		delete *it2;
		it2++;
	}
	it=(*BeadList).erase(first, last);
	return it;
}

Genome::i_bead Genome::ShuffleBead(i_bead it)
{
	i_bead tt, upstream;
	Bsite* bsite;
	bsite=dynamic_cast<Bsite*>(*it);

	is_mutated = true;
	Bsite* bsite_new = new Bsite(*bsite);	//Create copy and insert at random location.

	tt = FindRandomPosition(true);
	tt = (*BeadList).insert(tt, bsite_new);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	delete (bsite);	//Remove the old bead.
	it=(*BeadList).erase(it);

	return it;
}

Genome::i_bead Genome::FindFirstBsiteInFrontOfGene(i_bead it) const
{
	ri_bead rit(it);	//BEHIND in this case means an element to the left in the list (i.e if ii points to the 6th element, rii(ii) will make rii point to the 5th element). The function .base() as used below will make rii.base() point to the 6th element again.

	ri_bead rit2 = (*BeadList).rend();//search should be bounded
	while(rit != rit2)//begin not yet reached
	{
		if(WhatBead(*rit)=='B')	rit++;
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
	char wb;

	if(include_houses)	add_houses=gnr_houses;
	if(include_end)	end=1;

	if(gnr_regulators+gnr_transporters+add_houses==0)	return (*BeadList).end();
	else
	{
		it=(*BeadList).begin();
		while(it != (*BeadList).end())
		{
			wb = WhatBead(*it);
			if((wb=='R' || wb=='T') || (wb=='H' && include_houses))	pos.push_back(it);
			it++;
		}
		pos.push_back(it);
		randpos=(int)(uniform()*(gnr_regulators+gnr_transporters+add_houses+end));	//If you found the first gene, randpos will be 0 (no need to advance to another gene); if you find the last gene, randpos will be gnr_genes-1 which is enough to reach the gnr_genes'th gene.
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
	BeadList = new list<Bead*>();
	CopyPartOfGenome(ImageG->BeadList->begin(),ImageG->BeadList->end());

	g_length = ImageG->g_length;
	gnr_regulators = ImageG->gnr_regulators;
	gnr_transporters = ImageG->gnr_transporters;
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
	BeadList = new list<Bead*>();

	//Here the content of ReadBeadsFromString...

	fork_position = 0;
	terminus_position = g_length;
}

void Genome::ReadExpression(string expression)
{

}
