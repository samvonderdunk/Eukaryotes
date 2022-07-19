#include "Genome.hh"

Genome::Genome()
{
	BeadList=NULL;
	BeadList = new list<Bead*>();
	RegTypeList = {NULL};
	ExpressedGenes = NULL;
	ExpressedGenes = new list<Regulator*>();
	ExpressedGenesN = NULL;
	ExpressedGenesN = new list<Regulator*>();
	g_length=0;
	for (int i=0; i<4; i++)	gnr[i]=0;
	fork_position=0;
	terminus_position=0;
	nr_native_expressed=0;
	is_mutated=false;
	organelle=0;
}

Genome::~Genome()
{
	i_bead it;
	int i;

	it=BeadList->begin();
	while(it!=BeadList->end())
	{
		delete (*it);
		it++;
	}

	BeadList->clear();
	delete BeadList;
	BeadList=NULL;

	for (i=0; i<5; i++)
	{
		delete RegTypeList[i];
	}

	if (ExpressedGenes != NULL)	ExpressedGenes->clear();
	delete ExpressedGenes;
	ExpressedGenes = NULL;

	if (ExpressedGenesN != NULL)	ExpressedGenesN->clear();
	delete ExpressedGenesN;
	ExpressedGenesN = NULL;
}



void Genome::UpdateGeneExpression()
{
	i_bead it;
	int it_cntr, kind;
	double cum_effects = 0.;
	Bsite* bs;
	Regulator* reg;
	std::list<Regulator*>* PointerEG; //Temporary pointer for swapping old and new expression lists.

	//Determine regulatory dynamics and put result into the express variable of each gene (Gene or otherwise...).
	it = BeadList->begin();
	it_cntr = 0;
	while (it != BeadList->end())
	{
		kind = (*it)->kind;
		if (kind==REGULATOR)	//i.e. some kind of gene (REGULATOR or EFFECTOR).
		{
			reg = dynamic_cast<Regulator*>(*it);
			cum_effects -= (double)reg->threshold;
			reg->express = max(min((int)cum_effects+1,1),0);	//For the +1, see explanation of the gene threshold in Gene.hh
			reg->expression = reg->express;	//Store the expression.
			while (reg->express > 0)
			{
				ExpressedGenesN->push_back(reg);
				reg->express--;
			}
			cum_effects = 0.;
		}
		else if (kind==BSITE)
		{
			bs = dynamic_cast<Bsite*>(*it);
			reg = RegulatorCompetition(bs);
			if (reg != NULL)		cum_effects += bs->activity * reg->activity;
		}
		it++;
		it_cntr++;
		if (it_cntr == terminus_position)	cum_effects = 0.;
	}

	//Realise expression states by swapping pointers of old and new expression lists.
	ExpressedGenes->clear();
	PointerEG = ExpressedGenes;	//Pointer swapping using temporary pointer.
	ExpressedGenes = ExpressedGenesN;
	ExpressedGenesN = PointerEG;
}



void Genome::NativeExpression()
{
	i_bead it;
	Regulator* reg;

	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if ((*it)->kind==REGULATOR)
		{
			//See similar potential issue in UpdateExpression() and in BindingAffinity().
			reg = dynamic_cast<Regulator*>(*it);
			reg->express = reg->expression;	//Only difference with last paragraph of UpdateGeneExpression is that we here first have to read off gene->expression (which has been read from input files or which is stored in the just-divided genomes)
			while (reg->express > 0){
				ExpressedGenes->push_back(reg);
				reg->express--;
			}
		}
		it++;
	}
}



Regulator* Genome::RegulatorCompetition(Bsite* bsite)
{
	i_reg ir;
	double affinity, p_bind, z_partition = 1.;

	//Calculate partition function.
	ir = ExpressedGenes->begin();
	while (ir != ExpressedGenes->end())
	{
		affinity = (double)(*ir)->BindingAffinity(bsite->sequence, (*ir)->sequence);
		z_partition += 1 * k_zero * exp(affinity * epsilon);	//Exp===1, higher expression means that more copies will be present in ExpressedGenes.
		ir++;
	}

	//Pick one of the probabilities. BindingAffinity has been optimised so much, that this algorithm outperforms a simple one where we store the affinities in a list in the above iteration, and look them up here.
	double die_roll = uniform();
	die_roll -= 1 / z_partition;
	if (die_roll <= 0.)	return NULL;

	ir = ExpressedGenes->begin();
	while (ir != ExpressedGenes->end())
	{
		affinity = (double)(*ir)->BindingAffinity(bsite->sequence, (*ir)->sequence);
		p_bind = ( 1 * k_zero * exp(affinity * epsilon) ) / z_partition;	//Exp===1, see above.
		die_roll -= p_bind;
		if (die_roll <= 0.)	return *ir;
		ir++;
	}

	printf("Error: partition function out of bounds.\n");
	exit(1);
}



void Genome::ReplicateStep(double resource)
{
	i_bead it, start, end;
	Bead* bead;
	Regulator* reg;
	int replicate_beads = 0, repl_remaining_steps;
	double res_int, res_fract, fract_repl_remaining;
	bool last_round = false;

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
	while (fork_position+replicate_beads < terminus_position && repl_remaining_steps > 0)	//The maximal position of end is defined by pos_anti_ori - 1 (pointing to the last bead of the parental genome). pos_anti_ori holds the number of genes in the parental genome, so if your distance to the first bead is pos_anti_ori, you are actually one past the last bead of the parental genome. The first replication step, this will point to NULL, but in consecutive steps it will point to a child bead; we want to point to an end point that does not change, hence the last bead of the parental genome.
	{
		if ((*end)->kind>=REGULATOR)		repl_remaining_steps--;	//Genes for sure count for the repl_step_size.
		else if(!gene_replication)			repl_remaining_steps--;	//TFBSs count if replicate_entire_genes is set to false.
		replicate_beads++;
		end++;
	}
	end--;

	it = start;	//Loop over a number of beads (how many we will replicate in one step).
	if(start == end)	last_round = true;	//We go straight into the last round.
	while (it != BeadList->end())
	{
		//Copy bead.
		bead=(*it)->Clone();
		(*BeadList).push_back(bead);
		g_length++;
		gnr[(*it)->kind]++;

		//Copy expression.
		if (bead->kind==REGULATOR)	//If expressed gene is replicated, the new copy is also expressed.
		{
			reg = dynamic_cast<Regulator*>(bead);
			if (reg->expression > 0)
			{
				reg->express = reg->expression;
				while (reg->express > 0)
				{
					ExpressedGenes->push_back(reg);
					reg->express--;
				}
			}
		}

		it++;
		if (last_round)	break;
		if (it == end)	last_round = true;	//We have apparently hit the last bead of the parental genome, so time for one final replication step.
	}


	if (g_length > max_g_length)	//We only check after the full replication step, not after each replicated bead.
	{
		printf("Warning: genome sizes reached extravagant size (%d) during replication.\nExiting just to be safe...\n", g_length);
		cout << Show(NULL, true, false) << endl;
		exit(1);
	}

	fork_position += replicate_beads;	//Move the fork to the right.
	if (fork_position >= terminus_position)	//If this is the final replication step.
	{
		fork_position = terminus_position;
		assert(g_length == 2*terminus_position);
	}

}



void Genome::SplitGenome(Genome* parentG)	//Used to split a genome upon division
{
	int i;
	i_bead it, i_split;

	//Find the fork with i_split.
	i_split = parentG->BeadList->begin();
	advance(i_split, parentG->terminus_position);	//terminus_position points to the end of the parental genome, which is now the first bead of the child genome.

	BeadList->splice(BeadList->begin(), *parentG->BeadList, i_split, parentG->BeadList->end());

	for(i=0; i<5; i++)		RegTypeList[i] = new Regulator(*parentG->RegTypeList[i]);

	DevelopChildrenGenomes(parentG);	//including mutations.
}



void Genome::AbortChildGenome()
{
	i_bead it = BeadList->begin();
	advance(it, terminus_position);

	while (it != BeadList->end())
	{
		gnr[(*it)->kind]--;
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
	int i, del_length, dup_length, index, g_length_before_mut;
	int* pdup_length, * pdel_length;

	g_length = BeadList->size();
	g_length_before_mut = g_length;

	//Clean up the variables of the parental genome left behind.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		parentG->gnr[(*it)->kind]--;
		parentG->g_length--;
		it++;
	}

	parentG->fork_position = 0;
	//Copy gene numbers, mutations happen next.
	for (i=0; i<4; i++)	gnr[i] = parentG->gnr[i];

	if (mutations_on)	//START mutations.
	{
		del_length = 0;
		dup_length = 0;
		pdup_length = &dup_length;
		pdel_length = &del_length;
		index = 0;

		//Mutate beads and set gnr_genes.
		MutationList = new vector<bool>(g_length, false);
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			MutationList->at(index) = true;
			it = Mutation(it, pdel_length);
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
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			if ((*it)->duplicate)
			{
				it = Duplication(it, pdup_length);
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

		//Shuffling mutations.
		it = BeadList->begin();
		while(it != BeadList->end())
		{
			if (uniform() < mu[SHUFFLE][(*it)->kind])	it = Shuffle(it);
			else	it++;
		}

		//Innovations.
		Inventions(pdup_length);

		//WGDs.
		if (uniform() < muWGD)		WholeGenomeDuplication(pdup_length);

		assert(g_length == g_length_before_mut + (*pdup_length) - (*pdel_length));
	}	//END of mutations.

	assert((size_t)g_length == BeadList->size());
	CheckBeadCounts();
	terminus_position = g_length;
}


void Genome::PotentialTypeChange(i_bead it)	//Meant purely for regulatory genes
{
	i_bead it2;
	Regulator* reg, *reg2;
	list<int> UsedTypes;
	list<int>::iterator ix;
	int type_abundance, x, i;

	reg = dynamic_cast<Regulator*>(*it);
	type_abundance = CountTypeAbundance(reg->type);

	//Check for convergence to types 1-5, as they are defined before the start of mutations and redefined as we go on mutating.
	for (i=0; i<5; i++)
	{
		if ( reg->BindingAffinity(reg->sequence, RegTypeList[i]->sequence) + abs(reg->activity-RegTypeList[i]->activity) <= seq_hdist )
		{
			reg->type = i+1;
			return;
		}
	}

	//If there is at least one more copy of the gene, assign a new unique type to the gene, 6 at minimum.
	if(type_abundance > 1 || reg->type==0)	//We haven't been able to convert it to an existing type, so let's define it as a new type.
	{
		//Check which gene types we already have, if you want gene types >5 to be at least evolutionarily meaningful.
		it2 = BeadList->begin();
		while (it2 != BeadList->end())
		{
			if( (*it2)->kind==REGULATOR)
			{
				reg2 = dynamic_cast<Regulator*>(*it2);
				UsedTypes.push_back(reg2->type);
			}
			it2++;
		}

		x=6;
		while (ix != UsedTypes.end())
		{
			x++;
			ix = find(UsedTypes.begin(), UsedTypes.end(), x);
		}
		reg->type = x;
	}

	//If there exists only copy of this type, then type is not changed, but the type definition has to be updated (i.e. replaced).
	else if(reg->type < 6)	//We are only concerned with the defs of types 1-5.
	{
		delete RegTypeList[reg->type-1];
		RegTypeList[reg->type-1] = new Regulator(*reg);
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
		if ((*it)->kind==REGULATOR)
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


Genome::i_bead Genome::Mutation(i_bead it, int* pdel_length)
{
	is_mutated = false;

	double uu = uniform();
	if (uu < mu[DUPLICATION][(*it)->kind])			//Duplication.
	{
		(*it)->duplicate = true;
		is_mutated = true;
		it++;
	}

	else if (uu < mu[DUPLICATION][(*it)->kind] + mu[DELETION][(*it)->kind])		//Deletion.
	{
		it = Deletion(it, pdel_length);
		is_mutated = true;
	}

	else																						//Property mutations.
	{
		if ( (*it)->Mutate(organelle) )
		{
			if ( (*it)->kind == REGULATOR && mu[TYPE][(*it)->kind]==0.0 )	PotentialTypeChange(it);	//Check type change, maybe activity or sequence has been mutated.	Note that for effectors, the type change is immediately checked by the Mutate function, as it does not need to know other variables of the Genome (for regulatory types, we do need those).
			is_mutated = true;
		}
		it++;
	}
	return it;
}

Genome::i_bead Genome::Deletion(i_bead it, int* pdel_length)
{
	i_bead first, last, it2;
	int del_length;

	//Define the segment to be deleted.
	last = it;
	last++;			//Move last to one bead beyond the focal bead.

	switch ( (*it)->kind )
	{
		case HOUSE:
			first = it;
			gnr[HOUSE]--;
			break;
		case BSITE:
			first = it;
			gnr[BSITE]--;
			break;
		case REGULATOR:
			first = FindFirstBsiteInFrontOfGene(it);
			gnr[REGULATOR]--;
			break;
		case EFFECTOR:
			first = FindFirstBsiteInFrontOfGene(it);
			gnr[EFFECTOR]--;
			break;
	}

	//Adjust genome counters.
	del_length = distance(first, last);		//Store it in a variable so that we don't have to call the function multiple times.
	(*pdel_length) += del_length;
	g_length -= del_length;
	gnr[BSITE] -= del_length-1;		//If multiple beads are deleted, our focal bead (it) is a regulator or effector, and any surplus beads must be its binding sites.

	//Delete the selected genome segment.
	it2=first;	//Erase all the beads selected for the deletion (potentially multiple if focal bead was a regulator or effector).
	while( it2 != last )
	{
		delete *it2;
		it2++;
	}
	it=(*BeadList).erase(first, last);
	return it;
}

Genome::i_bead Genome::Duplication(i_bead it, int* pdup_length)
{
	int copy_length;
	i_bead insertsite, first, last;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	(*it)->duplicate = false;	//Remove the duplication flag.

	//Define the segment to be duplicated.
	last = it;
	last++;

	switch ( (*it)->kind )	//Whether we select a single bead or a chunk, and where we insert the single bead or chunk depends on the kind of bead.
	{
		case HOUSE:
			first = it;
			insertsite = FindRandomPosition(true);
			gnr[HOUSE]++;
			break;
		case BSITE:
			first = it;
			insertsite = FindRandomPosition(true);
			gnr[BSITE]++;
			break;
		case REGULATOR:
			first = FindFirstBsiteInFrontOfGene(it);
			insertsite = FindRandomGenePosition(true,true);	//Find position to insert gene (may be the end of the genome).
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
			gnr[REGULATOR]++;
			break;
		case EFFECTOR:
			first = FindFirstBsiteInFrontOfGene(it);
			insertsite = FindRandomGenePosition(true,true);
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);
			gnr[EFFECTOR]++;
			break;
	}

	//Adjust genome counters.
	copy_length = distance(first, last);
	g_length += copy_length;
	gnr[BSITE] += copy_length-1;
	(*pdup_length) += copy_length;

	//Duplicate the selected genome segment.
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a new piece of DNA with the duplicated beads on it.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.
	it = last;	//Make sure ii points to one position further than just duplicated gene.
	return it;
}

void Genome::Inventions(int* pdup_length)
{
	i_bead focus, insertsite;
	int k;
	double uu;
	House* house;
	Bsite* bsite;
	Regulator* reg;
	Effector* eff;

	for (k=0; k<4; k++)	//One potential invention per bead kind.
	{
		uu = uniform();

		if (k == HOUSE && uu < mu[INVENTION][HOUSE])
		{
			house = new House();
			house->Randomize(organelle);

			insertsite = FindRandomPosition(true);
			insertsite = BeadList->insert(insertsite, house);

			gnr[HOUSE]++;
			g_length++;
			(*pdup_length)++;
		}
		else if (k == BSITE && uu < mu[INVENTION][BSITE])
		{
			bsite = new Bsite();
			bsite->Randomize(organelle);

			insertsite = FindRandomPosition(true);
			insertsite = BeadList->insert(insertsite, bsite);

			gnr[BSITE]++;
			g_length++;
			(*pdup_length)++;
		}
		else if (k == REGULATOR && uu < mu[INVENTION][REGULATOR])
		{
			reg = new Regulator();
			reg->Randomize(organelle);

			insertsite = FindRandomGenePosition(true,true);
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);
			insertsite = BeadList->insert(insertsite, reg);
			if (mu[TYPE][REGULATOR]==0.0)	PotentialTypeChange(insertsite);
			else																			reg->MutateType(&reg->type, 1.0);	//m=1, so that "mutation" always takes place.

			gnr[REGULATOR]++;
			g_length++;
			(*pdup_length)++;
		}
		else if (k == EFFECTOR && uu < mu[INVENTION][EFFECTOR])
		{
			eff = new Effector();
			eff->Randomize(organelle);

			insertsite = FindRandomGenePosition(true,true);
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);
			insertsite = BeadList->insert(insertsite, eff);

			eff->DefineTypeFromSeq();

			gnr[EFFECTOR]++;
			g_length++;
			(*pdup_length)++;
		}
	}
}

Genome::i_bead Genome::Shuffle(i_bead it)
{
	i_bead insertsite, first, last;

	//Define the segment to be shuffled.
	last = it;
	last++;

	switch ( (*it)->kind )
	{
		case HOUSE:
			first = it;
			insertsite = FindRandomPosition(true);
			break;
		case BSITE:
			first = it;
			insertsite = FindRandomPosition(true);
			break;
		case REGULATOR:
			first = FindFirstBsiteInFrontOfGene(it);
			insertsite = FindRandomGenePosition(true,true);	//Find position to insert gene (may be the end of the genome).
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
			break;
		case EFFECTOR:
			first = FindFirstBsiteInFrontOfGene(it);
			insertsite = FindRandomGenePosition(true,true);
			insertsite = FindFirstBsiteInFrontOfGene(insertsite);
			break;
	}

	//Shuffle the selected genome segment.
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

void Genome::WholeGenomeDuplication(int* pdup_length)
{
	int k;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	//Update counters.
	for (k=0; k<4; k++){	gnr[k] *= 2;	}
	(*pdup_length) += g_length;
	g_length *= 2;

	//Do the actual copying.
	CopyPartOfGenomeToTemplate(BeadList->begin(), BeadList->end(), &BeadListTemp); //Makes a new piece of DNA with the duplicated beads on it.
	BeadList->splice(BeadList->end(), BeadListTemp);	//Splice temporary list into chromosome.
}

void Genome::CheckBeadCounts()	//Useful for debugging.
{
	int k;
	i_bead it;
	std::array<int,4> counts = {0};

	it = BeadList->begin();
	while (it != BeadList->end())
	{
		counts[(*it)->kind]++;
		it++;
	}

	for (k=0; k<4; k++)	assert(gnr[k] == counts[k]);
	assert(gnr[0]+gnr[1]+gnr[2]+gnr[3] == g_length);
}

Genome::i_bead Genome::FindFirstBsiteInFrontOfGene(i_bead it, bool ignore_houses) const
{
	ri_bead rit(it);	//BEHIND in this case means an element to the left in the list (i.e if ii points to the 6th element, rii(ii) will make rii point to the 5th element). The function .base() as used below will make rii.base() point to the 6th element again.

	ri_bead rit2 = (*BeadList).rend();//search should be bounded
	while(rit != rit2)//begin not yet reached
	{
		if((*rit)->kind==BSITE || ((*rit)->kind==HOUSE && ignore_houses))	rit++;
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

	if(include_houses)	add_houses=gnr[HOUSE];
	if(include_end)	end=1;

	if(gnr[REGULATOR]+gnr[EFFECTOR]+add_houses==0)	return (*BeadList).end();
	else
	{
		it=(*BeadList).begin();
		while(it != (*BeadList).end())
		{
			if((*it)->kind==EFFECTOR || (*it)->kind==REGULATOR || ((*it)->kind==HOUSE && include_houses))	pos.push_back(it);
			it++;
		}
		pos.push_back(it);
		randpos=(int)(uniform()*(gnr[REGULATOR]+gnr[EFFECTOR]+add_houses+end));	//If you found the first gene, randpos will be 0 (no need to advance to another gene); if you find the last gene, randpos will be gnr_genes-1 which is enough to reach the gnr_genes'th gene.
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
	int i;

	CopyPartOfGenome(ImageG->BeadList->begin(),ImageG->BeadList->end());
	for (i=0; i<5; i++)	RegTypeList[i] = new Regulator(*ImageG->RegTypeList[i]);

	g_length = ImageG->g_length;
	for (i=0; i<4; i++)	gnr[i] = ImageG->gnr[i];
	fork_position = ImageG->fork_position;
	terminus_position = ImageG->terminus_position;
	is_mutated = ImageG->is_mutated;
	organelle = ImageG->organelle;
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

	char* bead, *rseq_string, *eseq_string;
	int success, type, threshold, activity, i;
	std::bitset<regulator_length> RsequenceBITS;
	std::bitset<effector_length> EsequenceBITS;
	House* house;
	Bsite* bsite;
	Regulator* reg;
	Effector* eff;

	bead = strtok((char*)genome.c_str(),".");
	while (bead != NULL)
	{
		if(bead[1] == 'R')
		{
			rseq_string = new char[regulator_length+1];
			success = sscanf(bead, "(R%d:%d:%d:%[^)])", &type, &threshold, &activity, rseq_string);
			if(success != 4) cerr << "Could not read regulatory gene. Input seems corrupt. \n" << endl;

			for (i=0; i<regulator_length; i++)	RsequenceBITS[i] = (rseq_string[i]=='1');
			delete [] rseq_string;
			rseq_string = NULL;

			reg = new Regulator(type, threshold, activity, RsequenceBITS, 0);
			(*BeadList).push_back(reg);
			gnr[REGULATOR]++;
			g_length++;
		}
		else if (bead[1] == 'E')
		{
			eseq_string = new char[effector_length+1];
			success = sscanf(bead, "(E%d:%d:%[^)])", &type, &threshold, eseq_string);
			if(success != 4) cerr << "Could not read effector gene. Input seems corrupt. \n" << endl;

			for (i=0; i<effector_length; i++)		EsequenceBITS[i] = (eseq_string[i]=='1');
			delete [] eseq_string;
			eseq_string = NULL;

			eff = new Effector(type, threshold, EsequenceBITS, 0);
			(*BeadList).push_back(eff);
			gnr[EFFECTOR]++;
			g_length++;
		}
		else if(bead[1] == 'H')
		{
			house = new House();
			(*BeadList).push_back(house);
			gnr[HOUSE]++;
			g_length++;
		}
		else
		{
			rseq_string = new char[regulator_length+1];
			success = sscanf(bead, "(%d:%[^)])", &activity, rseq_string);
			if(success != 2) cerr << "Could not read binding site. Input seems corrupt. \n" << endl;

			for (i=0; i<regulator_length; i++)	RsequenceBITS[i] = (rseq_string[i]=='1');

			delete [] rseq_string;
			rseq_string = NULL;

			bsite = new Bsite(activity, RsequenceBITS);
			(*BeadList).push_back(bsite);
			gnr[BSITE]++;
			g_length++;
		}

		bead = strtok(NULL, ".");
	}

	fork_position = 0;
	terminus_position = g_length;
}

void Genome::ReadExpression(string expression)
{
	string data;
	int i;
	i_bead it;
	Gene* gene;

	data = expression.substr(1,expression.size()-2);	//Strip accolades or square brackets.

	it = BeadList->begin();
	for(i=0; i<gnr[REGULATOR]+gnr[EFFECTOR]; i++)
	{
		while (it != BeadList->end() && ((*it)->kind != REGULATOR && (*it)->kind != EFFECTOR)) it++;

		if ((size_t) i > data.length())
		{
			cerr << "Could not read expression states. Input seems corrupt.\n" << endl;
			exit(1);
		}

		gene = dynamic_cast<Gene*>(*it);
		gene->expression = (data[i]=='1');
		it++;
	}
}

void Genome::ReadDefinition(string definition)
{
	char* bead, *rseq_string;
	int i, success, activity, type;
	std::bitset<regulator_length> sequenceBITS;

	bead = strtok((char*)definition.c_str(),";");
	while (bead != NULL)
	{
		rseq_string = new char[regulator_length+1];
		success = sscanf(bead, "(R%d:%d:%[^)])", &type, &activity, rseq_string);
		if(success != 3) cerr << "Could not read regulatory type definitions. Input seems corrupt. \n" << endl;

		for (i=0; i<regulator_length; i++)	sequenceBITS[i] = (rseq_string[i]=='1');

		delete [] rseq_string;
		rseq_string = NULL;

		RegTypeList[type-1] = new Regulator(type, 0, activity, sequenceBITS, 0);

		bead = strtok(NULL, ";");
	}
}

string Genome::Show(list<Bead*>* chromosome, bool terminal, bool only_parent)
{
	string GenomeContent="";
	i_bead it, end;

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
		GenomeContent += (*it)->Show(terminal);
		it++;
	}
	return GenomeContent;
}

string Genome::ShowExpression(list<Bead*>* chromosome, bool only_parent)
{
	i_bead it, end;
	Gene* gene;
	string Content="{";

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
		if ((*it)->kind==REGULATOR || (*it)->kind==EFFECTOR)
		{
			std::stringstream ss;
			gene = dynamic_cast<Gene*>(*it);
			ss << gene->expression;
			Content += ss.str();
			ss.clear();
		}
		it++;
	}

	Content += "}";
	return Content;
}

string Genome::ShowDefinition(bool terminal)
{
	int i;
	string Content="";

	for (i=0; i<5; i++)
	{
		Content += RegTypeList[i]->Show(terminal, true);
		if(i<4)	Content += ";";
	}

	return Content;
}
