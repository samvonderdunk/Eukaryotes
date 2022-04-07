#include "Cell.hh"

Cell::Cell()
{
	kind = -1;
	barcode = -1;
}

Cell::Cell(int k)
{
	kind = k;
	barcode = -1;
}

Cell::~Cell()
{
}

Cell::i_bead Cell::TransferGene(i_bead it, Organelle* Source, Organelle* Target, bool include_distal, bool cut_and_paste)
{
	int copy_length;
	i_bead insertsite, first, last, ii;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	//If an expressed gene gets tranferred via cut-and-paste, we will have to flag the organelle, so that we can reset its gene expression (inside Population.cc).
	if (cut_and_paste)	//This has to happen before "it" is transferred and not reachable on source anymore...
	{
		Gene* gene = dynamic_cast<Gene*>(*it);
		if (gene->expression > 0)	Source->exp_gene_transfer = true;	//An expressed gene has been removed from the source genome.
	}

	last = it;		//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=Source->G->FindFirstBsiteInFrontOfGene(it, include_distal);	//First tfbs in front of gene.
	copy_length=distance(first, last);

	Source->G->CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	if (include_distal)
	{
		ii = BeadListTemp.begin();	//Don't transfer the household genes, because these will result in immediate fitness penalties.
		while (ii != BeadListTemp.end())
		{
			if ((*ii)->kind == HOUSE)
			{
				delete (*ii);
				ii = BeadListTemp.erase(ii);
				copy_length--;
			}
			else
			{
				ii++;
			}
		}
	}

	if (cut_and_paste)
	{
		ii = first;
		while (ii != last)
		{
			if ((*ii)->kind != HOUSE)	//Houses are never transferred, also not with include_distal.
			{
				Source->G->g_length--;
				Source->G->terminus_position--;
				Source->G->gnr[(*ii)->kind]--;
				delete (*ii);
				ii = Source->G->BeadList->erase(ii);
			}
			else
			{
				ii++;
			}
		}
		if (!Source->mutant)	//Cut-and-paste results in mutations during the lifetime of an organelle; these organelles (if not yet denoted as mutants) are flagged and added to the fossil record inside Population.cc.
		{
			Source->lifetime_mutant = true;
		}
	}

	//Splice the temporary chromosome into the full genome.
	insertsite=Target->G->FindRandomGenePosition(true,true);			//Find position to insert gene (may be the end of the genome).
	insertsite=Target->G->FindFirstBsiteInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	Target->G->BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	insertsite--;	//Go to the just inserted gene.
	if ((*insertsite)->kind == REGULATOR)	Target->G->PotentialTypeChange(insertsite);	//Regulator will be redefined based on the new genomic context, using the new organelle's RegTypeList.

	//Increment the number of beads and the number of genes.
	Target->G->g_length+=copy_length;
	if ((*insertsite)->kind == REGULATOR)			Target->G->gnr[REGULATOR]++;
	else if((*insertsite)->kind == EFFECTOR)	Target->G->gnr[EFFECTOR]++;
	Target->G->gnr[BSITE]+=copy_length-1;	//The length of the whole transferred piece except for the gene (i.e. you will always transfer 1 gene with x bsites and nothing else).
	Target->G->is_mutated = true;
	Target->mutant = true;
	Target->G->terminus_position = Target->G->g_length;
	return last;
}

void Cell::TransferBead(i_bead it, Organelle* Target)
{
	Bead* bead = (*it)->Clone();
	i_bead insertsite = Target->G->FindRandomPosition(true);
	Target->G->BeadList->insert(insertsite, bead);
	Target->G->g_length++;
	Target->G->terminus_position = Target->G->g_length;
	switch ( bead->kind )
	{
		case BSITE:
			Target->G->gnr[BSITE]++;
			break;
		case HOUSE:
			Target->G->gnr[HOUSE]++;
			break;
	}
	Target->G->is_mutated = true;
	Target->mutant = true;
}

string Cell::Show()
{
	string Content="";
	std::stringstream ss;

	ss << "<" << kind << "|" << barcode << ">" << endl;

	Content += ss.str();
	ss.clear();

	return Content;
}
