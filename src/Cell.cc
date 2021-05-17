#include "Cell.hh"

Cell::Cell()
{
	Host=NULL;
	for(int i=0;i<HS;i++)
	{
		Symbionts[i]=NULL;
	}
}

Cell::~Cell()
{
	delete Host;
	Host=NULL;
	for(int i=0;i<HS;i++)
	{
		if (Symbionts[i]!=NULL)
		{
			delete (Symbionts[i]);
			Symbionts[i]=NULL;
		}
	}
}

void Cell::UpdateOrganelles()
{
	int i;

	Host->UpdateExpression();
	for (i=0;i<HS;i++)	Symbionts[i]->UpdateExpression();

	RegulatorTransport();

	Host->UpdateState();
	for (i=0;i<HS;i++)	Symbionts[i]->UpdateState();

	//RegulatorTransport();	//You can do it here if leaking proteins should not affect cell-cycle stage (or before UpdateExpression(), which would be identical).
}

void Cell::RegulatorTransport()
{
	i_bead it;
	int nr_native_host_genes;

	nr_native_host_genes = (int) Host->ExpressedGenes->size();
	for (i=0;i<HS;i++)
	{
		//Movement of expressed regulators from symbionts to host.
		it = Symbionts[i]->ExpressedGenes->begin();
		while (it != Symbionts[i]->ExpressedGenes->end())
		{
			if (ActiveTransport(it, Symbionts[i]->ExpressedGenes, Host->ExpressedGenes) || uniform() < leakage_to_host)
			{
				Host->ExpressedGenes->push_back(*it);
			}
			it++;
		}

		//Movement of expressed regulators from host to symbionts.
		it = Host->ExpressedGenes->begin();
		while (distance(Host->ExpressedGenes->begin(), it) < nr_native_host_genes)
		{
			if (ActiveTransport(it, Host->ExpressedGenes, Symbionts[i]->ExpressedGenes) || uniform() < leakage_to_symbiont)
			{
				Symbionts[i]->ExpressedGenes->push_back(*it);
			}
			it++;
		}
	}
}

bool Cell::ActiveTransport(i_bead it, list<Bead*>* SourceCompartment, list<Bead*>* TargetCompartment)
{
	//For future use: if transporters are expressed, they can here lead to very likely transport of genes.
	//Maybe passing BeadLists is a bit cumbersome as well...
	return false;	//For now, there will be no active transport.
}
