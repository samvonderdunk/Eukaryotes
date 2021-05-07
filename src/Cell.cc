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
