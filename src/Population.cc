#include "Population.hh"

Population::Population()
{
	int i,j;
	cell_count=0;
	organelle_count=0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		Space[i][j]=NULL;
	}
}

Population::~Population()
{
	int i,j;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if((Space[i][j])!=NULL)
		{
			delete (Space[i][j]);
			Space[i][j]=NULL;
		}
	}
}
