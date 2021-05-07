#include "Organelle.hh"

Organelle::Organelle()
{
	G=NULL;
}

Organelle::~Organelle()
{
	if (G!=NULL)
	{
		delete (G);
		G=NULL;
	}
}
