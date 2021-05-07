#include "Fossils.hh"

Fossils::Fossils()
{
}

Fossils::~Fossils()
{
	i_fos it;
	it = FossilRecord.begin();
	while (it != FossilRecord.end())
	{
		if (*it != NULL)
		{
			delete (*it);
			it = FossilRecord.erase(it);
		}
		else	it++;
	}
}
