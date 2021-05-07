#ifndef FossilHeader
#define FossilHeader

#include "Header.hh"
#include "Cell.hh"

class Fossils
{
	public:
	std::list<Organelle*> FossilRecord;

	typedef std::list<Organelle*>::iterator i_fos;
	Fossils();
	~Fossils();

};
#endif
