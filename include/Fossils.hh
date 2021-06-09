#ifndef FossilHeader
#define FossilHeader

#include "Header.hh"
#include "Organelle.hh"

class Fossils
{
	public:
	std::list<Organelle*> FossilRecord;

	typedef std::list<Organelle*>::iterator i_fos;

	Fossils();
	~Fossils();

	void EraseFossil(unsigned long long fossilID);
	void BuryFossil(Organelle* O);
	void ExhibitFossils();

};
#endif
