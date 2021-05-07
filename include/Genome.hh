#ifndef GenomeHeader
#define GenomeHeader

#include "Header.hh"
#include "Bead.hh"
#include "Regulator.hh"
#include "Bsite.hh"
#include "Transporter.hh"
#include "House.hh"
#include <typeinfo>

class Genome {
 public:
  std::list<Bead*>* BeadList;
  typedef std::list<Bead*>::iterator i_bead;
  typedef std::list<Bead*>::reverse_iterator ri_bead;

  int g_length;
  int gnr_regulators;
	int gnr_bsites;
	int gnr_transporters;
  int gnr_houses;
  int fork_position, terminus_position;	//position of the replication fork and the terminus, where it stops.
  bool is_mutated;

  Genome();
  ~Genome();

	void UpdateExpression(list<Bead*>* ExpressedGenes);
	void EraseExpressedGenes(list<Bead*>* ExpressedGenes);
	i_bead RegulatorCompetition(i_bead i_bsite, list<Bead*>* ExpressedGenes);




	inline char WhatBead(Bead* bead) const
	//Testing; if this is not working, go back to seperate check functions and/or without inline.
	//Advantage of this is that I now have a single function, possible disadvantage may be slower speed.
	{
		switch ( (int)(typeid(*bead) )
		{
			case typeid(Regulator):
				return 'R'
			case typeid(Transporter):
				return 'T'
			case typeid(Bsite):
				return 'B'
			case typeid(House):
				return 'H'
		}
	}

	inline int BindingAffinity(Bead* i_bsite, Bead* i_reg)
	{
	//Testing; same as above. In addition, look at faster ways to do bitstring comparisons (or storage of bitstrings).
		Bsite* bsite = dynamic_cast<Bsite*>(*i_bsite);
		Regulator* reg = dynamic_cast<Regulator*>(*i_reg);

		int affinity = 0;
		for (int i=0;i<sequence_length;i++)
		{
			if (bsite->sequence[i] != reg->sequence[i])	affinity++;
		}
		return affinity;
	}

};
#endif
