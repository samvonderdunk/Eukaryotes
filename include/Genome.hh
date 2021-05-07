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
  typedef std::list<Bead*>::iterator iter;
  typedef std::list<Bead*>::reverse_iterator reviter;

  int g_length;
  int gnr_regulators;
	int gnr_bsites;
	int gnr_transporters;
  int gnr_houses;
  int fork_position, terminus_position;	//position of the replication fork and the terminus, where it stops.
  bool is_mutated;

  Genome();
  ~Genome();

};
#endif
