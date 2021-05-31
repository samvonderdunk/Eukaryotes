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

  int g_length, gnr_regulators, gnr_bsites, gnr_houses;
  int fork_position, terminus_position;	//position of the replication fork and the terminus, where it stops.
  bool is_mutated;

  Genome();
  ~Genome();

	void UpdateExpression(list<Bead*>* ExpressedGenes);
	void EraseExpression(list<Bead*>* ExpressedGenes);
	void SetExpression(list<Bead*>* ExpressedGenes, bool Updating);
	i_bead RegulatorCompetition(i_bead i_bsite, list<Bead*>* ExpressedGenes);

	void ReplicateStep(double resource);

	void SplitGenome(Genome* parentG);
	void DevelopChildrenGenomes(Genome* parentG);
	void PotentialTypeChange(i_bead ii);

	//Mutation functions.
	i_bead MutateRegulator(i_bead it, int* pdel_length);
	i_bead MutateBsite(i_bead it, int* pdel_length);
	i_bead MutateHouse(i_bead it, int* pdel_length);

	int ChangeParameter(int value);

	i_bead DeleteGene(i_bead it, int* pdel_length);
	i_bead DeleteBead(i_bead it);

	i_bead DuplicateGene(i_bead it, int* pdup_length);
	i_bead DuplicateBsite(i_bead it);
	i_bead DuplicateHouse(i_bead it);

	i_bead InventRegulator();
	void InventBsite();
	void InventHouse();

	i_bead ShuffleGene(i_bead it);
	i_bead ShuffleBead(i_bead it);

	i_bead FindFirstBsiteInFrontOfGene(i_bead it) const;
	i_bead FindRandomGenePosition(bool include_houses, bool include_end) const;
	i_bead FindRandomPosition(bool include_end) const;
	void CopyPartOfGenomeToTemplate(i_bead begin, i_bead end, list<Bead*>* template_beadlist);
	void CloneGenome(const Genome* ImageG);
	void CopyPartOfGenome(i_bead begin, i_bead end);

	void ReadGenome(string genome);
	void ReadExpression(string expression);
	void ReadBuffer(string buffer, bool* array, char stop_sign);


	inline char WhatBead(Bead* bead) const
	//Testing; if this is not working, go back to seperate check functions and/or without inline.
	//Advantage of this is that I now have a single function, possible disadvantage may be slower speed.
	{
		switch ( (int)(typeid(*bead) )
		{
			case typeid(Regulator):
				return 'R'
			case typeid(Bsite):
				return 'B'
			case typeid(House):
				return 'H'
		}
	}

	inline int BindingAffinity(Bead* bead1, Bead* bead2)
	{
	//Testing; same as above. In addition, look at faster ways to do bitstring comparisons (or storage of bitstrings).
	//If this flexible b1/b2 calling does not work, perhaps try with passing the sequence directly? If it is inline, maybe this does not cost anything extra.
	//In addition, perhaps it could be useful to only store sequences as integers, but unpack them as their true binary strings only in this function...
		switch (WhatBead(bead1))
		{
			case 'R':
				Regulator* b1 = dynamic_cast<Regulator*>(bead1);
			case 'B':
				Bsite* b1 = dynamic_cast<Bsite*>(bead1);
		}
		switch (WhatBead(bead2))
		{
			case 'R':
				Regulator* b2 = dynamic_cast<Regulator*>(bead2);
			case 'B':
				Bsite* b2 = dynamic_cast<Bsite*>(bead2);
		}

		int affinity = 0;
		for (int i=0;i<sequence_length;i++)
		{
			if (b1->sequence[i] != b2->sequence[i])	affinity++;
		}
		return affinity;
	}

};
#endif
