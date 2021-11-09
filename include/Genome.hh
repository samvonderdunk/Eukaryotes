#ifndef GenomeHeader
#define GenomeHeader

#include "Header.hh"
#include "Bead.hh"
#include "Gene.hh"
#include "Bsite.hh"
#include "House.hh"

class Genome {
 public:
  std::list<Bead*>* 		BeadList;
	std::array<Gene*,5>		RegTypeList;

  int g_length, gnr_genes, gnr_bsites, gnr_houses;
  int fork_position, terminus_position;	//position of the replication fork and the terminus, where it stops.
  bool is_mutated;
	bool is_symbiont;

	typedef std::list<Bead*>::iterator i_bead;
	typedef std::list<Bead*>::reverse_iterator ri_bead;

  Genome();
  ~Genome();

	void UpdateGeneExpression(list<Bead*>* ExpressedGenes);
	void NativeExpression(list<Bead*>* ExpressedGenes);
	i_bead GeneCompetition(i_bead i_bsite, list<Bead*>* ExpressedGenes);

	void ReplicateStep(double resource);

	void SplitGenome(Genome* parentG);
	void AbortChildGenome();
	void DevelopChildrenGenomes(Genome* parentG);
	void DetermineRegType(i_bead it);
	void PotentialTypeChange(i_bead ii);
	int CountTypeAbundance(int type);

	//Mutation functions.
	i_bead MutateGene(i_bead it, int* pdel_length, double muf);
	i_bead MutateBsite(i_bead it, int* pdel_length, double muf);
	i_bead MutateHouse(i_bead it, int* pdel_length, double muf);

	int ChangeParameter(int value);

	i_bead DeleteGene(i_bead it, int* pdel_length);
	i_bead DeleteBead(i_bead it);

	i_bead DuplicateGene(i_bead it, int* pdup_length);
	i_bead DuplicateBsite(i_bead it);
	i_bead DuplicateHouse(i_bead it);

	i_bead InventGene();
	void InventBsite();
	void InventHouse();

	i_bead ShuffleGene(i_bead it);
	i_bead ShuffleBead(i_bead it);

	i_bead FindFirstBsiteInFrontOfGene(i_bead it, bool ignore_houses=false) const;
	i_bead FindRandomGenePosition(bool include_houses, bool include_end) const;
	i_bead FindRandomPosition(bool include_end) const;
	void CopyPartOfGenomeToTemplate(i_bead begin, i_bead end, list<Bead*>* template_beadlist);
	void CloneGenome(const Genome* ImageG);
	void CopyPartOfGenome(i_bead begin, i_bead end);

	void ReadGenome(string genome);
	void ReadExpression(string expression);
	void ReadDefinition(string definition);
	void ReadBuffer(string buffer, bool* array, char start_sign, char stop_sign, int ith_start_sign = 1, int ith_stop_sign = 1);

	int WhatBead(Bead* bead) const;
	int BindingAffinity(const bool* sequenceA, const bool* sequenceB, int seqlen = sequence_length) const;

	string Show(list<Bead*>* chromosome, bool terminal, bool only_parent);
	string ShowExpression(list<Bead*>* chromosome, bool only_parent);
	string ShowDefinition(bool terminal);

};


inline int Genome::WhatBead(Bead* bead) const
//Ordered by expected number of beads (more houses than regulators), to increase speed.
{
	if ( typeid(*bead) == typeid(House) )						return HOUSE;
	else if( typeid(*bead) == typeid(Bsite) )				return BSITE;
	else if( typeid(*bead) == typeid(Gene) )		return REGULATOR;
	else
	{
		cerr << "Bead type error.\n" << endl;
		exit(1);
	}
}

inline int Genome::BindingAffinity(const bool* sequenceA, const bool* sequenceB, int seqlen) const
{
	int affinity = 0;
	for (int i=0; i<seqlen; i++)
	{
		affinity += (int)(sequenceA[i]^sequenceB[i]);
	}
	return affinity;
}

#endif
