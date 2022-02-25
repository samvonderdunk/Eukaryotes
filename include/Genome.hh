#ifndef GenomeHeader
#define GenomeHeader

#include "Header.hh"
#include "Bead.hh"
#include "House.hh"
#include "Bsite.hh"
#include "Gene.hh"
#include "Regulator.hh"
#include "Effector.hh"

class Genome {
 public:
  std::list<Bead*>* 					BeadList;
	std::array<Regulator*,5>		RegTypeList;	//Used to define 5 main regulatory types, not actually used when we have effectors.
	std::list<Bead*>*						ExpressedGenes;

  int g_length, gnr[4];
  int fork_position, terminus_position;	//position of the replication fork and the terminus, where it stops.
  bool is_mutated;
	int organelle;

	typedef std::list<Bead*>::iterator i_bead;
	typedef std::list<Bead*>::reverse_iterator ri_bead;

  Genome();
  ~Genome();

	void UpdateGeneExpression();
	void NativeExpression();
	i_bead RegulatorCompetition(i_bead i_bsite);

	void ReplicateStep(double resource);

	void SplitGenome(Genome* parentG);
	void AbortChildGenome();
	void DevelopChildrenGenomes(Genome* parentG);
	void PotentialTypeChange(i_bead ii);
	int CountTypeAbundance(int type);

	//Mutation functions.
	i_bead Mutation(i_bead it, int* pdel_length);
	i_bead Deletion(i_bead it, int* pdel_length);
	i_bead Duplication(i_bead it, int* pdup_length);
	void Inventions(int* pdup_length);
	i_bead Shuffle(i_bead it);
	void WholeGenomeDuplication(int* pdup_length);

	int CountBeads(int kind);

	i_bead FindFirstBsiteInFrontOfGene(i_bead it, bool ignore_houses=false) const;
	i_bead FindRandomGenePosition(bool include_houses, bool include_end) const;
	i_bead FindRandomPosition(bool include_end) const;
	void CopyPartOfGenomeToTemplate(i_bead begin, i_bead end, list<Bead*>* template_beadlist);
	void CloneGenome(const Genome* ImageG);
	void CopyPartOfGenome(i_bead begin, i_bead end);

	void ReadGenome(string genome);
	void ReadExpression(string expression);
	void ReadDefinition(string definition);

	string Show(list<Bead*>* chromosome, bool terminal, bool only_parent);
	string ShowExpression(list<Bead*>* chromosome, bool only_parent);
	string ShowDefinition(bool terminal);

};

#endif
