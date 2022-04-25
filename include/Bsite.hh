#ifndef BsiteHeader
#define BsiteHeader

#include "Bead.hh"
#include "Header.hh"

// Binding sites for transcription factors, i.e. regulatory elements.

class Bsite : public Bead {
	public:
	 	int activity;
		std::bitset<regulator_length> sequence;

		Bsite();
		Bsite(int act, std::bitset<regulator_length>& seq);
		explicit Bsite(const Bsite &bsite);
		virtual ~Bsite();

		virtual Bead* Clone() const;

		void Randomize(int organelle);
		bool Mutate(int organelle);
		string Show(bool terminal, bool type_only=false) const;
};

#endif
