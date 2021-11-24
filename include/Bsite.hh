#ifndef BsiteHeader
#define BsiteHeader

#include "Bead.hh"
#include "Header.hh"

// Binding sites for transcription factors, i.e. regulatory elements.

class Bsite : public Bead {
	public:
	 	int activity;
		bool sequence[regulator_length];

		Bsite();
		Bsite(int act, bool seq[]);
		explicit Bsite(const Bsite &bsite);
		virtual ~Bsite();

		virtual Bead* Clone() const;
		void Randomize();

		bool Mutate(double mut_factor);
		string Show(bool terminal, bool type_only=false) const;
};

#endif
