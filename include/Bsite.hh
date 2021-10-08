#ifndef BsiteHeader
#define BsiteHeader

#include "Bead.hh"
#include "Header.hh"

// Binding sites for transcription factors, i.e. regulatory elements.

class Bsite : public Bead {
	public:
	 	int activity;
		bool sequence[sequence_length];

		Bsite();
		Bsite(int act, bool seq[]);
		explicit Bsite(const Bsite &bsite);
		virtual ~Bsite();

		virtual Bead* Clone() const;
		void RandomBsite();

		string Show(bool terminal) const;
};

#endif
