#ifndef HouseHeader
#define HouseHeader

#include "Bead.hh"
#include "Header.hh"

// Household genes, required for normal fitness.

class House : public Bead {
	public:

		using Bead::Bead;	//All the same constructors as Bead class.
		House();	//Only this is slightly changed compared to the Bead con/destructors, bc we want to tell use the Bead constructor that sets kind to HOUSE.
		// House();
		// ~House();
		// explicit House(const House &house);
		virtual Bead* Clone() const;
		void Randomize();

		bool Mutate(int organelle);
		string Show(bool terminal, bool type_only=false) const;

};

#endif
