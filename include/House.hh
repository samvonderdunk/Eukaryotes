#ifndef HouseHeader
#define HouseHeader

#include "Bead.hh"
#include "Header.hh"

// Household genes, required for normal fitness.

class House : public Bead {
	public:

		House();
		~House();
		explicit House(const House &house);
		virtual Bead* Clone() const;

		string Show(bool terminal) const;

};

#endif
