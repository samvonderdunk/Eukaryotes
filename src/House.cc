#include "House.hh"

// House::House() : Bead()
// {
// }
//
// House::~House() {
// }
//
// House::House(const House &house) : Bead(house)
// {
// }

House::House() : Bead(HOUSE)
{
}

Bead* House::Clone() const
{
  return new House(*this);
}

void House::Randomize()
{
	cout << "House randomize" << endl;
}

bool House::Mutate(double mut_factor)	//Maybe don't need to define this? See what works for Randomize().
{
	return false;	//Nothing to mutate for now.
}

string House::Show(bool terminal, bool type_only) const
{
	string Content, color_prefix, color_suffix;
	std::stringstream ss;

	if(terminal){
		color_prefix = "\033[95m";
		color_suffix = "\033[0m";
	}
	else
	{
		color_prefix = "";
		color_suffix = "";
	}

	ss << "(" << color_prefix << "H" << color_suffix << ")";
	Content = ss.str();
	ss.clear();

	return Content;
}
