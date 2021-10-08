#include "House.hh"

House::House() : Bead()
{
}

House::~House() {
}

House::House(const House &house) : Bead(house)
{
	duplicate = false;
}

Bead* House::Clone() const
{
  return new House(*this);
}

string House::Show(bool terminal) const
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
