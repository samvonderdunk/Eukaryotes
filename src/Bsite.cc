#include "Bsite.hh"

Bsite::Bsite() : Bead()
{
	int i;

	duplicate = false;
  activity = 0;
  for(i=0; i<sequence_length; i++) sequence[i] = 0;
}

Bsite::Bsite(int act, bool seq[]) : Bead()
{
	int i;

	duplicate = false;
  activity = act;
  for(i=0; i<sequence_length; i++) sequence[i] = seq[i];
}

Bsite::Bsite(const Bsite &bsite) : Bead(bsite)
{
	int i;

	duplicate = bsite.duplicate;
  activity = bsite.activity;
  for(i=0; i<sequence_length; i++) sequence[i] = bsite.sequence[i];
}

Bsite::~Bsite()
{
}

Bead* Bsite::Clone() const
{
  return new Bsite(*this);
}

void Bsite::RandomBsite()
{
	int i;

	duplicate = false;
  activity = (uniform()>0.5) ? -1 : 1;
  for (i=0; i<sequence_length; i++)	sequence[i] = (uniform()>0.5) ? true : false;
}

string Bsite::Show(bool terminal) const
{
	int i;
	string Content, color_prefix, color_suffix;
	std::stringstream ss;

	if(terminal){
		color_prefix = "\033[92m";
		color_suffix = "\033[0m";
	}
	else
	{
		color_prefix = "";
		color_suffix = "";
	}

	ss << "(" << color_prefix << activity << ":";
	for(i=0; i<sequence_length; i++)	ss << sequence[i];
	ss << color_suffix << ")";
	Content = ss.str();
	ss.clear();

	return Content;
}
