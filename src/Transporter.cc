#include "Transporter.hh"
#include "Header.hh"

Transporter::Transporter() : Bead()
{
	type = 0;
	threshold = 0;
	for(int i=0; i<sequence_length; i++) sequence[i] = 0;
	expression = 0;
	express = 0;
}

Transporter::Transporter(const Transporter &tp) : Bead(tp)
{
	type=tp.type;
	threshold=tp.threshold;
	for(int i=0; i<sequence_length; i++) sequence[i] = tp.sequence[i];
	expression=tp.expression;
	express=tp.xpr;
}

Bead* Transporter::Clone() const
{
  return new Transporter(*this);
}

void Transporter::RandomTransporter()
{
  type = 0;	//This will be assigned somewhere else (in the genome)???
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Value between -WeightRange and +WeightRange (incl. borders).
  for (int i=0; i<sequence_length; i++)	sequence[k] = (uniform()>0.5) ? true : false;
  expression = 0;
	express = 0;
}

Transporter::~Transporter()
{
}
