#include "Transporter.hh"
#include "Header.hh"

Transporter::Transporter() : Bead()
{
}

Transporter::~Transporter()
{
}

Transporter::Transporter(const Transporter &tp) : Bead(tp)
{
}

Bead* Transporter::Clone() const
{
  return new Transporter(*this);
}
