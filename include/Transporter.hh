#ifndef TransporterHeader
#define TransporterHeader

#include "Bead.hh"

// Transporters (genes): regulated like regulatory genes, upon expression move to the membranes and allow expressed regulatory genes to move between compartments (i.e. expression of symbionts may be controlled by host and vice versa).
// For future use.

class Transporter : public Bead {
 public:

  Transporter();
  ~Transporter();
  explicit Transporter(const Transporter &tp);
  virtual Bead* Clone() const;

};

#endif
