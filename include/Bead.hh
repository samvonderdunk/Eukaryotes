#ifndef BeadHeader
#define BeadHeader

class Bead {
	public:
		bool duplicate;

		Bead();
		explicit Bead(const Bead &b);
		virtual ~Bead();

		virtual Bead* Clone() const=0;
};

#endif
