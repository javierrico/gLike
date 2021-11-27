//////////////////////////////////////////////////////////////////////
// Histogram containing dNdE information
//////////////////////////////////////////////////////////////////////

#ifndef HDNDE
#define HDNDE

class HdNdE : public TH1F
{
 public:
  
  // constructors
  HdNdE(TString inputString="");
  
  // destructor
  virtual ~HdNdE();

 protected:
  
 private:

  ClassDef(HdNdE,1) // Full Likelihood (unbinned)
};

#endif
