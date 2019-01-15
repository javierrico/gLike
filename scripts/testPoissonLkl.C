void testPoissonLkl()
{
  // create and configure the simplest possible PoissonLkl object
  const Int_t Non = 130;
  const Int_t Noff = 90;
  const Double_t tau = 1.;
  PoissonLkl* p = new PoissonLkl(Non,Noff,tau);
  
  // configure the minimization algorithm
  const Double_t errorDef = 4;
  const Double_t unitsOfG = 2;
  p->SetErrorDef(errorDef);
  p->SetUnitsOfG(unitsOfG);

  // call the minimization
  p->ComputeLklVsG();

  // print fit results
  cout << endl << "PrintOverview:" << endl;
  p->PrintOverview();     // print the details from the fit
  p->GetLklVsG()->Draw(); // plot the -2logL vs g curve
  
  // access to the gLike results
  const Double_t gmin = p->GetGLklMin();                // get the value of g that minimizes -2logL
  const Double_t gerr = p->GetGLklMinErr();             // get the value of gerr such that -2logL(gmin+/-gerr)=-2logLmin+fErrorDef
  const Double_t sig  = sqrt(p->GetLklVsG()->Eval(0)) ; // compute significance of detection

  // compute significance with LiMa formula for comparison
  const Double_t LiMasig = sqrt(2.*(Non*log((tau+1.)*(1.*Non/(Non+Noff)))+Noff*log((tau+1)/tau*(1.*Noff/(Non+Noff)))));

  // compare confidence interval with Rolke method for comparison
  TRolke* r = new TRolke(0.9544);
  r->SetPoissonBkgKnownEff(Non,Noff,tau,1);
  Double_t low,upp;
  r->GetLimits(low,upp);
  Double_t grolke    = unitsOfG*(upp+low)/2;
  Double_t grolkeerr = unitsOfG*(upp-low)/2;

  // print and compare results
  cout << endl;
  cout << "The max lkl  value of g = " << gmin << ",\t with "<< sqrt(errorDef) << "-sigma error bar = " << gerr << endl;
  cout << "The Rolke    value of g = " << grolke << ",\t with 2-sigma error bar = " << grolkeerr << endl;
  cout << "The gLike significance of the detection of signal is " << sig     << " sigma" << endl;
  cout << "The Li&Ma significande of the detection of signal is " << LiMasig << " sigma" << endl;


  Double_t Deff  = 0.1;
  Double_t Dtau  = 0.1;
  p->SetDEff(Deff);
  p->SetDTau(Dtau);
  p->ComputeLklVsG();
  p->PrintOverview();     // print the details from the fit
  p->GetLklVsG()->Draw("same"); // plot the -2logL vs g curve
}
