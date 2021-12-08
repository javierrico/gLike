/* ======================================================================== *\
!
!   Author(s): Javier Rico         01/2015 <mailto:jrico@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//						
// FermiTables2016Lkl
//
// This class takes as input a list of values of -2logL vs Flux in
// different bins of energy (in particular, those produced by Fermi)
// and computes likelihood for a given spectrum (fHdNdESignal).
// The spectrum has to be specified fully (i.e. with the correct units)
// and therefore SetUnitsOfG is deactivated for this class.
// 
// Usage example (for a complete example check exec/jointLklDM.cc):
// ----------------------------------------------------------------
//
// FermiTables2016Lkl** bfLkl = new FermiTables2016Lkl*[nsample];
// JointLkl*       jLkl = new JointLkl;
// 
// // Configure each FermiTables2016Lkl sample
// for(Int_t isample=0;isample<nsample;isample++)
//   {
//      // Configure the samples
//      bfLkl[isample] = new FermiTables2016Lkl(inputString[isample]);
//      bfLkl[isample]->ReaddNdESignal(dNdESignalFileName);
//fermiLkl->SetDMMass(mass);
//
//      // Set units for DM interpretation of results
//      bfLkl[isample]->SetDMMass(mass)
//      bfLkl[isample]->SetDUnitsOfG(Dlog10_J[isample],Lkl::invlog);    // if we know the error of log10(J)
//
//      // add the sample to the JointLkl object
//      jLkl->AddSample(bfLkl[isample]);
//   }
//
// // Configure the Joint Likelihood
// jLkl->SetErrorDef(deltaLogLkl); // set the error correponding to the required CL
// jLkl->SetGIsPositive();                       // only G values are considered (Fermi prescription)
//
// // Compute and profile -2logL curve
// jLkl->ComputeLklVsG();                        // compute joint -2logL vs g around minimum
//
// // get the results
// Double_t gminval     = jLkl->GetGLklMin();    // value of g minimizing joint -2logL
// Double_t glimval     = jLkl->GetGForLkl(2.7); // value of g for which -2logL = minimum+2.7
// Double_t gminvalErr  = glimval-gminval;       // determination of 95% CL error 
//
// jLkl->GetLklVsG()->Draw("al");                // draw joint -2logL vs g (in specified units) 
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TPRegexp.h"

#include "FermiTables2016Lkl.h"
#include "ParabolaLkl.h"

ClassImp(FermiTables2016Lkl);

using namespace std;

// static constants
static const TString  gName            = "FermiTables2016Lkl";
static const TString  gTitle           = "Flux vs E-Bin Likelihood";
static const Int_t    gNPars           = 1;                   // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g"};               // Name of parameters
static const Double_t gRefsv           = 3.e-26;              // [cm^3 s^-1] (reference sv for computing dPhi/dE)
static const Int_t    gNFineBins       = 5000;                // default number of fine bins for dNdESignal histo
static const Double_t gFineLEMin       = TMath::Log10(0.01);  // default minimum log(energy[GeV]) for dNdESignal hist
static const Double_t gFineLEMax       = TMath::Log10(1000);  // default maximum log(energy[GeV]) for dNdESignal histo
static const Int_t    gNMaxIndex       = 300;                 // number of empty-field-of-view files (for null-hypothesis pdf evaluation)


// -2logL function for minuit
void binnedFluxLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;


////////////////////////////////////////////////////////////////
// 
// Default constructor (see InterpretInputString for details)
//
FermiTables2016Lkl::FermiTables2016Lkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), JointLkl(inputString), HdNdE("","",gNFineBins,gFineLEMin,gFineLEMax),
  fNBins(0), fBinEMin(NULL), fBinEMax(NULL), 
  fEFluxInt(NULL), fMass(0), fLogJ(0)
{
  if(InterpretInputString(inputString))
    cout << "FermiTables2016Lkl::FermiTables2016Lkl Warning: there were problems interpreting the input string" << endl;
}
//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been passed to Lkl::InterpretInputString in constructor
// here it is searched for the following options:
//
// logJ=<val>:          mean value of the log10 of the J factor, in units of GeV^2 cm^-5 (annihilation) or GeV cm^-3 (decay)
// path=<val>:          path of the input file (will be appended to inputFileName)
// inputfile=<val>:     name of the input file, which contains the lkl vs flux tables in E bins
//                      in the format of the Fermi results published in
//                      Phys. Rev. D 89, 042001 (2014)
//                      available at http://www-glast.stanford.edu/pub_data/713/
// index=<val>:         number (-1, modulo gNMaxIndex) used to construct the input file for the case of empty fields of view (for null-hypothesis pdf evaluation)
//
Int_t FermiTables2016Lkl::InterpretInputString(TString inputString)
{
  // file name and path default values
  TString inputfileName = " (No file has been specified) ";
  TString path          = "";
  Int_t   index         = 0; // index for empty-field of view data file (for null-hypothesis pdf evaluation)

  // Prepeare to interpret inputString
  TPMERegexp re("\\s+");
  TPMERegexp fldre("=");
  
  // split the inputString into the different fields, and check the option and values specified in each of them
  UInt_t nfields = re.Split(inputString);
  for(UInt_t ifield=0;ifield<nfields;ifield++)
    {
      fldre.Split(re[ifield]);
      TString optname = fldre[0];
      if(optname.CompareTo("logJ",TString::kIgnoreCase)==0)
	fLogJ=fldre[1].Atof();
      if(optname.CompareTo("inputfile",TString::kIgnoreCase)==0)
	inputfileName=fldre[1];
      if(optname.CompareTo("index",TString::kIgnoreCase)==0)
	index=fldre[1].Atoi();
      else if(optname.CompareTo("path",TString::kIgnoreCase)==0)
	path=fldre[1];    
    }

  if(index>0)
    inputfileName.ReplaceAll(".txt",Form("_%05d.txt",(index-1)%gNMaxIndex));
  
  // open and read input file
  return ReadFermiInputData(path+"/"+inputfileName);
}
////////////////////////////////////////////////////////////////
// 
// Destructor
//
FermiTables2016Lkl::~FermiTables2016Lkl()
{
  if(fBinEMin)       delete [] fBinEMin;
  if(fBinEMax)       delete [] fBinEMax;
  if(fEFluxInt)      delete [] fEFluxInt;
  if(fHdNdESignal)   delete fHdNdESignal;
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of minuit subtleties)
//
void FermiTables2016Lkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN 
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
// ginit = initial value of g in the fit
//
void  FermiTables2016Lkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(binnedFluxLkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // Initialize minuit for added objects
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      sample->InitMinuit();
      sample->SetMinuitLink();
    }
  delete iter;

  Double_t pStart[gNPars] = {(ginit==0? 1 :ginit)};
  Double_t pDelta[gNPars] = {1};

  SetParameters(gParName,pStart,pDelta);
}  

////////////////////////////////////////////////////////////////
//
// Check that everything is correct before calling the minimization
//
// Return 0 in case of success, 1 otherwise
//
Int_t FermiTables2016Lkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // Look for dN/dE
  if(!fHdNdESignal)
    {
      cout << "FermiTables2016Lkl::MakeChecks (" << GetName() << ") Warning: missing fHdNdESignal histo" << endl;
      return 1;
    }

  // compute the E flux integral bin per bin (to be compared to Fermi results)
  if(ComputeEFluxIntegrals())
    {
      cout << "FermiTables2016Lkl::MakeChecks (" << GetName() << ") Warning: problems computing E-flux integrals" << endl;
      return 1;
    }

  // Set proper units
  Lkl::SetUnitsOfG(gRefsv);

  // Set weights of different samples
  TObjArrayIter* iter     = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl*          sample;
  Int_t cnt=0;
  while((sample=(Lkl*) iter->Next()))
    sample->SetUnitsOfG(fEFluxInt[cnt++]);
  
  // Check the JointLkl part
  ////////////////////////////
  if(ReorderSamples())
    return 1;

  // Check all parabolas
  iter->Reset();
  while((sample=(Lkl*) iter->Next()))
    if(!sample->IsChecked())
      sample->MakeChecks();
  delete iter;

  SetChecked();
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// compute the values of nuisance parameter of added samples
// This is an internal function called in Lkl::FixLklVsG
//
void FermiTables2016Lkl::SpreadFixLklVsG(Double_t g)
{
  // compute the values of nuisance parameter of added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();

  // loop over all samples and compute the -2logL curve in the corresponding range
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      Double_t   w = sample->GetUnitsOfG(); 
      if(w>0)
	sample->MinimizeLkl(g*w,kTRUE,kFALSE,kTRUE);
    }
  delete iter;
}


////////////////////////////////////////////////////////////////
//
// Disable SetUnitsOfG, since units are used internally by the 
// class and should not be changed.
// SetDMMass must be called instead!
//
void FermiTables2016Lkl::SetUnitsOfG(Double_t units)
{
  cout<< "FermiTables2016Lkl::SetUnitsOfG (" << GetName() << ") Warning: cannot set units for FermiTables2016Lkl. Use SetDMMass instead" << endl;
}

////////////////////////////////////////////////////////////////
//
// Read the Fermi input file and construct+add the 
// ParabolaLkl describing the -2logL function vs g 
//
// Return 0 in case of success, 1 otherwise
//
Int_t FermiTables2016Lkl::ReadFermiInputData(TString filename)
{
  const Double_t MeV2GeV = 1.e-3;

  // try to open file
  ifstream ff(filename.Data());
  if(!ff.is_open())
    {
      cout << "FermiTables2016Lkl::ReadFermiInputData (" << GetName() << ") Error: could not open file " << filename << endl;
      return 1;
    }
  else
    cout << "FermiTables2016Lkl::ReadFermiInputData (" << GetName() << ") Message: Reading likelihood values from file " << filename << endl;

  // skip first two lines
  ff.ignore(9999,'\n');
  ff.ignore(9999,'\n');

  // read file and create+add the ParabolaLogL objects
  Double_t emin,emax,eflux,logL;
  vector<Double_t> veflux;
  vector<Double_t> vlogL;
  veflux.reserve(100);
  vlogL.reserve(100);
  Double_t prevemin = -1;
  Double_t prevemax = -1;
  while(1)
    {
      ff >> emin >> emax >> eflux >> logL;
      if(ff.eof()) break;

      // a new bin starts
      if(emin!=prevemin)
	{
	  if(veflux.size()) // a bin is finished and ready to be saved
	    {
	      CreateAndAddNewParabola(prevemin*MeV2GeV,prevemax*MeV2GeV,veflux.size(),veflux.data(),vlogL.data());
	      veflux.clear();
	      vlogL.clear();
	    }

	  prevemin=emin;
	  prevemax=emax;
	}

      veflux.push_back(eflux*MeV2GeV);
      vlogL.push_back(-2*logL);
    }
  // include the last bin
  CreateAndAddNewParabola(prevemin*MeV2GeV,prevemax*MeV2GeV,veflux.size(),veflux.data(),vlogL.data());
  
  ff.close(); 

  SetChecked(kFALSE);

  cout << "FermiTables2016Lkl::ReadFermiInputData (" << GetName() << ") Message: Read " <<  GetSampleArray()->GetEntries() << " parabolas" << endl;
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Add a new ParabolaLkl object representing the -2logL values vs
// energy flux int npoints, log10-even spaced, for bin of energy
// between emin and emax
//
// Return 0 in case of success, 1 otherwise
//
Int_t FermiTables2016Lkl::CreateAndAddNewParabola(Double_t emin,Double_t emax,Int_t npoints,Double_t* eflux,Double_t* logL)
{
  fNBins++;

  // Add new bin bounds to the list
  Double_t* newEminList = new Double_t[fNBins];
  Double_t* newEmaxList = new Double_t[fNBins];

  // copy old lists into new ones
  if(fNBins>1)
    {
      for(Int_t ibin=0;ibin<fNBins-1;ibin++)
	{
	  newEminList[ibin] = fBinEMin[ibin];
	  newEmaxList[ibin] = fBinEMax[ibin];
	}
      delete [] fBinEMin;
      delete [] fBinEMax;
    }
  newEminList[fNBins-1] = emin;
  newEmaxList[fNBins-1] = emax;
  fBinEMin = newEminList;
  fBinEMax = newEmaxList;

  // Create and add the new ParabolaLkl object
  ParabolaLkl* newparabola = new ParabolaLkl(npoints,eflux,logL,Form("ParabolaLkl_%d",GetSampleArray()->GetEntries()));
  AddSample(newparabola);
  
  SetChecked(kFALSE);
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal and delete
// all histograms that depend on it
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t FermiTables2016Lkl::ResetdNdESignal()
{
  if(HdNdE::ResetdNdESignal())
    return 1;

  ResetGLklVsG();
  SetChecked(kFALSE);

  return 0;
}

//////////////////////////////////////////////////////////////
// 
// Read dN/dE signal histogram from Fermi-internal file (mostly for tests)
// Replacement of existing histogram is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t FermiTables2016Lkl::ReaddNdESignalFromFermi(TString filename)
{
  // open file and look for histo
 // try to open file
  ifstream ff(filename.Data());
  if(!ff.is_open())
    {
      cout << "FermiTables2016Lkl::ReaddNdESignalFromFermi (" << GetName() << ") Error: could not open file " << filename << endl;
      return 1;
    }
  else
    cout << "FermiTables2016Lkl::ReaddNdESignalFromFermi (" << GetName() << ") Message: Reading dN/dE values from file " << filename << endl;

  // skip first six lines
  for(Int_t i=0;i<6;i++)
    ff.ignore(9999,'\n');

  // read energy and dN/dE values - transform them into the right units
  Double_t logEMeV,dNdEMeV;
  vector<Double_t> vlogE;
  vector<Double_t> vdNdE;
  vlogE.reserve(100);
  vdNdE.reserve(100);
  while(1)
    {
      ff >> logEMeV >> dNdEMeV;
      if(ff.eof()) break;

      vlogE.push_back(logEMeV-3);
      vdNdE.push_back(dNdEMeV*1e3);
    }
  ff.close(); 

  // built histogram out of the read values
  Int_t nbins = vlogE.size();
  Float_t* binbound = new Float_t[nbins+1];
  for(Int_t ibin=1;ibin<nbins;ibin++)
    binbound[ibin] = (vlogE.data()[ibin-1]+vlogE.data()[ibin])/2.;
  binbound[0]     = vlogE.data()[0]-(binbound[1]-vlogE.data()[0]);
  binbound[nbins] = vlogE.data()[nbins-1]+(vlogE.data()[nbins-1]-binbound[nbins-1]);

  TH1F* hdNdE = new TH1F("fHdNdESignal","dN/dE for signal events",nbins,binbound);
  for(Int_t ibin=0;ibin<nbins;ibin++)
    fHdNdESignal->SetBinContent(ibin+1,vdNdE.data()[ibin]);

  // set this as the dNdE histogram
  SetdNdESignal(hdNdE);

  SetChecked(kFALSE);

  // clean and exit
  delete hdNdE;
  delete [] binbound;
  return 0;

}

//////////////////////////////////////////////////////////////
// 
// Integrate E*dPhi/dE over E for all bins
// Return 0 in case of success
//        1 if file is not found
//
Int_t FermiTables2016Lkl::ComputeEFluxIntegrals()
{  
  // basic checks
  if(fNBins<=0) return 0;
  if(!fHdNdESignal || !fBinEMin || !fBinEMax) return 1;

  // create the histogram E*dPhi/dE to be integrated
  Int_t nbins            = fHdNdESignal->GetNbinsX();
  TH1F* hEdPhidE         = new TH1F(*fHdNdESignal);
  hEdPhidE->SetDirectory(0);
  hEdPhidE->SetName("hEdPhidE");
  hEdPhidE->SetTitle("");
  for(Int_t ibin=0;ibin<nbins;ibin++)
    {
      Float_t yval = fHdNdESignal->GetBinContent(ibin+1);
      if(yval<0) yval=0;
      else
	{
	  Double_t E = TMath::Power(10,fHdNdESignal->GetBinCenter(ibin+1));
	  Double_t J = TMath::Power(10,fLogJ);
	  yval       = E*J*gRefsv/(8*TMath::Pi()*fMass*fMass)*yval;
	}
      hEdPhidE->SetBinContent(ibin+1,yval);
    }

  // create new array and fill it with the vales of the integrals
  if(fEFluxInt) delete [] fEFluxInt;
  fEFluxInt = new Double_t[fNBins];
  for(Int_t ibin=0;ibin<fNBins;ibin++)
    fEFluxInt[ibin] = IntegrateLogE(hEdPhidE,TMath::Log10(fBinEMin[ibin]),TMath::Log10(fBinEMax[ibin]));
  
  // clean and exit
  delete hEdPhidE;
  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Print the values of the parameters of the Lkl-based object
//
void FermiTables2016Lkl::PrintData(Int_t level)
{
  Lkl::PrintData(level);

  Margin(level); cout << "                # of bins = " << fNBins <<endl;
  Margin(level); cout << "                   E'_min = " << fBinEMin[0]        << " GeV" << endl;
  Margin(level); cout << "                   E'_max = " << fBinEMax[fNBins-1] << " GeV" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Plot a canvas with all material to compute the -2logL:
// - Array of graphs (1 per energy bin) with -2logL vs energy flux
// 
// return a pointer to the canvas, in case you want to modify it
//
TCanvas* FermiTables2016Lkl::PlotInputData() 
{
  const Int_t nlines = 5;  // number of lines in parabolas canvas
  Int_t ncols = TMath::Ceil(fNBins/Float_t(nlines)); // number of columns in parabolas canvas

  // create and divide canvas
  TCanvas* canvas = new TCanvas("inputDataCanvas","-2logL value per energy bin and energy flux",ncols*250,nlines*250);
  canvas->Divide(ncols,nlines);

  for(Int_t ibin=0;ibin<fNBins;ibin++)
    {
      ParabolaLkl* prbla   = (ParabolaLkl*)GetSample(ibin);
      TGraph*       grprbla = prbla->GetParabola(kFALSE);

      canvas->cd(ibin+1);
      TH1I *dummypr = new TH1I(Form("dummypr_%d",ibin),Form("-2logLkl vs E-flux for E=[%.1f,%.1f] GeV",fBinEMin[ibin],fBinEMax[ibin]),1,grprbla->GetX()[1]*0.1,grprbla->GetX()[grprbla->GetN()-1]);
      dummypr->SetDirectory(0);
      dummypr->SetStats(0);
      dummypr->SetXTitle("E-Flux GeV/cm^{2}/s");
      dummypr->SetYTitle("#Delta(-2logL)");
      dummypr->SetMinimum(0);
      dummypr->SetMaximum(10);
      dummypr->DrawCopy();
      gPad->SetLogx();
      delete dummypr;

      // plot -2logLkl vs E-flux
      grprbla->Draw("l");
      gPad->SetGrid();
      gPad->Modified();
      gPad->Update();
    }
  return canvas;
}


////////////////////////////////////////////////////////////////////////
// joint likelihood function (-2logL) for all bins
// To be minimized by TMinuit
// Free parameters:
// par[0] = g    (global scale for dPhidESignal we try)
//
void binnedFluxLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  f = 0;

  // get internal object
  FermiTables2016Lkl* mylkl   = dynamic_cast<FermiTables2016Lkl*>(minuit->GetObjectFit());
  
  // g is the global scale for dPhidESignal
  Double_t        g       = par[0];

  TObjArrayIter*  iter    = (TObjArrayIter*) mylkl->GetSampleArray()->MakeIterator();
  ParabolaLkl*   sample;

  while((sample=dynamic_cast<ParabolaLkl*>(iter->Next())))
    {
      Double_t   w_i = sample->GetUnitsOfG();
      Double_t   g_i = g*w_i;

      if(w_i>0)
	f+=sample->MinimizeLkl(g_i,kTRUE,kFALSE);      
    }   

  delete iter;
}
