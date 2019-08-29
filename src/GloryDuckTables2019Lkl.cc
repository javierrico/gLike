/* ======================================================================== *\
!
!   Author(s): Daniel Kerszberg         01/2019 <mailto:dkerszberg@ifae.es>
!   Author(s): Chiara Giuri             01/2019 <mailto:chiara.giuri@desy.de>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//						
// GloryDuckTables2019Lkl
//
// This class takes as input a list of values of -2logL vs <sv> for 
// different masses (in the framework of the Glory Duck project).
// 
// Usage example:
// --------------
//
// GloryDuckTables2019Lkl *Lkl = new GloryDuckTables2019Lkl(inputString);
// Lkl.PlotInputData();
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TMath.h"
#include "TPRegexp.h"

#include "GloryDuckTables2019Lkl.h"
#include "ParabolaLkl.h"

ClassImp(GloryDuckTables2019Lkl);

using namespace std;

// static constants
static const TString  gName            = "GloryDuckTables2019Lkl";
static const TString  gTitle           = "<sv> vs mass";
static const Int_t    gNPars           = 1;             // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g"};         // Name of parameters
static const Double_t gRefsv           = 1.e-26;        // [cm^3 s^-1] (reference sv)

// -2logL function for minuit
void gloryDuckTables2019Lkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

////////////////////////////////////////////////////////////////
// 
// Default constructor
// inputfilename is the name of a file with the input data
// (i.e. logL vs <sv> for several DM masses)
//
GloryDuckTables2019Lkl::GloryDuckTables2019Lkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), fMass(NULL), fNMasses(0),
  fActiveMass(0), fNsvVals(0), fSampleArray(NULL)
{
  if(InterpretInputString(inputString))
    cout << "GloryDuckTables2019Lkl::GloryDuckTables2019Lkl Warning: there were problems interpreting the input string" << endl;
}
//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been passed to Lkl::InterpretInputString in constructor
// here it is searched for the following options:
//
// inputfile=<val>:     name of the input file, which contains the lkl vs <sv> tables
//
Int_t GloryDuckTables2019Lkl::InterpretInputString(TString inputString)
{
  // file name and path default values
  TString inputfileName = " (No file has been specified) ";
  TString path          = "";

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
      else if(optname.CompareTo("inputfile",TString::kIgnoreCase)==0)
	inputfileName=fldre[1];
      else if(optname.CompareTo("path",TString::kIgnoreCase)==0)
	path=fldre[1];    
    }

  fSampleArray = new TObjArray();

  // open and read input file
  return ReadGloryDuckInputData(path+"/"+inputfileName);
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
GloryDuckTables2019Lkl::~GloryDuckTables2019Lkl()
{
  if(fSampleArray) delete fSampleArray;
  if(fMass)        delete [] fMass;
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of minuit subtleties)
//
void GloryDuckTables2019Lkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN 
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
//
void GloryDuckTables2019Lkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(gloryDuckTables2019Lkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // Initialize minuit for added objects
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      if(GetGIsPositive())
        sample->SetGIsPositive();
      sample->InitMinuit();
      sample->SetMinuitLink();
    }
  delete iter;

  Double_t pStart[gNPars] = {(ginit==0? gRefsv :ginit)};
  Double_t pDelta[gNPars] = {gRefsv}; // to be fine tuned when applied to real data

  SetParameters(gParName,pStart,pDelta);
}  

////////////////////////////////////////////////////////////////
//
// Check that everything is correct
//
// Return 0 in case of success, 1 otherwise
//
Int_t GloryDuckTables2019Lkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // Check that the number of mass is not 0
  if(!fNMasses)
    {
      cout << "GloryDuckTables2019Lkl::MakeChecks (" << GetName() << ") Warning: fNMasses = 0" << endl;
      return 1;
    }

  // Check that the array of masses is not NULL 
  if(!fMass)
    {
      cout << "GloryDuckTables2019Lkl::MakeChecks (" << GetName() << ") Warning: fMass is empty" << endl;
      return 1;
    }

  // Check that the number of <sv> values scanned is not 0
  if(!fNsvVals)
    {
      cout << "GloryDuckTables2019Lkl::MakeChecks (" << GetName() << ") Warning: fNsvVals = 0" << endl;
      return 1;
    }

  // Check all parabolas
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  iter->Reset();
  while((sample=(Lkl*) iter->Next()))
    if(!sample->IsChecked())
      sample->MakeChecks();
  delete iter;

  SetChecked();
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Read the Glory Duck input file and add the 
// ParabolaLkl describing the -2logL function vs <sv> 
//
// Return 0 in case of success, 1 otherwise
//
Int_t GloryDuckTables2019Lkl::ReadGloryDuckInputData(TString filename)
{
  // try to open file
  ifstream ff(filename.Data());
  if(!ff.is_open())
    {
      cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Error: could not open file " << filename << endl;
      return 1;
    }
  else
    cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Message: Reading masses and <sigmav> values from file " << filename << endl;

  // get the number of <sv> values (which is one less than the number of lines in the file)
  string line;
  while (getline(ff, line))
    fNsvVals++;
  fNsvVals -= 1;
  ff.clear();
  ff.seekg(0, ios::beg);

  vector<Double_t> sigmav ;                   // vector with <sv> values
  vector<Double_t> mass ;                     // vector with mass values
  vector<vector<Double_t> > vlkl2D(fNsvVals); // 2D vector with -2logL values
  vector<Double_t> vlkl;                      // vector with -2logL values

  // get first line of the file (with mass values)
  string massLine;
  getline(ff,massLine);

  // read first line of the file
  std::istringstream first_line(massLine);
  Bool_t read_logJ = kTRUE;
  Double_t field;
  while (first_line >> field)
    {
      if(read_logJ)
        {
          if(TMath::Abs(field-fLogJ) > 1.e-6)
            {
              cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Error: the logJ value stored in the file doesn't match the one given as an input in the configuration." << endl;
              return 1;
            }
          read_logJ = kFALSE;
        }
      else
        mass.push_back(field);
    }

  Double_t readingSigmav = 0.;
  Double_t readingLkl = 0.;
  Int_t row = 0;
  Int_t col = 0;

  // Initialisation of variables needed to symmetrise the readed parabola
  // The symmetrisation helps migrad to find the minimum
  Bool_t init_lklmin = kFALSE;
  vector<Double_t> lklmin = vector<Double_t>(mass.size()); // vector with minimum -2logL value for each mass
  vector<Double_t> imin = vector<Double_t>(mass.size()); // vector with position of value for which -2logL is minimal

  // get <sv> (first column)
  while(ff >> readingSigmav)
    {
      // store the <sv> values
      sigmav.push_back(readingSigmav);

      // get rest of the line
      string LklLine;
      getline(ff,LklLine);

      // read rest of the line
      vlkl2D[row] = vector<Double_t>(mass.size());
      col = 0;
      std::istringstream line(LklLine);
      while (line >> readingLkl)
        {
          vlkl2D[row][col] = readingLkl;
          // storing the minimum value and position
          if(!init_lklmin)
            {
              lklmin[col] = readingLkl;
              imin[col] = col;
            }
          else
            {
              if(readingLkl<lklmin[col])
                {
                  lklmin[col] = readingLkl;
                  imin[col] = col;
                }
            }
          col++;
        }
      if(!init_lklmin) init_lklmin = kTRUE;
      if(col != mass.size())
        {
          cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Error: You have a different line size length for line " << row+2 << " ! Number of mass values = " << mass.size() << " while number of likelihood values = " << col << endl;
          return 1;
        }
      row++;
    }

  reverse(sigmav.begin(),sigmav.end());

  vector<Double_t> sigmav_extended;    // vector with <sv> values after symmetrisation (if needed)
  vector<Double_t> vlkl_extended;      // vector with -2logL values after symmetrisation (if needed)

  // add new parabolas (after symmetrisation if needed)
  for (col = 0; col < mass.size(); col++) {
    for (row = fNsvVals-1; row >=0; row--) {
      vlkl.push_back(vlkl2D[row][col]);
    }
    if(sigmav[0] > 0. && sigmav[fNsvVals-1] > 0.) // symmetrise the parabola if the user only provided positive <sv> values
      {
        cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Message: Provided values for the parabola are only for positive <sv> values so the parabola will be extended by a simple symmetrisation to help migrad to find the minimum later on." << endl;
        imin[col] = fNsvVals-1-imin[col];
        Int_t isym = fNsvVals-1;
        while(vlkl[isym] > vlkl[0] && isym>0) isym--;
        for (int i=fNsvVals-1;i>isym;i--)
          {
            vlkl_extended.push_back(vlkl2D[fNsvVals-1-i][col]);
            sigmav_extended.push_back(-sigmav[i]+sigmav[isym]);
          }
      }
    for (int i=0;i<=fNsvVals-1;i++)
      {
        vlkl_extended.push_back(vlkl[i]);
        sigmav_extended.push_back(sigmav[i]);
      }
    CreateAndAddNewParabola(mass[col],sigmav_extended.size(),sigmav_extended.data(),vlkl_extended.data());
    vlkl.clear();
    sigmav.clear();
    vlkl_extended.clear();
    sigmav_extended.clear();
  }

  ff.close(); 

  SetChecked(kFALSE);

  cout << "GloryDuckTables2019Lkl::ReadGloryDuckInputData (" << GetName() << ") Message: Read " <<  GetSampleArray()->GetEntries() << " parabolas for a total of " << fNMasses << " DM masses." << endl;
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Add a new ParabolaLkl object representing the -2logL values 
// vs <sv> for one DM mass
//
// Parameters:
// mass    = mass of the DM particle to which this parabola corresponds
// npoints = number of likelihood values to store in the parabola
// *sv     = vector of <sv> values
// *logL   = vector of likelihood values
//
// Return 0 in case of success, 1 otherwise
//
Int_t GloryDuckTables2019Lkl::CreateAndAddNewParabola(Double_t mass,Int_t npoints,Double_t* sv,Double_t* logL) 
{
  fNMasses++;

  // Add new bin bounds to the list
  Double_t* newMassList = new Double_t[fNMasses];

  // copy old lists into new ones
  if(fNMasses>1)
    {
      for(Int_t ibin=0;ibin<fNMasses-1;ibin++)
	{
	  newMassList[ibin] = fMass[ibin];
	}
      delete [] fMass;
    }

  newMassList[fNMasses-1] = mass;
  fMass = newMassList;

  // Create and add the new ParabolaLkl object
  ParabolaLkl* newparabola = new ParabolaLkl(npoints,sv,logL,Form("ParabolaLkl_%d",GetSampleArray()->GetEntries()));
  AddSample(newparabola);
  
  SetChecked(kFALSE);
  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Function called by Lkl::PrintOverview
// Print info about the experimental data
//
void GloryDuckTables2019Lkl::PrintData(Int_t level)
{
  Lkl::PrintData(level);

  TGraph* grprbla = GetSample(0)->GetLklVsG(kFALSE);
  Double_t* sv = grprbla->GetX();

  Margin(level); cout << "         Number of masses = " << fNMasses << endl;
  Margin(level); cout << "                 <sv>_min = " << sv[0] << " [cm^3 s^-1]" << endl;
  Margin(level); cout << "                 <sv>_max = " << sv[grprbla->GetN()-1] << " [cm^3 s^-1]" << endl;
  Margin(level); cout << "                            " << endl;

  if(level>1)
    {
      Margin(level); cout << " Parabola's content:        " << endl;
      Margin(level); cout << " Masses [GeV]     <sv> [cm^3/s] -->  " ;
      for(Int_t ibin=0;ibin<grprbla->GetN();ibin++)
        {
          cout << sv[ibin] << "  ";
        }
      cout << endl;

      for(Int_t ibin=0;ibin<fNMasses;ibin++)
        {
          Margin(level); cout << "    " << fMass[ibin] << "                              ";
          for(Int_t jbin=0;jbin<grprbla->GetN();jbin++)
            {
              cout << GetSample(ibin)->GetLklVsG(kFALSE)->Eval(sv[jbin]) << "  ";
            }
          cout << endl;
        }
      cout << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Plot a canvas with all the -2logL, ie an array
// of graphs (1 per mass) with -2logL vs <sv>
// 
// Return a pointer to the canvas, in case you want to modify it
//
TCanvas* GloryDuckTables2019Lkl::PlotInputData() 
{
  const Int_t nlines = 5;  // number of lines in parabolas canvas
  Int_t ncols = TMath::Ceil(fNMasses/Float_t(nlines)); // number of columns in parabolas canvas

  // create and divide canvas
  TCanvas* canvas = new TCanvas("inputDataCanvas","-2logL value for range of <sv> values scanned",ncols*250,nlines*250);
  canvas->Divide(ncols,nlines);

  for(Int_t ibin=0;ibin<fNMasses;ibin++)
    {
      ParabolaLkl* prbla = (ParabolaLkl*)GetSample(ibin);
      TGraph*    grprbla = prbla->GetParabola(kFALSE);
      Double_t*       sv = grprbla->GetX();

      canvas->cd(ibin+1);
      TH1I *dummypr = new TH1I(Form("dummypr_%d",ibin),Form("-2logLkl vs <sv> for m_{DM}=%.1f \\ and <sv>=[%.1f,%.1f] [cm^3 s^-1]",fMass[ibin],sv[0],sv[grprbla->GetN()-1]),1,grprbla->GetX()[1]*0.1,grprbla->GetX()[grprbla->GetN()-1]);
      dummypr->SetDirectory(0);
      dummypr->SetStats(0);
      dummypr->SetXTitle("<sv> [cm^3 s^-1]");
      dummypr->SetYTitle("#Delta(-2logL)");
      dummypr->SetMinimum(0);
      dummypr->SetMaximum(10);
      dummypr->DrawCopy();
      gPad->SetLogx();
      delete dummypr;

      // plot -2logLkl vs <sv> 
      grprbla->Draw("l");
      gPad->SetGrid();
      gPad->Modified();
      gPad->Update();
    }
  return canvas;
}

////////////////////////////////////////////////////////////////////////
//
// Set mass index to the current active mass using the mass value
//
// Return 0 in case of success, 1 otherwise
//
Int_t GloryDuckTables2019Lkl::SetActiveMass(Double_t mass)
{
  // check if the mass exists in fMass[]
  for (Int_t counter = 0; counter < fNMasses; counter++)
    {
      if(TMath::Abs(fMass[counter]-mass) < 1e-3)
        {
          fActiveMass = counter;
          cout << "GloryDuckTables2019Lkl::SetActiveMass() Message: Mass m = " << mass << " GeV exists, setting the active mass to " << fMass[counter] << " GeV." << endl;
          return 0;
        }
    }

  cout << "GloryDuckTables2019Lkl::SetActiveMass() Message: Mass m = " << mass << " GeV doesn't exist." << endl;
  return 1;
}

////////////////////////////////////////////////////////////////////////
//
// Set mass index to the current active mass using the index
//
// Return 0 in case of success, 1 otherwise
//
Int_t GloryDuckTables2019Lkl::SetActiveMass(Int_t index)
{
  // check if the index is positive and smaller than the number of masses
  if(index >=0 && index < fNMasses)
    {
      fActiveMass = index;
      cout << "GloryDuckTables2019Lkl::SetActiveMass() Message: Index " << index << " is in the defined range, setting the active mass to " << fMass[index] << "." << endl;
      return 0;
     }

  cout << "GloryDuckTables2019Lkl::SetActiveMass() Message: Index " << index << " is outside the defined range." << endl;
  return 1;
}

////////////////////////////////////////////////////////////////////////
//
// Likelihood function (-2logL) for one mass
// To be minimized by TMinuit
// Free parameters:
// par[0] = g
//
void gloryDuckTables2019Lkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  f = 0;

  // get internal object
  GloryDuckTables2019Lkl* mylkl = dynamic_cast<GloryDuckTables2019Lkl*>(minuit->GetObjectFit());

  // g is the global scale for dPhidESignal
  Double_t g = par[0];

  // get likelihood for current active mass
  ParabolaLkl* sample = dynamic_cast<ParabolaLkl*>(mylkl->GetSample(mylkl->GetActiveMass()));

  // minimize
  f=sample->MinimizeLkl(g,kTRUE,kFALSE);
}
