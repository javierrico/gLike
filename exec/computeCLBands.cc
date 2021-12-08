//######################################################################
//##
//## computeCLBands.C
//##
//## Compute confidence level (CL) bands to Dark Matter annihilation
//## cross-section (<sv>) or decay lifetime (tau) using files obtained
//## with the jointLklDM.C script.
//##
//## You can specify two arguments:
//## 1. Name of the rc file (a TString): it should be the same as the 
//##    one used for jointLklDM.C
//## 2. The total number of simulation files (a Int_t) you computed 
//##    using jointLklDM.C. If you ran X simulations, it will assume
//##    that the files have seeds in the range [1,X]
//##
//## The script will automaticaaly stop if more than 5% of the simualtion
//## files are not found.
//##
//## The script produces (and saves) two files:
//## - One .root file containing the median curve, the right anf left 
//##   limits of the 68% CL intervall and the right and left limits of 
//##   the 95% CL intervall.
//## - One canavs limcanvas with the limits on <sv> or tau vs mass for
//##   the considered DM channel, range of masses and data. Overlaid on 
//##   this limit will be the 68% anf 95% CL bands as well as the median
//##   of the simulations.
//##
//## You need to run the compiled version of the script
//## by typing ".x computeCLBands.C+" in the ROOT command interpreter.
//##
//######################################################################

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include "TPRegexp.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"

#include "Lkl.h"

using namespace std;

void usage();
void setDefaultStyle();
Int_t GetNSkippedMasses(Int_t nm,const Double_t* vm,Double_t minm);
void decode_channel(TObjArray* coefficients, Int_t &nChannels, TString *channelval, Double_t *brval);
void compute_branonBR(Float_t &mass, Int_t &nChannels, TString *channelval, Double_t *brval, Double_t &translation_factor);

// number of bins in histograms
const Int_t nbins = 1000;

// show plots?
const Bool_t makePlots = kTRUE;

int main(int argc, char* argv[])
{
  TString configFileName="$GLIKESYS/rcfiles/jointLklDM.rc";
  Int_t nSimuFiles=300;
  // check input parameters
  for (int i = 1; i < argc; ++i) {
    TString arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      exit(-1);
    }
    if (arg == "--config") {
      configFileName = argv[i+1];
    } 
    if (arg == "--nsimufiles") {
      nSimuFiles = atof(argv[i+1]);
    }
  }
  setDefaultStyle();

  // Look for configuration file
  gSystem->ExpandPathName(configFileName);
  TPMERegexp re("\\s+");
  
  Lkl::PrintGLikeBanner();
  
  // Print-out configuration info (part 1)
  cout << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***                          RUNNING computeCLBands.C                          ***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***" << endl;
  cout << "*** CONFIGURATION FILE  : " << configFileName << endl;
  if (gSystem->AccessPathName(configFileName, kFileExists))
    {
      cout << endl << "    Oops! problems reading config file file " << configFileName << " <---------------- FATAL ERROR!!!"<< endl;
      exit(1);
    }

  // Read configuration file
  TEnv*  env = new TEnv(configFileName);
  TString  label             = env->GetValue("jointLklDM.Label","");
  TString  channel           = env->GetValue("jointLklDM.Channel","bb");
  TString  process           = env->GetValue("jointLklDM.Process","ann");
  TString  fInputDataPath    = env->GetValue("jointLklDM.path","");
  TString  fPlotsDir         = fInputDataPath+"/"+env->GetValue("jointLklDM.plotsDir","")+"/";
  TString  massList          = env->GetValue("jointLklDM.MassList","");
 
  // fill the list of masses to be studied
  UInt_t  nmass0  = re.Split(massList);
  Double_t* massval0 = new Double_t[nmass0];
  for(UInt_t imass=0;imass<nmass0;imass++)
    massval0[imass] = re[imass].Atof();

  // Set some flags
  Bool_t  isDecay            = !process.CompareTo("dec",TString::kIgnoreCase);
  
  // define some labels according to input data
  TString strprocess         = (isDecay?       "Decay" : "Annihilation");

  // Decode the channel and save the decoded channels and branching ratios in channelval and brval, respectively.
  TObjArray* coefficients = channel.Tokenize("+");
  Int_t nChannels = coefficients->GetEntries();
  TString* channelval = new TString[nChannels];
  Double_t* brval = new Double_t[nChannels];
  // call the function that decodes the coefficients of the channel
  decode_channel(coefficients,nChannels,channelval,brval);

  // annihilation/decay channel string
  TString strchannel;
  TString normchannel;
  Double_t minmass = 0;
  std::ostringstream buffer;
  for (Int_t iChannel = 0; iChannel < nChannels; iChannel++)
    {
      if(nChannels == 1) strchannel = normchannel = "";
      else
          if(iChannel == 0) 
            {
              Float_t rounded_brval = ((Int_t)(brval[iChannel] * 100 + .5) / 100.0);
              buffer << rounded_brval;
              strchannel  = (TString) buffer.str().substr(0, 4) + "#upoint";
              normchannel = (TString) buffer.str().substr(0, 4) + "*";
            }
          else
            {
              Float_t rounded_brval = ((Int_t)(brval[iChannel] * 100 + .5) / 100.0);
              buffer.clear();
              buffer.str("");
              buffer << rounded_brval;
              strchannel.Append("#plus" + buffer.str().substr(0, 4) + "#upoint");
              normchannel.Append("+" + buffer.str().substr(0, 4) + "*");
            }
      normchannel.Append(channelval[iChannel]);

      if     (!channelval[iChannel].CompareTo("bb",TString::kIgnoreCase))         {strchannel.Append("b#bar{b}");              minmass=5;}
      else if(!channelval[iChannel].CompareTo("cc",TString::kIgnoreCase))         {strchannel.Append("c#bar{c}");              minmass=1.3;}
      else if(!channelval[iChannel].CompareTo("tt",TString::kIgnoreCase))         {strchannel.Append("t#bar{t}");              minmass=173;}
      else if(!channelval[iChannel].CompareTo("tautau",TString::kIgnoreCase))     {strchannel.Append("#tau^{+}#tau^{-}");      minmass=1.8;}
      else if(!channelval[iChannel].CompareTo("mumu",TString::kIgnoreCase))       {strchannel.Append("#mu^{+}#mu^{-}");        minmass=0.106;}
      else if(!channelval[iChannel].CompareTo("WW",TString::kIgnoreCase))         {strchannel.Append("W^{+}W^{-}");            minmass=80.3;}
      else if(!channelval[iChannel].CompareTo("ZZ",TString::kIgnoreCase))         {strchannel.Append("ZZ");                    minmass=91.2;}
      else if(!channelval[iChannel].CompareTo("hh",TString::kIgnoreCase))         {strchannel.Append("HH");                    minmass=125;}
      else if(!channelval[iChannel].CompareTo("gammagamma",TString::kIgnoreCase)) {strchannel.Append("#gamma#gamma");          minmass=0;}
      else if(!channelval[iChannel].CompareTo("pi0pi0",TString::kIgnoreCase))     {strchannel.Append("#pi^{0}#pi^{0}");        minmass=0.135;}
      else if(!channelval[iChannel].CompareTo("gammapi0",TString::kIgnoreCase))   {strchannel.Append("#pi^{0}#gamma");         minmass=0.135/2.;}
      else if(!channelval[iChannel].CompareTo("pi0gamma",TString::kIgnoreCase))   {strchannel.Append("#pi^{0}#gamma");         minmass=0.135/2.;}
      else if(!channelval[iChannel].CompareTo("ee",TString::kIgnoreCase))         {strchannel.Append("e^{+}e^{-}");            minmass=0.511e-3;}
      else if(!channelval[iChannel].CompareTo("nuenue",TString::kIgnoreCase))     {strchannel.Append("#nu_{e}#nu_{e}");        minmass=120.0;}
      else if(!channelval[iChannel].CompareTo("numunumu",TString::kIgnoreCase))   {strchannel.Append("#nu_{#mu}#nu_{#mu}");    minmass=120.0;}
      else if(!channelval[iChannel].CompareTo("nutaunutau",TString::kIgnoreCase)) {strchannel.Append("#nu_{#tau}#nu_{#tau}");  minmass=120.0;}
      else if(!channelval[iChannel].CompareTo("VV-4e",TString::kIgnoreCase))      {strchannel.Append("VV-4e");                 minmass=2.044e-3;}
      else if(!channelval[iChannel].CompareTo("VV-4mu",TString::kIgnoreCase))     {strchannel.Append("VV-4#mu");               minmass=0.424;}
      else if(!channelval[iChannel].CompareTo("VV-4tau",TString::kIgnoreCase))    {strchannel.Append("VV-4#tau");              minmass=7.2;}
      else if(!channelval[iChannel].CompareTo("branon",TString::kIgnoreCase))     {strchannel.Append("branon");                minmass=0;}
      else strchannel = "";
    }
  if(isDecay) minmass*=2;
  
  // Remove mass values below kinematical threshold
  Int_t           nskippedmass = GetNSkippedMasses(nmass0,massval0,minmass);
  Int_t           nmass        = nmass0-nskippedmass; 
  const Double_t* massval      = massval0+nskippedmass;

  // Print-out configuration info (part 2)
  cout << "*** I/O PATH            : " << fInputDataPath << endl;
  cout << "***" << endl;
  cout << "*** LABEL               : " << label <<  endl;
  cout << "*** PROCESS             : " << strprocess << endl;
  cout << "*** CHANNEL             : " << channel << endl;
      
  // how many and which mass values are we considering?
  cout << "*** NUMBER OF MASSES    : " << nmass << endl;
  cout << " ** Mass values         : "; 
  for(Int_t imass=0;imass<nmass;imass++)
    cout << massval[imass] << ((imass<nmass-1)? ", " : "");
  cout << " GeV" << endl;
  
  cout << "***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << endl;
  // End of print-out configuration info

  // where the root files for the simulations are
  const TString simusPath     = fPlotsDir + "MC/root/";

  // where the root files for the data are
  const TString dataPath      = fPlotsDir + "Data/root/";

  // where the output files will be
  const TString resultsPath   = fPlotsDir + "BrasilianPlot/";
  gSystem->Exec(Form("mkdir -p %s",resultsPath.Data()));
  gSystem->Exec(Form("mkdir -p %s/root",resultsPath.Data()));
  gSystem->Exec(Form("mkdir -p %s/pdf",resultsPath.Data()));
  gSystem->Exec(Form("mkdir -p %s/txt",resultsPath.Data()));

  // File to store the results in txt format
  std::ofstream results2txt;
  const TString results2txtPath = resultsPath + "txt/" + label + ".txt";
  results2txt.open(results2txtPath);

  // open file and get the limits for the txt file
  TFile*   dataFile = TFile::Open(dataPath+label+"_Data_limits.root", "READ");
  TCanvas* canvas   = (TCanvas*)dataFile->Get("limcanvas");
  TGraph*  graph    = (TGraph*)canvas->GetPrimitive(isDecay? "grtaulim": "grsvlim");
  // open file and get the limits
  Char_t line[256];
  for(Int_t ival=0;ival<graph->GetN();ival++)
    {
      Double_t mass = graph->GetX()[ival];
      for(Int_t imass=0;imass<nmass;imass++)
        if(mass==massval[imass])
          {
            results2txt << massval[imass] << " ";
            break;
          }
    }
  results2txt << endl;
  for(Int_t ival=0;ival<graph->GetN();ival++)
    {
      Double_t mass = graph->GetX()[ival];
      for(Int_t imass=0;imass<nmass;imass++)
        if(mass==massval[imass])
          {
            Double_t lim  = graph->GetY()[ival];
            results2txt << lim << " ";
            break;
          }
    }
  results2txt << endl;
  delete dataFile;

  // Arrays for results
  Double_t sv1sigmaL[nmass];
  Double_t sv1sigmaR[nmass];
  Double_t sv0sigma[nmass];
  Double_t sv2sigmaL[nmass];
  Double_t sv2sigmaR[nmass];

  // vectors to read input data
  vector<Double_t>* svLimDist = new vector<Double_t>[nmass];
  for(Int_t imass=0;imass<nmass;imass++)
    svLimDist[imass].reserve(nSimuFiles);

  // reading input files
  cout << "Reading " << nSimuFiles << " input files" << endl;
  Int_t nFilesRead=0;
  for(Int_t ifile=1;ifile<nSimuFiles+1;ifile++)
    {
      // file name
      TString seedtag  = Form("_%d",ifile);
      TString filename = simusPath+label+"_MC_limits"+seedtag+".root";
      cout << "Opening file " << filename << endl;
      Int_t exist = gSystem->Exec("ls " + filename + " > /dev/null");
      if ( exist != 0 ){
        cout << "cannot find " << filename << endl << " ... skipping" << endl;
        continue;
      }
      else
        nFilesRead++;

      // open file and get the limits
      TFile*   infile  = TFile::Open(filename,"READ");
      TCanvas* canvas  = (TCanvas*) infile->Get("limcanvas");
      TGraph*  graph   = (TGraph*)  canvas->FindObject(isDecay? "grtaulim": "grsvlim");

      Char_t line[256];
      // First write the masses of the simulation
      for(Int_t ival=0;ival<graph->GetN();ival++)
        {
          Double_t mass = graph->GetX()[ival];
          for(Int_t imass=0;imass<nmass;imass++)
            if(mass==massval[imass])
              {
                results2txt << massval[imass] << " ";
                break;
              }
        }
      results2txt << endl;
      // Then write the limits of the simulation
      for(Int_t ival=0;ival<graph->GetN();ival++)
        {
          Double_t mass = graph->GetX()[ival];
          for(Int_t imass=0;imass<nmass;imass++)
            if(mass==massval[imass])
              {
                Double_t lim  = graph->GetY()[ival];
                svLimDist[imass].push_back(lim);
                results2txt << lim << " ";
                break;
              }
        }
      results2txt << endl;
      delete infile;
    }
  cout << "Could read " << nFilesRead << "/" << nSimuFiles << " of the input files." << endl << endl;

  // if we miss 5% or more of the simulation files
  if((Double_t)nFilesRead/(Double_t)nSimuFiles <0.95)
    {
      cout << "Too few MC files could be read compare to the expected number of files. Did you really run " << nSimuFiles << " simulations?" << endl;
      exit(1);
    }

  cout << "Computing the intervals for all masses..." << endl;
  // sort elements in the limit vector and put values in histos
  Int_t shift=0;
  TH1D* hsvLimDist[nmass];
  for(Int_t imass=0;imass<nmass;imass++)
    {
      hsvLimDist[imass]=NULL;
      // skip mass values with no associated data
      if(!svLimDist[imass].size())
        {
          shift = imass+1;
          continue;
        }

      // sort elements and fill them into histogram
      sort(svLimDist[imass].begin(),svLimDist[imass].end());
      hsvLimDist[imass] = new TH1D(Form("hsvLimDist%02d",imass),Form("hsvLimDist%02d",imass),nbins,svLimDist[imass][0]*0.99,svLimDist[imass][svLimDist[imass].size()-1]*1.01);
      hsvLimDist[imass]->SetStats(0);

      for(UInt_t ival=0;ival<svLimDist[imass].size();ival++)
        hsvLimDist[imass]->Fill(svLimDist[imass][ival]);

      // find the percentiles
      Double_t sv2sigmaLDist = 9e99;
      Double_t sv1sigmaLDist = 9e99;
      Double_t sv0sigmaDist  = 9e99;
      Double_t sv1sigmaRDist = 9e99;
      Double_t sv2sigmaRDist = 9e99;
      Double_t norm = hsvLimDist[imass]->Integral("width");
      for(Int_t ibin=0;ibin<nbins;ibin++)
        {
          Double_t perc = hsvLimDist[imass]->Integral(1,ibin+1,"width")/norm;
          if(TMath::Abs(perc-0.02275)<sv2sigmaLDist) {sv2sigmaLDist=TMath::Abs(perc-0.02275);sv2sigmaL[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.1587)< sv1sigmaLDist) {sv1sigmaLDist=TMath::Abs(perc-0.1587); sv1sigmaL[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.5)<    sv0sigmaDist)  {sv0sigmaDist= TMath::Abs(perc-0.5);    sv0sigma[imass]= hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.8414)< sv1sigmaRDist) {sv1sigmaRDist=TMath::Abs(perc-0.8414); sv1sigmaR[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.97725)<sv2sigmaRDist) {sv2sigmaRDist=TMath::Abs(perc-0.97725);sv2sigmaR[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
        }
      cout << "Mass = " << massval[imass] << ", read " << svLimDist[imass].size() << " values, limits: " <<sv2sigmaL[imass] << " " << sv1sigmaL[imass]<< " "<< sv0sigma[imass]<< " " <<  sv1sigmaR[imass] << " " << sv2sigmaR[imass] << endl;
    }

  // write results in a better format for using them later
  cout << endl;
  cout << "READY-TO-COPY results: " << endl;
  cout << endl;

  cout << "Double_t sv0sigma_"+channel+"[nmass]  = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv0sigma[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv1sigmaL_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv1sigmaL[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv2sigmaL_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv2sigmaL[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv1sigmaR_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv1sigmaR[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv2sigmaR_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv2sigmaR[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  // compute the branon translation factor for each DM mass for branon
  Double_t bt0sigma[nmass];
  Double_t bt1sigmaL[nmass];
  Double_t bt2sigmaL[nmass];
  Double_t bt1sigmaR[nmass];
  Double_t bt2sigmaR[nmass];
  if(!channel.CompareTo("branon",TString::kIgnoreCase))
    {
      // initialize the array of translation factors for the brane tension limits
      Double_t braneTensionVal[nmass];
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          Float_t mdm = massval[imass];
          // release the memory of channelval and brval
          delete [] channelval;
          delete [] brval;
          // number of considering channels
          nChannels = 9;
          channelval = new TString[nChannels];
          brval = new Double_t[nChannels];
          // call the function that computes the branching ratios and save them in channelval and brval
          compute_branonBR(mdm,nChannels,channelval,brval,braneTensionVal[imass]);
        }

      // To translate the sigmav values from [cm^3 s^-1] to [GeV^-2], the sigmav limits have to be divided by
      // h_bar^2 c^3 = (6.582119569^-25 [GeV s])^2 * (299792458x10^2 [cm s^-1])^3 = 1.167329990*10^-17 [GeV^2 cm^3 s^-1]
      // with following the values provided by the Particle Data Group in http://pdg.lbl.gov/2019/reviews/rpp2018-rev-phys-constants.pdf
      // h_bar = 6.582119569*10^-22 [MeV s] = 6.58211956910*^-25 [GeV s] and c = 299792458.0 [m s^−1] = 299792458.0*10^2 [cm s^−1].
      // The final conversion formula also contains:
      // - a factor (0.001) to translate the brane tension limit from [GeV] to [TeV]
      // - a power 1/8 to relate sigma to f (see equation 7 of https://arxiv.org/abs/hep-ph/0302041)
      Double_t unit_translation = 1.167329990*TMath::Power(10., -17.);

      cout << "Double_t bt0sigma[nmass]  = {";
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          bt0sigma[imass] = 0.001*TMath::Power((braneTensionVal[imass]*unit_translation)/sv0sigma[imass], 1./8.);
          cout << bt0sigma[imass] << (imass<nmass-1? "," : "");
        }
      cout << "};" << endl;

      cout << "Double_t bt1sigmaL[nmass] = {";
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          bt1sigmaL[imass] = 0.001*TMath::Power((braneTensionVal[imass]*unit_translation)/sv1sigmaL[imass], 1./8.);
          cout << bt1sigmaL[imass] << (imass<nmass-1? "," : "");
        }
      cout << "};" << endl;

      cout << "Double_t bt2sigmaL[nmass] = {";
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          bt2sigmaL[imass] = 0.001*TMath::Power((braneTensionVal[imass]*unit_translation)/sv2sigmaL[imass], 1./8.);
          cout << bt2sigmaL[imass] << (imass<nmass-1? "," : "");
        }
      cout << "};" << endl;

      cout << "Double_t bt1sigmaR[nmass] = {";
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          bt1sigmaR[imass] = 0.001*TMath::Power((braneTensionVal[imass]*unit_translation)/sv1sigmaR[imass], 1./8.);
          cout << bt1sigmaR[imass] << (imass<nmass-1? "," : "");
        }
      cout << "};" << endl;

      cout << "Double_t bt2sigmaR[nmass] = {";
      for(Int_t imass=shift;imass<nmass;imass++)
        {
          bt2sigmaR[imass] = 0.001*TMath::Power((braneTensionVal[imass]*unit_translation)/sv2sigmaR[imass], 1./8.);
          cout << bt2sigmaR[imass] << (imass<nmass-1? "," : "");
        }
      cout << "};" << endl;
    }

  // Save results as graphs in a root file
  TString outputFileName = resultsPath+"root/"+label+"_1and2sigmasBands.root";
  TFile*  outputFile     = new TFile(outputFileName,"RECREATE");

  TGraph* gsv0sigma  = new TGraph(nmass-shift,massval+shift,sv0sigma+shift);
  TGraph* gsv1sigmaL = new TGraph(nmass-shift,massval+shift,sv1sigmaL+shift);
  TGraph* gsv1sigmaR = new TGraph(nmass-shift,massval+shift,sv1sigmaR+shift);
  TGraph* gsv2sigmaL = new TGraph(nmass-shift,massval+shift,sv2sigmaL+shift);
  TGraph* gsv2sigmaR = new TGraph(nmass-shift,massval+shift,sv2sigmaR+shift);
  gsv0sigma->SetName("gsv0sigma");
  gsv1sigmaL->SetName("gsv1sigmaL");
  gsv1sigmaR->SetName("gsv1sigmaR");
  gsv2sigmaL->SetName("gsv2sigmaL");
  gsv2sigmaR->SetName("gsv2sigmaR");
  gsv0sigma->Write();
  gsv1sigmaL->Write();
  gsv1sigmaR->Write();
  gsv2sigmaL->Write();
  gsv2sigmaR->Write();

  cout << endl << "Bands saved in file: " << outputFileName << endl;
  delete outputFile;

  // Print statement that the results are saved in the txt file
  cout << "Results saved in file: " << results2txtPath << endl;

  // Save results as graphs in a root file for branon
  TGraph* gbt0sigma  = NULL;
  TGraph* gbt1sigmaL = NULL;
  TGraph* gbt1sigmaR = NULL;
  TGraph* gbt2sigmaL = NULL;
  TGraph* gbt2sigmaR = NULL;
  if(!channel.CompareTo("branon",TString::kIgnoreCase))
    {
      TString btoutputFileName = resultsPath+"root/"+label+"_branetension_1and2sigmasBands.root";
      TFile*  btoutputFile     = new TFile(btoutputFileName,"RECREATE");

      gbt0sigma  = new TGraph(nmass-shift,massval+shift,bt0sigma+shift);
      gbt1sigmaL = new TGraph(nmass-shift,massval+shift,bt1sigmaL+shift);
      gbt1sigmaR = new TGraph(nmass-shift,massval+shift,bt1sigmaR+shift);
      gbt2sigmaL = new TGraph(nmass-shift,massval+shift,bt2sigmaL+shift);
      gbt2sigmaR = new TGraph(nmass-shift,massval+shift,bt2sigmaR+shift);
      gbt0sigma->SetName("gbt0sigma");
      gbt1sigmaL->SetName("gbt1sigmaL");
      gbt1sigmaR->SetName("gbt1sigmaR");
      gbt2sigmaL->SetName("gbt2sigmaL");
      gbt2sigmaR->SetName("gbt2sigmaR");
      gbt0sigma->Write();
      gbt1sigmaL->Write();
      gbt1sigmaR->Write();
      gbt2sigmaL->Write();
      gbt2sigmaR->Write();

      cout << "Bands saved in file: " << btoutputFileName << endl;
      delete btoutputFile;
    }

  // PLOTS
  if(makePlots)
    {
      // style
      setDefaultStyle();
      gStyle->SetPadRightMargin(0.12);
      gStyle->SetPadLeftMargin(0.12);

      // canvas for distributions
      const Int_t nlines = 5;  // number of lines in parabolas canvas
      Int_t ncols = TMath::Ceil((nmass-shift)/Float_t(nlines)); // number of columns in parabolas canvas
      TCanvas* c = new TCanvas("canvas","limits distributions",ncols*250,nlines*250);
      c->Divide(ncols,nlines);

      // loop over masses and plot histos
      for(Int_t imass=shift;imass<nmass;imass++)
        {         
          c->cd(imass+1);
          hsvLimDist[imass]->DrawCopy("e");
          gPad->SetLogy();
          gPad->SetLogx();
        }
      c->Print(resultsPath+"root/"+label+"_distributionLimits.root");
      c->Print(resultsPath+"pdf/"+label+"_distributionLimits.pdf");

      // bands
      Int_t npnts = gsv0sigma->GetN();
      TGraph* band1sigma = new TGraph(2*npnts+1);
      band1sigma->SetNameTitle("band1sigma","band1sigma");
      TGraph* band2sigma = new TGraph(2*npnts+1);
      band2sigma->SetNameTitle("band2sigma","band2sigma");

      for(Int_t imass=0;imass<npnts;imass++)
        {
          band1sigma->SetPoint(imass,      gsv1sigmaL->GetX()[imass],gsv1sigmaL->GetY()[imass]);
          band1sigma->SetPoint(npnts+imass,gsv1sigmaR->GetX()[npnts-imass-1],gsv1sigmaR->GetY()[npnts-imass-1]);
          band2sigma->SetPoint(imass,      gsv2sigmaL->GetX()[imass],gsv2sigmaL->GetY()[imass]);
          band2sigma->SetPoint(npnts+imass,gsv2sigmaR->GetX()[npnts-imass-1],gsv2sigmaR->GetY()[npnts-imass-1]);
        }
      band1sigma->SetPoint(2*npnts,gsv1sigmaL->GetX()[0],gsv1sigmaL->GetY()[0]);
      band2sigma->SetPoint(2*npnts,gsv2sigmaL->GetX()[0],gsv2sigmaL->GetY()[0]);

      // result
      TFile* f = TFile::Open(dataPath+label+"_Data_limits.root", "READ");
      TCanvas* lim = (TCanvas*)f->Get("limcanvas");
      TGraph* grsvlim = (TGraph*)lim->GetPrimitive(isDecay? "grtaulim": "grsvlim");
      grsvlim->SetName(isDecay? "grtaulim": "grsvlim");
      grsvlim->SetLineColor(kBlack);
      grsvlim->SetLineWidth(2);

      band1sigma->SetFillColorAlpha(3,0.35);
      band2sigma->SetFillColorAlpha(5,0.35);
      gsv0sigma->SetLineColor(4);
      gsv0sigma->SetLineStyle(2);
      gsv0sigma->SetLineWidth(1);

      // canvas
      TCanvas* limcanvas = new TCanvas("limcanvas","1 and 2 sigma bands",1000,800);
      limcanvas->cd();

      TH1I *dummylim = new TH1I("dummylim","",1,massval[0],massval[nmass-1]);
      dummylim->SetStats(0);
      if(!channel.CompareTo("gammagamma",TString::kIgnoreCase))
        {
          dummylim->SetMinimum((isDecay? 1e23 : 1e-33));
          dummylim->SetMaximum((isDecay? 1e27 : 1e-22));
        }
      else
        {
          dummylim->SetMinimum((isDecay? 1e23 : 1e-28));
          dummylim->SetMaximum((isDecay? 1e27 : 1e-20));
        }
      dummylim->SetXTitle("m_{DM} [GeV]");
      dummylim->SetYTitle(Form("95%% %s [%s]",(isDecay? "#tau_{DM}^{LL}" : "<#sigma v>^{UL}"),(isDecay? "s" : "cm^{3}/s")));
      dummylim->DrawCopy();
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetGrid();

      band2sigma->Draw("f");
      band1sigma->Draw("f");
      gsv0sigma->Draw("c");
      grsvlim->Draw("c");

      // thermal relic density taken from Steigman G., Dasgupta B, and Beacom J. F., 
      // Precise relic WIMP abundance and its impact onsearches for dark matter annihilation, 
      // Phys.Rev. D86(2012) 023506, [arXiv:1204.3622]
      const Int_t nThermalRelicSv = 22;
      Double_t thermalRelicMass[nThermalRelicSv] = {1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e+00, 1.78e+00, 3.16e+00, 5.62e+00, 1.00e+01, 1.78e+01, 3.16e+01, 5.62e+01, 1.00e+02, 1.78e+02, 3.16e+02, 5.62e+02, 1.00e+03, 1.78e+03,3.16e+03, 5.62e+03, 1.00e+04, 1.00e+05};
      Double_t thermalRelicSigmav[nThermalRelicSv] = {4.8e-26, 4.9e-26, 5.1e-26, 5.0e-26, 4.7e-26, 4.5e-26, 3.9e-26, 2.8e-26, 2.5e-26, 2.3e-26, 2.2e-26, 2.2e-26, 2.2e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26,2.4e-26, 2.4e-26};

      if(!channel.CompareTo("gammagamma",TString::kIgnoreCase))
        {
          Double_t alpha = 0.0072973525664; //fine-structure constant
          for(int irelic=0;irelic<nThermalRelicSv;irelic++)
            thermalRelicSigmav[irelic]=thermalRelicSigmav[irelic]*alpha*alpha;
        }

      TGraph *relicDensity = new TGraph(nThermalRelicSv, thermalRelicMass, thermalRelicSigmav);
      relicDensity->SetTitle("Thermal relic cross section");
      relicDensity->SetLineColor(kRed);
      relicDensity->SetLineWidth(3);
      relicDensity->SetLineStyle(5);
      relicDensity->SetMarkerColor(kRed);
      if(!isDecay) relicDensity->Draw("same");

      // legend
      TLegend* limleg;
      limleg = new TLegend(0.35, 0.14, 0.85, 0.34);
      limleg->SetTextSize(0.032);
      limleg->SetMargin(0.40);
      limleg->SetBorderSize(0);

      if(grsvlim)    limleg->AddEntry(grsvlim,"Limit","l");
      if(gsv0sigma)  limleg->AddEntry(gsv0sigma,"H_{0} median","l");
      if(band1sigma) limleg->AddEntry(band1sigma,"H_{0} 68% containment","f");
      if(band2sigma) limleg->AddEntry(band2sigma,"H_{0} 95% containment","f");
      if(!isDecay)   limleg->AddEntry(relicDensity,"Thermal relic cross section","l");
      limleg->Draw();

      TLatex* txchannel = new TLatex(0.75,0.4,strchannel);
      txchannel->SetTextSize(0.05);
      txchannel->SetNDC();
      txchannel->Draw();

      limcanvas->Print(resultsPath+"root/"+label+"_bands.root");
      limcanvas->Print(resultsPath+"pdf/"+label+"_bands.pdf");

      if(!channel.CompareTo("branon",TString::kIgnoreCase))
        {
          // bands
          Int_t npnts = gbt0sigma->GetN();
          TGraph* btband1sigma = new TGraph(2*npnts+1);
          TGraph* btband2sigma = new TGraph(2*npnts+1);

          for(Int_t imass=0;imass<npnts;imass++)
            {
              btband1sigma->SetPoint(imass,      gbt1sigmaL->GetX()[imass],gbt1sigmaL->GetY()[imass]);
              btband1sigma->SetPoint(npnts+imass,gbt1sigmaR->GetX()[npnts-imass-1],gbt1sigmaR->GetY()[npnts-imass-1]);
              btband2sigma->SetPoint(imass,      gbt2sigmaL->GetX()[imass],gbt2sigmaL->GetY()[imass]);
              btband2sigma->SetPoint(npnts+imass,gbt2sigmaR->GetX()[npnts-imass-1],gbt2sigmaR->GetY()[npnts-imass-1]);
            }
          btband1sigma->SetPoint(2*npnts,gbt1sigmaL->GetX()[0],gbt1sigmaL->GetY()[0]);
          btband2sigma->SetPoint(2*npnts,gbt2sigmaL->GetX()[0],gbt2sigmaL->GetY()[0]);

          // result
          TFile* btf = TFile::Open(dataPath+label+"_Data_branetension_limits.root", "READ");
          TCanvas* btlim = (TCanvas*)btf->Get("branoncanvas");
          TGraph* grbtlim = (TGraph*)btlim->GetPrimitive("grbtlim");

          grbtlim->SetName("grbtlim");
          grbtlim->SetLineColor(kBlack);
          grbtlim->SetLineWidth(2);

          btband1sigma->SetFillColorAlpha(3,0.35);
          btband2sigma->SetFillColorAlpha(5,0.35);

          gbt0sigma->SetLineColor(4);
          gbt0sigma->SetLineStyle(2);
          gbt0sigma->SetLineWidth(1);

          // canvas
          TCanvas* branoncanvas = new TCanvas("branoncanvas","1 and 2 sigma bands",1000,800);
          branoncanvas->cd();

          TH1I *btdummylim = new TH1I("btdummylim","",1,massval[0],massval[nmass-1]);
          btdummylim->SetStats(0);
          btdummylim->SetMinimum(1e-1);
          btdummylim->SetMaximum(1e2);
          btdummylim->SetXTitle("m_{DM} [GeV]");
          btdummylim->SetYTitle(Form("f [TeV]"));
          btdummylim->DrawCopy();
          gPad->SetLogx();
          gPad->SetLogy();
          gPad->SetGrid();

          btband2sigma->Draw("f");
          btband1sigma->Draw("f");
          gbt0sigma->Draw("c");
          grbtlim->Draw("c");

          // thermal relic density taken from Steigman G., Dasgupta B, and Beacom J. F., 
          // Precise relic WIMP abundance and its impact onsearches for dark matter annihilation, 
          // Phys.Rev. D86(2012) 023506, [arXiv:1204.3622]. The thermal relic cross section was 
          // translated to the brane tension space using the formulas above.
          const Int_t nThermalRelicBt = 22;
          Double_t thermalRelicMass[nThermalRelicBt] = {1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e+00, 1.78e+00, 3.16e+00, 5.62e+00, 1.00e+01, 1.78e+01, 3.16e+01, 5.62e+01, 1.00e+02, 1.78e+02, 3.16e+02, 5.62e+02, 1.00e+03, 1.78e+03,3.16e+03, 5.62e+03, 1.00e+04, 1.00e+05};
          Double_t thermalRelicBt[nThermalRelicBt] = {0.00028197, 0.00131177, 0.00184612, 0.00250686, 0.00338571, 0.00742158, 0.0124374, 0.0199107, 0.0291449, 0.0400859, 0.054027, 0.0721836, 0.208693, 0.367468, 0.578445, 0.892658, 1.37518, 2.11902, 3.25889, 5.01881, 7.69107, 43.25};
          TGraph *btrelicDensity = new TGraph(nThermalRelicBt, thermalRelicMass, thermalRelicBt);
          btrelicDensity->SetTitle("Thermal relic cross section");
          btrelicDensity->SetLineColor(kRed);
          btrelicDensity->SetLineWidth(3);
          btrelicDensity->SetLineStyle(5);
          btrelicDensity->SetMarkerColor(kRed);
          btrelicDensity->Draw("same");

          // legend
          TLegend* btlimleg;
          btlimleg = new TLegend(0.2, 0.7, 0.45, 0.85);
          btlimleg->SetTextSize(0.032);
          btlimleg->SetMargin(0.40);
          btlimleg->SetBorderSize(0);

          if(grbtlim)    btlimleg->AddEntry(grbtlim,"Limit","l");
          if(gbt0sigma)  btlimleg->AddEntry(gbt0sigma,"H_{0} median","l");
          if(btband1sigma) btlimleg->AddEntry(btband1sigma,"H_{0} 68% containment","f");
          if(btband2sigma) btlimleg->AddEntry(btband2sigma,"H_{0} 95% containment","f");
          btlimleg->AddEntry(btrelicDensity,"Thermal relic cross section","l");
          btlimleg->Draw();

          branoncanvas->Print(resultsPath+"root/"+label+"_branetension_bands.root");
          branoncanvas->Print(resultsPath+"pdf/"+label+"_branetension_bands.pdf");
        }
    }
}

void usage(){
  cout << "computeCLBands usage:" << endl;
  cout << "\t -h or --help: display this message and quit" << endl;
  cout << "\t --config: configuration file, default jointLklDM.rc" << endl;
  cout << "\t --nsimufiles: number of simulations, default 300" << endl;
}

Int_t GetNSkippedMasses(Int_t nm,const Double_t* vm,Double_t minm)
{
  Int_t im=0;
  for(im=0;im<nm;im++)
    if(vm[im]>minm)
      break;
  return im;
}

void setDefaultStyle()
{
  // general settings
  //gStyle->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTextSize(1);
  gStyle->SetTitleSize(0.04, "xy");
  gStyle->SetTitleOffset(1.3, "x");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetOptStat(111110);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  //gStyle->SetColorModelPS(1);
  gStyle->SetPalette(1,0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetGridWidth(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTickLength(0.04, "xy");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetLegendBorderSize(4);
  gStyle->SetGridColor(13);
}

// Decode the channel and save the decoded channels and branching ratios in channelval and brval, respectively.
void decode_channel(TObjArray* coefficients, Int_t &nChannels, TString *channelval, Double_t *brval)
{
  Float_t sumBR = 0.0;
  TString factorstring, coefficientstring;
  for (Int_t iChannel = 0; iChannel < nChannels; iChannel++)
    {
      factorstring = (TString)((TObjString *)(coefficients->At(iChannel)))->String();
      TObjArray *factors = factorstring.Tokenize("*");
      // if there is only one channel selected
      if(nChannels == 1)
        {
          // save the decoded channel in channelval
          if(factors->GetEntries() == 1)       channelval[iChannel] = factorstring;
          else if (factors->GetEntries() == 2) channelval[iChannel] = (TString)((TObjString *)(factors->At(1)))->String();
          // check if the channel argument is in the right form
          else
            {
              cout << " ## Oops! Something went wrong with the parsing " << factorstring << " of the channel!  <---------------- FATAL ERROR!!!"<< endl;
              return;
            }
          // set the branching ratio to 1.0 (100%), since there is only one channel selected
          brval[iChannel] = sumBR = 1.0;
        }
      // if there is a linear combination of channels selected
      else
        {
          // check if the channel argument is in the right form
          if(factors->GetEntries() == 1)
            {
              channelval[iChannel] = factorstring;
              brval[iChannel] = 1.0;
              sumBR += brval[iChannel];
            }
          else if (factors->GetEntries() == 2)
            {
              // save the decoded channels and branching ratios in channelval and brval, respectively.
              for (Int_t i = 0; i < factors->GetEntries(); i++)
                {
                  coefficientstring = (TString)((TObjString *)(factors->At(i)))->String();
                  if (coefficientstring.IsFloat() && i == 0)
                    {
                      brval[iChannel] = coefficientstring.Atof();
                      sumBR += brval[iChannel];
                    }
                  else
                    channelval[iChannel] = coefficientstring;
                }
            }
          // check if the channel argument is in the right form
          else
            {
              cout << " ## Oops! Something went wrong with the parsing " << factorstring << " of the channel!  <---------------- FATAL ERROR!!!"<< endl;
              return;
            }
        }
    }
  // Normalize the branching ratios if they differ from 1.0 (100%)
  if (sumBR != 1.0)
      for (Int_t iChannel = 0; iChannel < nChannels; iChannel++) brval[iChannel] /= sumBR;
}

// Compute the branching ratios for the branon model for each DM mass and save the channels and branching ratios in channelval and brval, respectively.
void compute_branonBR(Float_t &mass, Int_t &nChannels, TString *channelval, Double_t *brval, Double_t &translation_factor)
  {
    // The included annihilation channels for the branon model
    TString particle_type[9]        = {"bb","cc","tt","ee","mumu","tautau","WW","ZZ","hh"};
    // The corresponding masses in GeV
    Float_t particle_mass[9]        = {4.18,1.28,173.1,0.511e-3,0.106,1.7768,80.39,91.19,125.0};
    // Classification of the elementary particles
    TString dirac_fermions          = "cc_tt_bb_ee_mumu_tautau";
    TString gauge_bosons            = "WW_ZZ";
    TString scalar_bosons           = "hh";
    // Loop over the channels and compute their branching ratios
    Double_t ann_crosssection[9]    = {0.0};
    Double_t total_ann_crosssection = 0.0;
    for (Int_t iChannel = 0; iChannel < nChannels; iChannel++)
      {
        // distinguish between the different elementary particles
        if(mass >= particle_mass[iChannel])
          {
            if (dirac_fermions.Contains(particle_type[iChannel]))
              ann_crosssection[iChannel] = (mass*mass * particle_mass[iChannel]*particle_mass[iChannel])/(16. * TMath::Pi()*TMath::Pi()) * (mass*mass - particle_mass[iChannel]*particle_mass[iChannel]) * TMath::Sqrt(1-((particle_mass[iChannel]*particle_mass[iChannel])/(mass*mass)));
            else if (gauge_bosons.Contains(particle_type[iChannel]))
              ann_crosssection[iChannel] = (mass*mass)/(64. * TMath::Pi()*TMath::Pi()) * (4. * TMath::Power(mass,4) - 4. * mass*mass * particle_mass[iChannel]*particle_mass[iChannel] + 3. * TMath::Power(particle_mass[iChannel],4)) * TMath::Sqrt(1-((particle_mass[iChannel]*particle_mass[iChannel])/(mass*mass)));
            else if (scalar_bosons.Contains(particle_type[iChannel]))
              ann_crosssection[iChannel] = (mass*mass)/(32. * TMath::Pi()*TMath::Pi()) * TMath::Power((2.* mass*mass + particle_mass[iChannel]*particle_mass[iChannel]),2) * TMath::Sqrt(1-((particle_mass[iChannel]*particle_mass[iChannel])/(mass*mass)));
            // WW with a factor 2 (because the W is complex)
            if(!particle_type[iChannel].CompareTo("WW",TString::kIgnoreCase)) ann_crosssection[iChannel] *= 2.;
            // hh with a factor 1/2 (because the Higgs is real)
            if(!particle_type[iChannel].CompareTo("hh",TString::kIgnoreCase)) ann_crosssection[iChannel] *= 0.5;
          }
        // add up all ann_crosssection for the normalization
        total_ann_crosssection += ann_crosssection[iChannel];
      }
    // Normalization of the branching ratios
    for (Int_t iChannel = 0; iChannel < nChannels; iChannel++)
      {
        ann_crosssection[iChannel] /= total_ann_crosssection;
        // save the computed branching ratios in brval
        brval[iChannel] = ann_crosssection[iChannel];
        // save the corresponding channels in channelval
        channelval[iChannel] = particle_type[iChannel];
      }
    // Computation of the translation factor for the tension of the brane
    translation_factor = total_ann_crosssection;
  }
