//######################################################################
//##
//##                         jointLklDM.C
//##
//##              AUTHOR: Javier Rico (jrico@ifae.es)
//##                        29/03/2017
//##
//##
//## IMPORTANT NOTE: THE USE OF THIS CODE TO PRODUCE PAPERS OF THE MAGIC
//## AND/OR CTA COLLABORATIONS IS ALLOWED FOLLOWING THEIR RESPECTIVE
//## PUBLICATION POLICIES FOR FULL-COLLABORATION PAPERS. FOR
//## PUBLICATIONS OUTSIDE THOSE FRAMEWORKS PLEASE CONTACT FIRST THE
//## AUTHORS (Jelena Aleksic <jelena@ifae.es> AND Javier Rico
//## <mailto:jrico@ifae.es>), WHO COULD CLAIM AUTHORSHIP OF THE
//## RESULTING PAPER.
//##
//## PLEASE CITE:
//## Aleksic, Rico & Martinez JCAP 10 (2012) 032
//##
//## jointLklDM.C
//##
//## Compute limits to Dark Matter annihilation cross-section (<sv>) or
//## decay lifetime (tau) using the full-likelihood approach as described
//## in Aleksic, Rico & Martinez JCAP 10 (2012) 032 and implemented in
//## gLike.
//##
//## You can specify two arguments:
//## 1. Name of the rc file (a TString)
//## 2. A random seed (a Int_t). If the seed is non-negative it is used 
//##    for the geenerator that simulates event energies according to the 
//##    background and signal pdf's, the observation time, signal 
//##    intensity, tau, dTau...
//##
//## The -2logLkl vs g0 curve is evaluated close to
//## its minimum and the limits are obtained from the point where the
//## curve crosses the minimum value plus deltaLogLkl (configurable).
//##
//## The macro produces (and saves) three different sets of plots:
//## - For each sample, a canvas (hadcanvas[isample]) containing the
//##   IRFs, dN_signal/dE and dN_signal/dE' and data plots. To obtain
//##   these plots set showSamplePlots==kTRUE.
//## - A canvas (lklcanvas) with one plot per considered
//##   mass value, where -2logLkl is plotted vs <sv> near the minimum.
//## - A canvas limcanvas with the limits on <sv> or tau vs mass for
//##   the considered DM channel, range of masses and data.
//##
//## You need to run the compiled version of the macro
//## by typing ".x jointLklDM.C+" in the ROOT command interpreter.
//##
//## Input files for MAGIC Stereo Segue observations can be downloaded
//## from 
//## https://dl.dropboxusercontent.com/u/14145932/MAGICStereoSegueIRFandData.tgz
//## and for CTA simulations from
//## http://www.mpi-hd.mpg.de/hfm/CTA/MC/performance-20deg/
//## BE AWARE THESE ARE MAGIC AND CTA INTERNAL FILES, RESPECTIVELY,
//## YOU CANNOT PUBLISH ON YOUR OWN IF YOU USE THEM!!
//##
//######################################################################

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "TPRegexp.h"
#include "TMath.h"
#include "TH1.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"

#include "Iact1dUnbinnedLkl.h"
#include "Iact1dBinnedLkl.h"
#include "FermiTables2016Lkl.h"
#include "GloryDuckTables2019Lkl.h"
#include "LineSearchLkl.h"
#include "JointLkl.h"

using namespace std;

void setDefaultStyle();
void decode_channel(TObjArray* coefficients, Int_t &nChannels, TString *channelval, Double_t *brval);
void compute_branonBR(Float_t &mass, Int_t &nChannels, TString *channelval, Double_t *brval);
Int_t GetNSkippedMasses(Int_t nm,const Double_t* vm,Double_t minm);
void builddNdESignal(Lkl* lkl,TString dNdEDir, Int_t nChannels,TString* channelval,Double_t* brval,Float_t mdm,Bool_t isDecay,Bool_t isSimulation);
TString mprecform(Double_t val);

const Int_t nMaxLkls = 1000;

void jointLklDM(TString configFileName="$GLIKESYS/rcfiles/jointLklDM.rc",Int_t seed=-1)
{
  setDefaultStyle();

  // Look for configuration file
  gSystem->ExpandPathName(configFileName);
  TPMERegexp re("\\s+");
  
  Lkl::PrintGLikeBanner();
  
  // Print-out configuration info (part 1)
  cout << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***                            RUNNING jointLklDM.C                             ***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***" << endl;
  cout << "*** CONFIGURATION FILE       : " << configFileName << endl;
  if (gSystem->AccessPathName(configFileName, kFileExists))
    {
      cout << endl << "    Oops! problems reading config file file " << configFileName << " <---------------- FATAL ERROR!!!"<< endl;
      return;
    }

  // Read configuration file
  TEnv*  env = new TEnv(configFileName);
  TString  label             = env->GetValue("jointLklDM.Label","");
  TString  channel           = env->GetValue("jointLklDM.Channel","bb");
  TString  process           = env->GetValue("jointLklDM.Process","ann");
  Bool_t   isGpositive       = env->GetValue("jointLklDM.isGpositive",kFALSE);
  Bool_t   showSamplePlots   = env->GetValue("jointLklDM.showSamplePlots",kTRUE);
  Bool_t   showParabolaPlots = env->GetValue("jointLklDM.showParabolaPlots",kTRUE);
  Bool_t   showLimitPlots    = env->GetValue("jointLklDM.showLimitPlots",kTRUE);
  Double_t plotmin           = env->GetValue("jointLklDM.plotmin",0.);
  Double_t plotmax           = env->GetValue("jointLklDM.plotmax",0.);
  Double_t plotScale         = env->GetValue("jointLklDM.plotScale",1.);
  Double_t deltaLogLkl       = env->GetValue("jointLklDM.deltaLogLkl",2.71);
  TString  fInputDataPath    = env->GetValue("jointLklDM.path","");
  TString  fdNdEDir          = fInputDataPath+"/"+env->GetValue("jointLklDM.dNdEDir","")+"/";
  TString  fPlotsDir         = fInputDataPath+"/"+env->GetValue("jointLklDM.plotsDir","")+"/";
  TString  provval           = env->GetValue("jointLklDM.dNdEpSignalDir","-");
  TString  fdNdEpSignalDir   = fInputDataPath+"/"+provval+"/";
  Float_t  mcsv              = env->GetValue("jointLklDM.mcsv",0.);    // value of <sv> assumed for simulations
  Float_t  mcmass            = env->GetValue("jointLklDM.mcmass",0.);  // value of DM mass assumed for simulations
  TString  mcChannel         = env->GetValue("jointLklDM.mcChannel","bb");
  TString  mcProcess         = env->GetValue("jointLklDM.mcProcess","ann");
  TString  massList          = env->GetValue("jointLklDM.MassList","");
  Bool_t   exportData        = env->GetValue("jointLklDM.exportData",kFALSE);
  Double_t svMin             = env->GetValue("jointLklDM.exportSvMin",0.);
  Double_t svMax             = env->GetValue("jointLklDM.exportSvMax",0.);
  Bool_t   svLogStep         = env->GetValue("jointLklDM.exportSvLogStep",kTRUE);
  Int_t    svNPoints         = env->GetValue("jointLklDM.exportSvNPoints",0.);  
  TString  fExportDataPath   = fInputDataPath+"/"+env->GetValue("jointLklDM.exportDataPath","")+"/";
      
  // fill the list of masses to be studied
  UInt_t  nmass0  = re.Split(massList);
  Double_t* massval0 = new Double_t[nmass0];
  for(UInt_t imass=0;imass<nmass0;imass++)
    massval0[imass] = re[imass].Atof();

  // Set some flags
  Bool_t  isSimulation       = seed>=0;
  Bool_t  ioHdNdEpSignal     = provval.CompareTo("-");
  Bool_t  isDecay            = !process.CompareTo("dec",TString::kIgnoreCase);
  
  // define some labels according to input data
  TString simulationlabel    = (isSimulation?  "MC"    : "Data");
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
  cout << "*** I/O PATH                 : " << fInputDataPath << endl;
  cout << "***" << endl;
  cout << "*** LABEL                    : " << label <<  endl;
  cout << "*** PROCESS                  : " << strprocess << endl;
  cout << "*** CHANNEL                  : " << channel << endl;
  // how many and which channels are we considering?
  cout << " ** Channel check            : " << normchannel << endl;
  cout << "*** DATA/MC                  : " << simulationlabel << endl;
  if(isSimulation)
    {
      cout << " ** Seed                     : " << seed << endl;
      cout << " ** <sv>_MC                  : " << mcsv << " cm^3/s" << endl;
      cout << " ** m_MC                     : " << mcmass << " GeV" << endl;
    }
  cout << "*** G IS POSITIVE            : " << (isGpositive? "YES" : "NO") << endl;
      
  // how many and which mass values are we considering?
  cout << "*** NUMBER OF MASSES         : " << nmass << endl;
  cout << " ** Mass values              : "; 
  for(Int_t imass=0;imass<nmass;imass++)
    cout << massval[imass] << ((imass<nmass-1)? ", " : "");
  cout << " GeV" << endl;
  
  cout << "*** IRF/DATA PLOTS      : " << (showSamplePlots?   "YES" : "NO") << endl;
  cout << "*** PARABOlA PLOTS      : " << (showParabolaPlots? "YES" : "NO") << endl;
  cout << "*** Plot y-axis range   : " << plotmin << " to " << plotmax << Form(" %s", (isDecay? "s" : "cm^3/s")) << endl;
  cout << "*** EXPORT DATA         : " << (exportData? "YES" : "NO") << endl;
  if(exportData)
    cout << " ** Format of export    : Glory Duck" << endl;

  cout << "***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << endl;
  // End of print-out configuration info

  // Arrays of Lkl parabolas (one per mass value) to be plotted is showParabolaPlots=kTRUE
  TGraph*  grLklParabola[nmass];
  for(Int_t imass=0;imass<nmass;imass++)
    grLklParabola[imass] = NULL;
  const Int_t nlines = 5;  // number of lines in parabolas canvas
  Int_t ncols = TMath::Ceil(nmass/Float_t(nlines)); // number of columns in parabolas canvas

  // Define logJ variable to check if the exportData option can be applied
  Double_t logJ = 0.;

  // The classes composing the joint likelihood
  Lkl** lkl    = new Lkl*[nMaxLkls];
  Lkl** sample = new Lkl*[nMaxLkls];
  
  // read rc file for likelihood construction
  cout << "****************************************************************" << endl;
  cout << "*** CONFIGURING Lkl OBJECTS " << endl;
  Int_t nLkls = 0;
  Int_t nsamples = 0;
  for(Int_t iLkl=0;iLkl<nMaxLkls;iLkl++)
    {
      // look for likelihood entries in the rc file
      TString lklString = env->GetValue(Form("jointLklDM.lklTerm%03d",iLkl),"");
      if(lklString.CompareTo("")==0) break; // if one index is missing, the search is stopped
      else nLkls++;
      
      // lklString should read: Type ParentTerm options
      UInt_t   nfields   = re.Split(lklString,3);
      if(nfields<2)
	{
	  cout << Form("jointLklDM Error: rc file entry for lklTerm%03d has %d entries (minimum of 2 is needed) <---------------- FATAL ERROR!!!",iLkl,nfields);
	  return;
	}
      
      TString  classType   = re[0];
      Int_t    parLkl      = (iLkl>0? re[1].Atoi() : -1);
      TString  inputString = (nfields>2 ? re[2] : "missingFileName");

      // Start reporting what we have read
      cout << "***" << endl;
      cout << "*** Found lkl Term #" << iLkl  << endl;
      cout << " ** Class type      : " << classType << endl;
      cout << " ** Parent term     : " << (parLkl<0? "NONE" : Form("%d",parLkl)) << endl;
      cout << " ** Input string    : " << inputString << endl;

      // add path
      inputString+=(" path="+fInputDataPath);

      // create the proper object according to classType 
      if(classType.CompareTo("Iact1dUnbinnedLkl")==0)
	{	  
	  lkl[iLkl] =  new Iact1dUnbinnedLkl(inputString);
	  lkl[iLkl]->SetName(Form("Iact1dUnbinnedLkl_%02d",iLkl));

	  // save as sample (as opposed to JointLkl)
	  sample[nsamples++] = lkl[iLkl];

          // if it's simulation of the null hypothesis, we can already simulate the data sample
          if(isSimulation && mcsv==0)
	    {
              cout << "  * Simulating observations for <sv> = " << mcsv << endl;
	      if(dynamic_cast<Iact1dUnbinnedLkl*>(lkl[iLkl])->SimulateDataSamples(seed,mcsv))
		{
		  cout << " ## Oops! Cannot simulate samples for " << lkl[iLkl]->GetName() << " <---------------- FATAL ERROR!!!"<< endl;
		  return;
		}
	    }
	  else
	    if(dynamic_cast<Iact1dUnbinnedLkl*>(lkl[iLkl])->GetNon()<1)
	      {
		  cout << " ## Oops! No data (from input or simulated) associated to " << lkl[iLkl]->GetName() << " <---------------- FATAL ERROR!!!"<< endl;
		  return;
	      }	 	      
	}
      else if(classType.CompareTo("Iact1dBinnedLkl")==0)
	{
	  // read input string
	  lkl[iLkl] =  new Iact1dBinnedLkl(inputString);

	  // save as sample (as opposed to JointLkl)
	  sample[nsamples++] = lkl[iLkl];

	  // configure
	  lkl[iLkl]->SetName(Form("Iact1dBinnedLkl_%02d",iLkl));	  

          // if it's simulation of the null hypothesis, we can already simulate the data sample
          if(isSimulation && mcsv==0)
	    {
              cout << "  * Simulating observations for <sv> = " << mcsv << endl;
	      if(dynamic_cast<Iact1dBinnedLkl*>(lkl[iLkl])->SimulateDataSamples(seed,mcsv))
		{
		  cout << " ## Oops! Cannot simulate samples for " << lkl[iLkl]->GetName() << " <---------------- FATAL ERROR!!!"<< endl;
		  return;
		}
	    }
	  else
	    if(dynamic_cast<Iact1dBinnedLkl*>(lkl[iLkl])->GetNon()<1)
	      {
		cout << " ## Oops! No data (from input or simulated) associated to " << lkl[iLkl]->GetName() << " <---------------- FATAL ERROR!!!"<< endl;
		return;
	      }	 
	}
      else if(classType.CompareTo("LineSearchLkl")==0)
	{
	  if(channel.CompareTo("gammagamma",TString::kIgnoreCase))
            {
              cout << " ## Oops! The LineSearchLkl class can only be used for the channel gammagamma and cannot be used for the channel " << channel << " <---------------- FATAL ERROR!!!"<< endl;
              return;
            }

	  // read input string
	  lkl[iLkl] =  new LineSearchLkl(inputString);

	  // save as sample (as opposed to JointLkl)
	  sample[nsamples++] = lkl[iLkl];

	  // configure
	  lkl[iLkl]->SetName(Form("LineSearchLkl_%02d",iLkl));
        }
      else if(classType.CompareTo("JointLkl")==0)
	{	  
	  lkl[iLkl] =  new JointLkl(inputString);
	  lkl[iLkl]->SetName(Form("JointLkl_%02d",iLkl));
	}
      else if(classType.CompareTo("FermiTables2016Lkl")==0)
	{
	  lkl[iLkl] =  new FermiTables2016Lkl(inputString);
	  lkl[iLkl]->SetName(Form("FermiTables2016Lkl_%02d",iLkl));
	}
      else if(classType.CompareTo("GloryDuckTables2019Lkl")==0)
	{
          if(isSimulation)
            {
              TPMERegexp re(".txt");
              UInt_t nfields = re.Split(inputString);
              inputString = re[0] + Form("_%d.txt",seed) + re[1];
              cout << "  * Modified        : " << inputString << endl;
            }
	  GloryDuckTables2019Lkl *tmpLkl =  new GloryDuckTables2019Lkl(inputString);
          Bool_t massAvailable = kFALSE;
          for(Int_t jmass=0;jmass<nmass;jmass++)
            {
              massAvailable = kFALSE;
              Double_t massToBeTested = massval[jmass];
              for(UInt_t kmass=0;kmass<tmpLkl->GetNMasses();kmass++)
                {
                  Double_t massInTheFile = tmpLkl->GetActiveMass(kmass);
                  if(TMath::Abs(massToBeTested-massInTheFile) < 1.e-6)
                    massAvailable = kTRUE;
                }
            }
          if(!massAvailable)
            {
              cout << " ## Oops! At least one mass to be tested is not in the file given as an input <---------------- FATAL ERROR!!!"<< endl;
              return;
            }
	  lkl[iLkl] = tmpLkl;
	}
      else
	{
	  cout << " ## Oops! Lkl class type " << classType << " unkonwn <---------------- FATAL ERROR!!!"<< endl;
	  return;
	}
      
      // Check if the input parameters allow to export data
      if(exportData)
        {
          if(iLkl == 0 && classType.CompareTo("JointLkl")==0)
            {
              if(lkl[iLkl]->GetDUnitsOfG() > 1e-5)
                {
                  cout << " ## Warning! This format of exported data implies that LogJ is fixed, your main JointLkl object has non-nul DlogJ error, by precaution the exported file will NOT be created." << endl;
                  exportData=kFALSE;
                }
              else
                continue;
            }
          else if(iLkl >=1 && classType.CompareTo("JointLkl")==0)
            {
              cout << " ## Warning! This format of exported data implies that one file corresponds to one source, your rc file contains more than one JointLkl object which contradicts this assumption, by precaution the exported file will NOT be created." << endl;
              exportData=kFALSE;
            }
          else if(iLkl==1)
            {
	      // casted pointer (for less messy code)
	      Iact1dUnbinnedLkl* fullLkl = NULL;
	      if(!strcmp(lkl[iLkl]->ClassName(),"Iact1dUnbinnedLkl")) fullLkl = dynamic_cast<Iact1dUnbinnedLkl*>(lkl[iLkl]);
	      if(!strcmp(lkl[iLkl]->ClassName(),"Iact1dBinnedLkl"))   fullLkl = dynamic_cast<Iact1dBinnedLkl*>(lkl[iLkl]);
              logJ = fullLkl->GetLogJ(); //initializing the logJ value with the first one found

              if (lkl[iLkl]->GetDUnitsOfG() > 1e-5)
                {
                  cout << " ## Warning! This format of exported data implies that LogJ is fixed, one of your objects has non-nul DlogJ error, by precaution the exported file will NOT be created." << lkl[iLkl]->GetDUnitsOfG() << endl;
                  exportData=kFALSE;
                }
            }
          else
            {
	      // casted pointer (for less messy code)
	      Iact1dUnbinnedLkl* fullLkl = NULL;
	      if(!strcmp(lkl[iLkl]->ClassName(),"Iact1dUnbinnedLkl")) fullLkl = dynamic_cast<Iact1dUnbinnedLkl*>(lkl[iLkl]);
	      if(!strcmp(lkl[iLkl]->ClassName(),"Iact1dBinnedLkl"))   fullLkl = dynamic_cast<Iact1dBinnedLkl*>(lkl[iLkl]);

              if(TMath::Abs(fullLkl->GetLogJ()-logJ) > 1e-5)
                {
                  cout << " ## Warning! This format of exported data implies that one file corresponds to one source, your rc file contains differents LogJ values which contradicts this assumption, by precaution the exported file will NOT be created." << endl;
                  exportData=kFALSE;
                }
              else if (lkl[iLkl]->GetDUnitsOfG() > 1e-5)
                {
                  cout << " ## Warning! This format of exported data implies that LogJ is fixed, one of your objects has non-nul DlogJ error, by precaution the exported file will NOT be created." << endl;
                  exportData=kFALSE;
                }
            }
        }

      // Construct the joint Lkl tree structure: link lkl[iLkl] to the proper JointLkl object with parLkl
      if(iLkl>0 && parLkl>=iLkl) // only lkls with lower indices (already processed) are accepted
	{
	  cout << Form(" ## Oops! Parent JointLkl for lklTerm%03d cannot be lklTerm%03d, should be smaller index",iLkl,parLkl) << " <---------------- FATAL ERROR!!!"<< endl;
	  return;	      
	}
      else if(parLkl>=0 && strcmp(lkl[parLkl]->ClassName(),"JointLkl")!=0) // only lkls of type JointLkl are accepted
	{
	  cout << Form(" ## Oops! Parent JointLkl for lklTerm%03d cannot be lklTerm%03d",iLkl,parLkl) << ", because it is of type" << lkl[parLkl]->ClassName() << " instead of JointLkl <---------------- FATAL ERROR!!!"<< endl;
	  return;
	}
      else if(parLkl>=0)
	dynamic_cast<JointLkl*>(lkl[parLkl])->AddSample(lkl[iLkl]);
    } // end of loop over likelihood terms in the rc file

  if(nLkls<2)
    {
      cout << " ## Oops! there must be at least 2 lkl terms... but there are " << nLkls << " <---------------- FATAL ERROR!!!"<< endl;
      return;      
    }

  lkl[0]->SetErrorDef(deltaLogLkl); // set the error correponding to the required CL
  if(isGpositive)
    lkl[0]->SetGIsPositive();
  cout << "***" << endl;
  cout << "****************************************************************" << endl;
  cout << "*** SUMMARY OF LKL TERMS: " << endl;
  cout << "****************************************************************" << endl;

  lkl[0]->PrintData();

  // Create stream to export data
  std::ofstream data;

  // Preparing data export to Glory Duck format
  Double_t svStep = 0.;
  Double_t isv = svMax;
  Double_t svScan[svNPoints+1];
  vector<vector<Double_t> > vlkl2D(svNPoints+1);
  if (exportData)
    {
      // Create directory and open file for data export
      TString exportDataDir = fExportDataPath+simulationlabel+"/";
      gSystem->Exec(Form("mkdir -p %s",exportDataDir.Data()));
      TString seedTag  = (seed<0? "" : Form("_%d",seed));
      TString dataFile = exportDataDir+label+seedTag+".txt";
      data.open(dataFile);

      // Write first line of the file
      data << left << setw(15) << logJ;
      for(Int_t imass=0;imass<nmass;imass++)
        data << left << setw(15) << massval[imass];
      data << endl;

      // Define svStep
      Int_t counter = 0;
      if (svLogStep) svStep = TMath::Exp(TMath::Log(TMath::Abs(svMax/svMin))/svNPoints);
      else           svStep = TMath::Abs(svMax-svMin)/svNPoints;

      // Initialise svScan with all <sv> values
      for(Int_t i=0; i<=svNPoints; i++)
        {
          svScan[counter] = isv;
          counter++;
          if (svLogStep) isv/=svStep;
          else           isv-=svStep;
          vlkl2D[i] = vector<Double_t>(nmass,0.);
        }
    }

  // Loop over masses and compute the limits
  ///////////////////////////////////////////
  Double_t svLimVal[nmass];
  Double_t svSenVal[nmass];
  TCanvas** hadcanvas = new TCanvas*[nsamples];
  for(Int_t isample=0;isample<nsamples;isample++)
    hadcanvas[isample] = NULL;
  TCanvas* lklcanvas = NULL;
  Bool_t Init_canvas_samples = kFALSE;
  Bool_t Init_canvas_parabolas = kFALSE;
  cout << endl;
  cout << "***********************************************************************************" << endl;
  cout << "**** LOOPING OVER MASSES, CONFIGURE MASS-DEPENDENT HISTOS AND CALL LKL MINIMIZATION" << endl;
  for(Int_t imass=0;imass<nmass;imass++)
    {
      // Configure
      const Double_t mass = massval[imass];
      Float_t mdm  = (isDecay? mass/2. : mass);
      
      cout << "****" << endl;
      cout << "**** DM mass = " << mass << " GeV  (" << imass+1 << "/" << nmass << ")" << endl;
      cout << "**************************************" << endl;
      // compute the branching ratios for the branon model for each DM mass
      if(!channel.CompareTo("branon",TString::kIgnoreCase))
        {
          // release the memory of channelval and brval
          delete [] channelval;
          delete [] brval;
          // number of considering channels
          nChannels = 9;
          channelval = new TString[nChannels];
          brval = new Double_t[nChannels];
          // call the function that computes the branching ratios and save them in channelval and brval
          compute_branonBR(mdm,nChannels,channelval,brval);
          // Print-Out the freshly computed branching ratios for validation
          cout << "**** Branon branching ratios: ";
          for(Int_t iChannel=0;iChannel<nChannels;iChannel++)
            cout << brval[iChannel] << "*" << channelval[iChannel] << ((iChannel<nChannels-1)? "+" : "");
          cout << endl;
        }

      // Loop over samples to read the dN/dE and dN/dE' histos for signal from files
      ///////////////////////////////////////////////////////////////////////////////////////////////
      Bool_t saveHistosInFile=kFALSE;
      cout << " *** Read or compute dN/dE and dN/dE' histos for each samples:" << endl;
      for(Int_t isample=0;isample<nLkls;isample++)
        {
	  if(!strcmp(lkl[isample]->ClassName(),"Iact1dUnbinnedLkl") || !strcmp(lkl[isample]->ClassName(),"Iact1dBinnedLkl") || !strcmp(lkl[isample]->ClassName(),"LineSearchLkl"))
	    {
	      builddNdESignal(lkl[isample],fdNdEDir,nChannels,channelval,brval,mdm,isDecay,isSimulation);

	      // Read or create dN/dE' for signal
	      Iact1dUnbinnedLkl* fullLkl = NULL;
	      if(!strcmp(lkl[isample]->ClassName(),"Iact1dUnbinnedLkl")) fullLkl = dynamic_cast<Iact1dUnbinnedLkl*>(lkl[isample]);
	      if(!strcmp(lkl[isample]->ClassName(),"Iact1dBinnedLkl"))   fullLkl = dynamic_cast<Iact1dBinnedLkl*>(lkl[isample]);
	      if(!strcmp(lkl[isample]->ClassName(),"LineSearchLkl"))     fullLkl = dynamic_cast<LineSearchLkl*>(lkl[isample]);

	      TString dNdEpSignalFileName;
	      if(ioHdNdEpSignal)
		{
		  TString dNdEpSignalFileNameForm = fdNdEpSignalDir+"dNdEpSignal_"+fullLkl->GetName()+"_"+normchannel+Form("_m%s",mprecform(mdm).Data())+".root";
		  dNdEpSignalFileName = Form(dNdEpSignalFileNameForm,mdm);
		  cout << "   * Reading dN/dE' for signal from file " << dNdEpSignalFileName << "... " << flush;
		  
		  
		  // if file does not exist, create dNdEpSignal histo (it will be saved in file later)
		  if(fullLkl->ReaddNdEpSignal(dNdEpSignalFileName))
		    {
		      cout << endl << "     Not found! dNdEpSignal histo will be calculated and saved in file " << dNdEpSignalFileName << endl;
		      saveHistosInFile=kTRUE;
		    }
		  // if file exists, read histo and use it
		  else
		    cout << "Ok!" << endl;	   	      
	        }
	      
	      // Set the units of g for the different samples (this also creates the fHdNdEpSignal histo
	      if(isDecay)
	        fullLkl->SetDMDecayUnitsForG(mass);
	      else	      
	        fullLkl->SetDMAnnihilationUnitsForG(mass);

              // in case simulations were requested with DM cont, do them now
	      if(isSimulation && mcsv>0 && mcmass>0)
		{
		  cout << "  * Simulating observations for <sv> = " << mcsv << " and m_DM = " << mcmass << endl;

		  if(fullLkl->SimulateDataSamples(seed,mcsv))
		    {
		      cout << " ## Oops! Cannot simulate samples for " << fullLkl->GetName() << " <---------------- FATAL ERROR!!!"<< endl;
		      exit(1);
		    }
		}

	      // Save the dN_signal/dE' created in the SetDMAnnihilationUnitsForG call
	      if(ioHdNdEpSignal && saveHistosInFile)
	        {
		  // saving fHdNdEpSignal histo
		  gSystem->Exec(Form("mkdir -p %s",fdNdEpSignalDir.Data()));
		  
		  TFile* dNdEpSignalFile = new TFile(dNdEpSignalFileName,"RECREATE");
		  TH1F*  hdNdEpSignal    = new TH1F(*fullLkl->GetHdNdEpSignal());
		  hdNdEpSignal->SetDirectory(0);
		  cout << "     Saving dN/dE' signal histo to " << dNdEpSignalFileName << "... " << flush;
		  hdNdEpSignal->Write();
		  dNdEpSignalFile->Close();
		  delete dNdEpSignalFile;
		  delete hdNdEpSignal;
		  cout << "Ok!" << endl;
	        } 
	    }
          else if(!strcmp(lkl[isample]->ClassName(),"GloryDuckTables2019Lkl"))
            {
              GloryDuckTables2019Lkl* gdLkl = dynamic_cast<GloryDuckTables2019Lkl*>(lkl[isample]);
              if(gdLkl->SetActiveMass(mass))
                {
                  cout << " ## Oops! The DM mass selected (" << mass << " GeV) doesn't exist in the input file <---------------- FATAL ERROR!!!" << endl;
                  return;
                } 
            }
        } // end of loop over samples
      cout << " *** End of reading dN_signal/dE and dN_signal/dE' histograms" << endl;  

      // Plot IRF and data
      /////////////////////
      if(showSamplePlots)
	// loop over samples
	for(Int_t isample=0;isample<nsamples;isample++)
	  {
	    Iact1dUnbinnedLkl* fullLkl = dynamic_cast<Iact1dUnbinnedLkl*>(sample[isample]);
 
	    if(!Init_canvas_samples)
	      {
		hadcanvas[isample] = fullLkl->PlotHistosAndData();
		hadcanvas[isample]->SetName(Form("hadcanvas_%d",isample));
		hadcanvas[isample]->SetTitle(Form("IRFs and data for sample %d",isample));
		hadcanvas[isample]->cd(5);
		TLatex* ltchannel;
		if(nChannels == 1) ltchannel = new TLatex(0.8,0.8,strchannel);
		else ltchannel = new TLatex(0.7,0.8,strchannel);
		ltchannel->SetTextSize(0.055);
		ltchannel->SetNDC();
		ltchannel->Draw();
	        if (!Init_canvas_samples && isample==nsamples-1) Init_canvas_samples = kTRUE;
	      }
	    else
	      {
		hadcanvas[isample]->cd(5);
		TH1F* hdNdESignal = new TH1F(*fullLkl->GetHdNdESignal());
		hdNdESignal->SetLineStyle((imass+1)%8);
		hdNdESignal->DrawCopy("same");
		delete hdNdESignal;

		hadcanvas[isample]->cd(6);
		TH1F* hdNdEpSignal = new TH1F(*fullLkl->GetHdNdEpSignal());
		hdNdEpSignal->SetLineStyle((imass+1)%8);
		hdNdEpSignal->Scale(hdNdEpSignal->GetBinContent(0));				
		hdNdEpSignal->DrawCopy("same");
		delete hdNdEpSignal;

		if(fullLkl->GetHdNdEpSignalOff())
		  {
		    TH1F* hdNdEpSignalOff = new TH1F(*fullLkl->GetHdNdEpSignalOff());
		    hdNdEpSignalOff->SetLineStyle((imass+1)%8);
		    hdNdEpSignalOff->Scale(hdNdEpSignalOff->GetBinContent(0));				
		    hdNdEpSignalOff->DrawCopy("same");
		    delete hdNdEpSignalOff;
		  }
	      }	    
	    gPad->Modified();
	    gPad->Update();
	  } // end of loop over samples
      
      // compute -2logLkl vs g for precise limit computation
      cout << " *** Computing -2logL (parabola) vs g:" << endl;
      if(!lkl[0]->ComputeLklVsG())
        {
          cout << " *** Skipping DM mass = " << mass << " GeV because checks were not successfull (maybe none of the samples will produce any signal event?)" << endl;
          svLimVal[imass] = 0.;
          svSenVal[imass] = 0.;
          continue;
        }

      cout << endl;
      cout << " *** Overview of likelihood maximization results:" << endl;
      cout << " ************************************************" << endl;
      lkl[0]->PrintOverview();
      cout << endl;
      
      // Compute the limits (<sv> or 1/tauDM) vs DM mass
      if(isGpositive) // use Fermi criterium
	{
	  grLklParabola[imass] = lkl[0]->GetLklVsG();
	  Double_t lklat0      = grLklParabola[imass]->Eval(0);
	  Double_t svminval    = lkl[0]->GetGLklMin();
	  Double_t svcutval    = svminval+lkl[0]->GetGLklMinErr();
	  Double_t svcutvalpos = (svminval>0? svcutval : lkl[0]->GetGForLkl(lklat0+deltaLogLkl));
	  svSenVal[imass]      = svcutval-svminval; // sensitivity
	  svLimVal[imass]      = (svminval>0?  svcutval : svcutvalpos); // convention used in the Fermi paper
	}
      else // use the Segue Stereo paper criterium
	{	  
	  Double_t svminval = lkl[0]->GetGLklMin();
	  Double_t svcutval = svminval+lkl[0]->GetGLklMinErr();
	  svSenVal[imass]  = svcutval-svminval; // sensitivity
	  svLimVal[imass]  = (svminval<0?  svcutval-svminval : svcutval); // convention used in the Segue paper
	}
      if(isDecay)
	{
	  svSenVal[imass]=1./svSenVal[imass];
	  svLimVal[imass]=1./svLimVal[imass];
	}

      // Plot -2logLkl vs <sv> (parabolas) if computed and requested
      //////////////////////////////////////////////////////////////
      if(showParabolaPlots)
	{
	  // get the graph of -2logLkl vs g*unitsOfG
	  grLklParabola[imass] = lkl[0]->GetLklVsG();
	  grLklParabola[imass]->SetName(Form("grLklParabola_%02d",imass));

	  if(!Init_canvas_parabolas)
	    {
	      gStyle->SetPadRightMargin(0.1);
	      lklcanvas = new TCanvas("lklcanvas","-2logLkl vs g curves",ncols*250,nlines*250);
	      lklcanvas->Divide(ncols,nlines);
	      Init_canvas_parabolas = kTRUE;
	    }
	  lklcanvas->cd(imass+1);

	  TString parabolaplotform = Form("-2logLkl vs %s for mass %s GeV",(isDecay? "1/#tau_{DM}":"<sv>"),mprecform(mass).Data());
	  
	  // plot empty histo with nice settings to hold the -2logLkl parabolas
	  TString dummytit = Form(parabolaplotform,mass);
	  TH1I *dymmyparabola = new TH1I(Form("dummyparabola_%d",imass),dummytit,1,grLklParabola[imass]->GetX()[0],grLklParabola[imass]->GetX()[grLklParabola[imass]->GetN()-1]);
	  dymmyparabola->SetDirectory(0);
	  dymmyparabola->SetStats(0);
	  dymmyparabola->SetXTitle((isDecay?"1/#tau_{DM} [s^{-1}]" : "<#sigma v> [cm^{3}/s]"));
	  dymmyparabola->SetYTitle("#Delta(-2logL)");
	  dymmyparabola->SetMinimum(0);
	  dymmyparabola->SetMaximum(10);
	  dymmyparabola->DrawCopy();
	  delete dymmyparabola;
	  
	  // plot -2logLkl vs <sv>
	  grLklParabola[imass]->Draw("l");
	  gPad->SetGrid();
	  gPad->Modified();
	  gPad->Update();
	}

      // Save -2logLkl vs <sv> in file
      //////////////////////////////////////////////////////////////
      if (exportData)
        for (Int_t isv=0; isv<=svNPoints; isv++)
          vlkl2D[isv][imass] = grLklParabola[imass]->Eval(svScan[isv]);

    } // end of loop over DM masses

  if (exportData)
    for (Int_t isv=0; isv<=svNPoints; isv++)
      {
        data << left << setw(15) << svScan[isv];

        for(Int_t imass=0;imass<nmass;imass++)
          data << left << setw(15) << vlkl2D[isv][imass];

        // go to next line in the file
        data << endl;
      }

  // Close file for exporting data
  if (exportData) data.close();

  //################
  // FINAL RESULTS
  //################

  // Scale limits
  //////////////////
  Float_t minparval =  9e99;
  Float_t maxparval = -9e99;
  cout << "Limit/sensitivity values are scaled by a factor " << plotScale << endl;
  for(Int_t imass=0;imass<nmass;imass++)
    {
      svLimVal[imass]*=plotScale;
      svSenVal[imass]*=plotScale;
      if(svLimVal[imass]<minparval) minparval=svLimVal[imass];
      if(svSenVal[imass]<minparval) minparval=svSenVal[imass];
      if(svLimVal[imass]>maxparval) maxparval=svLimVal[imass];
      if(svSenVal[imass]>maxparval) maxparval=svSenVal[imass];
    }


  // Report limits
  //////////////////
  cout << endl << endl;
  cout << "READY-TO-COPY results: " << endl;
  cout << endl;

  cout << "Double_t mass[nmass]  = {";
  for(Int_t imass=0;imass<nmass;imass++)
    cout << massval[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;
  cout << "Double_t limit[nmass]  = {";
  for(Int_t imass=0;imass<nmass;imass++)
    cout << svLimVal[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;
  cout << "Double_t snstvt[nmass]  = {";
  for(Int_t imass=0;imass<nmass;imass++)
    cout << svSenVal[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;
 
  cout << endl;
  cout << "**********************************" << endl;
  cout << "Joint likelihood results " << endl;
  cout << "**********************************" << endl;
  cout << endl;
  cout << Form("%s limit vs mass",(isDecay? "tauDM" : "<sv>")) << endl;
  cout << "*****************************************" << endl;
  for(Int_t imass=0;imass<nmass;imass++)
    cout << "mass = " << massval[imass] << " GeV, " << Form("%s^UL = ",(isDecay? "tauDM" : "<sv>")) << svLimVal[imass]<< Form(", %s_snstvty = ",(isDecay? "tauDM" : "<sv>")) << svSenVal[imass] << (isDecay? "s-1" : " cm^3s-1") << endl;
  cout << endl;
 
    
  // Create and plot the global limit and sensitivity curves   
  //////////////////////////////////////////////////////////
  if(plotmin==0 && plotmax==0)
    {
      plotmin=minparval*0.5;
      plotmax=maxparval*2;
    }

  // graph for <sv> upper limits
  TGraph* grsvlim  = new TGraph(nmass,massval,svLimVal);
  grsvlim->SetName(isDecay? "grtaulim": "grsvlim");
  grsvlim->SetLineColor(1);
  
  // graph for <sv> sensitivity
  TGraph* grsvsen  = new TGraph(nmass,massval,svSenVal);
  grsvsen->SetName(isDecay? "grtausen" : "grsvsen");
  grsvsen->SetLineColor(1);
  if(showLimitPlots)
    grsvsen->SetLineStyle(2);
  else
    grsvsen->SetLineStyle(1);
  

  // canvas for plots
  TCanvas* limcanvas  = new TCanvas("limcanvas",Form("Dark matter %s limits",(isDecay? "tauDM" : "<sv>")),800,800);

  TH1I *dummylim = new TH1I("dummylim",Form("%s ULs vs mass",(isDecay? "#tau_{DM}" : "<#sigma v>")),1,massval[0],massval[nmass-1]);
  dummylim->SetStats(0);
  dummylim->SetMinimum(plotmin);
  dummylim->SetMaximum(plotmax);
  dummylim->SetXTitle("m_{DM} [GeV]");
  dummylim->SetYTitle(Form("95%% CL %s^{UL} [%s]",(isDecay? "#tau_{DM}" : "<#sigma v>"),(isDecay? "s" : "cm^{3}/s")));
  dummylim->DrawCopy();
  grsvsen->Draw("l");
  if(showLimitPlots)
    grsvlim->Draw("l");  

  TLatex* txchannel;
  if(nChannels == 1) txchannel = new TLatex(0.8,0.2,strchannel);
  else txchannel = new TLatex(0.6,0.2,strchannel);
  txchannel->SetTextSize(0.055);
  txchannel->SetNDC();
  txchannel->Draw();

  TLegend* limleg = new TLegend(0.6, 0.7, 0.85, 0.85);
  limleg->SetFillColor(0);
  limleg->SetMargin(0.40);
  limleg->SetBorderSize(0);
  if(showLimitPlots)
    limleg->AddEntry(grsvlim,"Limit","L");
  limleg->AddEntry(grsvsen,"Sensitivity","L");
  limleg->Draw();

  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid();


  // save plots
  TString realPlotDir = fPlotsDir+simulationlabel+"/";
  gSystem->Exec(Form("mkdir -p %s/root",realPlotDir.Data()));
  gSystem->Exec(Form("mkdir -p %s/pdf",realPlotDir.Data()));
  
  TString seedTag  = (seed<0? "" : Form("_%d",seed));
  
  limcanvas->Print(realPlotDir+"root/"+label+"_"+simulationlabel+"_limits"+seedTag+".root");
  limcanvas->Print(realPlotDir+"pdf/" +label+"_"+simulationlabel+"_limits"+seedTag+".pdf");
  if(showParabolaPlots)
    {
      lklcanvas->Print(realPlotDir+"root/"+label+"_"+simulationlabel+"_2logLVsG"+seedTag+".root");
      lklcanvas->Print(realPlotDir+"pdf/" +label+"_"+simulationlabel+"_2logLVsG"+seedTag+".pdf");
    }
  if(showSamplePlots) 
    for(Int_t isample=0;isample<nsamples;isample++)
      {	
	hadcanvas[isample]->Print(realPlotDir+"root/"+label+"_"+simulationlabel+Form("_histos_sample%02d",isample)+seedTag+".root");
	hadcanvas[isample]->Print(realPlotDir+"pdf/" +label+"_"+simulationlabel+Form("_histos_sample%02d",isample)+seedTag+".pdf");
      }
    
  // Clean up and close 
  /////////////////////
  delete [] lkl;
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
  //  gStyle->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  //  gStyle->SetOptTitle(0);
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
  gStyle->SetColorModelPS(1);
  gStyle->SetPalette(1,0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetGridWidth(1);
  gStyle->SetHistLineWidth(1);
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
void compute_branonBR(Float_t &mass, Int_t &nChannels, TString *channelval, Double_t *brval)
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
  }

void builddNdESignal(Lkl* lkl,TString dNdEDir, Int_t nChannels,TString* channelval,Double_t* brval,Float_t mdm,Bool_t isDecay,Bool_t isSimulation)
{
  // casted pointer (for less messy code)
  Iact1dUnbinnedLkl* fullLkl = NULL;
  if(!strcmp(lkl->ClassName(),"Iact1dUnbinnedLkl")) fullLkl = dynamic_cast<Iact1dUnbinnedLkl*>(lkl);
  if(!strcmp(lkl->ClassName(),"Iact1dBinnedLkl"))   fullLkl = dynamic_cast<Iact1dBinnedLkl*>(lkl);

  cout << "  ** Building histos for sample " << fullLkl->GetName() << ":" << endl;

  // Delete existing fHdNdESignal and create empty one
  fullLkl->ResetdNdESignal();

  for(Int_t iChannel=0;iChannel<nChannels;iChannel++)
    if(brval[iChannel] != 0.0)
      {
	// first check if channel associated dN/dE can be represented by a simple formula
	if(!channelval[iChannel].CompareTo("gammagamma",TString::kIgnoreCase))
	  {
	    cout << "   * Setting dN/dE for a monochromatic line at energy " << mdm  << " GeV with BR = " << brval[iChannel] << " ... " << flush;
	    if(fullLkl->AdddNdESignalFunction("line",mdm,2,brval[iChannel]))
	      {
		cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
		exit(1);
	      }
	    else
	      cout << "Ok!" << endl;
	  }
	else if(!channelval[iChannel].CompareTo("pi0pi0",TString::kIgnoreCase))
	  {
	    const Float_t mpi = 0.135; // pi0 mass in GeV
	    Float_t emin =  mdm/2.*(1-TMath::Sqrt(1-mpi*mpi/(mdm*mdm)));
	    Float_t emax =  mdm/2.*(1+TMath::Sqrt(1-mpi*mpi/(mdm*mdm)));
	    cout << "   * Setting dN/dE for XX->pi0pi0 for mass = " << mdm  << " GeV (i.e. box between " << emin << " and " << emax << " GeV) with BR = " << brval[iChannel] << "... " << flush;

	    if(fullLkl->AdddNdESignalFunction("box",emin,emax,4,brval[iChannel]))
	      {
		cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
		exit(1);
	      }
	    else
	      cout << "Ok!" << endl;
	  }
	else if(!channelval[iChannel].CompareTo("pi0gamma",TString::kIgnoreCase) || !channelval[iChannel].CompareTo("gammapi0",TString::kIgnoreCase))
	  {
	    const Float_t mpi = 0.135; // pi0 mass in GeV
	    Float_t e0   = mdm-mpi*mpi/(4*mdm);
	    Float_t emin = mdm/2.*((1+mpi*mpi/(4*mdm*mdm))-(1-mpi*mpi/(4*mdm*mdm)));
	    Float_t emax = mdm/2.*((1+mpi*mpi/(4*mdm*mdm))+(1-mpi*mpi/(4*mdm*mdm)));
	    cout << "   * Setting dN/dE for XX->pi0gamma for mass = " << mdm  << " GeV (i.e. a line at E0=" << e0<< ", plus a box between " << emin << " and " << emax << " GeV) with BR = " << brval[iChannel] << "... " << flush;

	    if(fullLkl->AdddNdESignalFunction("line",e0,1,brval[iChannel]))
	      {
		cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
		exit(1);
	      }
	    if(fullLkl->AdddNdESignalFunction("box",emin,emax,2,brval[iChannel]))
	      {
		cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
		exit(1);
	      }
	    else
	      cout << "Ok!" << endl;
	  }
	// commented out because we read it from histograms like the rest of the usual channels
	// else if(!channelval[iChannel].CompareTo("ee",TString::kIgnoreCase))
	//   {
	//     const Float_t me = 0.511e-3; // e mass in GeV

	//     TF1* fee = new TF1("fee","1./137./TMath::Pi()/x*(TMath::Log(4*[0]*([0]-x)/([1]*[1]))-1)*(1+TMath::Power(4*[0]*([0]-x)/(4*[0]*[0]),2))",1e-4,mdm);
	//     fee->SetParameter(0,mdm);
	//     fee->SetParameter(1,me);
	//     Float_t emax = mdm*(1-TMath::Exp(1)/4.*me*me/(mdm*mdm));
	//     cout << "   * Setting dN/dE for XX->ee for mass = " << mdm  << " GeV (Emax = " << emax << ") with BR = " << brval[iChannel] << "... " << flush;

	//     if(fullLkl->AdddNdESignalFunction(fee,0,emax,brval[iChannel]))
	//       {
	//          cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
	//          return;
	//       }
	//     else
	//       cout << "Ok!" << endl;
	//     delete fee;
	//   }
	else
	  {
	    TString dNdESignalFileNameForm   = dNdEDir+"dNdESignal_"+channelval[iChannel]+Form("_%smass.root",mprecform(mdm).Data());
	    const TString dNdESignalFileName = Form(dNdESignalFileNameForm,mdm);
	    cout << "   * Reading dN/dE for signal from file " << dNdESignalFileName  << " with BR = " << brval[iChannel] << "... " << flush;

	    // Try to read dNdE from existing file
	    if(fullLkl->AdddNdESignal(dNdESignalFileName,brval[iChannel]))
	      {
		cout << "Failed! <---------------- FATAL ERROR!!!" << endl;
		exit(1);
	      }
	    else
	      cout << "Ok!" << endl;
	  }
      } // end of loop over channels
}

TString mprecform(Double_t val)
{
  // compute the minimum precision the value (e.g. mass) needs to be reported with
  Int_t mprec;
  for(mprec=10;mprec>=0;mprec--)
    if(val<TMath::Power(10,-mprec))
      break;
  return Form("%%.%df",mprec+2);
}
