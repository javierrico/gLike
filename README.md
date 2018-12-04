![gLike logo](https://github.com/javierrico/gLike/raw/master/logo/gLike_logo_small.png "gLike logo")
------------------------------------------------------------------------------------------------

gLike is a general-purpose [ROOT](root.cern.ch)-based code framework for the numerical maximization of joint likelihood functions.

The joint likelihood function has **one** *free* parameter and as many *nuisance* parameters as wanted, which will be profiled in the  maximization process. 

A (non-exhaustive) list of examples where gLike is useful are (in order of increasingly complexity):

 - Estimating the number of signal events (with uncertainties) in a dataset whose background content is in turn estimated from an independent measurement in a signal-free control-region.
 - Same as before, but considering in addition a systematic uncertainty in the estimation of the background. 
 - Estimating the intensity of a steady source of signal particles in the presence of background particles, merging datasets obtained under different experimental conditions.
 - Same as before, but each dataset actually comes from a different instrument and in different data format.
 - Estimating the dark matter annihilation cross-section combining observations of dwarf spheroidal galaxies by different ground-based gamma-ray telescopes, satellite gamma-ray detectors, neutrino telescopes, ....
 - Estimating the energy scale of quantum gravity by combining observations of fast gamma-ray flares observed by different ground-based gamma-ray telescopes.
 - ...

### INSTALLATION:

Linux & MAC OS:
Go to your main gLike directory and type 'make'.


### RUNNING gLike:
An example to load the gLike environment is shown in the rootlogon.C file under ./scripts. Then you can execute a script such as 'jointLklDM.C' that will read a config file from ./rcfiles such as 'jointLklDM.rc'.


### CAUTION:
On MAC OS, one might have to add the path to the gLike library (/path/to/gLike/lib) to 'DYLD_LIBRARY_PATH' to succesfully run scripts using gLike.


### DOCUMENTATION:
Go to your main gLike directory and type 'make doc'. The documentation will be generated under ./htmldoc. You can navigate in it openning the file './htmldoc/index.html'.
