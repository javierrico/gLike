![gLike logo](https://github.com/javierrico/gLike/raw/master/logo/gLike_logo_small.png "gLike logo")
------------------------------------------------------------------------------------------------

gLike is a general-purpose [ROOT](root.cern.ch)-based code framework for the numerical maximization of joint likelihood functions.

The joint likelihood function has **one** *free* parameter (named _g_) and as many *nuisance* parameters as wanted, which will be profiled in the  maximization process. 

Follows a non-exhaustive list of examples where gLike is useful (in order of increasingly complexity):

 - Estimating the number of signal events (with uncertainties) in a dataset whose background content is in turn estimated from an independent measurement in a signal-free control-region.
 - Same as before, but considering in addition a systematic uncertainty in the estimation of the background. 
 - Estimating the intensity of a steady source of signal particles in the presence of background particles from datasets obtained under different experimental conditions.
 - Same as before, but each dataset actually comes from a different instrument and in different data format.
 - Estimating the dark matter annihilation cross-section combining observations of dwarf spheroidal galaxies by different ground-based gamma-ray telescopes, satellite gamma-ray detectors, neutrino telescopes, ....
 - Estimating the energy scale of quantum gravity by combining observations of fast gamma-ray flares observed by different ground-based gamma-ray telescopes.
 - ...
 
### Prerequisites
`gLike`'s sole dependency is [ROOT](https://root.cern.ch). 

#### Note on IACT analyses with gLike 
To input [IACT data in FITS format](https://github.com/javierrico/gLike/wiki/Data-format-for-analysis-of-IACT-data#fits-format) `gLike` uses sepcial `TFITSHDU` methods available only in very recent ROOT versions (> 6.20.04).     
If you build `gLike` using an older ROOT version, the interface to FITS data for IACT analysis will not be available.

### Installation
1. Get the code from the GitHub [gLike repository](https://github.com/javierrico/gLike)
2. At user level you are recommended to download, compile and run the latest stable release. Once you become an expert you will want to develop your own classes, for which you will need to check out the master (for new developments) or release (for bug fixes) branches. Check the [release wiki entry](https://github.com/javierrico/gLike/wiki/Branch-releases-log) for more information about the gLike repository branch and release structure.
3. Define the environment variable GLIKESYS pointing to your main gLike directory, e.g., for bash:

    `export GLIKESYS="/Users/rico/gLike/"`

    and include it in the library path, eg. for Mac:

      `export DYLD_LIBRARY_PATH="${GLIKESYS}:${DYLD_LIBRARY_PATH}"`
4. Create the gLike library by typing _make_ in the main gLike directory
5. [Optional] Create the gLike online documentation with _make doc_

### gLike distribution structure
In gLike you find the following directories:
1. [`src`](https://github.com/javierrico/gLike/tree/master/src): source files (*.cc) with definition of every class,
2. [`include`](https://github.com/javierrico/gLike/tree/master/include): include files (*.h) with declaration of every class,
3. [`scripts`](https://github.com/javierrico/gLike/tree/master/scripts): root macros and scripts with some gLike example applications,
4. [`rcfiles`](https://github.com/javierrico/gLike/tree/master/rcfiles): examples of rcfiles (in principle the only thing a regular user should edit and modify),
5. [`data`](https://github.com/javierrico/gLike/tree/master/data): examples of input data files, e.g. events and corresponding IRFs mimicking a generic IACT telescope of the 2nd generation like MAGIC,
6. [`DM`](https://github.com/javierrico/gLike/tree/master/DM): files for DM-related analysis, e.g. the dN/dE functions for different masses and annihilation channels,
7. [`logo`](https://github.com/javierrico/gLike/tree/master/logo): the gLike logo

### Fast description of the code
 gLike is a general-purpose collection of root-based classes for maximum likelihood analysis. gLike provides a framework for producing an arbitrarily complicated joint likelihood as the product of likelihood functions of any kind. The structure is modular, so if the likelihood function you need to use is not included yet, you can program it and make a pull request. The likelihood function has one free parameter (_g_) and as many nuisance parameters as needed.
 
gLike consists of two main basic classes, [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)
      and [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h):
-  [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) contains all the machinery related to finding the minimum of the -2logL function and scanning it in the relevant range to be able to compute the confidence intervals with the desired confidence level. [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) is an abstract class and almost all the rest of gLike classes inherit from it. It does not implement any particular likelihood function, which is precisely what the daughter classes do.
- [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) inherits from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), it holds a list of [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)-based objects and implements a particular likelihood function that is simply the product of all the likelihood functions in the list. Because of its inheritance from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h)can include other [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) terms as part of the list. This is where one of the main strengths of gLike resides, because it allows to build joint likelihood functions of any level of complexity. 
 
 Almost all the other gLike classes ([`Iact1dUnbinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dUnbinnedLkl.h),[`Iact1dBinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dBinnedLkl.h), [`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/FermiTables2016Lkl.h),  [`Parabola`](https://github.com/javierrico/gLike/blob/master/include/Parabola.h), [`Poisson`](https://github.com/javierrico/gLike/blob/master/include/PoissonLkl.h)) just implement a particular likelihood function. Adding other likelihood functions can be easily done starting from the basic skeleton provided at [`TemplateLkl`](https://github.com/javierrico/gLike/blob/master/include/TemplateLkl.h) and following one of the previously listed classes as example.

### Basic Usage for dark matter searches
Within the gLike library, there is no assumption about the physical meaning of the free parameter _g_. That needs to be assigned externally in the macros or executables using the gLike library.
The gLike distribution provides, at the directory [`scripts`](https://github.com/javierrico/gLike/tree/master/scripts), a few example macros/executables for the most common use cases. In particular [`jointLklDM.C`](https://github.com/javierrico/gLike/blob/master/scripts/jointLklDM.C) can be used to search for dark matter, set limits in case of no detection and combine results from different targets and or instruments. 

The macro itself does not need to be edited. You can find some documentation in its [wiki entry](https://github.com/javierrico/gLike/wiki/jointLklDM.C). It is important to run it in compiled mode, i.e. with the ROOT command:
 `.x jointLklDM.C+(<rcfilename>)`

Example configuration files can be found in the [`rcfiles`](https://github.com/javierrico/gLike/tree/master/rcfiles) directory. There is a detailed documentation on how to configure it correctly both in the [wiki](https://github.com/javierrico/gLike/wiki/rc-files-for-jointLklDM.C), and the jointLklDM.rc rcfile provided with the distribution. 


### Transition from old cvs mdm distribution to the git gLike standalone distribution

If you are a pre-gitHub gLike user, and you whish to keep using new gLike distribututions with some old input files, you may need to adapt them following the advice of this [wiki entry](https://github.com/javierrico/gLike/wiki/Transition-from-old-cvs-mdm-distribution-to-the-git-gLike-standalone-distribution).
