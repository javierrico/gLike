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
1. `src`: source files (*.cc) with definition of every class
2. `include`: include files (*.h) with declaration of every class 
3. `scripts`: root macros and scripts with some gLike example applications
4. `rcfiles`: examples of rcfiles (in principle the only thing a regular user should edit and modify)
5. `log`: the gLike logo

### Fast description of the code
 gLike is a general-purpose collection of root-based classes for maximum likelihood analysis. gLike provides a framework for producing an arbitrarily complicated joint likelihood as the product of likelihood functions of any kind. The structure is modular, so if the likelihood function you need to use is not included yet, you can program it and make a pull request. The likelihood function has one free parameter (_g_) and as many nuisance parameters as needed.
 
gLike consists of two main basic classes, [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)
      and [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h):
-  [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) contains all the machinery related to finding the minimum of the -2logL function and scanning it in the relevant range to be able to compute the confidence intervals with the desired confidence level. [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) is an abstract class and almost all the rest of gLike classes inherit from it. It does not implement any particular likelihood function, which is precisely what the daughter classes do.
- [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) inherits from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), it holds a list of [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)-based objects and implements a particular likelihood function that is simply the product of the particular likelihood functions in the list. Because of its inheritance from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h)can include other [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) terms as part of the list. This is where one of the main strengths of gLike resides, because it allows to build joint likelihood functions of any level of complexity. 
 
 Almost all the other gLike classes ([`Iact1dUnbinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dUnbinnedLkl.h),[`Iact1dBinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dBinnedLkl.h), [`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/FermiTables2016Lkl.h),  [`Parabola`](https://github.com/javierrico/gLike/blob/master/include/Parabola.h), [`Poisson`](https://github.com/javierrico/gLike/blob/master/include/PoissonLkl.h)) just implement a particular likelihood function. Adding other likelihood functions can be easily done starting from the basic skeleton provided at [`TemplateLkl`](https://github.com/javierrico/gLike/blob/master/include/TemplateLkl.h) and following one of the previously listed classes as example.

### Basic Usage for dark matter searches
Within the gLike library, there is not any particular assumption about the physical meaning of the free parameter g. That needs to be assigned externally in the macros or executables using the gLike library.
The gLike distribution provides a few example macros/executables for the most common use cases. In particular `scripts/jointLklDM.C` can be used to search for dark matter, set limits in case of no detection and combine results from different targets and or instruments. 

The macro itself does not need to be edited. You can find some documentation in the [wiki](https://github.com/javierrico/gLike/wiki/jointLklDM.C). It is important to run it in compiled mode, i.e. with the ROOT command:
 `.x jointLklDM.C+(<rcfilename>)`

Example configuration files can be found in the `rcfiles` directory. There is a detailed documentation on how to configure it correctly both in the [wiki](https://github.com/javierrico/gLike/wiki/rc-files-for-jointLklDM.C), and the jointLklDM.rc rcfile provided with the distribution. 


### Transition from old cvs mdm distribution to the git gLike standalone distribution

In order to prepare gLike to grow, several actions have been recently taken:
- It has been removed from the dark matter oriented library mdm and to the gitHub repository
- More importantly, the naming convention of the different classes have changed. The following table provides the correspondence between old and new class names

<center>

 |Class old name| Class new name|
  |----:|-----:|
  |`MLkl`  | [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) |
  |`MJointLkl`  | [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) |
  |`MFullLkl`  | [`Iact1dUnbinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dUnbinnedLkl.h) |
  |`MBinnedFullLkl`  | [`Iact1dBinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dBinnedLkl.h) |
  |`MBinnedFluxLkl`  | [`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/FermiTables2016Lkl.h) |
  |`MPoissonLkl`  | [`PoissonLkl`](https://github.com/javierrico/gLike/blob/master/include/PoissonLkl.h) |
  |`MParabolaLkl`  | [`Parabola`](https://github.com/javierrico/gLike/blob/master/include/Parabola.h) |
  |[`MIACTEventListIRF`](https://github.com/javierrico/gLike/blob/master/include/MIACTEventListIRF.h)|[`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/IactEventListIrf.h)| 

</center>

- These changes should be transparent to the user except for two aspects:
	- The class names are used in the rcfile, and the new naming convention should be used from now on. This means that if you used to run gLike in the old mdm times, you need to change your rcfiles according to the previous table.
	- Data+IRF input files produced in the format defined in [`MIACTEventListIRF`](https://github.com/javierrico/gLike/blob/master/include/MIACTEventListIRF.h)will not work anymore. The new files should store date in the format defined in [`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/IactEventListIrf.h). At the moment both formats are identical, and in order to facilitate the conversion, the macro `convertSegueIntoIACTFormat.C` is provided. For the same reason, the old class [`MIACTEventListIRF`](https://github.com/javierrico/gLike/blob/master/include/MIACTEventListIRF.h) has not been yet removed. Both the conversion macro and the class will be removed from the gLike distribution in future releases.
