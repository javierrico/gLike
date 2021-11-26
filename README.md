![gLike logo](https://github.com/javierrico/gLike/raw/master/logo/gLike_logo_small.png "gLike logo")

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4601451.svg)](https://doi.org/10.5281/zenodo.4601451)

Numerical maximization of heterogeneous joint likelihood functions of a common free parameter plus nuisance parameters
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

## Citing gLike
If you are using `gLike` for your work you can acknowledge and reference it through its [zenodo record](https://zenodo.org/record/3967386).

## Prerequisites
`gLike`'s sole dependency is [ROOT](https://root.cern.ch). 

### Note on IACT analyses with gLike 
To input [IACT data in FITS format](https://github.com/javierrico/gLike/wiki/Data-format-for-analysis-of-IACT-data#fits-format) `gLike` uses special `TFITSHDU` methods available only in very recent ROOT versions (> 6.20.04).     
If you build `gLike` using an older ROOT version, the interface to FITS data for IACT analysis will not be available.

## Installation
1. Get the code from the GitHub [gLike repository](https://github.com/javierrico/gLike)

2. At user level you are recommended to download, compile and run the latest stable release. Once you become an expert you will want to develop your own classes, for which you will need to check out the master (for new developments) or release (for bug fixes) branches. Check the [release wiki entry](https://github.com/javierrico/gLike/wiki/Branch-releases-log) for more information about the gLike repository branch and release structure.

3. Since release 0.10, the sole installation method is using `cmake`. For  `Makefile`-based installation please refer to release 0.9 or earlier 

4. `cmake` installation
    - Define the environment variable `GLIKESYS` pointing to the directory _where you want to build_ gLike,    
      you can create the build directory within the original `gLike` repository itself    
      `cd gLike`   
      `mkdir build`     
      or outside in another repository    
      `mkdir ~/gLike_build`.    
      Just export this build directory as an environment variable in your `.bashrc`:     
      `export GLIKESYS="/home/user/gLike/build"`    
      or    
      `export GLIKESYS="/home/user/gLike_build"`   
      if the build directory is outside the one containing the source code.   

    - go in the build directory      
      `cd $GLIKESYS`   
    - run `cmake`, the first argument being the directory in which the `CMakeLists.txt` for building the project is defined 
      (the `gLike` directory we downloaded from git)     
      `cmake /path/to/the/source/gLike/dir/`    
      if you want to **unlock FITSIO support**, use the cache variable `USE_FITSIO` (can be set to `ON` or `True`, default is `OFF`)    
      `cmake /path/to/the/source/gLike/dir/ -DUSE_FITSIO=ON`  
      if you want to **generate the html documentation**, use the cache variable `MAKE_DOCS` (can be set to `ON` or `True`, default is `OFF`)    
      `cmake /path/to/the/source/gLike/dir/ -DMAKE_DOCS=ON`    
    - at this point `cmake` will have automatically generated the `Makefile` for us, so we just    
      `make`         
      and install    
      `make install`    
      this last step will make gLike libraries available system-wise, i.e. libraries will be saved in `/usr/local/lib` and executables in `/usr/local/bin`. If you do not have root privileges and cannot write into `/usr`, you can specify the path where you want the gLike libraries to be compiled and installed. Alternatives might be `$HOME/.local` or `$HOME/opt`. Pass the their path via the `CMAKE_INSTALL_PREFIX` option, e.g.       
      `cmake /path/to/the/source/gLike/dir/ -DCMAKE_INSTALL_PREFIX=$HOME/opt # install cmake without root privileges`      
      **NOTE** if you do not have root privileges and do not have installed `gLike` in `/usr/local`, but in `$HOME/.local` or `$HOME/opt` remember to add these paths to your `$PATH` to be able to use the executables without giving their full path.

### load gLike libraries in ROOT
If you want to load the gLike libraries each time you open ROOT, modify (or create) in your home a `.rootrc` file, inserting (editing) the line
```
Rint.Logon: $GLIKESYS/scripts/rootlogon.C
```
this should execute the `rootlogon.C` of gLike each time you launch root. 
An example `.rootrc` is available in `gLike` distribution repository.


## gLike distribution structure
In gLike you find the following directories:
1. [`src`](https://github.com/javierrico/gLike/tree/master/src): source files (*.cc) with definition of every class,
2. [`include`](https://github.com/javierrico/gLike/tree/master/include): include files (*.h) with declaration of every class,
3. [`exec`](https://github.com/javierrico/gLike/tree/master/exec): containing the source files (*.cc) for the executables generated by `cmake`,
4. [`scripts`](https://github.com/javierrico/gLike/tree/master/scripts): root macros and scripts with some gLike example applications,
5. [`rcfiles`](https://github.com/javierrico/gLike/tree/master/rcfiles): examples of rcfiles (in principle the only thing a regular user should edit and modify),
6. [`data`](https://github.com/javierrico/gLike/tree/master/data): examples of input data files, e.g. events and corresponding IRFs mimicking a generic IACT telescope of the 2nd generation like MAGIC,
7. [`DM`](https://github.com/javierrico/gLike/tree/master/DM): files for DM-related analysis, e.g. the dN/dE functions for different masses and annihilation channels,
8. [`logo`](https://github.com/javierrico/gLike/tree/master/logo): the gLike logo

#### Build directory structure
When compiling with `cmake` the following structure will be created in the build directory:
```shell
gLike_build
├── CMakeCache.txt
├── CMakeFiles
├── DM
├── Makefile
├── bin
├── cmake_install.cmake
├── data
├── gLikeDict.cxx
├── htmldoc
├── include
├── lib
├── rcfiles
└── scripts
```
directories of the same name of the ones found in the original repository are copied for
building purposes.     
The two directories created by `cmake` are
```shell
gLike_build/lib
├── libgLike.dylib
├── libgLikeDict.rootmap
└── libgLikeDict_rdict.pcm
```
containing the compiled library and the `ROOT` dictionary files;
```
gLike_build/bin
├── computeCLBands
└── jointLklDM
````
containing the `gLike` executables.     
`htmldoc`, containing the `html` documentation, will be generated only if the 
corresponding cache variable has been specified at building. 

## Fast description of the code
gLike is a general-purpose collection of root-based classes for maximum likelihood analysis. gLike provides a framework for producing an arbitrarily complicated joint likelihood as the product of likelihood functions of any kind. The structure is modular, so if the likelihood function you need to use is not included yet, you can program it and make a pull request. The likelihood function has one free parameter (_g_) and as many nuisance parameters as needed.
 
gLike consists of two main basic classes, [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)
      and [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h):
-  [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) contains all the machinery related to finding the minimum of the -2logL function and scanning it in the relevant range to be able to compute the confidence intervals with the desired confidence level. [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h) is an abstract class and almost all the rest of gLike classes inherit from it. It does not implement any particular likelihood function, which is precisely what the daughter classes do.
- [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) inherits from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), it holds a list of [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h)-based objects and implements a particular likelihood function that is simply the product of all the likelihood functions in the list. Because of its inheritance from [`Lkl`](https://github.com/javierrico/gLike/blob/master/include/Lkl.h), [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h)can include other [`JointLkl`](https://github.com/javierrico/gLike/blob/master/include/JointLkl.h) terms as part of the list. This is where one of the main strengths of gLike resides, because it allows to build joint likelihood functions of any level of complexity. 
 
 Almost all the other gLike classes ([`Iact1dUnbinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dUnbinnedLkl.h),[`Iact1dBinnedLkl`](https://github.com/javierrico/gLike/blob/master/include/Iact1dBinnedLkl.h), [`FermiTables2016Lkl`](https://github.com/javierrico/gLike/blob/master/include/FermiTables2016Lkl.h),  [`Parabola`](https://github.com/javierrico/gLike/blob/master/include/Parabola.h), [`Poisson`](https://github.com/javierrico/gLike/blob/master/include/PoissonLkl.h)) just implement a particular likelihood function. Adding other likelihood functions can be easily done starting from the basic skeleton provided at [`TemplateLkl`](https://github.com/javierrico/gLike/blob/master/include/TemplateLkl.h) and following one of the previously listed classes as example.

## Basic Usage for dark matter searches
Within the gLike library, there is no assumption about the physical meaning of the free parameter _g_. That needs to be assigned externally in the macros or executables using the gLike library.
The gLike distribution provides, at the directory [`scripts`](https://github.com/javierrico/gLike/tree/master/scripts), a few tutorial macros/executables for the most common use cases (e.g. minimisation of a Poisson likelihood, definition of a joint likelihood).    
They can be executed from the root terminal
```
root[0] .X tutorialPoissonLkl.C
```
or can be loaded and then executed
```
root[0] .L tutorialPoissonLkl.C
root[1] tutorialPoissonLkl()
```

After building the project the executable [`jointLklDM`](https://github.com/javierrico/gLike/blob/master/scripts/jointLklDM.C) can be used to search for dark matter, set limits in case of no detection and combine results from different targets and or instruments. 
The macros defining the executables do not need to be edited. 

You can find some documentation in its [wiki entry](https://github.com/javierrico/gLike/wiki/jointLklDM.C).
Example configuration files can be found in the [`rcfiles`](https://github.com/javierrico/gLike/tree/master/rcfiles) directory. 
There is a detailed documentation on how to configure it correctly both in the [wiki](https://github.com/javierrico/gLike/wiki/rc-files-for-jointLklDM.C), and the jointLklDM.rc rcfile provided with the distribution. 

## Transition from old cvs mdm distribution to the git gLike standalone distribution

If you are a pre-gitHub gLike user, and you whish to keep using new gLike distribututions with some old input files, you may need to adapt them following the advice of this [wiki entry](https://github.com/javierrico/gLike/wiki/Transition-from-old-cvs-mdm-distribution-to-the-git-gLike-standalone-distribution).
