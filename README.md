# g$\mathcal{L}$ike 

gLike is a general-purpose [ROOT](root.cern.ch)-based code framework for the numerical maximization of joint likelihood functions.

The joint likelihood function has **one** *free* parameter and as many *nuisance* parameters as wanted, which will be profiled in the  maximization process. 

A (non-exhaustive) list of examples where g$\mathcal{L}$ike is useful are (in order of increasingly complexity):

 - Estimating the number of signal events (with uncertainties) in a dataset whose background content is in turn estimated from an independent measurement in a signal-free control-region.
 - Same as before, but considering a systematic uncertainty in the estimation of the background. 
 - Estimating the intensity of a steady source of signal particles in the presence of background particles, merging datasets obtained under different experimental conditions.
 - Same as before, but each dataset actually comes from a different instrument and in different data format.
 - Estimating the dark matter annihilation cross-section combining observations of dwarf spheroidal galaxies by different ground-based gamma-ray telescopes, satellite gamma-ray detectors, neutrino telescopes, ....
 - Estimating the energy scale of quantum gravity combining observations of fast gamma-ray flares by different ground-based gamma-ray telescope
