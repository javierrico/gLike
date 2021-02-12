# Recipes for building gLike singularity containers
The following folder contains recipes for building singularity containers
with gLike. We plan to provide support for two different ROOT versions.

- root v5-34-36: for MARS users
- root v6-22-06: for users that want to use the CFITSIO extension to read gammapy data.

At the moment only the first is available.
Note 1: some of the tutorial scripts cannot properly be compiled in the singularity environment.
Note 2: the executables (`jointLlkDM`) are not properly installed globally in the container.
 
