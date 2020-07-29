# mocks_high_fidelity
Creation of high fidelity eROSITA mock catalogs (AGN, clusters and galaxies). Pipeline based on N-body simulations.

## Dependencies

### Code 

Dependencies:
 * python3: time, os, sys, numpy, scipy, astropy, h5py, extinction, dustmaps, matplotlib, astropy_healpix, healpy, sklearn
 * stilts/topcat
 * gawk 
 
directory `DM_LC' 
 * from halo lists, creates the light cone shells over the full sky 
 
directory `python' contains a pipeline split in four main steps: 
 * 001: coordinates
 * 002: galaxy properties
 * 003: AGN properties
 * 004: cluster properties

Follow readme files to understand the order of execution. 
 
### Files 

erosita flux limit map available here: http://www.mpe.mpg.de/~comparat/eROSITA_AGN_mock/catalogs/erosita/flux_limits.fits
to be placed or linked in data/erosita/flux_limits.fits

### N-body simulations

MultiDark :
 - https://www.cosmosim.org/cms/simulations/smdpl/
 - https://www.cosmosim.org/cms/simulations/mdpl2/

UNIT :
 - http://www.unitsims.org/

### Mock catalogs

The generated mock catalogs are available here :

http://www.mpe.mpg.de/~comparat/MultiDark/MDPL2/

