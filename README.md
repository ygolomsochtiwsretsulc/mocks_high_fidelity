# mocks_high_fidelity
Creation of high fidelity eROSITA mock catalogs (AGN, clusters and galaxies). Pipeline based on N-body simulations.

## Dependencies

### Code 

python3: time, os, sys, numpy, scipy, astropy, h5py, extinction, dustmaps, matplotlib, astropy_healpix, healpy, sklearn

stilts/topcat

### Files 

sdss depth maps available here: http://risa.stanford.edu/redmapper/
to be placed or linked in 
data/sdss_depth/
sdss_dr8_nodered_nside2048_g_model_10sigma.fits.gz  
sdss_dr8_nodered_nside2048_i_model_10sigma.fits.gz  
sdss_dr8_nodered_nside2048_r_model_10sigma.fits.gz  
sdss_dr8_nodered_nside2048_z_model_10sigma.fits.gz

erosita flux limit map available here: http://www.mpe.mpg.de/~comparat/eROSITA_AGN_mock/catalogs/erosita/flux_limits.fits
to be placed or linked in data/erosita/flux_limits.fits

### N-body simulations

MultiDark :
 - https://www.cosmosim.org/cms/simulations/smdpl/
 - https://www.cosmosim.org/cms/simulations/mdpl2/

UNIT :
 - http://www.unitsims.org/
