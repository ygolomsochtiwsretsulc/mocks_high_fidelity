# writes lines in 
# healpix_radius_nested.dat

import healpy

import numpy as n

for HEALPIX_8_id in n.arange(healpy.nside2npix(8)):
 out = healpy.pix2ang(8, HEALPIX_8_id, nest=True, lonlat=True)
 print(str(HEALPIX_8_id).zfill(3), out[0], out[1])


#HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - DEC*n.pi/180. , RA * n.pi/180. , nest=True)


