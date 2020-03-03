"""
generic functions to create the cluster mock

"""
from astropy.cosmology import FlatLambdaCDM
import astropy.io.fits as fits
import numpy as n
import os

# Create a circular image
# normalize the flux to the catalogs
# truncate to 2xr500c


def write_img(matrix, out='spherical_cc', n_pixel = 120, angularSize_per_pixel = 0.033):
    prihdr = fits.Header()
    prihdr['HDUCLASS'] = 'HEASARC/SIMPUT'
    prihdr['HDUCLAS1'] = 'IMAGE'
    prihdr['HDUVERS'] = '1.1.0'
    prihdr['EXTNAME'] = 'IMAGE'
    prihdr['CTYPE1'] = ('RA---TAN', 'first axis (column) is Right Ascension')
    prihdr['CRPIX1'] = ((n_pixel+1) / 2., 'middle pixel of array in col direction')
    prihdr['CRVAL1'] = (0, 'Dec of this middle pixel, in degrees')
    prihdr['CDELT1'] = (-angularSize_per_pixel/60., 'move 1column forward,decrease RA by CDELT1/deg')
    prihdr['CROTA1'] = 0
    prihdr['CUNIT1'] = 'deg'
    prihdr['CTYPE2'] = ('DEC--TAN', 'first axis (column) is Declination')
    prihdr['CRPIX2'] = ((n_pixel+1) / 2., 'middle pixel of array in row direction')
    prihdr['CRVAL2'] = (0, 'RA of this middle pixel, in degrees')
    prihdr['CDELT2'] = (angularSize_per_pixel/60., 'move 1column forward,increase Dec by CDELT1/deg')
    prihdr['CROTA2'] = 0
    prihdr['CUNIT2'] = 'deg'
    prihdr['EQUINOX'] = 2000
    prihdr['RADECSYS'] = 'FK5'
    prihdr['aMinpPix'] = angularSize_per_pixel
    prihdu = fits.PrimaryHDU(matrix, header=prihdr)
    #if os.path.isfile(out):
        #os.remove(out)
    prihdu.writeto(out, overwrite=True)


# b_a=0.71
def create_matrix(profile, n_pixel = 121, b_a = 0.7, truncation_radius = 12.):
	x_max = truncation_radius
	#sel = (profile.y < profile.y.max()/20)
	#x_max = n.min([n.min(profile.x[sel]), profile.x[-2]])  
	angularSize_per_pixel = x_max/(n_pixel/2.)
	matrix = n.zeros((n_pixel, n_pixel))
	xxx = (n.arange(n_pixel) - (n_pixel-1) / 2.) * angularSize_per_pixel
	x_matrix, y_matrix = n.meshgrid(xxx, xxx)
	r_matrix = ((x_matrix / b_a)**2 + y_matrix**2)**0.5
	matrix = profile(r_matrix)
	matrix = matrix / n.sum(matrix)
	return matrix, angularSize_per_pixel
