"""
from astropy_healpix import healpy
import numpy as n

N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels):
	print("nohup python 003_6_agn_magnitudes.py MD10 all "+str(HEALPIX_id).zfill(3)+" > 003_6_log/log.MD10.all."+str(HEALPIX_id).zfill(3)+" & ")
	print("nohup python 003_6_agn_magnitudes.py MD10 sat "+str(HEALPIX_id).zfill(3)+" > 003_6_log/log.MD10.sat."+str(HEALPIX_id).zfill(3)+" & ")
	print("nohup python 003_6_agn_magnitudes.py MD04 all "+str(HEALPIX_id).zfill(3)+" > 003_6_log/log.MD04.all."+str(HEALPIX_id).zfill(3)+" & ")
	print("nohup python 003_6_agn_magnitudes.py MD04 sat "+str(HEALPIX_id).zfill(3)+" > 003_6_log/log.MD04.sat."+str(HEALPIX_id).zfill(3)+" & ")

"""

import matplotlib.pyplot as p
from scipy.stats import norm
from scipy.interpolate import interp1d
from lib_magnitudes import *
import time
t0 = time.time()


#import matplotlib
# matplotlib.use('Agg')

env = sys.argv[1]
ftyp = sys.argv[2]
HEALPIX_ID = sys.argv[3]

catalog_input = os.path.join(
    os.environ[env],
    'cat_AGN_' + ftyp,
    HEALPIX_ID.zfill(6) + '.fit')
# catalog_inputs.sort()
# for catalog_input in catalog_inputs:
print(catalog_input)
agn_hdus_i = fits.open(catalog_input)
agn_hdus = Table(agn_hdus_i[1].data)
agn_hdus_i.close()

# retrieve the templates
all_tpl = sorted(
    n.array(
        glob.glob(
            os.path.join(
                template_cigale_dir,
                '*.fits'))))
# 2117 has the highest AGN fraction, z=1.749. Optical type 1
# 1460 is in the middle, z=0.257. Optical type 2
# 913 has the lowest AGN fraction, z=0.174. Optical elliptical

tpl_t1 = os.path.join(template_cigale_dir, '2117_best_model.fits')
tpl_t2 = os.path.join(template_cigale_dir, '1460_best_model.fits')
tpl_t3 = os.path.join(template_cigale_dir, '913_best_model.fits')
# then each template individually
flambda_t1, ll_t1, hdu_t1 = get_cigale_template(tpl_t1, redshift=1.749)
flambda_t2, ll_t2, hdu_t2 = get_cigale_template(tpl_t2, redshift=0.257)
flambda_t3, ll_t3, hdu_t3 = get_cigale_template(tpl_t3, redshift=0.174)


def get_rescaling_values(flambda, ll, r_mag, redshift):
    # interpolate the r magnitude with redshift
    zs = n.arange(0., 6.2, 0.1)
    ABmag_sdss = n.array([sdss_r_filter.get_ab_magnitudes(
        flambda, ll * (1 + z_i))['sdss2010-r'][0] for z_i in zs])
    ABmag_sdss_z = interp1d(zs, ABmag_sdss)
    # function to rescale to the right magnitude for a given redshift
    def rescale_by(r_mag_out, redshift): return 10**((r_mag_out +
                                                      48.6) / -2.5) / 10**((ABmag_sdss_z(redshift) + 48.6) / -2.5)
    # rescaling values
    rsbs = rescale_by(r_mag, redshift)
    return rsbs


r_mag = agn_hdus['AGN_SDSS_r_magnitude']
redshift = agn_hdus['redshift_R']
AGN_random_number = agn_hdus['AGN_random_number']
AGN_type = agn_hdus['AGN_type']

type_1 = (AGN_type == 11) | (AGN_type == 12)
type_2 = (AGN_type == 22) | (AGN_type == 21)
type_3 = (type_2) & (AGN_random_number < 0.2)

rsbs = n.zeros_like(redshift)
rsbs[type_1] = get_rescaling_values(
    flambda_t1, ll_t1, r_mag[type_1], redshift[type_1])
rsbs[type_2] = get_rescaling_values(
    flambda_t2, ll_t2, r_mag[type_2], redshift[type_2])
rsbs[type_3] = get_rescaling_values(
    flambda_t3, ll_t3, r_mag[type_3], redshift[type_3])

# computes the magnitudes 1 by 1
ABmag_sdss_out = all_filters.get_ab_magnitudes(
    flambda_t1 * rsbs[0], ll_t1 * (1 + redshift[0]))

for rsb, zz, t1, t2, t3 in zip(
        rsbs[1:], redshift[1:], type_1[1:], type_2[1:], type_3[1:]):
    if t1:
        ABmag_sdss_out.add_row(all_filters.get_ab_magnitudes(
            flambda_t1 * rsb, ll_t1 * (1 + zz))[0])
    elif t2:
        ABmag_sdss_out.add_row(all_filters.get_ab_magnitudes(
            flambda_t2 * rsb, ll_t2 * (1 + zz))[0])
    else:
        ABmag_sdss_out.add_row(all_filters.get_ab_magnitudes(
            flambda_t3 * rsb, ll_t3 * (1 + zz))[0])


def assign_mag(mag_name=ABmag_sdss_out.colnames[0]):
    # assign magnitude uncertainties
    # random variable in a norm distribution with sigma^2 = 1/80
    # scatter around the relation
    rds = norm.rvs(loc=0, scale=40**(-0.5), size=len(ABmag_sdss_out[mag_name]))
    # parameters of the mag - meg err relation
    pars = log_error_scaling_parameters[mag_name]
    # mean log error
    log_error = log_error_scaling(pars, ABmag_sdss_out[mag_name])
    # realization of the errors and magnitude
    error = 10**(log_error) + rds
    mag = ABmag_sdss_out[mag_name] + error
    mag_err = 10**(log_error)
    return mag, mag_err


for mag_name in ABmag_sdss_out.colnames:
    mag, err = assign_mag(mag_name)
    ABmag_sdss_out.replace_column(name=mag_name, col=mag)
    ABmag_sdss_out.add_column(
        col=Column(
            name=mag_name +
            '_err',
            data=err),
        index=None,
        name=mag_name +
        '_err')

for mag_name in ABmag_sdss_out.colnames:
    agn_hdus.add_column(ABmag_sdss_out[mag_name], name=mag_name)

os.system('rm ' + catalog_input)
agn_hdus.write(catalog_input, format='fits')
