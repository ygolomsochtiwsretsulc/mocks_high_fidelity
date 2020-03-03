import os
import sys
import glob
import numpy as n

# astropy
import astropy.units as u
from astropy.table import Table, Column
import astropy.io.fits as fits
import astropy.constants as cc

# speclite
import speclite.filters
from speclite.filters import FilterConvolution
from speclite.filters import ab_reference_flux

# where the filters are
path = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'data',
    'photometry',
    'filters')

# where the templates are
template_cigale_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'data',
    'templates',
    'template_cigale')
template_lephare_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'data',
    'templates',
    'template_lephare')


############################################################
############################################################
############################################################
# GALEX filters
############################################################
############################################################
############################################################

Tab_NUV = Table.read(os.path.join(path, 'galex', 'NUV.pb'), format='ascii')
Tab_FUV = Table.read(os.path.join(path, 'galex', 'FUV.pb'), format='ascii')

Pass_NUV = speclite.filters.FilterResponse(
    wavelength=Tab_NUV['col1'] * u.AA,
    response=Tab_NUV['col2'] / 100.,
    meta=dict(
        group_name='GALEX',
        band_name='NUV'))

Pass_FUV = speclite.filters.FilterResponse(
    wavelength=Tab_FUV['col1'] * u.AA,
    response=Tab_FUV['col2'] / 100.,
    meta=dict(
        group_name='GALEX',
        band_name='FUV'))

############################################################
############################################################
############################################################
# VISTA FILTERS
############################################################
############################################################
############################################################
Tab_Z = Table.read(
    os.path.join(
        path,
        'vista',
        'VISTA_Filters_at80K_forETC_Z.dat'),
    format='ascii')
ok_Z = (Tab_Z['col1'] > 700) & (Tab_Z['col1'] < 1000) & (Tab_Z['col2'] > 0)
Tab_Z['col2'][(ok_Z == False)] = 0.
ok_Z = (Tab_Z['col1'] > 700) & (Tab_Z['col1'] < 1000)

Tab_Y = Table.read(
    os.path.join(
        path,
        'vista',
        'VISTA_Filters_at80K_forETC_Y.dat'),
    format='ascii')
ok_Y = (Tab_Y['col1'] > 800) & (Tab_Y['col1'] < 1300) & (Tab_Y['col2'] > 0)
Tab_Y['col2'][(ok_Y == False)] = 0.
ok_Y = (Tab_Y['col1'] > 800) & (Tab_Y['col1'] < 1300)

Tab_J = Table.read(
    os.path.join(
        path,
        'vista',
        'VISTA_Filters_at80K_forETC_J.dat'),
    format='ascii')
ok_J = (Tab_J['col1'] > 1000) & (Tab_J['col1'] < 1500) & (Tab_J['col2'] > 0)
Tab_J['col2'][(ok_J == False)] = 0.
ok_J = (Tab_J['col1'] > 1000) & (Tab_J['col1'] < 1500)

Tab_H = Table.read(
    os.path.join(
        path,
        'vista',
        'VISTA_Filters_at80K_forETC_H.dat'),
    format='ascii')
ok_H = (Tab_H['col1'] > 1200) & (Tab_H['col1'] < 2000) & (Tab_H['col2'] > 0)
Tab_H['col2'][(ok_H == False)] = 0.
ok_H = (Tab_H['col1'] > 1200) & (Tab_H['col1'] < 2000)

Tab_K = Table.read(
    os.path.join(
        path,
        'vista',
        'VISTA_Filters_at80K_forETC_Ks.dat'),
    format='ascii')
ok_K = (Tab_K['col1'] > 1700) & (Tab_K['col1'] < 2700) & (Tab_K['col2'] > 0)
Tab_K['col2'][(ok_K == False)] = 0.
ok_K = (Tab_K['col1'] > 1700) & (Tab_K['col1'] < 2700)

Pass_Z = speclite.filters.FilterResponse(
    wavelength=Tab_Z['col1'] * u.nm,
    response=Tab_Z['col2'] / 100.,
    meta=dict(
        group_name='VISTA',
        band_name='Z'))

Pass_Y = speclite.filters.FilterResponse(
    wavelength=Tab_Y['col1'] * u.nm,
    response=Tab_Y['col2'] / 100.,
    meta=dict(
        group_name='VISTA',
        band_name='Y'))

Pass_J = speclite.filters.FilterResponse(
    wavelength=Tab_J['col1'] * u.nm,
    response=Tab_J['col2'] / 100.,
    meta=dict(
        group_name='VISTA',
        band_name='J'))

Pass_H = speclite.filters.FilterResponse(
    wavelength=Tab_H['col1'] * u.nm,
    response=Tab_H['col2'] / 100.,
    meta=dict(
        group_name='VISTA',
        band_name='H'))

Pass_K = speclite.filters.FilterResponse(
    wavelength=Tab_K['col1'] * u.nm,
    response=Tab_K['col2'] / 100.,
    meta=dict(
        group_name='VISTA',
        band_name='K'))


############################################################
############################################################
############################################################
# Herschel filters
############################################################
############################################################
############################################################

Tab_Hs100 = Table.read(
    os.path.join(
        path,
        'herschel',
        'PACS_filter_100.txt'),
    format='ascii')
ok = (
    Tab_Hs100['col1'] > 7.5e5) & (
        Tab_Hs100['col1'] < 1.4e6) & (
            Tab_Hs100['col2'] > 0)
Tab_Hs100['col2'][(ok == False)] = 0.

Tab_Hs160 = Table.read(
    os.path.join(
        path,
        'herschel',
        'PACS_filter_160.txt'),
    format='ascii')
ok = (
    Tab_Hs160['col1'] > 1.12e6) & (
        Tab_Hs160['col1'] < 2.6e6) & (
            Tab_Hs160['col2'] > 0)
Tab_Hs160['col2'][(ok == False)] = 0.

Tab_Hs250 = Table.read(
    os.path.join(
        path,
        'herschel',
        'SPIRE_filter_250.txt'),
    format='ascii')
ok = (
    Tab_Hs250['col1'] > 1.84e6) & (
        Tab_Hs250['col1'] < 3.1e6) & (
            Tab_Hs250['col2'] > 0)
Tab_Hs250['col2'][(ok == False)] = 0.

Tab_Hs350 = Table.read(
    os.path.join(
        path,
        'herschel',
        'SPIRE_filter_350.txt'),
    format='ascii')
ok = (
    Tab_Hs350['col1'] > 2.66e6) & (
        Tab_Hs350['col1'] < 4.48e6) & (
            Tab_Hs350['col2'] > 0)
Tab_Hs350['col2'][(ok == False)] = 0.

Tab_Hs500 = Table.read(
    os.path.join(
        path,
        'herschel',
        'SPIRE_filter_500.txt'),
    format='ascii')
ok = (
    Tab_Hs500['col1'] > 3.74e6) & (
        Tab_Hs500['col1'] < 7.56e6) & (
            Tab_Hs500['col2'] > 0)
Tab_Hs500['col2'][(ok == False)] = 0.

Pass_Hs100 = speclite.filters.FilterResponse(
    wavelength=Tab_Hs100['col1'] * u.AA,
    response=Tab_Hs100['col2'],
    meta=dict(
        group_name='HERSCHEL',
        band_name='h100'))

Pass_Hs160 = speclite.filters.FilterResponse(wavelength=Tab_Hs160['col1'][n.argsort(
    Tab_Hs160['col1'])] * u.AA, response=Tab_Hs160['col2'][n.argsort(Tab_Hs160['col1'])], meta=dict(group_name='HERSCHEL', band_name='h160'))

Pass_Hs250 = speclite.filters.FilterResponse(wavelength=Tab_Hs250['col1'][n.argsort(
    Tab_Hs250['col1'])] * u.AA, response=Tab_Hs250['col2'][n.argsort(Tab_Hs250['col1'])], meta=dict(group_name='HERSCHEL', band_name='h250'))

Pass_Hs350 = speclite.filters.FilterResponse(wavelength=Tab_Hs350['col1'][n.argsort(
    Tab_Hs350['col1'])] * u.AA, response=Tab_Hs350['col2'][n.argsort(Tab_Hs350['col1'])], meta=dict(group_name='HERSCHEL', band_name='h350'))

Pass_Hs500 = speclite.filters.FilterResponse(wavelength=Tab_Hs500['col1'][n.argsort(
    Tab_Hs500['col1'])] * u.AA, response=Tab_Hs500['col2'][n.argsort(Tab_Hs500['col1'])], meta=dict(group_name='HERSCHEL', band_name='h500'))


############################################################
############################################################
############################################################
# HSC filters
############################################################
############################################################
############################################################

Tab_HSC_g = Table.read(os.path.join(path, 'hsc', 'g_HSC.txt'), format='ascii')
ok = (
    Tab_HSC_g['col1'] > 3800) & (
        Tab_HSC_g['col1'] < 5677) & (
            Tab_HSC_g['col2'] > 0)
Tab_HSC_g['col2'][(ok == False)] = 0.

Tab_HSC_r = Table.read(os.path.join(path, 'hsc', 'r_HSC.txt'), format='ascii')
ok = (
    Tab_HSC_r['col1'] > 5190) & (
        Tab_HSC_r['col1'] < 7300) & (
            Tab_HSC_r['col2'] > 0)
Tab_HSC_r['col2'][(ok == False)] = 0.

Tab_HSC_i = Table.read(os.path.join(path, 'hsc', 'i_HSC.txt'), format='ascii')
ok = (
    Tab_HSC_i['col1'] > 6780) & (
        Tab_HSC_i['col1'] < 8800) & (
            Tab_HSC_i['col2'] > 0)
Tab_HSC_i['col2'][(ok == False)] = 0.

Tab_HSC_z = Table.read(os.path.join(path, 'hsc', 'z_HSC.txt'), format='ascii')
ok = (
    Tab_HSC_z['col1'] > 8288) & (
        Tab_HSC_z['col1'] < 9580) & (
            Tab_HSC_z['col2'] > 0)
Tab_HSC_z['col2'][(ok == False)] = 0.

Tab_HSC_y = Table.read(os.path.join(path, 'hsc', 'y_HSC.txt'), format='ascii')
ok = (
    Tab_HSC_y['col1'] > 9130) & (
        Tab_HSC_y['col1'] < 10800) & (
            Tab_HSC_y['col2'] > 0)
Tab_HSC_y['col2'][(ok == False)] = 0.

Pass_HSC_g = speclite.filters.FilterResponse(
    wavelength=Tab_HSC_g['col1'] * u.AA,
    response=Tab_HSC_g['col2'],
    meta=dict(
        group_name='HSC',
        band_name='g'))

Pass_HSC_r = speclite.filters.FilterResponse(
    wavelength=Tab_HSC_r['col1'] * u.AA,
    response=Tab_HSC_r['col2'],
    meta=dict(
        group_name='HSC',
        band_name='r'))

Pass_HSC_i = speclite.filters.FilterResponse(
    wavelength=Tab_HSC_i['col1'] * u.AA,
    response=Tab_HSC_i['col2'],
    meta=dict(
        group_name='HSC',
        band_name='i'))

Pass_HSC_z = speclite.filters.FilterResponse(
    wavelength=Tab_HSC_z['col1'] * u.AA,
    response=Tab_HSC_z['col2'],
    meta=dict(
        group_name='HSC',
        band_name='z'))

Pass_HSC_y = speclite.filters.FilterResponse(
    wavelength=Tab_HSC_y['col1'] * u.AA,
    response=Tab_HSC_y['col2'],
    meta=dict(
        group_name='HSC',
        band_name='y'))

#p.plot(Tab_HSC_g['col1'], Tab_HSC_g['col2'] )
#p.plot(Tab_HSC_r['col1'], Tab_HSC_r['col2'] )
#p.plot(Tab_HSC_i['col1'], Tab_HSC_i['col2'] )
#p.plot(Tab_HSC_z['col1'], Tab_HSC_z['col2'] )
#p.plot(Tab_HSC_y['col1'], Tab_HSC_y['col2'] )
# p.xscale('log')
# p.show()

############################################################
############################################################
############################################################
# WISE
############################################################
############################################################
############################################################
# RSR-W1.txt  RSR-W2.txt  RSR-W3.txt  RSR-W4.txt

Tab_WISE_1 = Table.read(
    os.path.join(
        path,
        'wise',
        'RSR-W1.txt'),
    format='ascii')
ok = (
    Tab_WISE_1['col1'] > 2.6) & (
        Tab_WISE_1['col1'] < 3.93) & (
            Tab_WISE_1['col2'] > 0)
Tab_WISE_1['col2'][(ok == False)] = 0.

Tab_WISE_2 = Table.read(
    os.path.join(
        path,
        'wise',
        'RSR-W2.txt'),
    format='ascii')
ok = (
    Tab_WISE_2['col1'] > 3.93) & (
        Tab_WISE_2['col1'] < 5.5) & (
            Tab_WISE_2['col2'] > 0)
Tab_WISE_2['col2'][(ok == False)] = 0.

Tab_WISE_3 = Table.read(
    os.path.join(
        path,
        'wise',
        'RSR-W3.txt'),
    format='ascii')
ok = (
    Tab_WISE_3['col1'] > 7.24) & (
        Tab_WISE_3['col1'] < 17.6) & (
            Tab_WISE_3['col2'] > 0)
Tab_WISE_3['col2'][(ok == False)] = 0.

Tab_WISE_4 = Table.read(
    os.path.join(
        path,
        'wise',
        'RSR-W4.txt'),
    format='ascii')
ok = (
    Tab_WISE_4['col1'] > 19.) & (
        Tab_WISE_4['col1'] < 28) & (
            Tab_WISE_4['col2'] > 0)
Tab_WISE_4['col2'][(ok == False)] = 0.

#p.plot(Tab_WISE_1['col1'], Tab_WISE_1['col2'] )
#p.plot(Tab_WISE_2['col1'], Tab_WISE_2['col2'] )
#p.plot(Tab_WISE_3['col1'], Tab_WISE_3['col2'] )
#p.plot(Tab_WISE_4['col1'], Tab_WISE_4['col2'] )
# p.xscale('log')
# p.show()

Pass_WISE_1 = speclite.filters.FilterResponse(
    wavelength=Tab_WISE_1['col1'] *
    u.micron,
    response=Tab_WISE_1['col2'],
    meta=dict(
        group_name='WISE',
        band_name='W1'))
Pass_WISE_2 = speclite.filters.FilterResponse(
    wavelength=Tab_WISE_2['col1'] *
    u.micron,
    response=Tab_WISE_2['col2'],
    meta=dict(
        group_name='WISE',
        band_name='W2'))
Pass_WISE_3 = speclite.filters.FilterResponse(
    wavelength=Tab_WISE_3['col1'] *
    u.micron,
    response=Tab_WISE_3['col2'],
    meta=dict(
        group_name='WISE',
        band_name='W3'))
Pass_WISE_4 = speclite.filters.FilterResponse(
    wavelength=Tab_WISE_4['col1'] *
    u.micron,
    response=Tab_WISE_4['col2'],
    meta=dict(
        group_name='WISE',
        band_name='W4'))


#sdss_filters = speclite.filters.load_filters('sdss2010-g','sdss2010-r','sdss2010-i', 'bessell-V')
sdss_r_filter = speclite.filters.load_filters('sdss2010-r')
all_filters = speclite.filters.load_filters(
    'GALEX-NUV',
    'GALEX-FUV',
    'HSC-g',
    'HSC-r',
    'HSC-i',
    'HSC-z',
    'HSC-y',
    #'VISTA-Z',
    #'VISTA-Y',
    'VISTA-J',
    'VISTA-H',
    'VISTA-K',
    'WISE-W1',
    'WISE-W2', )
    #'sdss2010-u', 
    #'sdss2010-g',
    #'sdss2010-r',
    #'sdss2010-i', 
    #'sdss2010-z')


def get_cigale_template(path_2_template, redshift=0.317):
    #all_tpl = n.array(glob.glob(os.path.join(template_cigale_dir,'*.fits')))
    # all_tpl.sort()
    #filename = all_tpl[0]
    hdu = fits.open(path_2_template)
    dat = hdu[1].data
    nu = (cc.c / (dat['wavelength'] * u.nm.to(u.m))).value  # Hz
    ll = dat['wavelength'] / 10 * u.AA
    flambda = dat['Fnu'] * 10**(-26) * nu / ll * u.erg * u.cm**(-2) * u.s**(-1)
    #dat['Fnu']*u.mJy (erg/cm2/s/Hz)
    ll_rf = ll / (1 + redshift)
    return flambda, ll_rf, hdu


def get_lephare_template(path_2_template):
    all_tpl = sorted(
        n.array(
            glob.glob(
                os.path.join(
                    template_lephare_dir,
                    '*.fits'))))
    filename = all_tpl[0]
    dat = fits.open(path_2_template)[1].data
    nu = (cc.c / (dat['wavelength'] * u.nm.to(u.m))).value  # Hz
    ll = dat['wavelength'] / 10 * u.AA
    flambda = dat['Fnu'] * 10**(-26) * nu / ll * u.erg * u.cm**(-2) * u.s**(-1)
    #dat['Fnu']*u.mJy (erg/cm2/s/Hz)
    ll_rf = ll / (1 + redshift)
    return flambda, ll_rf


log_error_scaling_parameters = {
	'GALEX-NUV': n.array([7.24492805e-03, -4.28564399e-04, -5.13126940e+00]), 
	'GALEX-FUV': n.array([6.39978442e-03, -2.43684633e-04, -4.67703435e+00]), 
	'HSC-g': n.array([0.33766153, -9.72731761]), 
	'HSC-r': n.array([0.34643741, -9.78832194]), 
	'HSC-i': n.array([0.35353582, -9.90507919]), 
	'HSC-z': n.array([0.35376108, -9.64655357]), 
	'HSC-y': n.array([0.3410503, -9.07559621]), 
	'VISTA-J': n.array([0.00979343, -0.0105481, -4.61121865]), 
	'VISTA-H': n.array([0.01291702, -0.09107075, -3.91056997]), 
	'VISTA-K': n.array([0.0141926, -0.12466629, -3.52130742]), 
	'WISE-W1': n.array([0.39463068, -9.12916917]),
	'WISE-W2': n.array([0.39637782, -8.7972894])}


def log_error_scaling(pars, magnitude): return n.polyval(pars, magnitude)


# log_error_scaling ={
# 'HSC-g'   : lambda mag : 0.344*mag - 9.49,
# 'HSC-r'   : lambda mag : 0.344*mag - 9.49,
# 'HSC-i'   : lambda mag : 0.344*mag - 9.49,
# 'HSC-z'   : lambda mag : 0.344*mag - 9.49,
# 'HSC-y'   : lambda mag : 0.344*mag - 9.49,
# 'VISTA-Z': lambda mag : 0.344*mag - 9.49,
# 'VISTA-Y': lambda mag : 0.344*mag - 9.49,
# 'VISTA-J': lambda mag : 0.344*mag - 9.49,
# 'VISTA-H': lambda mag : 0.344*mag - 9.49,
# 'VISTA-K': lambda mag : 0.344*mag - 9.49,
# 'WISE-W1' : lambda mag : 0.344*mag - 9.49,
# 'WISE-W2' : lambda mag : 0.344*mag - 9.49,
# 'WISE-W3' : lambda mag : 0.344*mag - 9.49,
# 'WISE-W4' : lambda mag : 0.344*mag - 9.49,
# 'HERSCHEL-h100':  lambda mag : 0.344*mag - 9.49,
# 'HERSCHEL-h160':  lambda mag : 0.344*mag - 9.49,
# 'HERSCHEL-h250':  lambda mag : 0.344*mag - 9.49,
# 'HERSCHEL-h350':  lambda mag : 0.344*mag - 9.49,
# 'HERSCHEL-h500':  lambda mag : 0.344*mag - 9.49
# }

#from scipy.stats import norm
# random variable in a norm distribution with sigma^2 = 1/80
#rds = norm.rvs( loc=0, scale = 80**(-0.5), size = len( ABmag_sdss_out[ ABmag_sdss_out.colnames[0] ] ) )
#log_error_scatter = log_error_scaling[ABmag_sdss_out.colnames[0]](ABmag_sdss_out[ABmag_sdss_out.colnames[0]])

#error = 10**(log_error_scatter + rds)
#mag = ABmag_sdss_out[ABmag_sdss_out.colnames[0]] + error
#mag_err = 10**(log_error_scatter)

# hdu_cols  = fits.ColDefs([
# fits.Column( name=ccn , unit='' , format='D',
# array=ABmag_sdss_out[ccn] ) for ccn in ABmag_sdss_out.colnames ])


# hdu_cols_e  = fits.ColDefs([
# fits.Column( name=ccn+"_err" , unit='' , format='D',
# array=error_scaling[ccn](ABmag_sdss_out[ccn]) ) for ccn in
# ABmag_sdss_out.colnames ])


#tb_hdu = fits.BinTableHDU.from_columns( agn_hdus[1].data.columns + hdu_cols + hdu_cols_e)
#prihdr = fits.Header()
#prihdr['author'] = 'JC'
#prihdu = fits.PrimaryHDU(header=prihdr)
#thdulist = fits.HDUList([prihdu, tb_hdu])
# if os.path.isfile(catalog_name):
#os.system("rm "+catalog_name)
# thdulist.writeto(catalog_name)


# add scatter ?
# add offsets ?

# sys.exit()


#p.figure(0, (10,5))
#p.plot(ll/(1+zz), flambda, ls='solid', lw=1.5, marker=',', label=os.path.basename(filename))
# p.title(os.path.basename(tp))
# p.xlim((100.,40000.))
##p.ylim((8e-17, 3e-14))
# p.yscale('log')
# p.xscale('log')
#p.axvline(1215, color='b', ls='dashed', label='1215 Lya')
#p.axvline(1546, color='c', ls='dashed', label='1546 CIV')
#p.axvline(2800, color='m', ls='dashed', label='2800 MgII')
#p.axvline(3727, color='g', ls='dashed', label='3727 [OII]')
#p.axvline(5007, color='r', ls='dashed', label='5007 [OIII]')
#p.axvline(6565, color='k', ls='dashed', label='6565 Ha')
#p.xlabel('wavelength A')
#p.ylabel('flux erg/cm2/s/A')
#p.legend(frameon=False, loc=0)
# p.savefig(os.path.join(template_cigale_dir,'images',os.path.basename(filename)[:-5]+".png"))
# p.clf()


# sys.exit()

# def read_template_image_s6(filename):
#h = fits.open(filename)
#FLUX = h[0].data
#LAMBDA = h[0].header['CRVAL1'] + h[0].header['CDELT1'] * n.arange(1, len(FLUX)+1, 1) - h[0].header['CRPIX1']
#rmag = float(os.path.basename(filename).split('_')[3][2:])/10.
# return os.path.basename(filename), rmag, LAMBDA, FLUX


#p.figure(0, (10,5))
# for ii, tp in enumerate(all_tpl):
#dat = fits.open(tp)[1].data
# print(tp)
#bns = os.path.basename(tp).split('_')
# if bns[1]=="AGN" or bns[1]=="qso":
# zz=(float(bns[4])*0.1+float(bns[6])*0.1)*0.5
# else:
# zz=(float(bns[3])*0.1+float(bns[5])*0.1)*0.5
# print(zz)
# if bns[1]=="LRG":
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z185_mr165_ellgl.fits"
#nnn, mmm, w185, f185 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z155_mr165_ellgl.fits"
#nnn, mmm, w155, f155 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z130_mr165_ellgl.fits"
#nnn, mmm, w130, f130 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z090_mr165_ellgl.fits"
#nnn, mmm, w090, f090 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z070_mr165_ellgl.fits"
#nnn, mmm, w070, f070 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z030_mr165_ellgl.fits"
#nnn, mmm, w030, f030 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v2.0_z010_mr165_ellgl.fits"
#nnn, mmm, w010, f010 = read_template_image_s6(file_s6)

# if zz<=1.9 and zz >1.7 :
#p.plot(w185,f185,label='v1 z185')
# elif zz<=1.7 and zz >1.5 :
#p.plot(w155,f155,label='v1 z155')
# elif zz<=1.5 and zz >1.1 :
#p.plot(w130,f130,label='v1 z130')
# elif zz<=1.1 and zz >0.8 :
#p.plot(w090,f090,label='v1 z090')
# elif zz<=0.8 and zz >0.5 :
#p.plot(w070,f070,label='v1 z070')
# elif zz<=0.5 and zz > 0.3  :
#p.plot(w030,f030,label='v1 z030')
# elif zz<=0.3 :
#p.plot(w010,f010,label='v1 z010')

# if bns[1]=="qso":
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z600_mr165_type1.fits"
#nnn, mmm, w600, f600 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z450_mr165_type1.fits"
#nnn, mmm, w450, f450 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z350_mr165_type1.fits"
#nnn, mmm, w350, f350 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z275_mr165_type1.fits"
#nnn, mmm, w275, f275 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z225_mr165_type1.fits"
#nnn, mmm, w225, f225 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z185_mr165_type1.fits"
#nnn, mmm, w185, f185 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z155_mr165_type1.fits"
#nnn, mmm, w155, f155 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z130_mr165_type1.fits"
#nnn, mmm, w130, f130 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z090_mr165_type1.fits"
#nnn, mmm, w090, f090 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z070_mr165_type1.fits"
#nnn, mmm, w070, f070 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z030_mr165_type1.fits"
#nnn, mmm, w030, f030 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z010_mr165_type1.fits"
#nnn, mmm, w010, f010 = read_template_image_s6(file_s6)

# if zz > 4.4 and zz<=4.6 :
##p.plot(w600,f600,label='v1 z600')
#p.plot(w450,f450,label='v1 z450')
##p.plot(w350,f350,label='v1 z350')
# elif zz<=3 and zz >2.5 :
#p.plot(w275,f275,label='v1 z275')
# elif zz<=2.5 and zz >1.9 :
#p.plot(w225,f225,label='v1 z225')
# elif zz<=1.9 and zz >1.7 :
#p.plot(w185,f185,label='v1 z185')
# elif zz<=1.7 and zz >1.5 :
#p.plot(w155,f155,label='v1 z155')
# elif zz<=1.5 and zz >1.1 :
#p.plot(w130,f130,label='v1 z130')
# elif zz<=1.1 and zz >0.8 :
#p.plot(w090,f090,label='v1 z090')
# elif zz<=0.8 and zz >0.5 :
#p.plot(w070,f070,label='v1 z070')
# elif zz<=0.5 and zz > 0.3  :
#p.plot(w030,f030,label='v1 z030')
# elif zz<=0.3 :
#p.plot(w010,f010,label='v1 z010')

# if bns[1]=="AGN":
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z600_mr165_type2.fits"
#nnn, mmm, w600, f600 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z450_mr165_type2.fits"
#nnn, mmm, w450, f450 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z350_mr165_type2.fits"
#nnn, mmm, w350, f350 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z275_mr165_type2.fits"
#nnn, mmm, w275, f275 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z225_mr165_type2.fits"
#nnn, mmm, w225, f225 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z185_mr165_type2.fits"
#nnn, mmm, w185, f185 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z155_mr165_type2.fits"
#nnn, mmm, w155, f155 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z130_mr165_type2.fits"
#nnn, mmm, w130, f130 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z090_mr165_type2.fits"
#nnn, mmm, w090, f090 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z070_mr165_type2.fits"
#nnn, mmm, w070, f070 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z030_mr165_type2.fits"
#nnn, mmm, w030, f030 = read_template_image_s6(file_s6)
#file_s6 = "/home/comparat/data-4most/opsim_data/survey_team_inputs/S6/info_2016-03-23_eROSITA-DE/templates/AGN_v1.0_z010_mr165_type2.fits"
#nnn, mmm, w010, f010 = read_template_image_s6(file_s6)
# if zz > 4.4 and zz<=4.6 :
##p.plot(w600,f600,label='v1 z600')
#p.plot(w450,f450,label='v1 z450')
##p.plot(w350,f350,label='v1 z350')
# elif zz<=3 and zz >2.5 :
#p.plot(w275,f275,label='v1 z275')
# elif zz<=2.5 and zz >1.9 :
#p.plot(w225,f225,label='v1 z225')
# elif zz<=1.9 and zz >1.7 :
#p.plot(w185,f185,label='v1 z185')
# elif zz<=1.7 and zz >1.5 :
#p.plot(w155,f155,label='v1 z155')
# elif zz<=1.5 and zz >1.1 :
#p.plot(w130,f130,label='v1 z130')
# elif zz<=1.1 and zz >0.8 :
#p.plot(w090,f090,label='v1 z090')
# elif zz<=0.8 and zz >0.5 :
#p.plot(w070,f070,label='v1 z070')
# elif zz<=0.5 and zz > 0.3  :
#p.plot(w030,f030,label='v1 z030')
# elif zz<=0.3 :
#p.plot(w010,f010,label='v1 z010')

#p.plot(dat['LAMBDA'], dat['FLUX'], ls='solid', lw=1.5, marker=',', label=os.path.basename(tp))
# p.title(os.path.basename(tp))
# p.xlim((3800.,9500.))
#p.ylim((8e-17, 3e-14))
# p.yscale('log')
#p.axvline(1215*(1+zz), color='b', ls='dashed', label='1215 Lya')
#p.axvline(1546*(1+zz), color='c', ls='dashed', label='1546 CIV')
#p.axvline(2800*(1+zz), color='m', ls='dashed', label='2800 MgII')
#p.axvline(3727*(1+zz), color='g', ls='dashed', label='3727 [OII]')
#p.axvline(5007*(1+zz), color='r', ls='dashed', label='5007 [OIII]')
#p.axvline(6565*(1+zz), color='k', ls='dashed', label='6565 Ha')
#p.xlabel('wavelength A')
#p.ylabel('flux erg/cm2/s/A')
#p.legend(frameon=False, loc=0)
# p.savefig(os.path.join(template_cigale_dir,'images',os.path.basename(tp)[:-5]+".png"))
# p.clf()
