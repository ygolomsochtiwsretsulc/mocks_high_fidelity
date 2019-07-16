"""
What it does
------------

Plots the tabulated the logNlogS,

More plots can be added.

References
----------

Command to run
--------------

python3 003_3_agn_plot_logNlogS.py environmentVAR ftyp

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

ftyp: 'all' or 'sat', type of AGN populating a central halo or satellite haloes.

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, matplotlib

"""
import glob
import sys
from astropy_healpix import healpy
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as p
import matplotlib
import astropy.io.fits as fits
import h5py
import numpy as n
print('PLOT XLF and logNlogS')
print('------------------------------------------------')
print('------------------------------------------------')

t0 = time.time()
#import astropy.io.fits as fits

area = healpy.nside2pixarea(8, degrees=True)

matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
# import all pathes

env = sys.argv[1]  # 'MD04'
ftyp = sys.argv[2]  # "sat"
print(env, ftyp)


if env == "MD10" or env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    z_maximum = 1.4
if env == "MD04":
    z_maximum = 0.45

DZ = 0.2
all_zs = n.arange(0., z_maximum, DZ)

# simulation setup

if env == "MD10" or env == "MD04":
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

root_dir = os.path.join(os.environ[env], 'cat_AGN_' + ftyp)

fig_dir = os.path.join(root_dir, 'figures')
logNlogS_dir = os.path.join(root_dir, 'logNlogS')
logNlogR_dir = os.path.join(root_dir, 'logNlogR')
XLF_dir = os.path.join(root_dir, 'XLF')

#DATA_XLF = n.transpose([xf, nhar, (nhar)*(1-(nharN)**(-0.5)), (nhar)*(1+(nharN)**(-0.5)), N_nh20, N_nh22, N_nh24, phi_h(10**xf,z_mean)])

fx_bins = n.arange(-20, -8., 0.2)
x_fx = fx_bins[:-1] + 0.1


def get_lognlogs(hpx_id='000061'):
    ff = os.path.join(logNlogS_dir, 'logNlogS_soft_' + hpx_id + '.ascii')
    # print(ff)
    outout = n.loadtxt(ff, unpack=True)
    #print(outout, outout.shape)
    # data=n.array(data)
    #N_pixels = len(data)
    NN = outout / area
    itp = interp1d(x_fx, n.log10(NN))
    return itp

#DATA_X = n.transpose(out)[(out[1]>0) & (out[2]!=n.inf)]
#n.savetxt(os.path.join(logNlogS_dir, 'logNlogS_soft_'+baseName+'.ascii'), DATA_X  )


hpx_ids = n.array(['000310', '000311', '000312', '000313',
                   '000314', '000315', '000316'])

p.figure(1, (6, 6))
#p.axhline(n.log10(300), ls='dashed')

for hpx_id in hpx_ids:
    itp = get_lognlogs(hpx_id)
    p.plot(itp.x, itp.y, rasterized=True, label=hpx_id, lw=1, ls='dashed')

#x, y_1 = get_lognlogs('C1')
#x, y_1_t1 = get_lognlogs_t1('C1')
#p.plot(x, y_1, rasterized = True, label = 'C1'  , lw=2, ls='dashed')

#x, y_2 = get_lognlogs('C2')
#x, y_2_t1 = get_lognlogs_t1('C2')
#p.plot(x, y_2, rasterized = True, label = 'C2' , lw=2, ls='dashed' )

#x, y_3 = get_lognlogs('C3')
#x, y_3_t1 = get_lognlogs_t1('C3')
#p.plot(x, y_3, rasterized = True, label = 'C3'  , lw=2, ls='dashed')

#x, y_4 = get_lognlogs('C4')
#x, y_4_t1 = get_lognlogs_t1('C4')
#p.plot(x, y_4, rasterized = True, label = 'C4'  , lw=2, ls='dashed')

#x, y_5 = get_lognlogs('C5')
#x, y_5_t1 = get_lognlogs_t1('C5')
#p.plot(x, y_5, rasterized = True, label = 'C5'  , lw=2, ls='dashed')

# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_VS"],
    'data',
    'logNlogS',
    'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='g')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_VS"],
    'data',
    'logNlogS',
    'logNlogS_Merloni_12_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='r')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_VS"],
    'data/logNlogS/logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(x_data, n.log10(y_data), lw=3, ls='dotted', color='b')

p.xlabel('log(F_X[0.5-2 keV])')
p.ylabel('log(>F_X) [/deg2]')
p.legend(frameon=False, loc=0)
# p.yscale('log')
p.xlim((-17, -11.5))
p.ylim((-2, 4.2))
# p.title('Mocks')
p.grid()
p.savefig(os.path.join(logNlogS_dir, "logN_logS_AGN.png"))
p.clf()


sys.exit()

for id_z in n.arange(len(all_zs)):
    zmin = all_zs[id_z]
    zmax = all_zs[id_z] + DZ
    z_mean = 0.5 * (zmin + zmax)
    file_list = sorted(n.array(glob.glob(os.path.join(XLF_dir,
                                                      'XLF_soft_*_' + str(n.round(zmin,
                                                                                  2)) + '_z_' + str(n.round(zmax,
                                                                                                            2)) + '.ascii'))))
    data = []
    for ff in file_list:
        data.append(n.loadtxt(ff, unpack=True))
        # DATA_XLF = n.transpose([xf, nhar, (nhar)*(1-(nharN)**(-0.5)), (nhar)*(1+(nharN)**(-0.5)), N_nh20, N_nh22, N_nh24, phi_h(10**xf,z_mean)])

p.figure(1, (6, 6))
p.axes([0.2, 0.18, 0.75, 0.75])
# Aird 2015
z_mean = (zmin + zmax) * 0.5 * n.ones_like(xf)
# mock
p.plot(xf, phi_h(10**xf, z_mean), c='cyan', ls='dashed',
       lw=2, label='Ai15')  # Aird 2-10 keV LADE')
p.fill_between(xf,
               y1=(nhar) * (1 - (nharN)**(-0.5)),
               y2=(nhar) * (1 + (nharN)**(-0.5)),
               color='g',
               alpha=0.7,
               label='Mock',
               lw=2)  # 2-10  keV')

p.plot(xf, N_nh20, label='nH<22', ls='dotted', lw=2)
p.plot(xf, N_nh22, label='22<nH<24', ls='dotted', lw=2)
p.plot(xf, N_nh24, label='24<nH', ls='dotted', lw=2)


p.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
p.legend(loc=3, fontsize=14)
p.yscale('log')
p.xlim((37., 46.5))
p.ylim((1 / (2 * vol), 1e-2))
p.title(str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
p.grid()
p.savefig(os.path.join(XLF_dir, "XLF_soft_" + baseName + ".png"))
p.clf()

DATA_XLF = n.transpose([xf,
                        nhar,
                        (nhar) * (1 - (nharN)**(-0.5)),
                        (nhar) * (1 + (nharN)**(-0.5)),
                        N_nh20,
                        N_nh22,
                        N_nh24,
                        phi_h(10**xf,
                              z_mean)])
n.savetxt(os.path.join(XLF_dir, 'XLF_soft_' + baseName + '.ascii'), DATA_XLF)


# XLF_ratio_
p.figure(1, (6, 6))
p.axes([0.18, 0.18, 0.75, 0.75])
p.plot(xf, nhar / phi_h(10**xf, z_mean), label='2-10 keV')
p.fill_between(xf, y1=1 - (nharN)**(-0.5), y2=1 + (nharN) **
               (-0.5), color='green', alpha=0.3, label='mock ERR')
p.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
p.ylabel(r'mock/model')
p.legend(frameon=False, loc=0)
p.xlim((37., 46.5))
p.ylim((0.7, 1.3))
p.title(str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
p.grid()
p.savefig(os.path.join(XLF_dir, "XLF_ratio_" + baseName + ".png"))
p.clf()


# LSAR binning
dlogf = 0.1
fbins = n.arange(30, 38, dlogf)
xf = n.arange(30, 38, dlogf)[:-1] + dlogf / 2.

nall = n.histogram(lsar, fbins)[0] / vol / dlogf
nallN = n.histogram(lsar, fbins)[0]

zsel = (logm >= 12)
nall_12 = n.histogram(lsar, fbins)[0] / vol / dlogf
nallN_12 = n.histogram(lsar, fbins)[0]

zsel = (logm >= 11) & (logm < 12)
nall_11 = n.histogram(lsar, fbins)[0] / vol / dlogf
nallN_11 = n.histogram(lsar, fbins)[0]

zsel = (logm >= 10) & (logm < 11)
nall_10 = n.histogram(lsar, fbins)[0] / vol / dlogf
nallN_10 = n.histogram(lsar, fbins)[0]

zsel = (logm >= 9) & (logm < 10)
nall_9 = n.histogram(lsar, fbins)[0] / vol / dlogf
nallN_9 = n.histogram(lsar, fbins)[0]

p.figure(1, (6, 6))
p.axes([0.16, 0.15, 0.8, 0.8])

#nall, nallN, nall_9, nallN_9, nall_10, nallN_10, nall_11, nallN_11, nall_12, nallN_12 = get_LSAR_hist(fbins = fbins, zmin=zmin, zmax=zmax)

fun = interp1d(xf, nall)
nrm = quad(fun, xf.min(), xf.max())[0]
p.plot(xf, nall / nrm, 'k', lw=3)  # , label='mock all'

fun = interp1d(xf, nall_9)
nrm = quad(fun, xf.min(), xf.max())[0]
p.plot(xf, nall_9 / nrm, 'g', lw=2)  # , label='9-10'

fun = interp1d(xf, nall_10)
nrm = quad(fun, xf.min(), xf.max())[0]
p.plot(xf, nall_10 / nrm, 'r', lw=2)  # , label='10-11'

p.xlabel(r'$\log_{10}(\lambda_{SAR})$')
p.ylabel(r'probability distribution function')
#p.legend(frameon=False, loc=3)
p.yscale('log')
p.xlim((30., 35.5))
p.ylim((1e-4, 4))
p.title('Specific accretion rate, ' +
        str(n.round(zmin, 2)) +
        r"<z<" +
        str(n.round(zmax, 2)))
p.grid()
p.savefig(os.path.join(XLF_dir, "LSAR_hist_" + baseName + ".png"))
p.clf()


# TABULATE LOG N LOG S, fx

def get_lognlogs_replicas(fx, area):
    log_f_05_20 = n.log10(fx[fx > 0])
    out = n.histogram(log_f_05_20, bins=n.arange(-18, -8., 0.2))
    # cumulative number density per square degrees
    x_out = 0.5 * (out[1][1:] + out[1][:-1])
    N_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0]))])
    c_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0]))]) / area
    c_out_up = (1 + N_out**(-0.5)) * c_out
    c_out_low = (1 - N_out**(-0.5)) * c_out
    c_err = (n.log10(c_out_up) - n.log10(c_out_low)) / 2.
    return x_out, c_out, c_err


Xgal = (abs(g_lat) > 20)
out = get_lognlogs_replicas(fx[Xgal], area=area)
DATA_X = n.transpose(out)[(out[1] > 0) & (out[2] != n.inf)]
n.savetxt(
    os.path.join(
        logNlogS_dir,
        'logNlogS_soft_' +
        baseName +
        '.ascii'),
    DATA_X)


# TABULATE LOG N LOG R, mag_r


def get_lognlogr(mag_r, area):
    # cumulative number density per square degrees
    dlogf = 0.1
    fbins = n.arange(8, 40 + 2 * dlogf, dlogf)
    #xf = fbins[:-1]+dlogf/2.
    out = n.histogram(mag_r, bins=fbins)
    x_out = 0.5 * (out[1][1:] + out[1][:-1])
    # n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])
    N_out = n.cumsum(out[0])
    # n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ]) / area
    c_out = N_out / area
    c_out_up = (1 + N_out**(-0.5)) * c_out
    c_out_low = (1 - N_out**(-0.5)) * c_out
    c_err = (n.log10(c_out_up) - n.log10(c_out_low)) / 2.
    return x_out, c_out, c_err


out = get_lognlogr(mag_r[Xgal], area=area)
DATA_R = n.transpose(out)[(out[1] > 0) & (out[2] != n.inf)]
n.savetxt(
    os.path.join(
        logNlogR_dir,
        'logNlogR_optical_' +
        baseName +
        '.ascii'),
    DATA_R)


sys.exit()

print('opens galaxy file ', time.time() - t0)
f3 = h5py.File(path_2_galaxy_file, 'r')
mass = f3['/galaxy/SMHMR_mass'][:]  # log of the stellar mass
f3.close()

f_duty = interp1d(
    n.array([0., 0.75, 2., 3.5, 6.1]),
    ##n.array([0.1, 0.2, 0.25, 0.25, 0.25])
    n.array([0.1, 0.2, 0.3, 0.3, 0.3])
)

dlogM = 0.1
bins = n.arange(8, 13, dlogM)
x_SMF = (bins[1:] + bins[:-1]) * 0.5
N_gal = n.histogram(mass, bins=bins)[0]

p.figure(2, (6, 6))
p.axes([0.16, 0.15, 0.8, 0.8])

p.axhline(f_duty(z_mean[0]), ls='dashed', lw=2)
p.axhline(f_duty(z_mean[0])**2., ls='dotted', lw=2)

N_agn_a = n.histogram(logm, bins=bins)[0]
y = N_agn_a * 1. / N_gal  # *dc_val
yerr = y * N_agn_a**(-0.5)  # *dc_val
p.errorbar(x_SMF, y, yerr=yerr, color='grey', label='all')

tsel = (lx > 41)
N_agn_41 = n.histogram(logm[tsel], bins=bins)[0]
y = N_agn_41 * 1. / N_gal  # *dc_val
yerr = y * N_agn_41**(-0.5)  # *dc_val
p.errorbar(x_SMF, y, yerr=yerr, color='black', label=r'L$_X>10^{41}$')

tsel = (lx > 42)
N_agn_42 = n.histogram(logm[tsel], bins=bins)[0]
y = N_agn_42 * 1. / N_gal  # *dc_val
yerr = y * N_agn_42**(-0.5)  # *dc_val
p.errorbar(x_SMF, y, yerr=yerr, color='red', label=r'L$_X>10^{42}$')

tsel = (lx > 43)
N_agn_43 = n.histogram(logm[tsel], bins=bins)[0]
y = N_agn_43 * 1. / N_gal  # *dc_val
yerr = y * N_agn_43**(-0.5)  # *dc_val
p.errorbar(x_SMF, y, yerr=yerr, color='blue', label=r'L$_X>10^{43}$')

tsel = (lx > 44)
N_agn_44 = n.histogram(logm[tsel], bins=bins)[0]
y = N_agn_44 * 1. / N_gal  # *dc_val
yerr = y * N_agn_44**(-0.5)  # *dc_val
p.errorbar(x_SMF, y, yerr=yerr, color='magenta', label=r'L$_X>10^{44}$')

p.xlabel(r'$\log_{10}(M^*/M_\odot)$')
p.ylabel(r'$f_{AGN}(M^*, ' + str(n.round(zmin, 2)) +
         r"<z<" + str(n.round(zmax, 2)) + r')$')
p.yscale('log')
p.ylim((5e-5, 0.4))
p.xlim((9.5, 12.))
p.grid()
p.title('Duty cycle')  # , '+str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)))
#p.legend( loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
p.legend(frameon=False, loc=0)
#p.legend( loc=8, ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
p.savefig(os.path.join(fig_dir, "duty_cycle_AGN_" + baseName + ".png"))
p.clf()


p.figure(1, (6, 6))
p.axes([0.17, 0.15, 0.73, 0.73])

p.plot(x_SMF, N_gal / (vol * dlogM), color='green', label='all galaxies')
p.plot(x_SMF, N_agn_a / (vol * dlogM), color='grey', label='all AGN')
p.plot(x_SMF, N_agn_41 / (vol * dlogM), color='black', label=r'L$_X>10^{41}$')
p.plot(x_SMF, N_agn_42 / (vol * dlogM), color='red', label=r'L$_X>10^{42}$')
p.plot(x_SMF, N_agn_43 / (vol * dlogM), color='blue', label=r'L$_X>10^{43}$')
p.plot(x_SMF, N_agn_44 / (vol * dlogM),
       color='magenta', label=r'L$_X>10^{44}$')

p.xlabel(r'$\log_{10}(M^*/[M_\odot])$')
p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
p.legend(frameon=False, loc=0, fontsize=12)
p.yscale('log')
p.xlim((9, 12.5))
p.ylim((1e-8, 1e-1))
p.title('Stellar mass function, ' +
        str(n.round(zmin, 2)) +
        "<z<" +
        str(n.round(zmax, 2)))
p.grid()
p.savefig(os.path.join(fig_dir, "SMF_AGN_" + baseName + ".png"))
p.clf()


p.figure(1, (6, 6))
p.axes([0.17, 0.15, 0.73, 0.73])

p.plot(logm, lx, 'k,', label='all AGNs')

p.xlabel(r'$\log_{10}(M^*/[M_\odot])$')
p.ylabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
p.legend(frameon=False, loc=0, fontsize=12)
# p.yscale('log')
# p.xlim((9,12.5))
# p.ylim((1e-8,1e-1))
p.title('LX-stellar mass, ' +
        str(n.round(zmin, 2)) +
        "<z<" +
        str(n.round(zmax, 2)))
p.grid()
p.savefig(os.path.join(fig_dir, "LX_MASS_" + baseName + ".png"))
p.clf()
