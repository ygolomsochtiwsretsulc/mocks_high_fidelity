"""
What it does
------------

Computes and tabulates the logNlogS, logNlogR, XLF

References
----------

Command to run
--------------

python3 003_2_agn_compute_XLF_logNlogS_R.py environmentVAR ftyp

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
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
# import all pathes

env = sys.argv[1]  # 'MD04'
ftyp = sys.argv[2]  # "sat"
print(env, ftyp)


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

if env == "MD10" or env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    z_maximum = 1.4
if env == "MD04":
    z_maximum = 0.45

DZ = 0.2
all_zs = n.arange(0., z_maximum, DZ)

root_dir = os.path.join(os.environ[env], 'cat_AGN_' + ftyp)
agn_catalog_list = sorted(
    n.array(
        glob.glob(
            os.path.join(
                os.environ[env],
                'cat_AGN_' +
                ftyp,
                '*.fit'))))
galaxy_catalog_list = sorted(n.array(
    glob.glob(
        os.path.join(
            os.environ[env],
            'cat_GALAXY_' +
            ftyp,
            '*.fit'))))

fig_dir = os.path.join(root_dir, 'figures')
logNlogS_dir = os.path.join(root_dir, 'logNlogS')
logNlogR_dir = os.path.join(root_dir, 'logNlogR')
XLF_dir = os.path.join(root_dir, 'XLF')
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)
if os.path.isdir(logNlogS_dir) == False:
    os.system('mkdir -p ' + logNlogS_dir)
if os.path.isdir(logNlogR_dir) == False:
    os.system('mkdir -p ' + logNlogR_dir)
if os.path.isdir(XLF_dir) == False:
    os.system('mkdir -p ' + XLF_dir)

area = healpy.nside2pixarea(8, degrees=True)

for path_2_eRO_catalog, path_2_GAL_catalog in zip(
        agn_catalog_list, galaxy_catalog_list):
    str_healpix_id = os.path.basename(path_2_eRO_catalog)[:-4]
    hd = fits.open(path_2_eRO_catalog)
    hg = fits.open(path_2_GAL_catalog)
    #In [24]: hd[1].data.columns
    # Out[24]:
    # ColDefs(
    #name = 'ra'; format = 'D'; unit = 'degree'
    #name = 'dec'; format = 'D'; unit = 'degree'
    #name = 'g_lat'; format = 'D'; unit = 'degree'
    #name = 'g_lon'; format = 'D'; unit = 'degree'
    #name = 'ecl_lat'; format = 'D'; unit = 'degree'
    #name = 'ecl_lon'; format = 'D'; unit = 'degree'
    #name = 'redshift_R'; format = 'D'; unit = 'real space'
    #name = 'redshift_S'; format = 'D'; unit = 'redshift space'
    #name = 'dL_cm'; format = 'D'; unit = 'cm'
    #name = 'galactic_NH'; format = 'D'; unit = 'cm'
    #name = 'galactic_ebv'; format = 'D'; unit = 'cm'
    #name = 'galaxy_stellar_mass'; format = 'D'; unit = 'M_sun'
    #name = 'galaxy_star_formation_rate'; format = 'D'; unit = 'M_sun'
    #name = 'galaxy_LX_hard'; format = 'D'; unit = 'M_sun'
    #name = 'AGN_LX_soft'; format = 'D'; unit = 'Luminosity/[erg/s] 0.5-2 keV'
    #name = 'AGN_FX_soft'; format = 'D'; unit = 'Flux/[erg/cm2/s] 0.5-2 keV'
    #name = 'AGN_LX_hard'; format = 'D'; unit = 'Luminosity/[erg/s] 2-10 keV'
    #name = 'AGN_FX_hard'; format = 'D'; unit = 'Flux/[erg/cm2/s] 2-10 keV'
    #name = 'AGN_SDSS_r_magnitude'; format = 'D'; unit = 'mag'
    #name = 'AGN_Nh'; format = 'D'; unit = 'log10(Nh/[cm-2])'
    #name = 'AGN_random_number'; format = 'D'
    #name = 'AGN_type'; format = 'D'; unit = 'X/opt type: 11, 12, 21, 22'
    #name = 'HALO_M500c'; format = 'D'; unit = 'log10(M/[M_sun])'
    #name = 'HALO_Mvir'; format = 'D'; unit = 'log10(M/[M_sun])'
    #name = 'HALO_halo_id'; format = 'K'
    #name = 'HALO_pid'; format = 'K'
    # )
    #
    # hg[1].data.columns
    #
    # ColDefs(
    #name = 'ra'; format = 'D'; unit = 'degree'
    #name = 'dec'; format = 'D'; unit = 'degree'
    #name = 'g_lat'; format = 'D'; unit = 'degree'
    #name = 'g_lon'; format = 'D'; unit = 'degree'
    #name = 'ecl_lat'; format = 'D'; unit = 'degree'
    #name = 'ecl_lon'; format = 'D'; unit = 'degree'
    #name = 'redshift_R'; format = 'D'; unit = 'real space'
    #name = 'redshift_S'; format = 'D'; unit = 'redshift space'
    #name = 'dL_cm'; format = 'D'; unit = 'cm'
    #name = 'galactic_NH'; format = 'D'; unit = 'cm-2'
    #name = 'galactic_ebv'; format = 'D'; unit = 'mag'
    #name = 'galaxy_stellar_mass'; format = 'D'; unit = 'log10(M/[M_sun])'
    #name = 'galaxy_star_formation_rate'; format = 'D'; unit = 'log10(SFR/[M_sun/year])'
    #name = 'galaxy_LX_hard'; format = 'D'; unit = 'log10(LX (2-10keV)/[erg/s])'
    #name = 'HALO_M500c'; format = 'D'; unit = 'log10(M/[M_sun])'
    #name = 'HALO_Mvir'; format = 'D'; unit = 'log10(M/[M_sun])'
    #name = 'HALO_halo_id'; format = 'K'
    #name = 'HALO_pid'; format = 'K'
    # )

    z = hd[1].data['redshift_R']
    lx = hd[1].data['AGN_LX_hard']
    logNH = hd[1].data['AGN_Nh']
    fx = hd[1].data['AGN_FX_soft']
    lx_0520 = hd[1].data['AGN_LX_soft']
    logm = hd[1].data['galaxy_stellar_mass']
    lsar = lx - logm
    mag_r = hd[1].data['AGN_SDSS_r_magnitude']
    g_lat = hd[1].data['g_lat']

    logm_gal = hg[1].data['galaxy_stellar_mass']
    z_gal = hg[1].data['redshift_R']

    n_agn = len(z)
    indexes = n.arange(n_agn)

    # TABULATE LOG N LOG S, fx

    def get_lognlogs_replicas(fx, area):
        log_f_05_20 = n.log10(fx[fx > 0])
        out = n.histogram(log_f_05_20, bins=n.arange(-20, -8., 0.2))
        # cumulative number density per square degrees
        #x_out = 0.5*(out[1][1:] + out[1][:-1])
        # n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])
        N_out = n.cumsum(out[0][::-1])[::-1]
        # n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ]) / area
        c_out = N_out / area
        c_out_up = (1 + N_out**(-0.5)) * c_out
        c_out_low = (1 - N_out**(-0.5)) * c_out
        c_err = (n.log10(c_out_up) - n.log10(c_out_low)) / 2.
        return N_out  # , c_err

    #Xgal = (abs(g_lat)>20)
    N_out = get_lognlogs_replicas(fx, area=area)
    # DATA_X = n.transpose(out)#[(out[1]>0) & (out[2]!=n.inf)]
    n.savetxt(
        os.path.join(
            logNlogS_dir,
            'logNlogS_soft_' +
            str_healpix_id +
            '.ascii'),
        N_out)

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
        return N_out  # , c_out#, c_err

    N_out = get_lognlogr(mag_r, area=area)
    # DATA_R = n.transpose(out)#[(out[1]>0) & (out[2]!=n.inf)]
    n.savetxt(
        os.path.join(
            logNlogR_dir,
            'logNlogR_optical_' +
            str_healpix_id +
            '.ascii'),
        N_out)

    for id_z in n.arange(len(all_zs)):
        zmin = all_zs[id_z]
        zmax = all_zs[id_z] + DZ
        z_mean = 0.5 * (zmin + zmax)
        z_selection = (z >= zmin) & (z < zmax)
        z_selection_GAL = (z_gal >= zmin) & (z_gal < zmax)
        baseName = str_healpix_id + '_' + \
            str(n.round(zmin, 2)) + '_z_' + str(n.round(zmax, 2))

        print(zmin, '<z<', zmax)
        vol = (cosmo.comoving_volume(zmax).value -
               cosmo.comoving_volume(zmin).value) * n.pi * area / 129600.
        DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
        print('volume', vol, 'Mpc3')

        #dataDir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data')
        #xlf_dir = os.path.join(dataDir, 'LF_SMF/LXFunction/aird_2015/')

        # LF binning
        dlogf = 0.05
        Lbin_min = 36
        fbins = n.arange(Lbin_min, 48, dlogf)
        xf = fbins[:-1] + dlogf / 2.

        #LXMIN = get_LX_min(lc_id)
        #pta = lambda name : os.path.join(xlf_dir, name)
        #z0, z1, z_center, Lxmin, Lxmax, Lx_c,   Nobj, phi, phierr = n.loadtxt(os.path.join(os.environ['GIT_VS'], 'data/LF_SMF', 'LXFunction', 'miyaji_2015.ascii'), unpack=True)
        #z0, z1, Lxmin, Lxmax, Nobj, Nmdl, phi_a = n.loadtxt(os.path.join(os.environ['GIT_VS'], 'data/LF_SMF', 'LXFunction', 'hasinger_2005_soft_XLF.ascii'), unpack=True)
        #z_center = ( z0 + z1 ) * 0.5
        #Lx_c = ( Lxmin + Lxmax ) * 0.5 - 2*n.log10(0.7)
        #phi = phi_a * 0.7**3.
        #x_aj12, y_aj12 = n.loadtxt(os.path.join(dataDir, 'LF_SMF', 'LXFunction', 'ajello_12.csv'), delimiter=',', unpack=True)

        # hard X-ray luminosity function Aird 2015)
        def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))
        def Ls_h(z): return 10**(44.84 - n.log10(((1 + 2.0) / (1 + z))
                                                 ** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))
        def phi_h(L, z): return kz_h(z) / \
            ((L / Ls_h(z))**0.48 + (L / Ls_h(z))**2.27)

        # Selections for the histogram
        NH20 = (z_selection) & (logNH < 22)
        NH22 = (z_selection) & (logNH >= 22) & (logNH < 24)
        NH24 = (z_selection) & (logNH >= 24)

        N_nh20 = n.histogram(lx[NH20], fbins)[0] / vol / dlogf
        N_nh22 = n.histogram(lx[NH22], fbins)[0] / vol / dlogf
        N_nh24 = n.histogram(lx[NH24], fbins)[0] / vol / dlogf

        nhar = n.histogram(lx[z_selection], fbins)[0] / vol / dlogf
        nharN = n.histogram(lx[z_selection], fbins)[0]  # /vol/dlogf

        p.figure(1, (6, 6))
        p.axes([0.2, 0.18, 0.75, 0.75])
        #sel = ((zmin+zmax)*0.5>z0)&((zmin+zmax)*0.5<z1)
        # if len(sel.nonzero()[0])>0:
        # if len(list(set(z_center[sel])))>=2:
        #sel = (sel) & (z_center==z_center[sel][0])
        #p.plot(Lx_c[sel], phi[sel], label='Ha05 '+str(n.round(z_center[sel][0],2)), ls='dashed')

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

        #DATA_XLF = n.transpose([xf, nhar, (nhar)*(1-(nharN)**(-0.5)), (nhar)*(1+(nharN)**(-0.5)), N_nh20, N_nh22, N_nh24, phi_h(10**xf,z_mean)])
        #n.savetxt(os.path.join(XLF_dir, 'XLF_soft_'+baseName+'.ascii'), DATA_XLF  )

        # XLF_ratio_
        p.figure(1, (6, 6))
        p.axes([0.18, 0.18, 0.75, 0.75])
        p.plot(xf, nhar / phi_h(10**xf, z_mean), label='2-10 keV')
        p.fill_between(xf,
                       y1=1 - (nharN)**(-0.5),
                       y2=1 + (nharN)**(-0.5),
                       color='green',
                       alpha=0.3,
                       label='mock ERR')
        p.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
        p.ylabel(r'mock/model')
        p.legend(frameon=False, loc=0)
        p.xlim((37., 46.5))
        p.ylim((0.7, 1.3))
        p.title(str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
        p.grid()
        p.savefig(os.path.join(XLF_dir, "XLF_ratio_" + baseName + ".png"))
        p.clf()

        # LSAR histogram
        dlogf = 0.1
        fbins = n.arange(30, 38, dlogf)
        xf = n.arange(30, 38, dlogf)[:-1] + dlogf / 2.

        nall = n.histogram(lsar[z_selection], fbins)[0] / vol / dlogf
        nallN = n.histogram(lsar[z_selection], fbins)[0]

        zsel = (logm >= 12) & (z_selection)
        nall_12 = n.histogram(lsar[zsel], fbins)[0] / vol / dlogf
        nallN_12 = n.histogram(lsar[zsel], fbins)[0]

        zsel = (logm >= 11) & (logm < 12) & (z_selection)
        nall_11 = n.histogram(lsar[zsel], fbins)[0] / vol / dlogf
        nallN_11 = n.histogram(lsar[zsel], fbins)[0]

        zsel = (logm >= 10) & (logm < 11) & (z_selection)
        nall_10 = n.histogram(lsar[zsel], fbins)[0] / vol / dlogf
        nallN_10 = n.histogram(lsar[zsel], fbins)[0]

        zsel = (logm >= 9) & (logm < 10) & (z_selection)
        nall_9 = n.histogram(lsar[zsel], fbins)[0] / vol / dlogf
        nallN_9 = n.histogram(lsar[zsel], fbins)[0]

        p.figure(1, (6, 6))
        p.axes([0.16, 0.15, 0.8, 0.8])

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
                str(n.round(zmin, 2)) + r"<z<" + str(n.round(zmax, 2)))
        p.grid()
        p.savefig(os.path.join(XLF_dir, "LSAR_hist_" + baseName + ".png"))
        p.clf()

        f_duty = interp1d(
            n.array([0., 0.75, 2., 3.5, 6.1]),
            ##n.array([0.1, 0.2, 0.25, 0.25, 0.25])
            n.array([0.1, 0.2, 0.3, 0.3, 0.3])
        )

        dlogM = 0.1
        bins = n.arange(8, 13, dlogM)
        x_SMF = (bins[1:] + bins[:-1]) * 0.5
        z_selection_gal = (z_gal >= zmin) & (z_gal < zmax)
        N_gal = n.histogram(logm_gal[z_selection_gal], bins=bins)[0]

        p.figure(2, (6, 6))
        p.axes([0.16, 0.15, 0.8, 0.8])

        p.axhline(f_duty(z_mean[0]), ls='dashed', lw=2)
        p.axhline(f_duty(z_mean[0])**2., ls='dotted', lw=2)

        N_agn_a = n.histogram(logm, bins=bins)[0]
        y = N_agn_a * 1. / N_gal  # *dc_val
        yerr = y * N_agn_a**(-0.5)  # *dc_val
        p.errorbar(x_SMF, y, yerr=yerr, color='grey', label='all')

        tsel = (lx > 41) & (z_selection)
        N_agn_41 = n.histogram(logm[tsel], bins=bins)[0]
        y = N_agn_41 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_41**(-0.5)  # *dc_val
        p.errorbar(x_SMF, y, yerr=yerr, color='black', label=r'L$_X>10^{41}$')

        tsel = (lx > 42) & (z_selection)
        N_agn_42 = n.histogram(logm[tsel], bins=bins)[0]
        y = N_agn_42 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_42**(-0.5)  # *dc_val
        p.errorbar(x_SMF, y, yerr=yerr, color='red', label=r'L$_X>10^{42}$')

        tsel = (lx > 43) & (z_selection)
        N_agn_43 = n.histogram(logm[tsel], bins=bins)[0]
        y = N_agn_43 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_43**(-0.5)  # *dc_val
        p.errorbar(x_SMF, y, yerr=yerr, color='blue', label=r'L$_X>10^{43}$')

        tsel = (lx > 44) & (z_selection)
        N_agn_44 = n.histogram(logm[tsel], bins=bins)[0]
        y = N_agn_44 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_44**(-0.5)  # *dc_val
        p.errorbar(
            x_SMF,
            y,
            yerr=yerr,
            color='magenta',
            label=r'L$_X>10^{44}$')

        p.xlabel(r'$\log_{10}(M^*/M_\odot)$')
        p.ylabel(r'$f_{AGN}(M^*, ' + str(n.round(zmin, 2)) +
                 r"<z<" + str(n.round(zmax, 2)) + r')$')
        p.yscale('log')
        p.ylim((5e-5, 0.4))
        p.xlim((9.5, 12.))
        p.grid()
        # , '+str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)))
        p.title('Duty cycle')
        #p.legend( loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
        p.legend(frameon=False, loc=0)
        #p.legend( loc=8, ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
        p.savefig(os.path.join(fig_dir, "duty_cycle_AGN_" + baseName + ".png"))
        p.clf()

        p.figure(1, (6, 6))
        p.axes([0.17, 0.15, 0.73, 0.73])

        p.plot(x_SMF, N_gal / (vol * dlogM),
               color='green', label='all galaxies')
        p.plot(x_SMF, N_agn_a / (vol * dlogM), color='grey', label='all AGN')
        p.plot(x_SMF, N_agn_41 / (vol * dlogM),
               color='black', label=r'L$_X>10^{41}$')
        p.plot(x_SMF, N_agn_42 / (vol * dlogM),
               color='red', label=r'L$_X>10^{42}$')
        p.plot(x_SMF, N_agn_43 / (vol * dlogM),
               color='blue', label=r'L$_X>10^{43}$')
        p.plot(x_SMF, N_agn_44 / (vol * dlogM),
               color='magenta', label=r'L$_X>10^{44}$')

        p.xlabel(r'$\log_{10}(M^*/[M_\odot])$')
        p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
        p.legend(frameon=False, loc=0, fontsize=12)
        p.yscale('log')
        p.xlim((9, 12.5))
        p.ylim((1e-8, 1e-1))
        p.title('Stellar mass function, ' + str(n.round(zmin, 2)) +
                "<z<" + str(n.round(zmax, 2)))
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
        p.title('LX-stellar mass, ' + str(n.round(zmin, 2)) +
                "<z<" + str(n.round(zmax, 2)))
        p.grid()
        p.savefig(os.path.join(fig_dir, "LX_MASS_" + baseName + ".png"))
        p.clf()
