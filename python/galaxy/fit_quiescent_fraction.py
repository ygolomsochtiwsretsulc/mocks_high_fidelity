"""
Adds the empirical sdss r band magnitude using the empirical relations derived from COSMOS

python3 /home/comparat/software/lss_mock_dev/python/galaxy/fit_quiescent_fraction.py

"""
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
import matplotlib.pyplot as p
import matplotlib
from scipy.stats import norm
import linmix
from scipy.special import erf
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import linear_model
from sklearn.model_selection import cross_val_predict
from sklearn import datasets
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import sys
from scipy.interpolate import interp1d
import numpy as n
import glob
import os
print('Fit quiescent fraction on COSMOS')
print('------------------------------------------------')
print('------------------------------------------------')
cosmoMD = FlatLambdaCDM(
    H0=67.77 * u.km / u.s / u.Mpc,
    Om0=0.307115)  # , Ob0=0.048206)

plotDir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'figures',
    'galaxies',
    'SFR')


def DL_z(redshift): return cosmoMD.luminosity_distance(redshift).to(u.cm).value


def get_lx(z, fx): return n.log10(
    4. * n.pi * (DL_z(z)**2) * 0.6777 * 0.6777 * fx)


matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})


# IMPORTING COSMOS, Marchesi 2016

data_file = os.path.join(
    os.environ['OBS_REPO'],
    'COSMOS',
    'photoz_vers2.0_010312.fits')
# J_ApJ_817_34_catalog.dat.gz.fits.gz
# J_ApJ_817_34_ccosmos.dat.fits.gz
# "zbest"
#
# "R"
#Clph ==1, 2, 3


dd = fits.open(data_file)[1].data
U = dd['U']
V = dd['V']
J = dd['J']
R = dd['R']
K = dd['K']
NUV = dd['NUV']
photoz = dd['photoz']
z_sel = (
    photoz > 0) & (
        photoz < 1.3) & (
            dd['ssfr_med'] > -
            30) & (
                dd['logMS'] > 6) & (
                    R < 24.5)

qu_w12 = (z_sel) & (U > 0) & (V > 0) & (J > 0) & (
    U - V > 0.8 * (V - J) + 0.7) & (U - V > 1.3) & (V - J < 1.5)
sf_w12 = (z_sel) & (U > 0) & (V > 0) & (J > 0) & (qu_w12 == False)

qu_a13 = (z_sel) & (NUV > 0) & (R > 0) & (K > 0) & (
    (NUV - R) > 3.75) & ((NUV - R) > 1.37 * (R - K) + 3.2)
sf_a13 = (z_sel) & (NUV > 0) & (R > 0) & (K > 0) & (qu_a13 == False)

QU = (z_sel) & (dd['ssfr_med'] > -30) & (dd['ssfr_med'] < -11)
SF = (z_sel) & (dd['ssfr_med'] > -11)


def nl(sel): return len(sel.nonzero()[0])
# print(nl(qu_w12))
# print(nl(SF_w12))
# print(nl(qu_a13))
# print(nl(SF_a13))


# dd['logMS']
# dd['logMS']
# dd['sfr_med']

p.figure(figsize=(6, 6))
p.plot(dd['photoz'][SF], dd['logMS'][SF], 'b+', rasterized=True,
       alpha=0.5, label='SF I13, N=' + str(nl(SF)))
p.plot(dd['photoz'][QU], dd['logMS'][QU], 'rx', rasterized=True,
       alpha=0.5, label='QU I13, N=' + str(nl(QU)))
p.xlabel(r'$z$')
p.ylabel(r'$M_\odot$')
p.legend(frameon=False, loc=3)
p.grid()
p.title('I13')  # + r", $F_X>1\times10^{-17}$ ")
p.savefig(os.path.join(plotDir, "mass-redshift-SSFRm11cut.png"))
p.clf()

m_bins = n.arange(7, 12.5, 0.1)
x_m_bins = m_bins[:-1] + 0.1 / 2.
DZ = 0.2
z_vals = n.arange(0.2, 1.4, DZ)
N_SF = n.zeros((len(z_vals), len(m_bins) - 1))
N_QU = n.zeros((len(z_vals), len(m_bins) - 1))
N_zz = n.zeros((len(z_vals), len(m_bins) - 1))
N_mm = n.zeros((len(z_vals), len(m_bins) - 1))

p.figure(2, (6, 6))
for jj, z_i in enumerate(z_vals):
    z_bin = (photoz > z_i) & (photoz < z_i + DZ)
    N_SF[jj] = p.hist(dd['logMS'][SF & z_bin],
                      bins=m_bins,
                      histtype='step',
                      rasterized=True,
                      alpha=0.5,
                      label=str(n.round(z_i,
                                        1)) + ' SF I13, N=' + str(nl(SF & z_bin)))[0]
    N_QU[jj] = p.hist(dd['logMS'][QU & z_bin],
                      bins=m_bins,
                      histtype='step',
                      rasterized=True,
                      alpha=0.5,
                      label=str(n.round(z_i,
                                        1)) + ' QU I13, N=' + str(nl(QU & z_bin)))[0]
    N_zz[jj][:] = z_i + DZ / 2.
    N_mm[jj][:] = x_m_bins

p.ylabel('Counts')
p.yscale('log')
p.xlabel(r'$M_\odot$')
p.legend(frameon=False, loc=3)
p.grid()
p.title('I13')  # + r", $F_X>1\times10^{-17}$ ")
p.savefig(os.path.join(plotDir, "mass-redshift-SSFRm11cut_hist.png"))
p.clf()


def fun(mass, M0, width): return 0.5 + 0.5 * erf((mass - M0) / width)


M0_s, widths_s = n.meshgrid(
    n.arange(
        9.7, 12., 0.01), n.arange(
            0.3, 1.8, 0.005))
all_Ms = n.hstack((M0_s))
all_Ws = n.hstack((widths_s))


def curve_lib(mass): return n.array(
    [fun(mass, M0_i, width_i) for M0_i, width_i in zip(all_Ms, all_Ws)])


chis_all = []

for jj, z_i in enumerate(z_vals):
    p.figure(0, (6, 6))
    N_total = N_SF[jj] + N_QU[jj]
    y_data_i = n.log10(1. * N_QU[jj] / N_total)
    x_data_i, x_err_i = x_m_bins, n.ones_like(x_m_bins) * 0.05
    y_err_val_i = N_QU[jj]**(-0.5) * y_data_i
    y_err_i = 0.434 * N_QU[jj]**(-0.5)
    ok = (
        n.isnan(y_data_i) == False) & (
        n.isinf(y_data_i) == False) & (
            y_err_i > 0)
    x_data, y_data, x_err, y_err = x_data_i[ok], y_data_i[ok], x_err_i[ok], y_err_i[ok]
    y_err_val = y_err_val_i[ok]
    # linear fitting
    #lm = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
    #lm.run_mcmc(miniter=200, maxiter=600)
    #ys = n.median(lm.chain['alpha']) + x_m_bins * n.median(lm.chain['beta'])
    #p.plot(x_m_bins, ys, ls='dotted', label=r'SFR=Ms'+str(n.round(n.median(lm.chain['beta']),2))+'+'+str(n.round(n.median(lm.chain['alpha']),2)), lw=3)
    print(x_data)
    print(10**y_data)
    chi2s = n.sum(abs(10**y_data - curve_lib(x_data)), axis=1)
    chis_all.append(chi2s)
    idx_best = n.argmin(chi2s)
    M0_best, width_best = all_Ms[idx_best], all_Ws[idx_best]
    p.errorbar(x_data, 10**y_data, xerr=x_err, yerr=y_err_val,
               label=str(n.round(z_i, 1)) + '<z<' + str(n.round(z_i + 0.3, 1)))
    #pt2, ct2 = curve_fit(fun, x_data, 10**y_data, p0=(11.,0.6), sigma=y_err_val)
    p.plot(x_data, fun(x_data, M0_best, width_best), label='erf(' +
           str(n.round(M0_best, 2)) +
           ', ' +
           str(n.round(width_best, 2)) +
           ')', ls='dashed')
    print(jj, z_i)
    print(M0_best)
    print(width_best)
    print('==========================================')
    p.ylabel('Counts')
    # p.yscale('log')
    p.ylim((0, 1))
    p.xlim((6.5, 12.0))
    p.xlabel(r'$M_\odot$')
    p.ylabel('Quiescent fraction')
    p.legend(frameon=False, loc=0)
    p.grid()
    # + r", $F_X>1\times10^{-17}$ ")
    p.title(str(n.round(z_i, 1)) + '<z<' + str(n.round(z_i + DZ, 1)))
    p.savefig(
        os.path.join(
            plotDir,
            "mass-redshift-SSFRm11cut_hist_fractionQU_z_" +
            str(
                n.round(
                    z_i +
                    DZ /
                    2.,
                    1)) +
            ".png"))
    p.clf()
    # figure chi2 vs parameters
    p.figure(1, (6, 6))
    p.scatter(all_Ms, all_Ws, c=chi2s, s=4, marker='s',
              edgecolors='face')  # , vmin=0, vmax=3)
    p.colorbar(label='chi2')
    p.plot(all_Ms[n.argmin(chi2s)], all_Ws[n.argmin(chi2s)], 'k+')
    p.xlabel(r'$M_0$')
    p.ylabel(r'$\sigma$')
    # p.ylim((10,28))
    # p.xlim((-16,-11))
    #p.legend(frameon=False, loc=0)
    p.grid()
    # + r", $F_X>1\times10^{-17}$ ")
    p.title(str(n.round(z_i, 1)) + '<z<' + str(n.round(z_i + DZ, 1)))
    p.tight_layout()
    p.savefig(os.path.join(plotDir, "bayesian_M0_sigma_z_" +
                           str(n.round(z_i + DZ / 2., 1)) + ".png"))
    p.clf()
    n.savetxt(os.path.join(plotDir,
                           "mass-redshift-SSFRm11cut_hist_fractionQU_z_" + str(n.round(z_i + DZ / 2.,
                                                                                       1)) + ".txt"),
              n.transpose([x_data,
                           10**y_data,
                           x_err,
                           y_err_val,
                           z_i * n.ones_like(x_data),
                           (z_i + DZ) * n.ones_like(x_data),
                           fun(x_data,
                               M0_best,
                               width_best)]))
    n.savetxt(os.path.join(plotDir,
                           "bayesian_M0_sigma_z_" + str(n.round(z_i + DZ / 2.,
                                                                1)) + ".txt"),
              n.transpose([all_Ms,
                           all_Ws,
                           chi2s,
                           z_i * n.ones_like(chi2s),
                           (z_i + DZ) * n.ones_like(chi2s)]))

data_out_list = n.array(
    glob.glob(
        os.path.join(
            plotDir,
            "mass-redshift-SSFRm11cut_hist_fractionQU_z_*.txt")))
param_out_list = n.array(
    glob.glob(
        os.path.join(
            plotDir,
            "bayesian_M0_sigma_z_*.txt")))
data_out_list.sort()
param_out_list.sort()


p.figure(2, figsize=(6, 6))
for filename in data_out_list:
    str_z = os.path.basename(filename).split('_')[-1][:-4]
    x_data, y_data, x_err, y_err_val, z_min, z_max, y_model = n.loadtxt(
        filename, unpack=True)
    p.errorbar(
        x_data,
        y_data,
        xerr=x_err,
        yerr=y_err_val,
        label='data z=' +
        str_z,
        alpha=0.5)
    p.plot(x_data, y_model, label='model z=' + str_z, ls='solid', lw=3)

p.ylabel('Counts')
p.ylim((0, 1))
p.xlim((6.5, 12.0))
p.xlabel(r'$M_\odot$')
p.ylabel('Quiescent fraction')
p.legend(frameon=False, loc=0)
p.grid()
p.savefig(
    os.path.join(
        plotDir,
        "summary_mass-redshift-SSFRm11cut_hist_fractionQU.png"))
p.clf()

# figure parameter vs parameter vs parameters
p.figure(3, (6, 6))
for filename in param_out_list:
    str_z = os.path.basename(filename).split('_')[-1][:-4]
    all_Ms, all_Ws, chi2s, z_min, z_max = n.loadtxt(filename, unpack=True)
    ok = (chi2s < n.min(chi2s) * 1.2)
    p.plot(all_Ms[ok], all_Ws[ok], label=str_z, ls='', marker='+')

p.xlabel(r'$M_0$')
p.ylabel(r'$\sigma$')
p.legend(frameon=False, loc=0)
p.grid()
p.tight_layout()
p.savefig(os.path.join(plotDir, "summary_bayesian_M0_sigma.png"))
p.clf()

# figure parameters vs redshift
p.figure(4, (6, 6))
x_fit = []
y_fit = []
for filename in param_out_list:
    str_z = os.path.basename(filename).split('_')[-1][:-4]
    all_Ms, all_Ws, chi2s, z_min, z_max = n.loadtxt(filename, unpack=True)
    ok = (chi2s < n.min(chi2s) * 1.2)
    p.plot(1 + (z_min[ok] + z_max[ok]) * 0.5,
           all_Ws[ok], label=str_z, ls='', marker='+')
    x_fit.append(1 + (z_min[ok] + z_max[ok]) * 0.5)
    y_fit.append(all_Ws[ok])

x_fit = n.hstack((x_fit))
y_fit = n.hstack((y_fit))
p_out = n.polyfit(x_fit, y_fit, deg=1)
p.plot(x_fit, n.polyval(p_out, x_fit), 'k--', label=r'$\sigma=z$' +
       str(n.round(p_out[0], 2)) + '+' + str(n.round(p_out[1], 2)))

p.xlabel(r'$1+z$')
p.ylabel(r'$\sigma$')
p.legend(frameon=False, loc=0)
p.grid()
p.tight_layout()
p.savefig(os.path.join(plotDir, "summary_bayesian_redshift_sigma.png"))
p.clf()

# figure parameters vs redshift
p.figure(5, (6, 6))
x_fit = []
y_fit = []
for filename in param_out_list:
    str_z = os.path.basename(filename).split('_')[-1][:-4]
    all_Ms, all_Ws, chi2s, z_min, z_max = n.loadtxt(filename, unpack=True)
    ok = (chi2s < n.min(chi2s) * 1.2)
    p.plot(1 + (z_min[ok] + z_max[ok]) * 0.5,
           all_Ms[ok], label=str_z, ls='', marker='+')
    x_fit.append(1 + (z_min[ok] + z_max[ok]) * 0.5)
    y_fit.append(all_Ms[ok])

x_fit = n.hstack((x_fit))
y_fit = n.hstack((y_fit))
p_out = n.polyfit(x_fit, y_fit, deg=1)
p.plot(x_fit, n.polyval(p_out, x_fit), 'k--', label=r'$M_0=z$' +
       str(n.round(p_out[0], 2)) + '+' + str(n.round(p_out[1], 2)))

p.xlabel(r'$1+z$')
p.ylabel(r'$\log_{10}(M_0)$')
p.legend(frameon=False, loc=0)
p.grid()
p.tight_layout()
p.savefig(os.path.join(plotDir, "summary_bayesian_redshift_M0.png"))
p.clf()

# sys.exit()

# p.figure(figsize=(6,6))
#p.plot(dd['photoz'][sf_w12], dd['logMS'][sf_w12], 'b+', rasterized=True, alpha=0.5, label='SF W12, N='+str(nl(sf_w12)))
#p.plot(dd['photoz'][qu_w12], dd['logMS'][qu_w12], 'rx', rasterized=True, alpha=0.5, label='QU W12, N='+str(nl(qu_w12)))
# p.xlabel(r'$z$')
# p.ylabel(r'$M_\odot$')
#p.legend(frameon=False, loc=3)
# p.grid()
# p.title('W12')#+ r", $F_X>1\times10^{-17}$ ")
#p.savefig(os.path.join( plotDir, "mass-redshift-W12.png"))
# p.clf()


# p.figure(figsize=(6,6))
#p.plot(dd['photoz'][sf_a13], dd['logMS'][sf_a13], 'b+', rasterized=True, alpha=0.5, label='SF A13, N='+str(nl(sf_a13)))
#p.plot(dd['photoz'][qu_a13], dd['logMS'][qu_a13], 'rx', rasterized=True, alpha=0.5, label='QU A13, N='+str(nl(qu_a13)))
# p.xlabel(r'$z$')
# p.ylabel(r'$M_\odot$')
#p.legend(frameon=False, loc=3)
# p.grid()
# p.title('A13')#+ r", $F_X>1\times10^{-17}$ ")
#p.savefig(os.path.join( plotDir, "mass-redshift-A13.png"))
# p.clf()

# Now plot MASS -- SFR for the 2 population
# add main sequence equation from W12
# equation for quiescent population

def mean_SFR(zz, mass): return (0.70 - 0.13 * zz) * \
    (mass - 10.5) + 0.38 + 1.14 * zz - 0.19 * zz**2.


def fun(x, loc, scale): return norm.cdf(x, loc=loc, scale=scale)


n.random.seed(2)

for jj, z_i in enumerate(z_vals):
    p.figure(6, (6, 6))
    z_bin = (photoz > z_i) & (photoz < z_i + DZ)
    # SF
    p.plot(dd['logMS'][SF & z_bin], dd['sfr_med'][SF & z_bin], ls='', marker='+',
           rasterized=True, alpha=0.5, label=str(n.round(z_i + DZ / 2., 1)) + ' SF')
    x_data, y_data, x_err, y_err = dd['logMS'][SF & z_bin], dd['sfr_med'][SF & z_bin], (
        dd['mass_sup'][SF & z_bin] - dd['mass_inf'][SF & z_bin]) / 2., (dd['sfr_sup'][SF & z_bin] + dd['sfr_inf'][SF & z_bin]) / 2.
    lm = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
    lm.run_mcmc(miniter=200, maxiter=600)
    ys = n.median(lm.chain['alpha']) + x_m_bins * n.median(lm.chain['beta'])
    #p.plot(x_m_bins, ys, ls='dotted', label=r'SFR=Ms'+str(n.round(n.median(lm.chain['beta']),2)) +'+'+str(n.round(n.median(lm.chain['alpha']),2)), lw=3)
    # QU
    p.plot(dd['logMS'][QU & z_bin], dd['sfr_med'][QU & z_bin], ls='', marker='x',
           rasterized=True, alpha=0.5, label=str(n.round(z_i + DZ / 2., 1)) + ' QU')
    x_data_Q, y_data_Q, x_err_Q, y_err_Q = dd['logMS'][QU & z_bin], dd['sfr_med'][QU & z_bin], (
        dd['mass_sup'][QU & z_bin] - dd['mass_inf'][QU & z_bin]) / 2., (dd['sfr_sup'][QU & z_bin] + dd['sfr_inf'][QU & z_bin]) / 2.
    lmQ = linmix.LinMix(x_data_Q, y_data_Q, x_err_Q, y_err_Q, K=2)
    lmQ.run_mcmc(miniter=200, maxiter=600)
    ys = n.median(lmQ.chain['alpha']) + x_m_bins * n.median(lmQ.chain['beta'])
    p.plot(x_m_bins, ys, ls='dashed', label=r'SFR=Ms' +
           str(n.round(n.median(lmQ.chain['beta']), 2)) +
           '+' +
           str(n.round(n.median(lmQ.chain['alpha']), 2)), lw=3)
    # W12
    SFR = mean_SFR((z_i + DZ / 2.) * n.ones_like(x_m_bins), x_m_bins)
    p.plot(x_m_bins, SFR, ls='dashed', label='W12', lw=2)
    p.ylabel(r'$log10(SFR)$')
    p.xlabel(r'$log10(M_\odot)$')
    p.legend(frameon=False, loc=3)
    p.grid()
    p.xlim((6.5, 12.0))
    p.title(str(n.round(z_i, 1)) + '<z<' + str(n.round(z_i + DZ, 1)))
    p.savefig(os.path.join(plotDir, "mass-SFR-SSFRm11cut_z_" +
                           str(n.round(z_i + DZ / 2., 1)) + ".png"))
    p.clf()

    n.savetxt(os.path.join(plotDir, "mass-SFR-z_" + str(n.round(z_i + DZ / 2., 1)) +
                           ".txt"), n.transpose([n.median(lmQ.chain['beta']), n.median(lmQ.chain['alpha'])]))

    # data
    y_model_W12 = mean_SFR(dd['logMS'][SF & z_bin], dd['photoz'][SF & z_bin])
    dy_W12 = y_data - y_model_W12
    y_model_SF = n.median(lm.chain['alpha']) + \
        x_data * n.median(lm.chain['beta'])
    dy_SF = y_data - y_model_SF
    y_model_Q = n.median(lmQ.chain['alpha']) + \
        x_data_Q * n.median(lmQ.chain['beta'])
    dy_Q = y_data_Q - y_model_Q

    # figure
    bins = n.arange(-3, 3, 0.1)
    xbins = n.arange(-3, 3, 0.1)[1:] - 0.05
    p.figure(7, (6, 6))
    llxs = n.arange(38, 45, 0.1)
    n_t2 = p.hist(
        dy_W12,
        bins=bins,
        normed=True,
        cumulative=True,
        histtype='step',
        label='W12')[0]
    pt2, ct2 = curve_fit(fun, xbins, n_t2, p0=(0, 1))
    p.plot(xbins,
           norm.cdf(xbins,
                    loc=pt2[0],
                    scale=pt2[1]),
           label='W12 N(' + str(n.round(pt2[0],
                                        2)) + ', ' + str(n.round(pt2[1],
                                                                 2)) + ')',
           ls='dashed')

    n_t2 = p.hist(
        dy_SF,
        bins=bins,
        normed=True,
        cumulative=True,
        histtype='step',
        label='SF')[0]
    pt2, ct2 = curve_fit(fun, xbins, n_t2, p0=(0, 1))
    p.plot(xbins,
           norm.cdf(xbins,
                    loc=pt2[0],
                    scale=pt2[1]),
           label='SF N(' + str(n.round(pt2[0],
                                       2)) + ', ' + str(n.round(pt2[1],
                                                                2)) + ')',
           ls='dashed')

    n_t2 = p.hist(
        dy_Q,
        bins=bins,
        normed=True,
        cumulative=True,
        histtype='step',
        label='QU')[0]
    pt2, ct2 = curve_fit(fun, xbins, n_t2, p0=(0, 1))
    p.plot(xbins,
           norm.cdf(xbins,
                    loc=pt2[0],
                    scale=pt2[1]),
           label='QU N(' + str(n.round(pt2[0],
                                       2)) + ', ' + str(n.round(pt2[1],
                                                                2)) + ')',
           ls='dashed')

    p.ylabel('counts')
    p.xlabel(r'delta mag')
    p.grid()
    p.legend(frameon=False, loc=0)
    p.title(str(n.round(z_i, 1)) + '<z<' + str(n.round(z_i + DZ, 1)))
    p.savefig(os.path.join(plotDir, "residual_mass-SFR-SSFRm11cut_z_" +
                           str(n.round(z_i + DZ / 2., 1)) + ".png"))
    p.clf()

    n.savetxt(os.path.join(plotDir, "mass-SFR-z-residuals_" + \
              str(n.round(z_i + DZ / 2., 1)) + ".txt"), n.transpose([pt2[0], pt2[1]]))

    # figure
    # p.figure(figsize=(10,5))
    #p.plot(fx_mock, r_mock, 'yx', label='mock', rasterized=True, alpha=0.1)
    #p.plot(n.log10(fx_S), r_S, 'ko', rasterized=True, alpha=0.25, label='Stripe82X, N='+str(len(r_S)))
    #p.plot(n.log10(fx_C), r_C, 'cs', rasterized=True, alpha=0.25, label='CDFS, N='+str(len(r_C)))
    #p.plot(n.log10(fx_R), r_R, 'm^', rasterized=True, alpha=0.25, label='2RXS, N='+str(len(r_R)))
    #p.plot(n.log10(fx_M), r_M, 'r+', rasterized=True, alpha=0.25, label='COSMOS, N='+str(len(r_M)))

    #xs = n.arange(-17,-10,0.02)
    # for beta_val, alpha_val, redshift in zip(beta_s, alpha_s, (zmins+zmaxs)/2. ):
    #ys = beta_val+ xs *alpha_val
    #print(beta_val, alpha_val, redshift)
    #p.scatter(xs, ys, c=n.ones_like(xs)*redshift, s=2, marker='o', edgecolors='face', vmin=0, vmax=3)

    # p.colorbar(label='redshift')
    ##p.plot(xs, -2.*xs-6.5, color='k', lw=3, ls='dashed', label=r'$r_{AB}=-6.5-2\times\log_{10}(F_X)$')
    #p.plot(xs, -2.*xs-7., color='k', lw=3, ls='dashed', label=r'$r_{AB}=-7-2\times\log_{10}(F_X)$')
    #p.xlabel(r'$\log_{10}(F_X/$[erg cm$^{-2}$ s$^{-1}$]$)$')
    # p.ylabel(r'$r_{AB}$')
    # p.ylim((10,28))
    # p.xlim((-16,-11))
    #p.legend(frameon=False, loc=3)
    # p.grid()
    # p.title(title_summary)
    # p.tight_layout()
    #p.savefig(os.path.join( plotDir, "bayesian_mag_fx_all_fits.png"))
    # p.clf()


beta, alpha, loc, scale = n.zeros_like(z_vals), n.zeros_like(
    z_vals), n.zeros_like(z_vals), n.zeros_like(z_vals)
for jj, z_i in enumerate(z_vals):
    beta[jj], alpha[jj] = n.loadtxt(os.path.join(
        plotDir, "mass-SFR-z_" + str(n.round(z_i + DZ / 2., 1)) + ".txt"), unpack=True)
    loc[jj], scale[jj] = n.loadtxt(os.path.join(
        plotDir, "mass-SFR-z-residuals_" + str(n.round(z_i + DZ / 2., 1)) + ".txt"), unpack=True)


def fun(x, a, b): return a * x + b


p.figure(8, (6, 6))
z_data = z_vals + DZ / 2.
p.plot(z_data, beta, label='beta')
pt2b, ct2b = curve_fit(fun, z_data, beta, p0=(0, 1))
p.plot(z_data,
       fun(z_data,
           pt2b[0],
           pt2b[1]),
       label='betaM,' + str(n.round(pt2b[0],
                                    2)) + ', ' + str(n.round(pt2b[1],
                                                             2)),
       ls='dashed')

p.plot(z_data, alpha + 10, label='alpha')
pt2a, ct2a = curve_fit(fun, z_data, alpha, p0=(0, 1))
p.plot(z_data, 10 +
       fun(z_data, pt2a[0], pt2a[1]), label='alphaM,' +
       str(n.round(pt2a[0], 2)) +
       ', ' +
       str(n.round(pt2a[1], 2)), ls='dashed')

p.plot(z_data, loc, label='loc')
pt2l, ct2l = curve_fit(fun, z_data, loc, p0=(0, 1))
p.plot(z_data,
       fun(z_data,
           pt2l[0],
           pt2l[1]),
       label='locM,' + str(n.round(pt2l[0],
                                   2)) + ', ' + str(n.round(pt2l[1],
                                                            2)),
       ls='dashed')

p.plot(z_data, scale, label='scale')
pt2s, ct2s = curve_fit(fun, z_data, scale, p0=(0, 1))
p.plot(z_data,
       fun(z_data,
           pt2s[0],
           pt2s[1]),
       label='scaleM,' + str(n.round(pt2s[0],
                                     2)) + ', ' + str(n.round(pt2s[1],
                                                              2)),
       ls='dashed')
p.ylabel('parameter')
p.xlabel(r'redshift')
p.grid()
p.legend(frameon=False, loc=0)
# p.title(str(n.round(z_i,1))+'<z<'+str(n.round(z_i+DZ,1)))
p.savefig(os.path.join(plotDir, "mass-SFR-z-model.png"))
p.clf()

sys.exit()

# TREATMENT OF THE RESIDUALS
# data
y_model_B = n.median(lm.chain['alpha']) + x_data * n.median(lm.chain['beta'])
dy_B = y_data - y_model_B
# fit


def fun(x, loc, scale): return norm.cdf(x, loc=loc, scale=scale)


# figure
bins = n.arange(-3, 3, 0.1)
xbins = n.arange(-3, 3, 0.1)[1:] - 0.05
p.figure(2, (6, 6))
llxs = n.arange(38, 45, 0.1)
n_t2 = p.hist(
    dy_B,
    bins=bins,
    normed=True,
    cumulative=True,
    histtype='step',
    label='residuals')[0]
pt2, ct2 = curve_fit(fun, xbins, n_t2, p0=(0, 1))
p.plot(xbins,
       norm.cdf(xbins,
                loc=pt2[0],
                scale=pt2[1]),
       label='Gaussian N(' + str(n.round(pt2[0],
                                         2)) + ', ' + str(n.round(pt2[1],
                                                                  2)) + ')',
       ls='dashed')
p.ylabel('counts')
p.xlabel(r'delta mag')
p.grid()
p.legend(frameon=False, loc=0)
p.title(str(n.round(zmin, 2)) + r'$<z<$' + str(n.round(zmax, 2)))
p.savefig(os.path.join(plotDir, "hist_residual_mag_fx_" +
                       str(n.round(zbar, 2)) + ".png"))
p.clf()
residuals = [pt2, ct2]
# return popt, pcov, n.median(lm.chain['beta']),
# n.median(lm.chain['alpha']), n.std(lm.chain['beta']),
# n.std(lm.chain['alpha']), residuals

a_s = n.zeros_like(zmins)
b_s = n.zeros_like(zmins)
a_e = n.zeros_like(zmins)
b_e = n.zeros_like(zmins)

alpha_s = n.zeros_like(zmins)
beta_s = n.zeros_like(zmins)
alpha_e = n.zeros_like(zmins)
beta_e = n.zeros_like(zmins)

a_t2_s = n.zeros_like(zmins)
b_t2_s = n.zeros_like(zmins)
a_t2_e = n.zeros_like(zmins)
b_t2_e = n.zeros_like(zmins)


for index in n.arange(len(zmins)):
    pi, c, alpha, beta, alpha_err, beta_err, residuals = get_params(
        zmins[index], zmaxs[index])
    a_s[index] = pi[0]
    b_s[index] = pi[1]
    a_e[index] = c[0][0]**(-0.5)
    b_e[index] = c[1][1]**(-0.5)
    alpha_s[index] = alpha
    beta_s[index] = beta
    alpha_e[index] = alpha_err
    beta_e[index] = beta_err
    try:
        pt2, ct2 = residuals
        a_t2_s[index] = pt2[0]
        b_t2_s[index] = pt2[1]
        a_t2_e[index] = ct2[0][0]**(0.5)
        b_t2_e[index] = ct2[1][1]**(0.5)
    except(ValueError):
        print('error')

n.savetxt(os.path.join(plotDir,
                       "parameters.ascii"),
          n.transpose([zmins,
                       zmaxs,
                       a_s,
                       b_s,
                       a_e,
                       b_e,
                       alpha_s,
                       beta_s,
                       alpha_e,
                       beta_e,
                       a_t2_s,
                       b_t2_s,
                       a_t2_e,
                       b_t2_e]))

zmins, zmaxs, a_s, b_s, a_e, b_e, alpha_s, beta_s, alpha_e, beta_e, a_t2_s, b_t2_s, a_t2_e, b_t2_e = n.loadtxt(
    os.path.join(plotDir, "parameters.ascii"), unpack=True)

zmins, zmaxs, a_s, b_s, a_e, b_e, alpha_s, beta_s, alpha_e, beta_e, a_t2_s, b_t2_s, a_t2_e, b_t2_e = n.transpose(
    [zmins, zmaxs, a_s, b_s, a_e, b_e, alpha_s, beta_s, alpha_e, beta_e, a_t2_s, b_t2_s, a_t2_e, b_t2_e])[zmins > 0.1].T

# alpha vs. redshift
p.figure(1, (6, 6))
p.errorbar(
    (zmins + zmaxs) / 2.,
    alpha_s,
    yerr=alpha_e,
    xerr=(
        zmaxs - zmins) / 2.,
    color='g',
    label="data")
n.random.seed(2)
x_data = (zmins + zmaxs) / 2.
y_data = alpha_s
x_err = (zmaxs - zmins) / 2.
y_err = alpha_e
lm_1 = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
lm_1.run_mcmc(miniter=200, maxiter=2000)
for i in range(0, len(lm_1.chain), 5):
    xs = n.arange(
        n.min(x_data),
        n.max(x_data),
        (n.max(x_data) - n.min(x_data)) / 20.)
    ys = lm_1.chain[i]['alpha'] + xs * lm_1.chain[i]['beta']
    p.plot(xs, ys, color='g', alpha=0.02)
ys = n.median(lm_1.chain['alpha']) + xs * n.median(lm_1.chain['beta'])
p.plot(xs, ys, color='k', ls='dashed', lw=2, label=r'$y=x\times$' +
       str(n.round(n.median(lm_1.chain['beta']), 2)) +
       '+' +
       str(n.round(n.median(lm_1.chain['alpha']), 2)))
p.xlabel(r'z')
p.ylabel(r"$\alpha$")
p.grid()
# p.ylim((-3,-1))
p.legend(frameon=False, loc=0)
p.title(option)
p.savefig(os.path.join(plotDir, "mag_llx_a_b.png"))
p.clf()

# beta vs. redshift
p.figure(2, (6, 6))
p.errorbar(
    (zmins + zmaxs) / 2.,
    beta_s,
    yerr=beta_e,
    xerr=(
        zmaxs - zmins) / 2.,
    color='g',
    label="data")
n.random.seed(2)
x_data = (zmins + zmaxs) / 2.
y_data = beta_s
x_err = (zmaxs - zmins) / 2.
y_err = beta_e
lm_1 = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
lm_1.run_mcmc(miniter=200, maxiter=2000)
for i in range(0, len(lm_1.chain), 5):
    xs = n.arange(
        n.min(x_data),
        n.max(x_data),
        (n.max(x_data) - n.min(x_data)) / 20.)
    ys = lm_1.chain[i]['alpha'] + xs * lm_1.chain[i]['beta']
    p.plot(xs, ys, color='g', alpha=0.02)
ys = n.median(lm_1.chain['alpha']) + xs * n.median(lm_1.chain['beta'])
p.plot(xs, ys, color='k', ls='dashed', lw=2, label=r'$y=x\times$' +
       str(n.round(n.median(lm_1.chain['beta']), 2)) +
       '+' +
       str(n.round(n.median(lm_1.chain['alpha']), 2)))
p.xlabel(r'z')
p.ylabel(r"$\beta$")
p.grid()
# p.ylim((-20,0))
p.legend(frameon=False, loc=0)
p.title(option)
p.savefig(os.path.join(plotDir, "mag_llx_beta.png"))
p.clf()

# mean of the residuals vs. redshift
p.figure(3, (6, 6))
p.errorbar(
    (zmins + zmaxs) / 2.,
    a_t2_s,
    yerr=a_t2_e,
    xerr=(
        zmaxs - zmins) / 2.,
    color='g',
    label='data')
n.random.seed(2)
x_data = (zmins + zmaxs) / 2.
y_data = a_t2_s
x_err = (zmaxs - zmins) / 2.
y_err = a_t2_e
lm_1 = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
lm_1.run_mcmc(miniter=1000, maxiter=2000)
for i in range(0, len(lm_1.chain), 5):
    xs = n.arange(
        n.min(x_data),
        n.max(x_data),
        (n.max(x_data) - n.min(x_data)) / 20.)
    ys = lm_1.chain[i]['alpha'] + xs * lm_1.chain[i]['beta']
    p.plot(xs, ys, color='g', alpha=0.02)

ys = n.median(lm_1.chain['alpha']) + xs * n.median(lm_1.chain['beta'])
p.plot(xs, ys, color='k', ls='dashed', lw=2, label=r'$y=x\times$' +
       str(n.round(n.median(lm_1.chain['beta']), 2)) +
       '+' +
       str(n.round(n.median(lm_1.chain['alpha']), 2)))
p.xlabel(r'z')
p.ylabel(r'mean residual')
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(plotDir, "residuals_loc.png"))
p.clf()

# variance of the residuals vs. redshift
p.figure(4, (6, 6))
p.errorbar(
    (zmins + zmaxs) / 2.,
    b_t2_s,
    yerr=b_t2_e,
    xerr=(
        zmaxs - zmins) / 2.,
    color='g',
    label='data')
n.random.seed(2)
x_data = (zmins + zmaxs) / 2.
y_data = b_t2_s
x_err = (zmaxs - zmins) / 2.
y_err = b_t2_e
lm_1 = linmix.LinMix(x_data, y_data, x_err, y_err, K=2)
lm_1.run_mcmc(miniter=200, maxiter=2000)
for i in range(0, len(lm_1.chain), 5):
    xs = n.arange(
        n.min(x_data),
        n.max(x_data),
        (n.max(x_data) - n.min(x_data)) / 20.)
    ys = lm_1.chain[i]['alpha'] + xs * lm_1.chain[i]['beta']
    p.plot(xs, ys, color='g', alpha=0.02)
ys = n.median(lm_1.chain['alpha']) + xs * n.median(lm_1.chain['beta'])
p.plot(xs, ys, color='k', ls='dashed', lw=2, label=r'$y=x\times$' +
       str(n.round(n.median(lm_1.chain['beta']), 2)) +
       '+' +
       str(n.round(n.median(lm_1.chain['alpha']), 2)))
p.xlabel(r'z')
p.ylabel(r'scatter in the residuals')
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(plotDir, "residuals_scale.png"))
p.clf()

if option == 't1':
    title_summary = 'type 1'
if option == 't2':
    title_summary = 'type 2'
if option == 'tall':
    title_summary = 'type 1 or 2'

# figure
p.figure(figsize=(10, 5))
p.plot(fx_mock, r_mock, 'yx', label='mock', rasterized=True, alpha=0.1)
p.plot(n.log10(fx_S), r_S, 'ko', rasterized=True,
       alpha=0.25, label='Stripe82X, N=' + str(len(r_S)))
p.plot(n.log10(fx_C), r_C, 'cs', rasterized=True,
       alpha=0.25, label='CDFS, N=' + str(len(r_C)))
p.plot(n.log10(fx_R), r_R, 'm^', rasterized=True,
       alpha=0.25, label='2RXS, N=' + str(len(r_R)))
p.plot(n.log10(fx_M), r_M, 'r+', rasterized=True,
       alpha=0.25, label='COSMOS, N=' + str(len(r_M)))

xs = n.arange(-17, -10, 0.02)
for beta_val, alpha_val, redshift in zip(
        beta_s, alpha_s, (zmins + zmaxs) / 2.):
    ys = beta_val + xs * alpha_val
    print(beta_val, alpha_val, redshift)
    p.scatter(
        xs,
        ys,
        c=n.ones_like(xs) *
        redshift,
        s=2,
        marker='o',
        edgecolors='face',
        vmin=0,
        vmax=3)

p.colorbar(label='redshift')
#p.plot(xs, -2.*xs-6.5, color='k', lw=3, ls='dashed', label=r'$r_{AB}=-6.5-2\times\log_{10}(F_X)$')
p.plot(xs, -2. * xs - 7., color='k', lw=3, ls='dashed',
       label=r'$r_{AB}=-7-2\times\log_{10}(F_X)$')
p.xlabel(r'$\log_{10}(F_X/$[erg cm$^{-2}$ s$^{-1}$]$)$')
p.ylabel(r'$r_{AB}$')
p.ylim((10, 28))
p.xlim((-16, -11))
p.legend(frameon=False, loc=3)
p.grid()
p.title(title_summary)
p.tight_layout()
p.savefig(os.path.join(plotDir, "bayesian_mag_fx_all_fits.png"))
p.clf()
