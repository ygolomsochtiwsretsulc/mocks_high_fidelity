"""
Adds in the list of snapshots the distances of the borders between spherical shells and the number of times the box needs to be replicated along each axis for a shell.

"""
import numpy as n
import os
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.interpolate import interp1d

l_box = float(sys.argv[1]) # 1000.
env = sys.argv[2] # MD10
print('runs geometry.py with arguments ')
print(sys.argv)

if env == 'TEST_DIR' or env == 'MD04' or env == 'MD10' or env == 'MD40':
    cosmo = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Om0=0.307115)
    h=0.6777
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR" or env == "UNIT_fA2i_DIR":
    cosmo = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h=0.6774

# deduce the L_box without h
L_box = l_box / h

# read shell list
A_snap = n.loadtxt(
    os.path.join(
        os.environ[env], 'snap_list.txt'),
    unpack=True)
Z_snap = 1. / A_snap - 1.
N_snap = n.arange(len(A_snap))

z_array = n.arange(0, 20., 0.001)
dcs = cosmo.comoving_distance(z_array)
dc_to_z = interp1d(dcs, z_array)
z_to_dc = interp1d(z_array, dcs)

DC_snap = z_to_dc(Z_snap)
DC_borders = 0.5 * (DC_snap[:-1] + DC_snap[1:])

DC_max = n.hstack((DC_snap[0], DC_borders))
DC_min = n.hstack((DC_borders, DC_snap[-1]))
N_replicas = n.ceil(DC_max / L_box)
n.savetxt(
    os.path.join(
        os.environ[env], 'snap_list_with_border.txt'),
    n.transpose(
        [
            N_snap,
            Z_snap,
            A_snap,
            DC_max,
            DC_min,
            N_replicas]),
    header=' N_snap redshift aexp DC_min DC_max N_replicas')

print('finished')