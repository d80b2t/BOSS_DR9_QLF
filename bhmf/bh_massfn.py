#!/usr/bin/env python
#
# Numerically integrate the BH mass function forward in time, loosely
# following Shankar et al. (2009; ApJ, 690, 20)
#

import numpy as N
import sys
import os


OmM = 0.274
OmL = 0.726
hub = 0.700


def Ez(zz):
    """
    Ez(zz):
    The dimensionless Hubble parameter.
    """
    tmp = N.sqrt( OmM*(1+zz)**3+OmL )
    return(tmp)
    #



def Lbol2Mbj(lgLbol):
    """
    Lbol2Lbj(lgLbol):
    Convert from (log) bolometric luminosities in Watts to bJ magnitudes.
    Uses Shen++09 bolometric calibration and Richards++06 conversion.
    """
    lgL  = lgLbol + 7.		# Convert from W to erg/s.
    Miz2 = -2.5*lgL + 90.	# The i-band magnitude, k-corrected to z=2.
    Mbj  = Miz2 + 0.71		# Finally convert to bJ.
    return(Mbj)
    #




def dPhi_C09(mag, zz):
    """
    dPhi_C09(mag, zz):
    The derivative of the QSO luminosity function at redshift zz, i.e.
    d^2n/dM^2, with magnitudes in the bJ system.
    Volumes are per Mpc/h (cubed).
    Parameters for the LF fit from Croom++04, as modified by Croton09.
    This agrees with the C04 for z<3.
    """
    phi0,m0,b=1.67e-6,-21.61,-1.09
    if zz<3.0:
        k1,k2,a = 1.39,-0.29,-3.31
    else:
        k1,k2,a = 1.22,-0.23,-3.31+0.5*(zz-3.0)
    mstar = m0 - 2.5*(k1*zz + k2*zz*zz);
    dm    = mag - mstar;
    phi   = 1.0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
    dphi  = phi*phi*(0.4*N.log(10.0)*(a+1)*10.0**(0.4*(a+1)*dm)+\
                     0.4*N.log(10.0)*(b+1)*10.0**(0.4*(b+1)*dm));
    #print ""
    #print dphi, phi0, hub
    dphi *= phi0/hub**3  # Convert from Mpc to Mpc/h volumes.
    #print dphi, phi0, hub
    return(dphi)
    #






def evolve(f01=0.7,lam=0.50,zobs=0.0):
    """
    evolve(f01=0.7,lam=0.50,zobs=0.0):
    f01 is f=epsilon/(1-epsilon) divided by 0.1 where epsilon is the
    accretion efficiency (L=epsilon dot{M}_inflow c^2).
    lam is L/L_{Edd}=dot{m}f01, assumed constant.
    """
    ell = 1.26e31	# W/Msun, L_{Edd}/M_BH.
    mdot= lam/f01	# Mass Eddington ratio: L/Ledd=mdot.f01
    tS  = 4.5e7		# Salpeter time for f=0.1, in yr.
    #
    # Do a stupid brute-force integration of the continuity equation
    # in the form:
    #  dn/dt = -[dot{m}/t_s/M_BH] [2.5/ln(10)]^2  dPhi(Mag,z)/dMag
    # where Phi is n(x).x.ln(10) for x=Lbol or M_BH per dex (not ln)
    # or the equivalent in magnitudes.
    #
    # Start with drho/dt+div[rho.v]=0, but we're flowing in mass not
    # space so we have dn/dt+d[n Mdot]/dM.  Write Mdot/M=mdot/tS and
    # convert d/dM to (1/M)d/dlnM.  Note
    #   n_M(M)dM = N_L(L)dL ==> M.n_M = L.n_L
    # for L\propto M.  Convert from ln to lg to get the above.
    # We actually compute our LFs per magnitude, so there is an additional
    # factor of 2.5 converting lgL to magnitude.
    #
    mbh = N.logspace(6.0,10.0,25)
    nbh = N.zeros_like(mbh)
    zz  = 10.00 + zobs
    dz  =-0.001
    # Do the evolution.
    while zz>zobs:
        lgL  = N.log10(f01*mdot*ell*mbh)
        Mbj  = Lbol2Mbj(lgL)
        dPhi = (2.5/N.log(10))**2 * dPhi_C09(Mbj,zz) # d^2n/dlnL^2
        dndt = -mdot/tS/mbh * dPhi
        dtdz = -9.7778e9/hub/Ez(zz)/(1+zz)
        nbh += dndt*dtdz*dz
        zz  += dz
    # Convert to mass function.
    PhiM = nbh*mbh*N.log(10)	# Number density per dex per (Mpc/h)^3.
    return( (mbh,PhiM) )
    #


if __name__=="__main__":
    if len(sys.argv) != 2:
        raise RuntimeError,"Usage: bh_massfn.py <zobs>"
    else:
        zobs   = float(sys.argv[1])
        M,PhiM = evolve(zobs=zobs)
        print "# BH mass function from LF assuming constant Eddington"
        print "# ratio and constant efficiency accretion."
        print "# Output is for z=%.2f."%zobs
        print "# Masses are in Msun, volumes (Mpc/h)^3, mass function per dex."
        print "# %10s %12s"%("MBH","Phi")
        for i in range(0,PhiM.size):
            print "%12.4e %12.4e"%(M[i],PhiM[i])
