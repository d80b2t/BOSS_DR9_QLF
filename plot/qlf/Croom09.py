def lumfn_C09(mag, zz):
   """
   lumfn_C09(mag, zz):
   The QSO luminosity function at redshift zz: dn/dM with magnitudes
   in the bJ system.
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
   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
   print phi
   
   #phi  /= hub**3  # Convert from Mpc to Mpc/h volumes.
   return(phi)
   #
