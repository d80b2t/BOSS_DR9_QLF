

;; No. of mag bins
mag_bins = 60

;; No. of redshift bins
zbins = 100. 
;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
zz = findgen(zbins)/20.    ;;  Redshift bins..

;; Set up phi...
Phi = fltarr(zbins, mag_bins)


omega_m   = 0.30
omega_lam = 0.70
H0        = 70. 
little_h  = H0/100.


;;
;;  The two different models...
;;


;---------------------------------------------------------------
;
; P. L. E.   ``Full''
; 

;; PLE FULL

Mstar0_g        = -22.17  ;; not actually needed...
z_ref           =  2.00   ;; fixed
log_phi_star    = -4.825  ;; +/- 0.060 in Mpc^-3
log_L_star_zero = 13.036  ;; +/- 0.045 in L_star = 3.9 x10^33 erg s^-1
k_L_1           =  0.632  ;; +/- 0.077
k_L_2           = -11.76  ;; +/- 0.38
k_L_3           = -14.25  ;; +/- 0.80
gamma_one_zero  =  0.417  ;; 0.055
k_gamma_one     = -0.623  ;; 0.132 
gamma_two_zero  =  2.174  ;; 0.055
k_gamma_two_one =  1.460  ;; 0.096
k_gamma_two_two = -0.793  ;; 0.057

alpha = (-1.)*(gamma_one_zero +1.)
beta  = (-1.)*(gamma_two_zero +1.)
phi_star_prime       = log_phi_star / (alog(10))
phi_star_prime_prime = 0.4 * log_phi_star 

log_L_star = fltarr(zbins)


openw, 10, 'HRH07_PLE_full_temp.dat'
printf, 10, '# '
printf, 10, '# '
printf, 10, '# ii, z[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom '
printf, 10, '# '
for ii=0L, zbins-1 do begin
   Hz = H0 * sqrt(omega_m*(1.0+zz[ii])^3.0 + omega_lam)
   Mstar_g  = Mstar0_g  - 2.5*(k_L_1*zz[ii] + k_L_2*zz[ii]*zz[ii])

   xi = alog10( (1+zz[ii])/(1.+z_ref))
   log_L_star[ii] = log_L_star_zero + (k_L_1*xi) + (k_L_2*xi*xi) + (k_L_3*xi*xi*xi) 
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)
      ;; => if mag_bins = 60, -32.00 < M_g < -17.25



      dmag   = mag - Mstar_g 
      ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = log_Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

;      print, 'A', ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
 ;     printf, 10, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom, $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8,e18.8)'
   endfor
;   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
;   phi  /= hub**3  ;; Convert from Mpc to Mpc/h volumes.
;   print, ii, zz[ii], 
endfor
close, 10
close, /all



end

;;

 
