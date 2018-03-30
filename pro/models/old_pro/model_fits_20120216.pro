

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
;;  The three different models...
;;


;---------------------------------------------------------------
;
; P. L. E. 
; 

;;
;; Croom et al (2009) best-fit PLE params
;; e.g. Table 2. Row 4
;;
    alpha = -3.33
     beta = -1.41
 Mstar0_g = -22.17
       k1 = 1.46
       k2 = -0.328
 Phi_star = 10^(-5.84)

openw, 10, 'Croom09b_PLE_temp.dat'
printf, 10, '# '
printf, 10, '# alpha = -3.33,  beta = -1.41, k1 = 1.46,  k2==0.328, Phi_star = 10^(-5.84)'
printf, 10, '# ii, z[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom '
printf, 10, '# '

for ii=0L, zbins-1 do begin
   Hz = H0 * sqrt(omega_m*(1.0+zz[ii])^3.0 + omega_lam)
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)
      ;; => if mag_bins = 60, -32.00 < M_g < -17.25

      dmag   = mag - Mstar_g 
      ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

      print, 'A', ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 10, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom, $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8,e18.8)'
   endfor
;   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
;   phi  /= hub**3  ;; Convert from Mpc to Mpc/h volumes.
;   print, ii, zz[ii], 
endfor
Phi_h = Phi / (little_h^3)
print

close, 10




;---------------------------------------------------------------
; 
;  L  D  D  E
;   Luminosity-dependent density evolution 
;   Section 6.2 from Croom et al. (2009b)

   alpha =  -3.70
    beta =  -2.34
M_star_g = -26.69 
  M_g_c  = -23.90
  gamma  =   0.68 
zc_zero  =   2.47
      p1 =   6.28
      p2 =  -2.85
       A = 10^(-9.21)

e_d = fltarr(zbins, mag_bins)

openw, 11, 'Croom09b_mLDDE_temp.dat'

;
; "Modified Luminosity-dependent density evolution" (mLDDE)
;  where the modification is to the behaviour of the 
;  "e_d" parameter... 
;
for ii=0L, zbins-1 do begin
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)      

      dmag    = mag - M_star_g 
      ;; N.B. This is a different M_star_g than the Mstar_g from
      ;; above!!! Grrr.....
      dmag_c  = mag - M_g_c 

      ;; Equnation (18) of Croom et al. (2009b)
      zc = zc_zero / ( 1. + (10^(0.4*gamma*dmag_c)))

      ;; Equnation (17) of Croom et al. (2009b)
      top_ed    = 2.*(1+zc)^(p1)
      bottom_ed = (((1+zz[ii])/(1+zc))^(-1.*p1)) + (((1+zz[ii])/(1+zc))^(-1.*p2))
      e_d(ii,jj) = top_ed / bottom_ed 

;      Abs_Mag_model[jj] = mag

      ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = (A * e_d[ii,jj])  / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

      print, ii, jj, zz[ii], mag, top_ed, bottom_ed, A, e_d[ii, jj],  (A * e_d[ii,jj]), ple_denom
;      print,      ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj],  $ 
 ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.6)'

      printf, 11, ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj], $
              format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
   endfor
endfor
close, 11



;---------------------------------------------------------------
;  
;  "Q  D  E"
;  Quadratic Density Evolution
;   Section 6.3 from Croom et al. (2009b)


;---------------------------------------------------------------
;
; "Modified PLE"
;   Section 6.4, first...

;;   alpha_z    = alpha_ref * (( 1.+ zz[ii] / 1.+ z_ref)^(p_alpha)
;;   phi_star_z = phi_star_ref * (( 1.+ zz[ii] / 1.+ z_ref)^(p_phi)


;  
;  L E D E
;   Luminosity Evolution + Density Evolution 
; Table 4. 

    alpha_ref = -3.48
      p_alpha =  0.220
         beta = -1.38
     Mstar0_g = -22.24
           k1 =   1.23
           k2 =  -0.206
       k_phi1 =   0.430
       k_phi2 =   1.139
phi_star_zero = 10^(-5.79)

;; "where z_ref is fixed at z_ref =2.0"
z_ref = 2.0

openw, 12, 'Croom09b_LEDE_temp.dat'
for ii=0L, zbins-1 do begin
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   
   alpha_z    = alpha_ref * (( 1.+ zz[ii] / 1.+ z_ref)^(p_alpha))
;   phi_star_z = phi_star_ref * (( 1.+ zz[ii] / 1.+ z_ref)^(p_phi)
   
   log_phi_star = alog10(phi_star_zero) + ( k_phi1 * zz[ii]*(1.0-(0.5*zz[ii])/k_phi2))
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)      
      
      dmag   = mag - Mstar_g 
      Phi[ii, jj] = 10^(log_phi_star) / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      
;      print, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 12, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8)'
   endfor
endfor
close, 12
close, /all


end
