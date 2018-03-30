;+
;
;
;
;-


print
print

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

;;---------------------------------------------------------------
;;
;; P. L. E. 
;; 
;;
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

     ; print, 'A', ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 10, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom, $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8,e18.8)'
   endfor
;   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
;   phi  /= hub**3  ;; Convert from Mpc to Mpc/h volumes.
;   print, ii, zz[ii], 
endfor
Phi_h = Phi / (little_h^3)
close, 10
print, 'Croom09b_PLE_temp.dat made...  '
print
print




fltarr_bins = 150
;; No. of mag bins
;mag_bins = 60
mag_bins = 40

;; No. of redshift bins
;zbins = 100. 
zbins = 150

starting_redshift = 0.30

;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
;zz = findgen(zbins)/20.    ;;  Redshift bins..
zz = ((findgen(fltarr_bins))/25.)+starting_redshift


z_out = ((findgen(fltarr_bins))/25.)+starting_redshift
M_range_out   = (findgen(40)*0.30)-31.05


;; Set up phi...
;Phi = fltarr(zbins, mag_bins)
 Phi = fltarr(fltarr_bins,40)


;---------------------------------------------------------------
; 
;  (m)  L  D  D  E
;
;  (modified) Luminosity-dependent density evolution 
;   Section 6.2 from Croom et al. (2009b)
;   and Table 3

   alpha =  -3.70
    beta =  -2.34
M_star_g = -26.69 
  M_g_c  = -23.90
  gamma  =   0.68 
zc_zero  =   2.47
      p1 =   6.28
      p2 =  -2.85
       A = 10^(-9.21)


;; from *my* "best-fitting" models...
;   alpha =  -4.50
;    beta =  -2.50
;M_star_g = -27.75 
;  M_g_c  = -23.90
;  gamma  =   0.68 
;zc_zero  =   2.47
;      p1 =   5.90
 ;    p2 =   -4.20
  ;   A = 10^(-9.80)

;; Bongiorno et al. (2007) VVDS
;   alpha =  -3.29
;    beta =  -2.0
;M_star_g = -24.38 
;e  M_g_c  = -23.90
;  gamma  =   0.68 
;zc_zero  =   2.47
;      p1 =   5.90
 ;    p2 =   -4.20
  ;   A = 10^(-9.80)


e_d = fltarr(zbins, mag_bins)

openw, 11, 'Croom09b_mLDDE_temp.dat'

printf, 11, '#' 
printf, 11, '# alpha  beta  M_star_g  M_g_c  gamma  zc_zero  p1 p2  A '
printf, 11,    alpha, beta, M_star_g,  M_g_c,  gamma,  zc_zero,  p1, p2,  A ,  $
        format='(f8.3, f8.3, f8.3, f8.3, f8.3, f8.3, f8.3, f6.3, f6.3, f6.3, e14.7)' 
printf, 11, '# '
;
; "Modified Luminosity-dependent density evolution" (mLDDE)
;  where the modification is to the behaviour of the 
;  "e_d" parameter... 
;
for ii=0L, zbins-1 do begin
   for jj=0L, mag_bins-1 do begin
;      mag = -32.0 + (jj*0.25)      
      mag = M_range_out[jj]

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

      ;print, ii, jj, zz[ii], mag, top_ed, bottom_ed, A, e_d[ii, jj],  (A * e_d[ii,jj]), ple_denom
;      print,      ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj],  $ 
 ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.6)'

      printf, 11, zz[ii], mag, Phi[ii,jj];, $
;              format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
;      printf, 11, ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj], $
;              format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
   endfor
endfor
close, 11
print, 'Croom09b_mLDDE_temp.dat made...  '
print
print


;;
;; Bongiorno et al. (2007) VVDS
;;
  alpha =  -3.29
   beta =  -2.0
 M_star =   -24.38 
    M_c =  -27.36
  gamma =   0.21 
zc_zero =   2.08
     p1 =   6.54
     p2 =   -1.37
     A = 10^(-7.5544)

e_d = fltarr(zbins, mag_bins)

openw, 15, 'Bongiorno07_VVDS_LDDE_temp.dat'

printf, 15, '#' 
printf, 15, '# alpha  beta  M_star_g  M_g_c  gamma  zc_zero  p1 p2  A '
printf, 15,    alpha, beta, M_star_g,  M_g_c,  gamma,  zc_zero,  p1, p2,  A ,  $
        format='(f8.3, f8.3, f8.3, f8.3, f8.3, f8.3, f8.3, f6.3, f6.3, f6.3, e14.7)' 
printf, 15, '# '
;
; "Modified Luminosity-dependent density evolution" (mLDDE)
;  where the modification is to the behaviour of the 
;  "e_d" parameter... 
;
for ii=0L, zbins-1 do begin

   for jj=0L, mag_bins-1 do begin
;      mag = -32.0 + (jj*0.25)      
      mag = M_range_out[jj]

      dmag    = mag - M_star 
      ;; N.B. This is a different M_star_g than the Mstar_g from
      ;; above!!! Grrr.....
      dmag_c  = mag - M_c 

      ;; Equnation (11) of Bongiorno et al (2007)
      if mag ge M_c then begin
         zc = zc_zero * (10^ ( (-1.)*0.4*gamma*dmag_c) )
      endif else begin
         zc= zc_zero
      endelse

      ;; Equnation (10) of Bongiorno et al (2007)
      if zz[ii] le zc then begin
         e_d(ii,jj) = (1+zz[ii])^p1
      endif else begin
         e_d(ii,jj) = ((1+zc)^p1) * ( ((1+zz[ii])/(1+zc))^p2 ) 
      endelse

;      Abs_Mag_model[jj] = mag

      ple_denom   = ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = (A * e_d[ii,jj])  / ple_denom

      ;print, ii, jj, zz[ii], mag, zc, A, e_d[ii, jj],  (A * e_d[ii,jj]), ple_denom
;      print,      ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj],  $ 
 ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.6)'

      printf, 15, zz[ii], mag, Phi[ii,jj];, $
;              format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
;      printf, 11, ii, zz[ii], jj, mag, e_d[ii,jj], zc, Phi[ii,jj], $
;              format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
   endfor
endfor
close, 15
print, 'Bongiorno07_VVDS_LDDE_temp made... '
print
print



;;---------------------------------------------------------------
;;  
;;   L  E  D  E
;;
;;   Luminosity Evolution + Density Evolution 
;;   Table 4 from Croom et al. (2009b)

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
print, 'Croom09b_LEDE_temp.dat made.... '
print
print



fltarr_bins = 150
;; No. of mag bins
;mag_bins = 60
mag_bins = 40

;; No. of redshift bins
;zbins = 100. 
zbins = 150

starting_redshift = 0.30

;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
;zz = findgen(zbins)/20.    ;;  Redshift bins..
zz = ((findgen(fltarr_bins))/25.)+starting_redshift


z_out = ((findgen(fltarr_bins))/25.)+starting_redshift
M_range_out   = (findgen(40)*0.30)-31.05


;; Set up phi...
;Phi = fltarr(zbins, mag_bins)
 Phi = fltarr(fltarr_bins,40)

;;---------------------------------------------------------------
;;  
;;   L  A  D  E
;;
;;   Table 4 from Aird et al. (2010)

    alpha = -3.33
     beta = -1.41
 Mstar0_g = -22.17
       k1 = 1.46
       k2 = -0.328
 Phi_star = 10^(-5.84)

 openw, 14, 'Aird10_LADE_temp.dat'
printf, 14, '# '
printf, 14, '# alpha = -3.33,  beta = -1.41, k1 = 1.46,  k2==0.328, Phi_star = 10^(-5.84)'
printf, 14, '# ii, z[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom '
printf, 14, '# '

for ii=0L, zbins-1 do begin
   Hz = H0 * sqrt(omega_m*(1.0+zz[ii])^3.0 + omega_lam)
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)
      ;; => if mag_bins = 60, -32.00 < M_g < -17.25

      dmag   = mag - Mstar_g 
      ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

      ;print, 'A', ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 14, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom, $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8,e18.8)'
   endfor
;   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
;   phi  /= hub**3  ;; Convert from Mpc to Mpc/h volumes.
;   print, ii, zz[ii], 
endfor
close, 14
print, 'Aird10_LADE_temp.dat made ....'
print
print




;; No. of mag bins
mag_bins = 60

;; No. of redshift bins
zbins = 100. 
;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
zz = findgen(zbins)/20.    ;;  Redshift bins..

;; Set up phi...
Phi = fltarr(zbins, mag_bins)



;;------------------------------------------------------------------------------------------
;;
;;  F O R    T H E    B O S S    D R 9    Q L F    P A P E R !!    arXiv v1 and v2!!
;;
;;------------------------------------------------------------------------------------------
;;  
;;   P L E  
;;
;;   Pure Luminosity Evolution 
;;
;;
;;  From Table 8 from Ross et al. arXiv:1210.6389v1 (Version O-N-E... ;-)
;;
;;  0.3 < z < 2.2
   alpha =  -1.16
    beta =  -3.37
Mstar0_g = -22.85
      k1 =  1.241  
      k2 =  -0.249 
Phi_star = 10^(-5.96)

;;  2.2 < z < 3.5
   alpha =  -1.52
    beta =  -3.10
Mstar0_g = -24.29
      k1 =   1.134  
      k2 =  -0.273 
Phi_star = 10^(-6.37)


;;  From Table 8 from Ross et al. "rv1"   (aka the Version T-W-O... ;-)
;; 

;;  0.30 < z < 2.2
   alpha =  -1.16
    beta =  -3.37
Mstar0_g = -22.85
      k1 =   1.241  
      k2 =  -0.249 
Phi_star = 10^(-5.96)

;;  2.2 < z < 3.5
;   alpha =  -1.53
;    beta =  -3.26
;Mstar0_g = -24.45
;      k1 =   1.095  
;      k2 =  -0.263 
;Phi_star = 10^(-6.36)


openw, 20, 'chi_sq_PLE_model_McG_XzY_temp.dat'
printf, 20, '# '
printf, 20, '#      alpha,        beta         k1          k2    alog10(Phi_star)  ' 
printf, 20, '# ', alpha,    beta,      k1,     k2,    alog10(Phi_star)
printf, 20, '#  ii, z[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom '
printf, 20, '# '

for ii=0L, zbins-1 do begin
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)
      ;; => if mag_bins = 60, -32.00 < M_g < -17.25

      dmag   = mag - Mstar_g 
      ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      Phi[ii, jj] = Phi_star /  ple_denom 

  ;    print, 'A', ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 20, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj], ple_denom, $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8,e18.8)'
   endfor
endfor
close, 20
close, /all
print, 'chi_sq_PLE_model_McG_2.2z3.5_temp.dat  made....'
print
print





;;---------------------------------------------------------------
;;  
;;   L  E  D  E
;;
;;   Luminosity Evolution + Density Evolution 

;;
;;  From Table 8 from Ross et al. arXiv:1210.6389v1 (Version O-N-E... ;-)
;;
        alpha = -1.42
        beta  = -3.53
     Mstar0_g = -26.70
           c1 = -0.604
           c2 = -0.678
phi_star_zero = 10^(-6.08)

;;  From Table 8 from Ross et al. arXiv:1210.6389v1 (Version T-W-O... ;-)
        alpha = -1.42
        beta  = -3.53
     Mstar0_g = -26.70
           c1 = -0.604
           c2 = -0.678
phi_star_zero = 10^(-6.08)


;;  From Table 8 from Ross et al. "rv1"   (aka the Version T-W-O... ;-)
;; 
;;  2.2 < z < 3.5
        alpha = -1.30
        beta  = -3.44
     Mstar0_g = -26.51
           c1 = -0.684
           c2 = -0.843
phi_star_zero = 10^(-5.92)


;; "where z_ref is fixed at z_ref =2.0"
z_ref = 2.20

openw, 22, 'chi_sq_loglinear_LEDE_model_McG_2.2z3.5_temp.dat'
printf, 22, '## ' 
printf, 22, '##    alpha,  beta,   Mstar0_g,  c1, c2, alog10(phi_star_zero)'
printf, 22, '## ', alpha,  beta,   Mstar0_g,  c1, c2, alog10(phi_star_zero), format='(a3,6f12.6)'
;print, '## ', alpha,  beta,   Mstar0_g,  c1, c2, alog10(phi_star_zero), format='(a4,f12.6)'
printf, 22, '## ii, zz[ii],   jj,   mag,  c1, c2, alpha, Phi[ii,jj] '
printf, 22, '##' 

for ii=0L, zbins-1 do begin
   
;   log_phi_star = alog10(phi_star_zero) + ( k_phi1 * zz[ii]*(1.0-(0.5*zz[ii])/k_phi2))
   log_phi_star = alog10(phi_star_zero) + (c1 * (zz[ii] - z_ref))
   Mstar_g      = Mstar0_g              + (c2 * (zz[ii] - z_ref))
   
   for jj=0L, mag_bins-1 do begin
      mag = -32.0 + (jj*0.25)      
      
      dmag   = mag - Mstar_g 
      Phi[ii, jj] = 10^(log_phi_star) / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
      
;      print, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
      printf, 22, ii, zz[ii], jj, mag, c1, c2, alpha, Phi[ii,jj], $
              format='(i5, f9.5, i5, f12.5, f9.5, f9.5, f9.5, e18.8)'
   endfor
endfor
close, 22
close, /all
print, 'chi_sq_loglinear_LEDE_model_McG_2.2z3.5_temp.dat  made... '
print
print








end
