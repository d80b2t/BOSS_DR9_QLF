;+
; NAME:
;   chi_sq_PLE
;
; PURPOSE:
;   Okay, so this is just a fairly wee, fairly simple script that
;   takes the ("pure") PLE model from Croom et al. (2009b), 
;   and computes some chi_squared values for the BOSS DR9 dataset
;
; CALLING SEQUENCE:
;   .run chi_sq_PLE
;
; INPUTS:
;   A QLF datafile, with columns:
;      z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma  
;
; OUTPUTS:
;   Plot to screen, chi_sq values...
;
; NOTES:
;    I'm initially testing on the narrowZ, S82 dataset...
;
; REVISION HISTORY:
;   11-Jan-2011  v0.0.1     NPR
;-


;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent



;filename = 'My_QLF_iband_boss_narrowZ_S82_wExtinc_20120622.dat'
filename = 'My_QLF_iband_boss_narrowZ_wgt_20120619_4modelfits.dat'
dir_file = '../qlf/'+filename

;;
;;  THE   S T R I P E   8 2   DATA  - 6,250 Quasars strong
;;
;readcol, '../qlf/My_QLF_iband_boss_narrowZ_S82.dat', $  ;; USE THIS FOR THE S82 plots!!!
readcol, dir_file, $
         z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Phi_BOSS_S82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82


boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)

print
print, ' READ-IN::  ', filename, ' with ', n_elements(z_bin_boss_s82), ' n_elements ' 
print

;
; boss_ilimit = 21.80
;
;;red_bins = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]
;
;dlums = DLUMINOSITY(red_bins) / 1e6
;Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ;+ kcor(fix(red_bins/0.01))
red_bin_limit  = [   2.40,    2.25,    2.35,    2.45,    2.55,    2.65,    2.75,    2.90,    3.13,   3.35,     3.75,    4.00]
Abs_iMag_limit = [-24.640, -24.470, -24.585, -24.695, -24.800, -24.901, -24.998, -25.137, -25.337, -25.514, -25.806, -25.973]
;; Note, Abs_iMag_limit[0] never used...

z_bin_range_min = 2.50
z_bin_range_max = 2.60

print
print, 'z_bin_range_min  set as  ', z_bin_range_min 
print, 'z_bin_range_max  set as  ', z_bin_range_max
print 


w = where(red_bin_limit ge z_bin_range_min and red_bin_limit le z_bin_range_max, N)
Abs_iMag_limit_bin = Abs_iMag_limit[w]

;ii = 45  ;; picking out zz[ii]=2.25 
;ii = 47  ;; picking out zz[ii]=2.35 
;ii = 49  ;; picking out zz[ii]=2.45
ii = 51  ;; picking out zz[ii]=2.55
;ii = 53  ;; picking out zz[ii]=2.65
;ii = 55  ;; picking out zz[ii]=2.75
;ii = 58  ;; picking out zz[ii]=2.90
;ii = 63  ;; picking out zz[ii]=3.15
;ii = 67  ;; picking out zz[ii]=3.35
;ii = 75  ;; picking out zz[ii]=3.75

;w = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max and raw_N_boss_QSOs_s82 ge 1, N)

w = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin[0], N_qlf_bins)

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03

;;
;; Trying to match figs in the paper(s)...
;;
;; x-ranges
;xmin = -18.001
xmin = -22.001
xmax = -30.50

;; y-ranges
ymin = -9.20
ymax = -5.10


plotsym, 0, 1.6, /fill
plot, Abs_mag_bin_boss_s82[w] , log_Phi_BOSS_S82[w], $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4, $ 
      xtitle='!6Abs_mag_bin_boss_s82[w]', ytitle='!6log_Phi_BOSS_S82[w]'

Mi_boss_err_S82 = Abs_mag_bin_boss_s82[w] - Abs_mag_bin_boss_s82[w]
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_up_S82[w],  $
            /hibar, errthick=3., psym=8 ;, color=160, ERRcolor=160
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_down_S82[w], $   
            /lobar, errthick=3., psym=8 ;, color=160, errcolor=160


;;
;; No. of redshift bins
;;
;;zbins = 140. 
;;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
;;zz = findgen(zbins)/20.    ;;  Redshift bins..

;; 51 bins, /10 +offset gives  e.g. 0.90, 1.00, ... , 5.80, 5.90
;zbins = 51. 
;zz = findgen(zbins)/10. + 0.90   ;;  Redshift bins..

;; 3 bins, /10 +offset gives  e.g. 2.40, 2.50, 2.60
zbins = 3 
zz = findgen(zbins)/10. + 2.35   ;;  Redshift bins..


;; No. of mag bins
mag_bins = 30
mag_PLE = fltarr(mag_bins)

;; A, the normalization
phi_star_bins = 50
phi_star = fltarr(phi_star_bins)

;; alpha, the (bright-end in C09b) slope
;alpha_bins = 100   ;; for realz...
alpha_bins = 10     ;; for testing..
alpha = fltarr(alpha_bins)

;; A, the normalization
;beta_bins = 100    ;; for realz...
beta_bins = 10      ;; for testing...
beta = fltarr(beta_bins)

;; Set up phi...
;Phi_LDDE = fltarr(zbins, mag_bins)
;Phi_PLE = fltarr(mag_bins, phi_star_bins, alpha_bins, beta_bins)
Phi_PLE = fltarr(zbins,mag_bins, phi_star_bins, alpha_bins, beta_bins)

print
print
;print, 'Things to note: ii=', ii, '  implies you are picking out zbin[ii] =  ', zz[ii]
print 
print, 'z_bin_range_min , z_bin_range_max, Abs_iMag_limit_bin'
print,  z_bin_range_min , z_bin_range_max, Abs_iMag_limit_bin
print
print, 'Gives ',  N_qlf_bins, ' N_qlf_bins to perform a fit over...' 
print

;---------------------------------------------------------------
; 
;  P  L  E 
;
;   Good ol' PLE!! 
;   Section 6.1 from Croom et al. (2009b)
;   e.g. Table 2, Row 4, 

;    alpha = -3.33
;    beta = -1.41
Mstar0_g = -22.17
k1 = 1.46
k2 = -0.328
; Phi_star = 10^(-5.84)

Phi_star_norm = -7.27
alpha_norm = -5.50
beta_norm  = -3.00

alpha_step = 0.5      ;; 0.05 for "realz", 0.5 for testing..
beta_step  = 0.5      ;; 0.05 for "realz", 0.5 for testing..
phi_step    = 0.05

openw, 9, 'Stripe82_PLE_temp.dat'

;;
;; P A R T   A 
;;
;; First part of the program, is to set-up the
;; e.g. PLE model QLFs...
;;
;;

;Mstar_g = -22.0
for zbin_i = 0L, zbins-1 do begin
;;
;; Equation (11) of Croom et al. (2009b)
;;
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[zbin_i] + k2*zz[zbin_i]*zz[zbin_i])
   print, 'zbin_i, zz[zbin_i], Mstar_g  ', zbin_i, zz[zbin_i], Mstar_g 
   ;; => Mstar_g(z=2.25) = -26.2313
   ;; => Mstar_g(z=2.45) = -26.1905
   
   for alpha_i = 0LL, alpha_bins-1 do begin
      alpha[alpha_i] = alpha_norm + (alpha_i*alpha_step)
      print, 'alpha_i, ', alpha_i 
      
      for beta_i = 0LL, beta_bins-1 do begin
         beta[beta_i] = beta_norm + (beta_i*beta_step)
         
         for phi_i = 0L, phi_star_bins-1 do begin
            Phi_star[phi_i] = 10^(Phi_star_norm+(phi_i*phi_step))
            
            for jj=0L, mag_bins-1 do begin
               ;;mag_PLE[jj] = -32.0 + (jj*0.30)                       ;; To match the data's delta_mags
               mag_PLE[jj] = min(Abs_mag_bin_boss_s82[w]) + (jj*0.30)  ;; and put them on the same "absolute scale" 
               
               dmag    = mag_PLE[jj] - Mstar_g 
               
               ple_denom =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
               Phi_PLE[zbin_i, jj, phi_i, alpha_i, beta_i] = Phi_star[phi_i]  / ple_denom
               
               printf, 9, zz[zbin_i], mag_PLE[jj], alog10(Phi_PLE[zbin_i, jj, phi_i, alpha_i, beta_i]), alpha[alpha_i]  , beta[beta_i], $
                       format='(d,d,d, d,d)'
               
;            print, ii, jj, alpha_i, beta_i, phi_i, zz[ii], alog10(phi_star[phi_i]), mag_PLE[jj], ple_denom, Phi_PLE[jj, phi_i, alpha_i, beta_i] 
;      printf, 9, ii, zz[ii], jj, mag_LDDE[jj], e_d[ii,jj], zc, Phi_LDDE[ii,jj], $
                                ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
            endfor
;         if (phi_i mod 2) eq 1 then oplot, mag_PLE, alog10(Phi_PLE[*, phi_i, alpha_i, beta_i]), color=100
         endfor
      endfor ;; beta
   endfor    ;; alpha
endfor
close, 9

min = 0. 
mid = 5.
max = 9.
;; help, phi_ple
;; PHI_PLE         FLOAT     = Array[3, 30, 50, 10, 10]

oplot, mag_PLE, alog10(Phi_PLE[1,*, 0,  min, min]), color=40, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 10, mid, min]), color=60, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 10, max, min]), color=80, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, min, mid]), color=120, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, mid, mid]), color=160, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, max, mid]), color=200, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 24, mid, mid]), color=240, thick=2
oplot, mag_PLE, alog10(Phi_PLE[1,*, 30, mid, mid]), color=240, thick=2, linestyle=1
oplot, mag_PLE, alog10(Phi_PLE[1,*, 36, mid, mid]), color=240, thick=2, linestyle=2

openw, 11, 'chi_sq_values_PLE_inside_mag_loop_temp.dat'
printf, 11, '# Normalization, A  Abs i-mag     Phi_BOSS_S82, sigma_Phi_BOSS_S82,    Phi_LDDE_model,  chi_sq_bin, chi_sq '
printf, 11, '# (-9.57+(phi_i*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]), sigma_Phi_BOSS_S82[w[kk]], Phi_LDDE[ii,kk, phi_i], chi_sq_bin, chi_sq '

openw, 10, 'chi_sq_values_PLE_temp.dat'
printf, 10, '# Normalization A,     chi_sq_bin, chi_sq '
printf, 10, '# (-9.57+(phi_i*0.01)),   chi_sq_bin, chi_sq '


;; http://en.wikipedia.org/wiki/Goodness_of_fit

w = where(z_bin_boss_S82 ge min(zz) and z_bin_boss_S82 lt max(zz), N_qlf_bins)

chi_sq_min = 1000.
for zbin_i = 0LL, zbins - 1 do begin
   for alpha_i = 0LL, alpha_bins-1 do begin
      for beta_i = 0LL, beta_bins-1 do begin
         for phi_i = 0L, phi_star_bins-1 do begin
            chi_sq=0         
            
            for kk=0L, N_qlf_bins-1 do begin
; Phi_PLE[zbin_i, jj, phi_i, alpha_i, beta_i]
                                ;chi_sq_bin = ( log_Phi_BOSS_S82[w[kk]] -  alog10(Phi_LDDE[ii,kk]))^2 /  (alog10(sigma_Phi_BOSS_S82[w[kk]])^2)
               chi_sq_bin = ( 10^(log_Phi_BOSS_S82[w[kk]]) -  Phi_PLE[kk, phi_i, alpha_i, beta_i])^2 $
                            /  (sigma_Phi_BOSS_S82[w[kk]]^2)
               chi_sq = chi_sq + chi_sq_bin
               
                                ;print, phi_i, kk, Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]])*1e7, sigma_Phi_BOSS_S82[w[kk]]*1e7, $
                                ;       Phi_LDDE[ii,kk,phi_i]*1e7, chi_sq_bin, chi_sq
                                ;print, Abs_mag_bin_boss_s82[w[kk]], log_Phi_BOSS_S82[w[kk]], alog10(sigma_Phi_BOSS_S82[w[kk]]), Phi_LDDE[ii,kk], $
                                ;       chi_sq_bin, chi_sq[kk]
               
;            printf, 11, (Phi_star_norm+(phi_i*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]),  sigma_Phi_BOSS_S82[w[kk]], $
;;                    Phi_LDDE = fltarr(mag_bins, A_norm_bins, 1, beta_bins)
                                ;                   Phi_PLE[kk, phi_i, alpha_i, beta_i], $
                                ;                  chi_sq_bin, chi_sq, $
;               format='(f9.4,2x, f9.4,2x, e,e,e, f, f)'
            endfor
            
            if chi_sq lt chi_sq_min then begin
               print,  phi_i, alpha_i, beta_i, $
                       chi_sq_min, chi_sq,  $
                       (Phi_star_norm+(phi_i*phi_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step) 
               
               chi_sq_min  = chi_sq 
               alpha_i_min = alpha_i
               beta_i_min  = beta_i
               phi_i_min   = phi_i
            endif
            
         endfor
      endfor
   endfor
endfor
   

print
print
;print,  'alpha_i, beta_i, Ai, chi_sq_min, chi_sq,  (-9.57+(Ai*0.01)),  -3.00 +(beta_i*0.02), -4.00 +(alpha_i*0.02)'
print
print, 'phi_i_min,   ', phi_i_min
print, 'alpha_i_min, ', alpha_i_min
print, 'beta_i_min,  ', beta_i_min
print, 'chi_sq_min,  ', chi_sq_min
print, 'chi_sq,      ', chi_sq
print, '(Phi_star_norm+(phi_i*phi_step))       ', (Phi_star_norm+(phi_i*phi_step))
print, '(alpha_norm+(alpha_i_min*alpha_step)), ', (alpha_norm+(alpha_i_min*alpha_step))
print, '(beta_norm+(beta_i_min*beta_step))     ', (beta_norm+(beta_i_min*beta_step))  
print
;print,  phi_i_min, alpha_i_min, beta_i_min, chi_sq_min, chi_sq, (Phi_star_norm+(phi_i*phi_step)), (alpha_norm+(alpha_i_min*alpha_step)), (beta_norm+(beta_i_min*beta_step)) 
print
print

close, 10
close, 11
close, /all

;; Current best fit....
oplot, mag_PLE, alog10(Phi_PLE[*,phi_i_min,alpha_i_min,beta_i_min]), thick=4





end
