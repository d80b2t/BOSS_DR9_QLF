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






;; 
;; You *HAVE* the (model) values of redshift and magnitude from 
;; the **data** (!!), no need to do anything else fancy here!
;; 

;;
;; No. of redshift bins
;;
;;  Surely, this just needs to be the numnber of mag_bins
;;  in the data file...
zbins = 3 
zz = findgen(zbins)/10. + 2.35   ;;  Redshift bins..

;;
;; No. of mag bins
;;  Surely, this just needs to be the numnber of mag_bins
;;  in the data file...
mag_bins = n_elements(ABS_MAG_BIN_BOSS_S82)
mag_PLE = fltarr(mag_bins)


;; phi_star, the normalization
phi_star_bins = 50
phi_star = fltarr(phi_star_bins)

;; alpha, the (bright-end in C09b) slope
alpha_bins = 100   ;; for realz...
;alpha_bins = 10     ;; for testing..
alpha = fltarr(alpha_bins)

;; A, the normalization
beta_bins = 100    ;; for realz...
;beta_bins = 10      ;; for testing...
beta = fltarr(beta_bins)

k1_bins = 10 
;k1 = fltarr(k1_bins)

k2_bins = 10 


;; 
;; Setting up the double-power law denominator...
;; 
ple_denom= fltarr(alpha_bins,beta_bins)


;; Set up phi...
;Phi_LDDE = fltarr(zbins, mag_bins)
;
Phi_PLE = fltarr(phi_star_bins, alpha_bins, beta_bins)
;Phi_PLE = fltarr(n_elements(z_bin_boss_s82), phi_star_bins, alpha_bins, beta_bins)
;Phi_PLE = fltarr(mag_bins, phi_star_bins, alpha_bins, beta_bins)
;Phi_PLE = fltarr(zbins,mag_bins, phi_star_bins, alpha_bins, beta_bins)

print
;print, 'Gives ',  N_qlf_bins, ' N_qlf_bins to perform a fit over...' 
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

alpha_step = 0.05      ;; 0.05 for "realz", 0.5 for testing..
beta_step  = 0.05      ;; 0.05 for "realz", 0.5 for testing..
phi_step   = 0.05

openw, 9, 'Stripe82_PLE_temp.dat'

;;
;;
;;


w_zbin = where(Z_BIN_BOSS_S82 ge 2.80 and Z_BIN_BOSS_S82 le 3.00, N)

chi_sq_min= 100000.

ple_denom = findgen(N_elements(z_bin_boss_s82 ))

zz        = z_bin_boss_s82 
zz        = z_bin_boss_s82[w_zbin]

;;
;; Equation (11) of Croom et al. (2009b)
;;
;Mstar_g = -22.0
for k1_i = 0ll, k1_bins do begin
   k1 =  0.960000 + (k1_i*0.1)
   print, 'k1_i', k1_i
;   k1 = 1.46
   
  for k2_i = 0ll, k2_bins do begin
      k2 = -0.828000 + (k2_i*0.1)
      print, 'k2_i', k2_i
;      k2 = -0.328

      Mstar_g  = Mstar0_g  - 2.5*(k1*zz + k2*zz*zz)
      
      mag_PLE = Abs_mag_bin_boss_s82[w_zbin]
      dmag    = mag_PLE - Mstar_g 
      
      for alpha_i = 0LL, alpha_bins-1 do begin
         alpha[alpha_i] = alpha_norm + (alpha_i*alpha_step)
;         print, 'alpha_i, ', alpha_i 
         
         for beta_i = 0LL, beta_bins-1 do begin
            beta[beta_i] = beta_norm + (beta_i*beta_step)
            ple_denom =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
;      ple_denom[alpha_i,beta_i] =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
            
            for phi_i = 0L, phi_star_bins-1 do begin
               Phi_star[phi_i] = 10^(Phi_star_norm+(phi_i*phi_step))
               xxx = Phi_star[phi_i]  / ple_denom
;         Phi_PLE[alpha_i, beta_i, phi_i] = Phi_star[phi_i]  / ple_denom
               
;;;               chi_sq_bin = ( 10^(log_Phi_BOSS_S82[w[kk]]) -  Phi_PLE[kk, phi_i, alpha_i, beta_i])^2 $
;;;                            /  (sigma_Phi_BOSS_S82[w[kk]]^2)
;;;              chi_sq = chi_sq + chi_sq_bin
               
               
               chi_sq_bin = ( (10^(log_Phi_BOSS_S82) -  xxx)^2) $
                            /  (sigma_Phi_BOSS_S82^2)
               
               chi_sq = total(chi_sq_bin)
;         print,  phi_i, alpha_i, beta_i,  chi_sq,xxx[0:10]
;         print, phi_i, alpha_i, beta_i,  chi_sq
               if chi_sq lt chi_sq_min then begin
                  print,  phi_i, alpha_i, beta_i, $
                          chi_sq_min, chi_sq,  $
                          (Phi_star_norm+(phi_i*phi_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step), k1, k2
                  
                  chi_sq_min  = chi_sq 
                  alpha_i_min = alpha_i
                  beta_i_min  = beta_i
                  phi_i_min   = phi_i
;                  k1_i_min    = k1_i
                  xxx_min = xxx
               endif
            endfor ;; beta
         endfor    ;; alpha
      endfor
   endfor
endfor

close, 9

min = 0. 
mid = 5.
max = 9.
;; help, phi_ple
;; PHI_PLE         FLOAT     = Array[3, 30, 50, 10, 10]

;oplot, mag_PLE, alog10(Phi_PLE[1,*, 0,  min, min]), color=40, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 10, mid, min]), color=60, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 10, max, min]), color=80, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, min, mid]), color=120, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, mid, mid]), color=160, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 35, max, mid]), color=200, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 24, mid, mid]), color=240, thick=2
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 30, mid, mid]), color=240, thick=2, linestyle=1
;oplot, mag_PLE, alog10(Phi_PLE[1,*, 36, mid, mid]), color=240, thick=2, linestyle=2

openw, 11, 'chi_sq_values_PLE_inside_mag_loop_temp.dat'
printf, 11, '# Normalization, A  Abs i-mag     Phi_BOSS_S82, sigma_Phi_BOSS_S82,    Phi_LDDE_model,  chi_sq_bin, chi_sq '
printf, 11, '# (-9.57+(phi_i*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]), sigma_Phi_BOSS_S82[w[kk]], Phi_LDDE[ii,kk, phi_i], chi_sq_bin, chi_sq '

openw, 10, 'chi_sq_values_PLE_temp.dat'
printf, 10, '# Normalization A,     chi_sq_bin, chi_sq '
printf, 10, '# (-9.57+(phi_i*0.01)),   chi_sq_bin, chi_sq '


;; http://en.wikipedia.org/wiki/Goodness_of_fit

;               chi_sq_bin = ( 10^(log_Phi_BOSS_S82[w[kk]]) -  Phi_PLE[kk, phi_i, alpha_i, beta_i])^2 $
;                            /  (sigma_Phi_BOSS_S82[w[kk]]^2)
 ;              chi_sq = chi_sq + chi_sq_bin


;            if chi_sq lt chi_sq_min then begin
 ;              print,  phi_i, alpha_i, beta_i, $
  ;                     chi_sq_min, chi_sq,  $
   ;                    (Phi_star_norm+(phi_i*phi_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step) 
    ;           
     ;          chi_sq_min  = chi_sq 
      ;         alpha_i_min = alpha_i
       ;        beta_i_min  = beta_i
        ;       phi_i_min   = phi_i
         ;   endif
            
close, 10
close, 11
close, /all


;;
;; Current best fit....
;;

red_bin_limit  = [   2.40,    2.25,    2.35,    2.45,    2.55,    2.65,    2.75,    2.90,    3.13,   3.35,     3.75,    4.00]
Abs_iMag_limit = [-24.640, -24.470, -24.585, -24.695, -24.800, -24.901, -24.998, -25.137, -25.337, -25.514, -25.806, -25.973]
;; Note, Abs_iMag_limit[0] never used...


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

;;
;; !p.multi=[0,4,0]  ==> 1 row, 4 columns....
;;
!p.multi=[0,4,0]

z_bin_range_min = 2.20
z_bin_range_max = 2.30
w = where(red_bin_limit ge z_bin_range_min and red_bin_limit le z_bin_range_max, N)
Abs_iMag_limit_bin = Abs_iMag_limit[w]

w_plot = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin[0], N_qlf_bins)

plotsym, 0, 1.6, /fill
plot, Abs_mag_bin_boss_s82[w_plot] , log_Phi_BOSS_S82[w_plot], $
      position=[0.20, 0.20, 0.36, 0.96], $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 ;, $ 
;      xtitle='!6Abs_mag_bin_boss_s82[w]', ytitle='!6log_Phi_BOSS_S82[w]'

Mi_boss_err_S82 = Abs_mag_bin_boss_s82[w_plot] - Abs_mag_bin_boss_s82[w_plot]
oploterror, Abs_Mag_bin_boss_S82[w_plot], log_Phi_boss_S82[w_plot], Mi_boss_err_S82, boss_delta_up_S82[w_plot],  $
            /hibar, errthick=3., psym=8 ;, color=160, ERRcolor=160
oploterror, Abs_Mag_bin_boss_S82[w_plot], log_Phi_boss_S82[w_plot], Mi_boss_err_S82, boss_delta_down_S82[w_plot], $   
            /lobar, errthick=3., psym=8 ;, color=160, errcolor=160

;oplot, mag_PLE, alog10(Phi_PLE[*,phi_i_min,alpha_i_min,beta_i_min]), thick=4
oplot, ABS_MAG_BIN_BOSS_S82[w_plot], alog10(xxx_min[w_plot]), thick=4



z_bin_range_min = 2.40
z_bin_range_max = 2.50
w = where(red_bin_limit ge z_bin_range_min and red_bin_limit le z_bin_range_max, N)
Abs_iMag_limit_bin = Abs_iMag_limit[w]

w_plot = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin[0], N_qlf_bins)

plot, Abs_mag_bin_boss_s82[w_plot] , log_Phi_BOSS_S82[w_plot], $
      position=[0.40, 0.20, 0.56, 0.96], $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4

oplot, ABS_MAG_BIN_BOSS_S82[w_plot], alog10(xxx_min[w_plot]), thick=4


z_bin_range_min = 2.60
z_bin_range_max = 2.70
w = where(red_bin_limit ge z_bin_range_min and red_bin_limit le z_bin_range_max, N)
Abs_iMag_limit_bin = Abs_iMag_limit[w]

w_plot = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin[0], N_qlf_bins)

plot, Abs_mag_bin_boss_s82[w_plot] , log_Phi_BOSS_S82[w_plot], $
      position=[0.60, 0.20, 0.76, 0.96], $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4

oplot, ABS_MAG_BIN_BOSS_S82[w_plot], alog10(xxx_min[w_plot]), thick=4



z_bin_range_min = 2.80
z_bin_range_max = 3.00
w = where(red_bin_limit ge z_bin_range_min and red_bin_limit le z_bin_range_max, N)
Abs_iMag_limit_bin = Abs_iMag_limit[w]

w_plot = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin[0], N_qlf_bins)

plot, Abs_mag_bin_boss_s82[w_plot] , log_Phi_BOSS_S82[w_plot], $
      position=[0.80, 0.20, 0.96, 0.96], $
      psym=8, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4

oplot, ABS_MAG_BIN_BOSS_S82[w_plot], alog10(xxx_min[w_plot]), thick=4

!p.multi=0

end
