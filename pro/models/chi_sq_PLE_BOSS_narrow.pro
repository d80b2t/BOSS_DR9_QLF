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


;;
;;   B O S S  +  S D S S
;;

;readcol, '../../pro/qlf/temp_for_chi_sq_fits.dat', $
;readcol, 'QLF_iband_boss_narrowZ_wgt_formatted_20120731.dat', $
readcol, '../qlf/QLF_iband_boss_narrowZ_formatted4paper_5731_tozeq4.dat', $
         z_full,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, $ 
         NQ_full, ba,  log_Phi_full, bar, error_full, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)'

;;w_zbin = where(z_2SLAQ_full ge 0.40 and Z_2SLAQ_full le 2.60, N)
;w_zbin = where(z_full ge 0.20 and Z_full le 2.20, N)
;w_zbin = where(z_full ge 1.0 and Z_full le 2.20, N)
w_zbin = where(z_full ge 2.20 and Z_full le 2.30, N)


print
print
print, 'No of points that we will be using...', N
print

;; actually kinda have options whether to use the Mi "bin centre" 
;; or the actual mean of the Mi values... 
          Mg_2SLAQ = Mi_bin_full[w_zbin] 
           z_2SLAQ = z_full[w_zbin]
                NQ = NQ_full[w_zbin]
     log_phi_2SLAQ = log_phi_full[w_zbin]
         sigma_Phi = error_full[w_zbin]*1e-9

ple_denom = findgen(N_elements(z_2SLAQ ))
zz        = z_2SLAQ

    error_delta  = (log_phi_2slaq - alog10(sigma_Phi))

      delta_up   = alog10(10^(log_phi_2slaq) +  sigma_Phi)
error_delta_up   = delta_up - log_phi_2slaq
      delta_down = alog10(10^(log_phi_2slaq) -  sigma_Phi)
error_delta_down = abs(delta_down - log_phi_2slaq)



;; 
;; You *HAVE* the (model) values of redshift and magnitude from 
;; the **data** (!!), no need to do anything else fancy here!
;; 

;; phi_star, the normalization
phi_star_bins = 250
phi_star = fltarr(phi_star_bins)

;; alpha, the (bright-end in C09b) slope
alpha_bins = 250   ;; for realz...
alpha_bins = 100   ;; for realz...
alpha_bins = 10     ;; for testing..
alpha = fltarr(alpha_bins)

;; beta, the (faint-end slope in C09b)
beta_bins = 250    ;; for realz...
beta_bins = 100    ;; for realz...
beta_bins = 10      ;; for testing...
beta = fltarr(beta_bins)

;; 
;; Mstar_g  = Mstar0_g  - 2.5*(k1*zz + k2*zz*zz)
;;
;Mstar0_g_bins = 250  ;; this is gonna kill me...
;Mstar0_g_bins = 50  ;; this is gonna kill me...
Mstar0_g_bins = 25  ;; this is gonna kill me...
Mstar0_g  = fltarr(Mstar0_g_bins)


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
;Mstar0_g = -22.17
k1 = 1.46
k2 = -0.328
; Phi_star = 10^(-5.84)

Phi_star_norm = -7.27
alpha_norm    = -6.50
beta_norm     = -3.00

Mstar0_g_norm = -26.5    

alpha_step    = 0.50      ;; 0.01 or 0.05 for "realz", 0.5 for testing..
beta_step     = 0.50      ;; 0.01, 0.05 for "realz", 0.5 for testing..
phi_step      = 0.01
Mstar0_g_step = 0.25
k1_step       = 0.1
k2_step       = 0.1

openw, 9, 'BOSS_PLE_temp.dat'


chi_sq_min = 100000.
xxx_0min    = 100000.

;;
;; Equation (11) of Croom et al. (2009b)
;;
;Mstar_g = -22.0
;for k1_i = 0ll, k1_bins do begin
;   k1 =  0.960000 + (k1_i*k1_step)
;  print, 'k1_i  ', k1_i, ' k1 ', k1
;   k1 = 1.46  ;; C09b fit
    k1 = 1.16
   
;  for k2_i = 0ll, k2_bins do begin
 ;    k2 = -0.828000 + (k2_i*k2_step)
;      print, 'k2_i  ', k2_i, '  k2  ', k2
    ;  k2 = -0.328 ;; C09b fit
      k2= -0.228

      for Mstar0_g_i = 0ll, Mstar0_g_bins-1 do begin
;      Mstar_g  = Mstar0_g[Mstar0_g_i] - 2.5*(k1*zz + k2*zz*zz)
      Mstar0_g[Mstar0_g_i] = Mstar0_g_norm  + ( Mstar0_g_i * Mstar0_g_step)
      print,  Mstar0_g_i,  Mstar0_g[Mstar0_g_i] 
      Mstar_g  = Mstar0_g[Mstar0_g_i]  - 2.5*(k1*zz + k2*zz*zz)

      
;      mag_PLE = Abs_mag_bin_boss_s82[w_zbin]
      mag_PLE = Mg_2SLAQ
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
               
               if xxx[0] lt xxx_0min then  xxx_0min = xxx[0] 

           ;    print, Phi_star[phi_i], alpha[alpha_i], beta[beta_i], log_phi_2SLAQ[0], xxx[0], sigma_Phi_2SLAQ[0], xxx_0min 

               chi_sq_bin = ( (10^(log_phi_2SLAQ) -  xxx)^2) $
                            /  (sigma_Phi^2)

               chi_sq = total(chi_sq_bin)
;         print,  phi_i, alpha_i, beta_i,  chi_sq,xxx[0:10]
;         print, phi_i, alpha_i, beta_i,  chi_sq
               if chi_sq lt chi_sq_min then begin
                  print,  phi_i,      alpha_i,  beta_i, $
                          chi_sq_min, chi_sq,   $
                          (Phi_star_norm+(phi_i*phi_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step), $
                          k1, k2,  Mstar0_g[MSTAR0_G_I],  $
                          format='(i8,i8,i8,  f15.4,f15.4,  f9.4,f9.4,f9.4,  f9.4,f9.4, f9.4 )' 
                  
                  chi_sq_min  = chi_sq 
                  alpha_i_min = alpha_i
                  beta_i_min  = beta_i
                  phi_i_min   = phi_i
               MSTAR0_G_I_min = MSTAR0_G_I
 ;                 k1_i_min    = k1_i
  ;                k2_i_min    = k2_i
                  xxx_min = xxx
               endif

            endfor ;; beta
         endfor    ;; alpha
      endfor
   endfor ;;Mg_star
;   endfor ;; k2
;endfor   ;; k1


;alpha_test =  -3.33
; beta_test =  -1.41
alpha_test =    alpha[alpha_i_min]
 beta_test =  beta[beta_i_min]
;Mstar0_g_test = -22.17
Mstar0_g_test = Mstar0_g[ MSTAR0_G_I_min]
k1_test =  1.16   ; k1_i_min
k2_test = -0.228  ; k2_i_min
;Phi_star_test =  10^(-5.84)
Phi_star_test =  Phi_star[phi_i_min]  ;10^(-5.84)
z_test  = z_full
Mg_test = MI_BIN_FULL

Mstar_g_test  = Mstar0_g_test  - 2.5*(k1_test*z_test + k2_test*z_test*z_test)
dmag_test = Mg_test - Mstar_g_test
ple_denom_test =  ( (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) )

phi_best_fit = Phi_star_test /ple_denom_test

chi_sq_bin_test = ( (10^(log_phi_2SLAQ) -  phi_best_fit)^2)  /  (sigma_Phi^2)


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

close, 10
close, 11
close, /all


;;
;; Current best fit....
;;

red_bin_limit  = [   0.40,    0.68, 1.06, 1.44, 1.86, 2.20, 2.60]
;red_bin_limit  = [   2.40,    2.25,    2.35,    2.45,    2.55,    2.65,    2.75,    2.90,    3.13,   3.35,     3.75,    4.00]
Abs_iMag_limit = [-24.640, -24.470, -24.585, -24.695, -24.800, -24.901, -24.998, -25.137, -25.337, -25.514, -25.806, -25.973]
;; Note, Abs_iMag_limit[0] never used...

;; Colour Table
clr_table =13
loadct, clr_table

;; Colours for clr_table =13
black      =   0
purple     =  32
deep_blue  =  48
blue       =  64
light_blue =  80
turquiose  = 128
green      = 150
yellow     = 210
orange     = 232
red        = 254

charsize  = 3.6
charthick = 4.8
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.04
YTICKLEN  = 0.06


;;
;; Trying to match figs in the paper(s)...
;;
;; x-ranges
;xmin = -18.001
xmin = -22.001
xmax = -30.50

;; y-ranges
ymin = -9.80  ;; -9.20
ymax = -4.40  ;; -5.10


;;
;; !p.multi=[0,4,0]  ==> 1 row, 4 columns....
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_iband_PLE_narrow_temp.ps', $
        xsize=18.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated
color_2slaq = red
;;plot_sym_size_BOSS = 2.2
sym_size = 2.2
line_thick = 6.0

;;
;; Do a mini-reset with the read-in data, so that *all* the data gets plotted...
;;
w_zbin = where(z_full ge 0.20 and Z_full le 5.60, N)
;; actually kinda have options whether to use the Mi "bin centre" 
;; or the actual mean of the Mi values... 
          Mg_2SLAQ = Mi_bin_full[w_zbin] 
           z_2SLAQ = z_full[w_zbin]
                NQ = NQ_full[w_zbin]
     log_phi_2SLAQ = log_phi_full[w_zbin]
         sigma_Phi = error_full[w_zbin]*1e-9
        delta_up   = alog10(10^(log_phi_2slaq) +  sigma_Phi)
  error_delta_up   = delta_up - log_phi_2slaq
        delta_down = alog10(10^(log_phi_2slaq) -  sigma_Phi)
  error_delta_down = abs(delta_down - log_phi_2slaq)



readcol, '../../pro/qlf/temp_for_chi_sq_fits.dat', $
         z_firstbin,  fo,  Mi_mean_firstbin, foo, Mi_bin_firstbin, fooo, $ 
         NQ_first, ba,  log_Phi_firstbin, bar, error_firstbin, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)', /silent

;; z = 2.40
z_bin_range_min = 2.20
z_bin_range_max = 2.60
w_plot = where(z_firstbin ge z_bin_range_min and z_firstbin le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, Mi_bin_firstbin[w_plot] , LOG_PHI_firstbin[w_plot], psym=8, $
      position=[0.14, 0.70, 0.35, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN,  xtickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
;oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
;oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

;oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color



;;  z  =  0.87
z_bin_range_min = 2.20
z_bin_range_max = 2.30
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.35, 0.70, 0.56, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color



;;  z  =  1.25
z_bin_range_min = 2.30
z_bin_range_max = 2.40
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.70, 0.77, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


z_bin_range_min = 2.40
z_bin_range_max = 2.50
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
       position=[0.77, 0.70, 0.98, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;;      z = 2.55       
z_bin_range_min = 2.50
z_bin_range_max = 2.60
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.14, 0.42, 0.35, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;   z = 2.65
z_bin_range_min = 2.60
z_bin_range_max = 2.70
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.35, 0.42, 0.56, 0.70], $
;      position=[0.20, 0.20, 0.36, 0.96], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 ;, $ 
;      xtitle='!6Abs_mag_bin_boss_s82[w]', ytitle='!6log_Phi_BOSS_S82[w]'

;color_2slaq = 0
;plotsym, sym_2slaq, plot_sym_size_R06, /fill
M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]

oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;  z = 2.80
z_bin_range_min = 2.70
z_bin_range_max = 2.80
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.42, 0.77, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;  z = 3.25
z_bin_range_min = 2.80
z_bin_range_max = 3.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.77, 0.42, 0.98, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;  z = 3.
z_bin_range_min = 3.00
z_bin_range_max = 3.25
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.14, 0.14, 0.35, 0.42], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

;color_2slaq = orange
M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq
oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;  z = 4.25
z_bin_range_min = 3.25
z_bin_range_max = 3.50
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
       position=[0.35, 0.14, 0.56, 0.42], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq
oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;  z = 4.75
z_bin_range_min = 3.50
z_bin_range_max = 4.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

if N_qlf_bins gt 0 then begin
   plotsym, 0, sym_size, /fill
   plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
         position=[0.56, 0.14, 0.77, 0.42], $
         xrange=[xmin, xmax], yrange=[ymin, ymax], $
         xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
         XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
         charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 
   
   M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
   oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
   oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq
   oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color
endif

device, /close
set_plot, 'X'
close, /all

!p.multi=0

end
