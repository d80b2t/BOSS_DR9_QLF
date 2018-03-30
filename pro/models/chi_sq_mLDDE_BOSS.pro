;+
; NAME:
;   chi_sq_mLDDE
;
; PURPOSE:
;   Okay, so this is just a fairly wee, fairly simple script that
;   takes the ("modifed") mLDDE model from Croom et al. (2009b), 
;   and computes some chi_squared values for the BOSS DR9 dataset
;
; CALLING SEQUENCE:
;   .run chi_sq_mLDDE
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


coarse_or_fine  = 0
;read, coarse_or_fine, PROMPT=' (0) coarse,  (1) Medium, or (2) Fine model fitting?? '
;; ``fine'' still to be fully implemented!


;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent


;;
;;   B O S S  +  S D S S
;;

readcol, '../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_20120731.dat', $
         z_full,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, NQ_full, $
         ba,  log_Phi_full, bar, error_full, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)'

readcol, '../../pro/qlf/QLF_iband_boss_wgt_formatted4paper_S82.dat', $
;readcol, 'My_QLF_iband_boss_wgt_formatted4paper_temp.dat', $
         z_full,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, NQ_full, $
         format='(d,a, d,a, d,a,  d)'


;w_zbin = where(z_full ge 2.20 and Z_full le 3.50, N)  ;; N= 74
;w_zbin = where(z_full ge 0.20 and Z_full le 2.20, N)  ;; N= 74
;w_zbin = where(z_full ge 1.0 and Z_full le 2.20, N)    ;; N = 50
;w_zbin = where(z_full ge 0.2 and Z_full le 2.60, N)    ;; N = 89
;w_zbin = where(z_full ge 0.2 and Z_full le 3.50, N)    ;; N = 120
w_zbin = where(z_full ge 0.2 and Z_full le 3.80, N)    ;; N = 128

;w_zbin = where(z_full ge 2.20 and Z_full le 3.50, N)
;w_zbin = where(z_full ge 2.20 and Z_full le 2.30, N)
;w_zbin = where(z_full ge 1.0 and Z_full le 3.50, N)
;w_zbin = where(z_full ge 3.00 and Z_full le 3.50, N)


print
print
print, 'No of points that we will be using...', N
print
print, ' minmax(z_full[w_zbin]) ', minmax(z_full[w_zbin])
print
print

;; actually kinda have options whether to use 
;; the Mi "bin centre" or the actual mean of 
;; the Mi values... 
          Mg_2SLAQ = Mi_bin_full[w_zbin] 
           z_2SLAQ = z_full[w_zbin]
                NQ = NQ_full[w_zbin]
     log_phi_2SLAQ = log_phi_full[w_zbin]
         sigma_Phi = error_full[w_zbin]*1e-9

denom = findgen(N_elements(z_2SLAQ ))
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

;; A, the normalization
;A_bins = 250
;A_bins = 100
A_bins = 50
A = fltarr(A_bins)

;; alpha, the (bright-end in C09b) slope
alpha_bins = 250   ;; for realz...
alpha_bins = 100   ;; for realz...
alpha_bins = 20     ;; for testing..
;alpha_bins = 10     ;; for testing..
alpha = fltarr(alpha_bins)

;; beta, the (faint-end slope in C09b)
beta_bins = 250    ;; for realz...
beta_bins = 100    ;; for realz...
beta_bins = 20     ;; for testing..
;beta_bins = 10      ;; for testing...
beta = fltarr(beta_bins)

;; 

;;      top_ed    = 2.*(1+zc)^(p1)
;;      bottom_ed = (((1+zz[ii])/(1+zc))^(-1.*p1)) + (((1+zz[ii])/(1+zc))^(-1.*p2))
;;
;Mstar0_g_bins = 250  ;; this is gonna kill me...
;Mstar0_g_bins = 50  
Mstar0_g_bins = 25 
Mstar0_g  = fltarr(Mstar0_g_bins)


;; 
p1_bins = 10 
p1 = fltarr(p1_bins)

p2_bins = 10 
p2 = fltarr(p2_bins)


;;  M_g_c  staying fixed for the time being...
;M_g_c_bins = 100 
M_g_c_bins = 20
M_g_c  = -23.90

;; as is the thing "gamma"
gamma_bins = 20
;gamma  =   0.68 
gamma = fltarr(gamma_bins)

;; 
;; Setting up the double-power law denominator...
;; 
denom= fltarr(alpha_bins,beta_bins)

;;
;; Set up phi...
;;
Phi_mLDDE = fltarr(A_bins, alpha_bins, beta_bins)
;Phi_mLDDE = fltarr(n_elements(z_bin_boss_s82), phi_star_bins, alpha_bins, beta_bins)
;Phi_mLDDE = fltarr(mag_bins, phi_star_bins, alpha_bins, beta_bins)
;Phi_mLDDE = fltarr(zbins,mag_bins, phi_star_bins, alpha_bins, beta_bins)

print
;print, 'Gives ',  N_qlf_bins, ' N_qlf_bins to perform a fit over...' 
print

;---------------------------------------------------------------
; 
;  m  L  D  D  E 
;
;   The ``modifed'' m L D D E!! 
;   Section 6.2 from Croom et al. (2009b)
;   e.g. Table 3

;    alpha = -3.70    ;; varying...
;    beta = -2.34     ;; varying...
; Mstar0_g = -26.69   ;; varying...
;  M_g_c  = -23.90    
;  gamma  =   0.68 
zc_zero  =   2.47
;p1 = 6.28           ;;  varying...
;p2 = -2.85          ;;  varying...
;      A = 10^(-9.21)
log_A = -9.21

if coarse_or_fine  eq 0 then begin
   A_norm        = -10.50
   alpha_norm    = -5.50
   beta_norm     = -3.50
   
   Mstar0_g_norm = -28.5    
;   Mstar0_g_norm = -26.5    

   alpha_step    = 0.50      ;; 0.01 or 0.05 for "realz", 0.5 for testing..
   beta_step     = 0.50      ;; 0.01, 0.05 for "realz", 0.5 for testing..
   A_step        = 0.1
   Mstar0_g_step = 0.25

   p1_norm       = 5.80    ;; 
   p1_step       = 0.1       ;; p1_bins = 10 

   p2_norm       = -3.30  ;; (for 20120920v1-v3) 
   p2_step       = 0.1     ;; use to be 0.1

print, 'means max(Mstar0_g) = ', Mstar0_g_norm  + ( Mstar0_g_bins * Mstar0_g_step)

    gamma_norm  =  -0.30
    gamma_step  =  0.1
endif

if coarse_or_fine eq 1 then begin
   A_norm        = -10.50
   alpha_norm    =  -4.50
   beta_norm     =  -2.00
   
   Mstar0_g_norm = -24.5    

   alpha_step    = 0.05      ;; 0.01 or 0.05 for "realz", 0.5 for testing..
   beta_step     = 0.05      ;; 0.01, 0.05 for "realz", 0.5 for testing..
   A_step        = 0.01
   Mstar0_g_step = 0.25
   p1_step       = 0.05
   p1_norm       = 1.30      ;; 
   p2_norm       = -0.82800
   p2_step       = 0.01
end



;;openw, 9, 'BOSS_mLDDE_temp.dat'

chi_sq_min = 100000.
xxx_0min   = 100000.

mag = Mg_2SLAQ

;;
;; Equation (11) of Croom et al. (2009b)
;;
;Mstar_g = -22.0
for p1_i = 0ll, p1_bins do begin
   p1 = p1_norm + (p1_i*p1_step)       
   
   for p2_i = 0ll, p2_bins do begin
      p2 = p2_norm + (p2_i*p2_step)
      
      print, ' p1_i, p2_i',  p1_i, p2_i, '   p1, p2, ', p1, p2
      
      for Mstar0_g_i = 0ll, Mstar0_g_bins-1 do begin
         Mstar0_g[Mstar0_g_i] = Mstar0_g_norm  + ( Mstar0_g_i * Mstar0_g_step)
;         print,  ' Mstar0_g_i,  Mstar0_g[Mstar0_g_i] ', Mstar0_g_i,  Mstar0_g[Mstar0_g_i] 
         
         dmag    = mag - Mstar0_g[Mstar0_g_i]     
;         print, ' Mstar_g ', Mstar_g
         
         ;; N.B. This is a different Mstar_g than the Mstar_g from
         ;; above!!! Grrr.....
         dmag_c  = mag - M_g_c 
         
         for gamma_i = 0ll, gamma_bins-1 do begin
            gamma[gamma_i] = gamma_norm + (gamma_i*gamma_step)
            ;; Equnation (18) of Croom et al. (2009b)
;         zc = zc_zero / ( 1. + (10^(0.4*gamma*dmag_c)))
            zc = zc_zero / ( 1. + (10^(0.4*gamma[gamma_i]*dmag_c)))
           ; print, 'gamma_i, gamma[gamma_i]', gamma_i, gamma[gamma_i]
            
            ;; Equnation (17) of Croom et al. (2009b)
            top_ed    = 2.*(1+zc)^(p1)
            bottom_ed = (((1+zz)/(1+zc))^(-1.*p1)) + (((1+zz)/(1+zc))^(-1.*p2))
            e_d = top_ed / bottom_ed 
            
            for alpha_i = 0LL, alpha_bins-1 do begin
               alpha[alpha_i] = alpha_norm + (alpha_i*alpha_step)
;         print, 'alpha_i, ', alpha_i 
               
               for beta_i = 0LL, beta_bins-1 do begin
                  beta[beta_i] = beta_norm + (beta_i*beta_step)
                  denom =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
                  
                  for A_i = 0L, A_bins-1 do begin
                     A[A_i] = 10^(A_norm+(A_i*A_step))
                     xxx = (A[A_i] * e_d ) / denom
;         Phi_mLDDE[alpha_i, beta_i, phi_i] = A[A_i]  / denom
                     
                     if xxx[0] lt xxx_0min then  xxx_0min = xxx[0] 
                                ;print, Phi_star[phi_i], alpha[alpha_i], beta[beta_i], log_phi_2SLAQ[0], xxx[0], sigma_Phi_2SLAQ[0], xxx_0min 
                     ;print, 'log_phi_2SLAQ[0], xxx[0], sigma_phi[0]', log_phi_2SLAQ[0], xxx[0], sigma_phi[0]
                     chi_sq_bin = ( (10^(log_phi_2SLAQ) -  xxx)^2) $
                                  /  (sigma_Phi^2)
                     
                     chi_sq = total(chi_sq_bin)
;         print,  phi_i, alpha_i, beta_i,  chi_sq,xxx[0:10]
;         print, phi_i, alpha_i, beta_i,  chi_sq
                     if chi_sq lt chi_sq_min then begin
                        print,  A_i, alpha_i,  beta_i, chi_sq_min, chi_sq,   $
                                (A_norm+(A_i*A_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step), $
                                p1, p2,  Mstar0_g[MSTAR0_G_I],  gamma[gamma_i], $
                                format='(i8,i8,i8,  f15.4,f15.4,  f9.4,f9.4,f9.4,  f9.4,f9.4, f9.4, f9.4 )' 
                        print
                        print
                        print, '10^(log_phi_2SLAQ[0]), xxx[0], (10^(log_phi_2SLAQ[0]) -  xxx[0])^2,  ', $
                             10^(log_phi_2SLAQ[0]), xxx[0], (10^(log_phi_2SLAQ[0]) -  xxx[0])^2,   sigma_Phi[0]^2
                     
                        chi_sq_min  = chi_sq 
                        alpha_i_min = alpha_i
                        beta_i_min  = beta_i
                        A_i_min     = A_i
                        MSTAR0_G_I_min = MSTAR0_G_I
                        p1_min    = p1
                        p2_min    = p2
                        gamma_min = gamma
                        xxx_min = xxx
                     endif ;;chi_sq_min
                     
                  endfor ;; A
               endfor    ;; beta
            endfor       ;; alpha
         endfor                 ;gamma
      endfor                    ;;Mg_star
   endfor                       ;; p2
endfor                          ;; p1


print
print
print, '    alpha =  ',  alpha[alpha_i_min]
print, '     beta =  ',  beta[beta_i_min]
print, ' M_star_g =  ',  Mstar0_g[ MSTAR0_G_I_min]
print, '    M_g_c =  ', M_g_c
print, '    gamma =  ', gamma
print, '  zc_zero =  ', zc_zero
print, '       p1 =  ', p1_min
print, '       p2 =  ', p2_min
print, '       A  =   ', 10^A[A_i]
print
print

   
;; 
;; 
;; 
;alpha_test =  -3.33
; beta_test =  -1.41
alpha_test =    alpha[alpha_i_min]
 beta_test =  beta[beta_i_min]
;Mstar0_g_test = -22.17
Mstar0_g_test = Mstar0_g[ MSTAR0_G_I_min]
M_g_c  = -23.90
;; as is the thing "gamma"
gamma      =   0.68 
gamma_test =   gamma

p1_test =  1.160   ; p1_i_min
p2_test = -0.228   ; p2_i_min
;; as is the thing "A"
A_test = 10^(-9.21)
A_test = A[A_i_min]  ;10^(-5.84)
z_test  = z_full
Mg_test = MI_BIN_FULL


dmag_c_test = dmag_c

zc_zero_test = zc_zero
zc_test = zc_zero_test / ( 1. + (10^(0.4*gamma_test*dmag_c_test)))

dmag_test = Mg_test - Mstar_g_test
  top_ed_test = 2.*(1+zc_test)^(p1_test)
bottom_e_test = (((1+z_test)/(1+zc))^(-1.*p1)) + (((1+z_test)/(1+zc))^(-1.*p2))
      e_d(ii,jj) = top_ed / bottom_ed 

denom_test =  ( (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) )

phi_best_fit = A_test * e_d_test / denom_test

chi_sq_bin_test = ( (10^(log_phi_2SLAQ) -  phi_best_fit)^2)  /  (sigma_Phi^2)


;close, 9


openw, 11, 'chi_sq_values_mLDDE_inside_mag_loop_temp.dat'
printf, 11, '# Normalization, A  Abs i-mag     Phi_BOSS_S82, sigma_Phi_BOSS_S82,    Phi_LDDE_model,  chi_sq_bin, chi_sq '
printf, 11, '# (-9.57+(phi_i*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]), sigma_Phi_BOSS_S82[w[kk]], Phi_LDDE[ii,kk, phi_i], chi_sq_bin, chi_sq '

openw, 10, 'chi_sq_values_mLDDE_temp.dat'
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
device, filename='QLF_iband_mLDDE_temp.ps', $
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



;; z = 0.49
z_bin_range_min = 0.30
z_bin_range_max = 0.68
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.14, 0.70, 0.35, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN,  xtickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq
;;
;; oplot the model fits....
oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color



;;  z  =  0.87
z_bin_range_min = 0.68
z_bin_range_max = 1.06
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.35, 0.70, 0.56, 0.98], $
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



;;  z  =  1.25
z_bin_range_min = 1.06
z_bin_range_max = 1.44
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.70, 0.77, 0.98], $
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

;; z = 1.63
z_bin_range_min = 1.44
z_bin_range_max = 1.82
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


;;;;;;;;;      z = 2.01       
z_bin_range_min = 1.82
z_bin_range_max = 2.20
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.14, 0.42, 0.35, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_best_fit[w_plot]), thick=3, color=line_color


;;   z = 2.40     !!!  F i r s t   BOSS  BIN  !!!
z_bin_range_min = 2.20
z_bin_range_max = 2.60
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
z_bin_range_min = 2.60
z_bin_range_max = 3.00
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
z_bin_range_min = 3.00
z_bin_range_max = 3.50
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


;;  z = 3.75
z_bin_range_min = 3.50
z_bin_range_max = 4.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.14, 0.14, 0.35, 0.42], $
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


;;  z = 4.25
z_bin_range_min = 4.00
z_bin_range_max = 4.50
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
z_bin_range_min = 4.50
z_bin_range_max = 5.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.14, 0.77, 0.42], $
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


device, /close
set_plot, 'X'
close, /all

!p.multi=0

end
