;+
; NAME:
;   chi_sq_PLE_plot
;
; PURPOSE:
;   Okay, so this is just a fairly wee, fairly simple script that
;   takes the various "PLE models from Croom et al. (2009b), 
;   and overplots to the SDSS+boss21+BOSS dataset
;
; CALLING SEQUENCE:
;   .run chi_sq_PLE_plot
;
; INPUTS:
;
; OUTPUTS:
;   .ps 3x4 plot... 
;
; NOTES:
;
;-



;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent


;;
;;   B O S S  +  S D S S
;;

;readcol, '../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_20120731.dat', $
readcol, '../../pro/qlf/QLF_iband_sdss_S82_boss21_wSelfn_formatted_McGkcor.dat', $
         z_full,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, $ 
         NQ_full, ba,  log_Phi_full, bar, error_full, survey, $
         format='(d,a, d,a, d,a,  d,a, d,a, d, a)'

;w_zbin = where(z_2SLAQ_full ge 0.40 and Z_2SLAQ_full le 2.60, N)
w_zbin = where(z_full ge 0.20 and Z_full le 2.20, N)
;w_zbin = where(z_full ge 1.0 and Z_full le 2.20, N)
;w_zbin = where(z_full ge 2.20 and Z_full le 2.30, N)


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
;; For plotting in the 3x4 panel and the 
;;
;;  R E D S H I F T  -   L O G _ P H I  plane
;;
;; A list of given models...
;;
fltarr_bins = 150
z_range_a = fltarr(fltarr_bins)+0.49
z_range_b = fltarr(fltarr_bins)+0.89
z_range_c = fltarr(fltarr_bins)+1.25
z_range_d = fltarr(fltarr_bins)+1.63

z_range_e = fltarr(fltarr_bins)+2.01
z_range_f = fltarr(fltarr_bins)+2.35
z_range_g = fltarr(fltarr_bins)+2.70
z_range_h = fltarr(fltarr_bins)+3.25

z_range_i = fltarr(fltarr_bins)+3.75
z_range_j = fltarr(fltarr_bins)+4.25
z_range_k = fltarr(fltarr_bins)+4.75

z_range = [z_range_a, z_range_b, z_range_c, z_range_d, z_range_e, z_range_f, z_range_g, z_range_h, z_range_i, z_range_j, z_range_k]

M_range = (findgen(fltarr_bins)/10.)-33.                   
M_range_full = [M_range, M_range, M_range, M_range, M_range, M_range, M_range, M_range, M_range, M_range, M_range]


;;
;;   M O D E L     M c G R E E R    P L E       z = 0.30 - 2.20 
;; 
   alpha_test =  -1.44
    beta_test =  -3.83
Mstar0_g_test = -23.45
      k1_test =   1.259
      k2_test =  -0.273
Phi_star_test =  10^(-6.28)
       z_test = z_full
      Mg_test = MI_BIN_FULL

Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_test + k2_test*z_test*z_test)
dmag_test      = Mg_test - Mstar_g_test
ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
phi_fit_a      = Phi_star_test /ple_denom_test

Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_range  + k2_test*z_range*z_range)
dmag_test      = M_range_full - Mstar_g_test
ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
phi_fit_a_full = Phi_star_test /ple_denom_test

starting_redshift = 0.30

        z_out = ((findgen(fltarr_bins))/25.)+starting_redshift
  M_range_out = (findgen(40)*0.30)-31.05
;; 150 * 40 gives me my array of 6000 lines...
phi_fit_a_out = fltarr(fltarr_bins,40)
  Mstar_g_out = fltarr(fltarr_bins)

openw, 10, 'chi_sq_PLE_model_McG_0.3z2.2_temp.dat' 
for ii=0LL, n_elements(z_out)-1 do begin
      Mstar_g_out[ii]   = Mstar0_g_test  - 2.5*(k1_test*z_out[ii]  + k2_test*z_out[ii]*z_out[ii])
   for jj=0LL, n_elements(M_range_out)-1 do begin
      dmag_out      = M_range_out[jj] - Mstar_g_out[ii]
      ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_out))) + (10.0^(0.4*(beta_test+1)*(dmag_out))) 
      phi_fit_a_out[ii,jj] = Phi_star_test /ple_denom_test
      printf, 10, Z_out[ii], M_range_out[jj], phi_fit_a_out[ii,jj]
   endfor
endfor
close, 10


;;
;;   M O D E L     M c G R E E R    P L E       z = 2.20 - 3.50
;; 
   alpha_test =  -1.53
    beta_test =  -3.05
Mstar0_g_test = -24.03
      k1_test =  1.152   
      k2_test = -0.267  
Phi_star_test =  10^(-6.30)
       z_test = z_full
      Mg_test = MI_BIN_FULL

Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_test + k2_test*z_test*z_test)
dmag_test      = Mg_test - Mstar_g_test
ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
phi_fit_b   = Phi_star_test /ple_denom_test

Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_range + k2_test*z_range*z_range)
dmag_test      = M_range_full - Mstar_g_test
ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
phi_fit_b_full = Phi_star_test /ple_denom_test

        z_out = ((findgen(fltarr_bins))/25.)+starting_redshift
  M_range_out = (findgen(40)*0.30)-31.05
phi_fit_b_out = fltarr(fltarr_bins,40)
  Mstar_g_out = fltarr(fltarr_bins)

openw, 10, 'chi_sq_PLE_model_McG_2.2z3.5_temp.dat' 
for ii=0LL, n_elements(z_out)-1 do begin
      Mstar_g_out[ii]   = Mstar0_g_test  - 2.5*(k1_test*z_out[ii]  + k2_test*z_out[ii]*z_out[ii])
   for jj=0LL, n_elements(M_range_out)-1 do begin
      dmag_out      = M_range_out[jj] - Mstar_g_out[ii]
      ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_out))) + (10.0^(0.4*(beta_test+1)*(dmag_out))) 
      phi_fit_b_out[ii,jj] = Phi_star_test /ple_denom_test
      printf, 10, Z_out[ii], M_range_out[jj], phi_fit_b_out[ii,jj]
   endfor
endfor
close, 10



;;
;;   M O D E L     M c G R E E R    ``L I N E A R''    z = 2.20 - 3.50
;; 
    alpha_test =  -1.40
     beta_test =  -3.57
 Mstar0_g_test = -26.78
       c1_test =  -0.606
       c2_test =  -0.505
;Phi_star0_test =  10^(-6.05)
Phi_star0_test =  -6.05
        z_test = z_full
       Mg_test = MI_BIN_FULL
       z_pivot = 2.20

;;log(Phi*) = Phi*(0) + k1*(z-2.2)
;;M* = M*(0) + k2*(z-2.2)

Mstar_g_test   = Mstar0_g_test  + (c2_test*(z_test-z_pivot))
dmag_test      = Mg_test - Mstar_g_test
ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 

 Phi_star_test = Phi_star0_test + (c1_test * (z_test-z_pivot))
  phi_fit_out  = Phi_star_test /ple_denom_test

starting_redshift = 0.30
        z_out = ((findgen(fltarr_bins))/25.)+starting_redshift
  Mstar_g_out = fltarr(fltarr_bins)
 Phi_star_out = fltarr(fltarr_bins)

;; 150 * 40 gives me my array of 6000 lines...
  M_range_out = (findgen(40)*0.30)-31.05

  phi_fit_out = fltarr(fltarr_bins,40)
  
openw, 15, 'chi_sq_linearPhi_model_McG_2.2z3.5_temp.dat' 
for ii=0LL, n_elements(z_out)-1 do begin
   Phi_star_out[ii] = Phi_star0_test + (c1_test * (z_out[ii] - z_pivot))
    Mstar_g_out[ii] = Mstar0_g_test  + (c2_test * (z_out[ii] - z_pivot))

   for jj=0LL, n_elements(M_range_out)-1 do begin
      dmag_out      = M_range_out[jj] - Mstar_g_out[ii]
      ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_out))) + (10.0^(0.4*(beta_test+1)*(dmag_out))) 
      phi_fit_out[ii,jj] = (10^(Phi_star_out[ii])) /ple_denom_test
      printf, 15, Z_out[ii], M_range_out[jj], phi_fit_out[ii,jj], 10^Phi_star_out[ii]
   endfor
endfor
close, 15




;;
;;   M O D E L    m L D D E
;;
;;   Never implememnted; just went back to the model_fits.pro
;;   of the old Croom09b basis...
;; 
;   alpha_test =  -3.70
;    beta_test =  -2.34
;Mstar0_g_test = -22.75
;      k1_test =  1.60   
;Phi_star_test =  10^(-6.32)
;       z_test = z_full
;      Mg_test = MI_BIN_FULL

;;Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_test + k2_test*z_test*z_test)
;;dmag_test      = Mg_test - Mstar_g_test
;;ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
;;phi_fit_a      = Phi_star_test /ple_denom_test;
;
;;Mstar_g_test   = Mstar0_g_test  - 2.5*(k1_test*z_range  + k2_test*z_range*z_range)
;;dmag_test      = M_range_full - Mstar_g_test
;;ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_test))) + (10.0^(0.4*(beta_test+1)*(dmag_test))) 
;;phi_fit_a_full = Phi_star_test /ple_denom_test

;
;starting_redshift = 0.30
;z_out         = ((findgen(fltarr_bins))/25.)+starting_redshift
;M_range_out   = (findgen(40)*0.30)-31.05
;; 150 * 40 gives me my array of 6000 lines...
;phi_fit_a_out = fltarr(fltarr_bins,40)
;Mstar_g_out   = fltarr(fltarr_bins)

;openw, 11, 'chi_sq_mLDDE_model_temp.dat' 
;for ii=0LL, n_elements(z_out)-1 do begin
 ;     Mstar_g_out[ii]   = Mstar0_g_test  - 2.5*(k1_test*z_out[ii]  + k2_test*z_out[ii]*z_out[ii])
  ; for jj=0LL, n_elements(M_range_out)-1 do begin
   ;   dmag_out      = M_range_out[jj] - Mstar_g_out[ii]
    ;  ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_out))) + (10.0^(0.4*(beta_test+1)*(dmag_out))) 
     ; phi_fit_a_out[ii,jj] = Phi_star_test /ple_denom_test
     ; printf, 11, Z_out[ii], M_range_out[jj], phi_fit_a_out[ii,jj]
   ;endfor
;endfor
;close, 11





;;
;; 
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
light_blue =  88
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

;; xy_outs 
x_xyouts = -25.00
y_xyouts =  -5.05

;;
;; !p.multi=[0,4,0]  ==> 1 row, 4 columns....
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_iband_PLE_plot_temp.ps', $
        xsize=18.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated
color_2slaq = red
;;plot_sym_size_BOSS = 2.2
sym_size = 2.4
line_thick = 12.0
err_thick = 6.0 

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

;;
;;  z = 0.49
;;
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
            /hibar, errthick=err_thick, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=err_thick, psym=8, color=color_2slaq, errcolor=color_2slaq
;; oplot the model fits....
oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,   linestyle=0
oplot, M_range,          alog10(phi_fit_a_full),    thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range,          alog10(phi_fit_b_full), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range,          alog10(phi_fit_c_full), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_d_full), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full), thick=line_thick/2., color=orange,linestyle=1

xyouts, x_xyouts, y_xyouts, '!80.30<z<0.68!3', charsize=2.2, charthick=6



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
;; oplotting the models..
oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]),       thick=line_thick,    color=black, linestyle=0
oplot, M_range,          alog10(phi_fit_a_full[150:299]), thick=line_thick/2., color=black, linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]),       thick=line_thick,    color=purple, linestyle=0 
oplot, M_range,          alog10(phi_fit_b_full[150:299]), thick=line_thick/2., color=purple, linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]),       thick=line_thick,    color=light_blue, linestyle=1 
;oplot, M_range,          alog10(phi_fit_c_full[150:299]), thick=line_thick/2., color=light_blue, linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]),       thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_d_full[150:299]), thick=line_thick/2., color=orange, linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),       thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[150:299]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!80.68<z<1.06!3', charsize=2.2, charthick=6



;;  z  =  1.25
z_bin_range_min = 1.06
z_bin_range_max = 1.44
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.70, 0.77, 0.98], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq
;; oplotting the models...
oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]),       thick=line_thick, color=black,  linestyle=0
oplot, M_range,          alog10(phi_fit_a_full[300:449]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]),       thick=line_thick, color=purple, linestyle=0 
oplot, M_range,          alog10(phi_fit_b_full[300:449]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]),       thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range,          alog10(phi_fit_c_full[300:449]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]),       thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_d_full[300:449]), thick=line_thick/2., color=orange, linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),       thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[300:449]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!81.06<z<1.44!3', charsize=2.2, charthick=6



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
;; oplotting the models...
oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[450:599]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[450:599]), thick=line_thick/2., color=black,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[450:599]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[450:599]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_e_full[450:599]), thick=line_thick/2., color=orange,linestyle=1

xyouts, x_xyouts, y_xyouts, '!81.44<z<1.82!3', charsize=2.2, charthick=6



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
;;
oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]),       thick=line_thick,    color=black, linestyle=0
oplot, M_range,          alog10(phi_fit_a_full[600:749]), thick=line_thick/2., color=black, linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[600:749]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[600:749]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[600:749]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),       thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[600:749]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!81.82<z<2.20!3', charsize=2.2, charthick=6


;;
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

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[750:899]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[750:899]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[750:899]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[750:899]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_e_full[750:899]), thick=line_thick/2., color=orange,linestyle=1

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.60!3', charsize=2.2, charthick=6



;;  z = 2.80
z_bin_range_min = 2.60
z_bin_range_max = 3.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.42, 0.77, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[900:1049]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[900:1049]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[900:1049]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[900:1049]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),        thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[900:1049]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!82.60<z<3.00!3', charsize=2.2, charthick=6



;;  z = 3.25
z_bin_range_min = 3.00
z_bin_range_max = 3.50
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.77, 0.42, 0.98, 0.70], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, xtickformat='(a1)', ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[1050:1200]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[1050:1200]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[1050:1200]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[1050:1200]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),         thick=line_thick,    color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[1050:1200]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!83.00<z<3.50!3', charsize=2.2, charthick=6



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

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[1200:1349]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[1200:1349]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[1200:1349]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[1200:1349]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[1200:1349]), thick=line_thick/2., color=orange,linestyle=1

xyouts, x_xyouts, y_xyouts, '!83.50<z<4.00!3', charsize=2.2, charthick=6




;;  z = 4.25
z_bin_range_min = 4.00
z_bin_range_max = 4.50
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
       position=[0.35, 0.14, 0.56, 0.42], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[1350:1499]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[1350:1499]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1
;oplot, M_range, alog10(phi_fit_c_full[1350:1499]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1
;oplot, M_range, alog10(phi_fit_d_full[1350:1499]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]),         thick=line_thick,    color=orange, linestyle=1
;oplot, M_range,          alog10(phi_fit_e_full[1350:1499]), thick=line_thick/2., color=orange, linestyle=1

xyouts, x_xyouts, y_xyouts, '!84.00<z<4.50!3', charsize=2.2, charthick=6



;;  z = 4.75
z_bin_range_min = 4.50
z_bin_range_max = 5.00
w_plot = where(z_2SLAQ ge z_bin_range_min and z_2SLAQ le z_bin_range_max and Nq gt 1, N_qlf_bins)

plotsym, 0, sym_size, /fill
plot, MG_2SLAQ[w_plot] , LOG_PHI_2SLAQ[w_plot], psym=8, $
      position=[0.56, 0.14, 0.77, 0.42], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, xthick=xthick/1.4, ythick=ythick/1.4, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, ytickformat='(a1)', $
      charsize=charsize, charthick=charthick/1.4, thick=thick/1.4 

M_2SLAQ_err = MG_2SLAQ[w_plot] - MG_2SLAQ[w_plot]
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_up[w_plot], $
            /hibar, errthick=4, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, MG_2SLAQ[w_plot], log_phi_2SLAQ[w_plot], M_2SLAQ_err, error_delta_down[w_plot], $
            /lobar, errthick=4, psym=8, color=color_2slaq, errcolor=color_2slaq

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_a[w_plot]), thick=line_thick, color=black,  linestyle=0
oplot, M_range, alog10(phi_fit_a_full[1500:1649]), thick=line_thick/2., color=black,linestyle=0

oplot, MG_2SLAQ[w_plot], alog10(phi_fit_b[w_plot]), thick=line_thick, color=purple, linestyle=0 
oplot, M_range, alog10(phi_fit_b_full[1500:1649]), thick=line_thick/2., color=purple,linestyle=0

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_c[w_plot]), thick=line_thick, color=light_blue, linestyle=1 
;oplot, M_range, alog10(phi_fit_c_full[1500:1649]), thick=line_thick/2., color=light_blue,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_d[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range, alog10(phi_fit_d_full[1500:1649]), thick=line_thick/2., color=orange,linestyle=1

;oplot, MG_2SLAQ[w_plot], alog10(phi_fit_e[w_plot]), thick=line_thick, color=orange, linestyle=1 
;oplot, M_range,          alog10(phi_fit_e_full[1500:1649]), thick=line_thick/2., color=orange,linestyle=1

xyouts, x_xyouts, y_xyouts, '!84.50<z<5.00!3', charsize=2.2, charthick=6




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  G L O B A L   A X E S 
;;
;; mind that the positions are tied to the ranges 
;; of the final panel...

xyouts, -2.5, -7.0,  $
        'log!I10!N !7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
        ORIENTATION=90, charsize=charsize, $
       charthick=charthick*1.8


xyouts, -15., -11.25,  $
        'M!Ii!N[z=2]', $
       charthick=charthick*1.6, charsize=charsize


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick=4.8
charsize=2.0

xpos     = -30.60 
;xpos_off = 
ypos     =  -4.60
ypos_off =  -0.60

;purple     =  32
;deep_blue  =  48
;blue       =  64
;light_blue =  80
;turquiose  = 128
;green      = 150
;yellow     = 210
;orange     = 232
;red    

loadct, clr_table
legend, ' ' , $
        position=[xpos, ypos], box=0, linestyle=0, color=black, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(1.*ypos_off)], box=0, linestyle=2, color=purple-8, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(2.*ypos_off)], box=0, linestyle=1, color=deep_blue-8, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(3.*ypos_off)], box=0, linestyle=4, color=blue, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(4.*ypos_off)], box=0, linestyle=3, color=light_blue+8, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(5.*ypos_off)], box=0, linestyle=0, color=turquiose, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(6.*ypos_off)], box=0, linestyle=1, color=green, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(7.*ypos_off)], box=0, linestyle=6, color=orange, $
        charsize=charsize/1.4, thick=20.
legend, ' ' , $
        position=[xpos, ypos+(8.*ypos_off)], box=0, linestyle=5, color=red, $
        charsize=charsize/1.4, thick=20.



x_xyouts  = -34.40

y_xyouts  = -5.00
y_offset  = -0.60

print, 'charsize  ', charsize, '  charthick', charthick
print
charsize  = charsize / 1.4
charthick = charsize / 1.2
print
print, 'charsize  ', charsize, '  charthick', charthick
print

;xyouts, x_xyouts, y_xyouts, 'C09b',  color=blackl, $
xyouts, x_xyouts, y_xyouts, 'PLE 1.0<z<2.2',  color=blackl, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(1.*y_offset), 'PLE 0.3<z<2.2', color=purple-8, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(2.*y_offset), 'PLE 0.3<z<2.6', color=deep_blue-8, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(3.*y_offset), 'PLE sqrt', color=blue, $
        charsize=charsize*1.4, charthick =charthick*2.4
   
xyouts, x_xyouts, y_xyouts+(4.*y_offset), 'PLE sqrt', color=light_blue+8, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(5.*y_offset), 'PLE sqrt', color=turquiose, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(6.*y_offset), 'HRH07 PLE', color=green, $
        charsize=charsize*1.4, charthick =charthick*2.4

;xyouts, x_xyouts, y_xyouts+(7.*y_offset), 'HRH07 LDDE', color=yellow, $
 ;       charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(7.*y_offset), 'PLE sqrt', color=orange, $
        charsize=charsize*1.4, charthick =charthick*2.4

xyouts, x_xyouts, y_xyouts+(8.*y_offset), 'HRH07 ``full''', color=red, $
        charsize=charsize*1.4, charthick =charthick*2.4


loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all


;;
;; e.g. Line 1 from Table 2 of Boyle00
;;  http://www.astro.rug.nl/~onderwys/ACTUEELONDERZOEK/JAAR2003/college2/Cosmology.html
;;
;red, omega0=0.30, omegalambda=0.70, h100=0.70    
;red, omega0=1.00, omegalambda=0.00, h100=0.70    

;   alpha_test =  -3.00 ;-3.37
;    beta_test =  -1.50 ;-1.55
;Mstar0_r_test = -23.75 ;-24.00  ;; guessed conversion from Mstar_B = -21.16
;  k1_tau_test =  5.00 ; 7.11   
;Phi_star_test = 10^(-6.2700) ;0.47e-6
;       z_test = z_full
;      Mr_test = MI_BIN_FULL
;     tau_test = getage(z_test)






end
