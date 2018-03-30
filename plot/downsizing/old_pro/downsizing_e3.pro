;+
; NAME:
;       downsizing
; 
; PURPOSE:
;        To plot Quasar Luminosity Functions, as several M_i lines on 
;        the redshift-number density plane, a la e.g. Croom et
;        al. 2009b. This is trying to show you, visually the 
;        (so-called) ``AGN downsizing''.
;
; CALLING SEQUENCE:
;       .run downsizing
;
; INPUTS:
;        Richards06_Table06.dat   - QLF from Richards06, z=2.01 bin only
;
; OUTPUTS:
;       .ps file
;
; COMMENTS:
;       /usr/common/rsi/lib/general/LibAstro/ 
;         
;-

print
print
print, '=========================================================='
print, '== '
print, '==  red, omega0=0.30, omegalambda=0.70, h100=0.700 '
print, '=='
print, '=========================================================='
print
print

red, omega0=0.30, omegalambda=0.70, h100=0.700

;; choice of "every other" line...
choice_eo = 'y'
;; if choice_eo eq 'n' then 

choice_line = 'n' 
;read, choice_line, PROMPT=' - Plot a line for the data??  y/n  '


choice_models    = 'n'

choice_plot_PLE  = 'n'
;read, choice_plot_PLE, PROMPT=' - Plot PLE models?? y/n  '

choice_plot_LDDE = 'n'
;read, choice_plot_LDDE, PROMPT=' - Plot LDDE models?? y/n  '

choice_plot_LEDE = 'n'
;read, choice_plot_LEDE, PROMPT=' - Plot LEDE models?? y/n  '

;; want something a bit shorter...
plotPLE = 'y' 


plot_log_age = 0
read, plot_log_age, PROMPT='  Plot in (0) redshift  or (1) log10(age/Gyr)?? y/n  '



;;
;;   B O S S  +  S D S S
;;
;readcol, '../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_20120731.dat', $
;readcol, '../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_McGkcor.dat', $
readcol, '../../pro/qlf/QLF_iband_sdss_S82_narrow_boss21_wSelfn_formatted_McGkcor.dat', $
         z_input,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, $ 
         NQ_full, ba,  log_Phi_full, bar, error_full, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)', /silent

 sigma_Phi = error_full*(1e-9)
  delta_up = alog10 ((10^log_Phi_full + sigma_Phi))
        up = delta_up       - log_Phi_full
delta_down = alog10 ((10^log_Phi_full - sigma_Phi))
      down = abs(delta_down - log_Phi_full)

z_full = z_input
if plot_log_age eq 0 then z_full = z_input

zz = getage(z_full)
;z_z = zz
;z_full = alog(zz+1.)
;z_full = alog10(zz+1.)
;z_full = alog10(VCOMOVING(z_input)/1e27)
;z_full = alog10((1+z_input)^3.)

if plot_log_age eq 1 then z_full = alog10(zz)   ;; Works well....



;;
;;  M o d e l   f i t s...
;;
if (choice_models eq 'y') then begin
   readcol, '../../pro/models/Croom09b_PLE_toz5.dat', $
            ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   LOG_PHI_PLE = alog10(PHI_PLE)
      
   readcol, '../../pro/models/Croom09b_mLDDE_toz5.dat', $
            ii_LDDE, z_LDDE, jj_LDDE, mag_LDDE, e_d_LDDE, zc_LDDE, Phi_LDDE
   LOG_PHI_LDDE = alog10(PHI_LDDE)
   
   readcol, '../../pro/models/Croom09b_LEDE_toz5.dat', $
            ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
   LOG_PHI_LEDE = alog10(PHI_LEDE)

endif

;;
;;  Best-fit PLE model (aka "model A"...)
;;
;readcol, '../../pro/models/chi_sq_PLE_model_a_20120810.dat', $
;readcol, '../../pro/models/chi_sq_PLE_model_a_20120914.dat', $
;readcol, '../../pro/models/chi_sq_PLE_model_a_20120919.dat', $
readcol, '../../pro/models/chi_sq_PLE_model_McG_0.3z2.2_temp.dat', $
         z_PLE, Mi_PLE, Phi_PLE

w_PLE_range  = where(z_PLE ge 0.0 and  z_PLE le 2.05, N)
  z_PLE =  z_PLE[w_PLE_range]
 Mi_PLE = Mi_PLE[w_PLE_range]
Phi_PLE = Phi_PLE[w_PLE_range]

LOG_PHI_PLE = alog10(PHI_PLE)
zz_ple = getage(z_ple)
if plot_log_age eq 1 then z_ple = alog10(zz_ple)


;; 
;;   L  E  D  E 
;;
;readcol, '../../pro/models/chi_sq_linearPhi_model_McG_2.2z3.5_alpha1.38.dat', $
readcol, '../../pro/models/chi_sq_linearPhi_model_McG_2.2z3.5_temp.dat', $
         z_LEDE, Mi_LEDE, Phi_LEDE

w_LEDE_range  = where( z_LEDE ge 1.98 and  z_LEDE le 4.0, N)
  z_LEDE = z_LEDE[w_LEDE_range]
 Mi_LEDE = Mi_LEDE[w_LEDE_range]
Phi_LEDE = Phi_LEDE[w_LEDE_range]

LOG_PHI_LEDE = alog10(Phi_LEDE)
     zz_LEDE = getage(z_LEDE)
if plot_log_age eq 1 then  z_LEDE = alog10(zz_LEDE)

   
;;
;;  The modified LDDE model from 2SLAQ...
;;
;;  ( Croom09b_mLDDE_20120920.dat produced to match 
;;    the binning of chi_sq_PLE_model_a_20120919.dat...)
;;
;readcol, '../../pro/models/Croom09b_mLDDE_20120920.dat', $
;readcol, '../../pro/models/Croom09b_mLDDE_20120920v2.dat', $
;readcol, '../../pro/models/Croom09b_mLDDE_20120920v3.dat', $
;readcol, '../../pro/models/Croom09b_mLDDE_zc_zero4.47_temp.dat', $
;readcol, '../../pro/models/Croom09b_mLDDE_zc_zero2.47_gamma1.68.dat',$
;readcol, '../../pro/models/Croom09b_mLDDE_zc_zero2.47_gamma0.68_M_g_c_25.9.dat', $
;readcol, '../../pro/models/Croom09b_mLDDE_temp.dat', $
readcol, '../../pro/models/Bongiorno07_VVDS_LDDE_temp.dat', $
         z_mLDDE, Mi_mLDDE, Phi_mLDDE

LOG_PHI_MLDDE = alog10(Phi_mLDDE)
zz_mLDDE = getage(z_mLDDE)
if plot_log_age eq 1 then z_mLDDE = alog10(zz_mLDDE)



;;
;;  The Conroy & White (2012) model...
;;
readcol, '../../pro/models/CW12/CW12_qsolf_allz_temp.dat', $
         z_CW12, Mi_CW12, log_Phi_CW12

zz_CW12 = getage(z_CW12)
if plot_log_age eq 1 then z_CW12 = alog10(zz_CW12)




; print, minmax(MI_BIN_FULL)         
;      -29.550000      -22.950000

delta_bin = 0.30 
delta = delta_bin/2.
;binranges=[-29.85, -29.55,-29.25, -28.95, -28.65, -28.35, -28.05, -27.75, -27.45, -27.15, -26.85, -26.55, -26.25, -25.95, -25.65, -25.35, -25.05, -24.75, -24.45, -24.15, -23.85, -23.55, -23.25, -22.95,  -22.65, -22.35, -22.05]
binranges=[-29.55,-29.25, -28.95, -28.65, -28.35, -28.05, -27.75, -27.45, -27.15, -26.85, -26.55, -26.25, -25.95, -25.65, -25.35, -25.05, -24.75, -24.45, -24.15, -23.85, -23.55, -23.25, -22.95, -22.65]

binranges = binranges + delta


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

expand = 1.5

charsize  = 4.2*expand
charthick = 10.2*expand
thick     = 10.2 ;; doesn't do anything!!
xthick    = 20.0
ythick    = 20.0
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; x-ranges
;; for z, redshift...
xmin = 0.0
xmax = 5.0
xstyle = 1

;; For Age of the Universe... 
;xmin = 14.0
;xmax =  0.0

;;  in alog10 (age  +1)
if plot_log_age eq 1 then begin
   xmin =  1.25
   xmax =  0.0
   xstyle = 8+1
endif

;;  for = alog10(VCOMOVING(z_input)/1e27)
;xmin =  1.25
;xmax =  3.50

;; for (1+z)^3
;xmin =    0.00
;xmax =  200.

;; for alog10(1+z)^3
;xmin =  0.50
;xmax =  2.30


;; y-ranges
ymin = -9.0   ;-10.0
ymax =  -5.5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  
;;
set_plot, 'ps'
device, filename='downsizing_SDSS_BOSS_boss21_temp_e3.eps', $
;        xsize=(14.0*expand), ysize=16.0*(expand), $
        xsize=(22.0*expand), ysize=17.0*(expand), $
;        xsize=22.0, ysize=16.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated
print, 'charthick', charthick

plotsym, 0, 1.4, /fill
plot, z_full, log_phi_full, $
      psym=8, $
      position=[0.16, 0.12, 0.96, 0.88], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $;
;      xtitle=' redshift ', $
;      xtitle='log!I10!N(age / Gyr)', $
      ytitle='log!I10!N !7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
      color=black

thick=thick*2.0
clr_step  = 18
linethick = 30

;;
;;  M O D E L S
;;

;;
;;   C O N R OY    A N D    W H I T E    2 0 1 2 
;;
linestyle = 2
choice_plot_CW12 = 'n'

ep = 0.139

if choice_plot_CW12 eq 'y' then begin

   w_ple = where(Mi_CW12      le binranges[3]-ep and Mi_CW12      gt binranges[2]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*0), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[4]-ep and Mi_CW12      gt binranges[3]+.125, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*1), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[5]-ep and Mi_CW12      gt binranges[4]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*2), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[6]-ep and Mi_CW12      gt binranges[5]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*3), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[7]-ep and Mi_CW12      gt binranges[6]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*4), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[8]-ep and Mi_CW12      gt binranges[7]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*5), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[9]-ep and Mi_CW12      gt binranges[8]+.125, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*6), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[10]-ep and Mi_CW12      gt binranges[9]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*7), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[11]-ep and Mi_CW12      gt binranges[10]+.125, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*8), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[12]-ep and Mi_CW12      gt binranges[11]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*9), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[13]-ep and Mi_CW12      gt binranges[12]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*10), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[14]-ep and Mi_CW12      gt binranges[13]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*11), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[15]-ep and Mi_CW12      gt binranges[14]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*12), thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[16]-ep and Mi_CW12      gt binranges[15]+.125, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*13)+2, thick=linethick, linestyle=linestyle
   w_ple = where(Mi_CW12      le binranges[17]-ep and Mi_CW12      gt binranges[16]+ep, N_ple)
   oplot, z_CW12[w_ple], log_phi_CW12[w_ple], color=black+(clr_step*14)+4, thick=linethick, linestyle=linestyle
endif


;;
;;  P L E 
;;

;thick=thick*0.5
choice_plot_PLE = 'y'
ple_linestyle = 0

if choice_plot_PLE eq 'y' then begin
   w_ple = where(Mi_PLE      le binranges[4] and Mi_PLE      gt binranges[3], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*1), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[6] and Mi_PLE      gt binranges[5], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*3), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[8] and Mi_PLE      gt binranges[7], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*5), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[10] and Mi_PLE      gt binranges[9], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*7), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[12] and Mi_PLE      gt binranges[11], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*9), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[14] and Mi_PLE      gt binranges[13], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*11), thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[16] and Mi_PLE      gt binranges[15], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*13)+2, thick=linethick, linestyle=ple_linestyle
   w_ple = where(Mi_PLE      le binranges[18] and Mi_PLE      gt binranges[17], N_ple)
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*15)+2, thick=linethick, linestyle=ple_linestyle
endif


;;
;;  L E D E 
;;
;thick=thick*0.5
choice_plot_LEDE = 'y'
linestyle = 0

if choice_plot_LEDE eq 'y' then begin

   w_lede = where(Mi_LEDE le binranges[4] and Mi_LEDE      gt binranges[3], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*1), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[6] and Mi_LEDE      gt binranges[5], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*3), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[8] and Mi_LEDE      gt binranges[7], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*5), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[10] and Mi_LEDE      gt binranges[9], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*7), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[12] and Mi_LEDE      gt binranges[11], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*9), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[14] and Mi_LEDE      gt binranges[13], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*11), thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[16] and Mi_LEDE      gt binranges[15], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*13)+2, thick=linethick, linestyle=linestyle
   w_lede = where(Mi_LEDE      le binranges[18] and Mi_LEDE      gt binranges[17], N_lede)
   oplot, z_LEDE[w_lede], log_phi_LEDE[w_lede], color=black+(clr_step*15)+4, thick=linethick, linestyle=linestyle
endif



;;
;;  m  L D D E 
;;
;thick=thick*0.5
choice_plot_MLDDE = 'y'
linestyle = 5

if choice_plot_MLDDE eq 'y' then begin

   w_mldde = where(Mi_MLDDE le binranges[4] and Mi_MLDDE      gt binranges[3], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*1), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[6] and Mi_MLDDE      gt binranges[5], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*3), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[8] and Mi_MLDDE      gt binranges[7], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*5), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[10] and Mi_MLDDE      gt binranges[9], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*7), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[12] and Mi_MLDDE      gt binranges[11], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*9), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[14] and Mi_MLDDE      gt binranges[13], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*11), thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[16] and Mi_MLDDE      gt binranges[15], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*13)+2, thick=linethick, linestyle=linestyle
   w_mldde = where(Mi_MLDDE      le binranges[18] and Mi_MLDDE      gt binranges[17], N_mldde)
   oplot, z_MLDDE[w_mldde], log_phi_MLDDE[w_mldde], color=black+(clr_step*15)+4, thick=linethick, linestyle=linestyle
endif



;;
;;   D  A  T  A
;;
solid  = 0.8

plotsym, 0, 5.0, /fill
errthick = 36.


w = where(Mi_bin_full le binranges[2] and Mi_bin_full gt binranges[1], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*0), thick=thick
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)

   print, N, mean(Mi_bin_full[w]), 0.5*(binranges[2] +binranges[1] )
endif

w = where(Mi_bin_full le binranges[5] and Mi_bin_full gt binranges[4], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*2), thick=linethick*solid
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[8] and Mi_bin_full gt binranges[7], N)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*3), psym=8
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*3), ERRcolor=black+(clr_step*3)

   if choice_line eq 'y' then oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*3), thick=linethick*solid
   print, N, mean(Mi_bin_full[w]), ' le binranges[6] and Mi_bin_full gt binranges[5] '
endif

w  = where(Mi_bin_full le binranges[11] and Mi_bin_full gt binranges[10], N)
if N ge 5 then begin 
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*5), psym=8
   
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*5), ERRcolor=black+(clr_step*5)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*5), ERRcolor=black+(clr_step*5)
   
   if choice_line eq 'y' then oplot, z_full[w],    log_phi_full[w], color=black+(clr_step*5), thick=linethick*solid
   print, N, mean(Mi_bin_full[w])
endif

w  = where(Mi_bin_full le binranges[14] and Mi_bin_full gt binranges[13], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*7), psym=8
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*7), ERRcolor=black+(clr_step*7)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*7), ERRcolor=black+(clr_step*7)

   if choice_line eq 'y' then oplot, z_full[w], log_phi_full[w], color=black+(clr_step*7), thick=linethick*solid
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[17] and Mi_bin_full gt binranges[16], N)
if N ge 5 then begin
   oplot, z_full[w],  log_phi_full[w],    color=black+(clr_step*9), psym=8
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*9), ERRcolor=black+(clr_step*9)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*9), ERRcolor=black+(clr_step*9)

   if choice_line eq 'y' then oplot, z_full[w],  log_phi_full[w],    color=black+(clr_step*9), thick=linethick*solid
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[20] and Mi_bin_full gt binranges[19], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*11), psym=8
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*11), ERRcolor=black+(clr_step*11)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*11), ERRcolor=black+(clr_step*11)

   if choice_line eq 'y' then oplot, z_full[w], log_phi_full[w], color=black+(clr_step*11), thick=linethick*solid
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[23] and Mi_bin_full gt binranges[22], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*13), psym=8
   z_full_err = z_full[w] - z_full[w]
   oploterror, z_full[w], log_Phi_full[w], z_full_err, up[w], $
               /hibar, errthick=errthick, psym=8, color=black+(clr_step*13), ERRcolor=black+(clr_step*13)
   oploterror, z_full[w], log_Phi_full[w], z_full_err, down[w], $
               /lobar, errthick=errthick, psym=8, color=black+(clr_step*13), ERRcolor=black+(clr_step*13)

   if choice_line eq 'y' then oplot, z_full[w], log_phi_full[w], color=black+(clr_step*13), thick=linethick*solid
   print, N, mean(Mi_bin_full[w])
endif





;;
;; LABELS
;;
charsize=charsize*1.1
charthick=charthick*1.2
;xpos =  0.10
xpos =  1.20
ypos = -8.85
yoff =  0.20

;xyouts, xpos, ypos+(yoff*0),  '-28.65',  charsize=charsize, charthick=charthick, color=black
xyouts, xpos, ypos+(yoff*0),  '-28.35',  charsize=charsize, charthick=charthick, color=black+(clr_step*1)
;xyouts, xpos, ypos+(yoff*2),  '-28.05',  charsize=charsize, charthick=charthick, color=black+(clr_step*2)
xyouts, xpos, ypos+(yoff*2),  '-27.75',  charsize=charsize, charthick=charthick, color=black+(clr_step*3)
;xyouts, xpos, ypos+(yoff*4),  '-27.45',  charsize=charsize, charthick=charthick, color=black+(clr_step*4)
xyouts, xpos, ypos+(yoff*4),  '-27.15',  charsize=charsize, charthick=charthick, color=black+(clr_step*5)
;xyouts, xpos, ypos+(yoff*6),  '-26.85',  charsize=charsize, charthick=charthick, color=black+(clr_step*6)
xyouts, xpos, ypos+(yoff*6),  '-26.55',  charsize=charsize, charthick=charthick, color=black+(clr_step*7)
;xyouts, xpos, ypos+(yoff*8),  '-26.25',  charsize=charsize, charthick=charthick, color=black+(clr_step*8)
xyouts, xpos, ypos+(yoff*8),  '-25.95',  charsize=charsize, charthick=charthick, color=black+(clr_step*9)
;xyouts, xpos, ypos+(yoff*10), '-25.65',  charsize=charsize, charthick=charthick, color=black+(clr_step*10)
xyouts, xpos, ypos+(yoff*10), '-25.35',  charsize=charsize, charthick=charthick, color=black+(clr_step*11)
;xyouts, xpos, ypos+(yoff*12), '-25.05',  charsize=charsize, charthick=charthick, color=black+(clr_step*12)
xyouts, xpos, ypos+(yoff*12), '-24.75',  charsize=charsize, charthick=charthick, color=black+(clr_step*13)
;xyouts, xpos, ypos+(yoff*14), '-24.45',  charsize=charsize, charthick=charthick, color=black+(clr_step*14)
xyouts, xpos, ypos+(yoff*14), '-24.15',  charsize=charsize, charthick=charthick, color=black+(clr_step*15)

xyouts, xpos, ypos+(yoff*15.5), 'M!Ii!N[z=2]',  charsize=charsize, charthick=charthick, color=black

xyouts, 0.23, -5.68, 'PLE!Dz=0.3-2.2!N',  charsize=charsize/1.4, charthick=charthick, color=black
xyouts, 0.23, -5.82, 'LEDE!Dz=2.2-3.5!N',  charsize=charsize/1.4, charthick=charthick, color=black

legend, '' , $
        position=[0.34, -5.65], box=0, linestyle=1, color=black, $
        charsize=charsize/4., charthick=16., thick=48


;xyouts, -31.20, -0.75,'(Richards+ 2006)',  color=0, charsize=charsize, charthick =charthick



if plot_log_age eq 1 then begin 
   charthick=charthick/1.2
   charsize=charsize/1.1
   
;;red, omega0=0.30, omegalambda=0.70, h100=0.70
   zzz = [-0.01, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0]   
   nticks=n_elements(zzz)
   ttt = getage(zzz)   
   
;; zzz = -0.0100000     0.500000      1.00000      2.00000      3.00000      5.00000      10.0000
;; gives *where* on the top axis the labels should go. 
;; needs to be converted to age, and flipped around. 
   
   zzz_name = ['0.0', '0.5', '1.0', '2', '3', '5', '7']
   ttt_norm =  ttt / 13.608
   
   xtickv = alog10(ttt)
;xtickv = alog10(VCOMOVING(z_input)/1e27)
   
   print, 'charthick', charthick
   
   axis, xaxis=1,       xtickv=xtickv, $
         xtickname=zzz_name, $
         xticks=nticks-1, $
         xthick=xthick, charthick=charthick, charsize=charsize, $ 
;      XTICKFORMAT='(F5.2)',  $
         xtitle='!8z!3, redshift' 
;      xtitle='!8 Age of the Universe (Gyr)',charsize=2

endif
   
device, /close
set_plot, 'X'     

end
