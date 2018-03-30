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

;;
;;   B O S S  +  S D S S
;;

readcol, '../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_20120731.dat', $
         z_input,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, $ 
         NQ_full, ba,  log_Phi_full, bar, error_full, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)', /silent
z_full = z_input


zz = getage(z_full)
;z_full = zz
;z_full = alog(zz+1.)
z_full = alog10(zz)
;z_full = alog10(zz+1.)

zmax  = 2.60
;w_zbin = where(z_2SLAQ_full ge 0.40 and Z_2SLAQ_full le 2.60, N)
;w_zbin = where(z_full ge 0.20 and Z_full le 2.20, N)  ;; N= 74
;w_zbin = where(z_full ge 1.0 and Z_full le 2.20, N)    ;; N = 50
w_zbin = where(z_full ge 0.2 and Z_full le zmax, N)    ;; N = 89
;w_zbin = where(z_full ge 2.20 and Z_full le 3.50, N)
;w_zbin = where(z_full ge 2.20 and Z_full le 2.30, N)
;w_zbin = where(z_full ge 1.0 and Z_full le 3.50, N)
;w_zbin = where(z_full ge 3.00 and Z_full le 3.50, N)

print
print
print, 'No of points that we will be using...', N, '   with zmax =', zmax
print
print, ' minmax(z_full[w_zbin]) ', minmax(z_full[w_zbin])
print
print

;; actually kinda have options whether to use the Mi "bin centre" 
;; or the actual mean of the Mi values... 
          Mi_bin = Mi_bin_full[w_zbin] 
             z_z = z_full[w_zbin]
              NQ = NQ_full[w_zbin]
         log_phi = log_phi_full[w_zbin]
       sigma_Phi = error_full[w_zbin]*1e-9

       sigma_Phi_full = error_full*1e-9

      delta_up   = alog10(10^(log_phi_full) +  sigma_Phi_full)
error_delta_up   = delta_up - log_phi_full
      delta_down = alog10(10^(log_phi_full) -  sigma_Phi_full)
error_delta_down = abs(delta_down - log_phi_full)



choice_boss_line  = 'n'
;read, choice_BOSS_line, PROMPT=' - Plot BOSS line?? y/n  '


choice_models    = 'y'

choice_plot_PLE  = 'n'
;read, choice_plot_PLE, PROMPT=' - Plot PLE models?? y/n  '

choice_plot_LDDE = 'n'
;read, choice_plot_LDDE, PROMPT=' - Plot LDDE models?? y/n  '

choice_plot_LEDE = 'n'
;read, choice_plot_LEDE, PROMPT=' - Plot LEDE models?? y/n  '



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


;readcol, '../../pro/models/chi_sq_PLE_model_a_20120810.dat', $
readcol, '../../pro/models/chi_sq_PLE_model_a_temp.dat', $
         z_PLE, Mi_PLE, Phi_PLE
   LOG_PHI_PLE = alog10(PHI_PLE)

zz_ple = getage(z_ple)
z_ple = alog10(zz_ple)




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

charsize = 4.2*expand
charthick= 10.2*expand
thick    = 10.2 ;; doesn't do anything!!
xthick   = 20.0
ythick   = 20.0
XTICKLEN = 0.03
YTICKLEN = 0.03


;; x-ranges
;; for z, redshift...
xmin = 0.0
xmax = 7.5

;; For Age of the Universe... 
;xmin = 14.0
;xmax =  0.0

;;  in alog10 (age  +1)
xmin =  1.25
xmax =  0.0
;xmax =  -0.2


;; y-ranges
ymin = -9.0   ;-10.0
ymax =  -5.5


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  
;;
set_plot, 'ps'
device, filename='downsizing_SDSS_BOSS_boss21_temp.eps', $
        xsize=(14.0*expand), ysize=16.0*(expand), $
;        xsize=21.0, ysize=24.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated
print, 'charthick', charthick

plotsym, 0, 1.4, /fill
plot, z_full, log_phi_full, $
      psym=8, $
      position=[0.20, 0.12, 0.96, 0.86], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
;      xstyle=1, ystyle=1, $
      xstyle=8+1, ystyle=1, $
 ;     xstyle=8+1, ystyle=3, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $;
;      xtitle=' redshift ', $
      xtitle='log!I10!N(age / Gyr)', $
      ytitle='log!I10!N !7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
      color=black

plotsym, 0, 1.7, /fill

thick=thick*2.0
;;
;; SDSS DR3 from RIchards 06
;;
clr_step  = 18
linethick = 30
solid     = 1.4

ple_linestyle = 

w = where(Mi_bin_full le binranges[0], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black, thick=thick
   print, N, mean(Mi_bin_full[w])
endif

w = where(Mi_bin_full le binranges[1] and Mi_bin_full gt binranges[0], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*0), thick=thick
   print, N, mean(Mi_bin_full[w]), 0.5*(binranges[1] +binranges[0] )
endif

w = where(Mi_bin_full le binranges[2] and Mi_bin_full gt binranges[1], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*0), thick=thick
   print, N, mean(Mi_bin_full[w]), 0.5*(binranges[2] +binranges[1] )
endif

;;
;; first  N ge 5  bin
;;
w     = where(Mi_bin_full le binranges[3] and Mi_bin_full gt binranges[2], N)
w_ple = where(Mi_PLE      le binranges[3] and Mi_PLE      gt binranges[2], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*0), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*0), thick=linethick, linestyle=ple_linestyle

   print, N, mean(Mi_bin_full[w]), '0.5*(binranges[3] +binranges[2] )  ', 0.5*(binranges[3] +binranges[2] )
endif

w     = where(Mi_bin_full le binranges[4] and Mi_bin_full gt binranges[3], N)
w_ple = where(Mi_PLE      le binranges[4] and Mi_PLE      gt binranges[3], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*1), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*1), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w]), 0.5*(binranges[4] +binranges[3] )
endif

w =     where(Mi_bin_full le binranges[5] and Mi_bin_full gt binranges[4], N)
w_ple = where(Mi_PLE      le binranges[5] and Mi_PLE      gt binranges[4], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*2), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*2), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[6] and Mi_bin_full gt binranges[5], N)
w_ple = where(Mi_PLE      le binranges[6] and Mi_PLE      gt binranges[5], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*3), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*3), thick=linethick, linestyle=2

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[7] and Mi_bin_full gt binranges[6], N)
w_ple = where(Mi_PLE      le binranges[7] and Mi_PLE      gt binranges[6], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w], color=black+(clr_step*4), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*4), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[8] and Mi_bin_full gt binranges[7], N)
w_ple = where(Mi_PLE      le binranges[8] and Mi_PLE      gt binranges[7], N_ple)
if N ge 5 then begin 
   oplot, z_full[w],    log_phi_full[w], color=black+(clr_step*5), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*5), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[9] and Mi_bin_full gt binranges[8], N)
w_ple = where(Mi_PLE      le binranges[9] and Mi_PLE      gt binranges[8], N_ple)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*6), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*6), thick=linethick, linestyle=1

endif
print, N, mean(Mi_bin_full[w])

w     = where(Mi_bin_full le binranges[10] and Mi_bin_full gt binranges[9], N)
w_ple = where(Mi_PLE      le binranges[10] and Mi_PLE      gt binranges[9], N_ple)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*7), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*7), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[11] and Mi_bin_full gt binranges[10], N)
w_ple = where(Mi_PLE      le binranges[11] and Mi_PLE      gt binranges[10], N_ple)
if N ge 5 then begin 
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*8), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*8), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[12] and Mi_bin_full gt binranges[11], N)
w_ple = where(Mi_PLE      le binranges[12] and Mi_PLE      gt binranges[11], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*9), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*9), thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[13] and Mi_bin_full gt binranges[12], N)
w_ple = where(Mi_PLE      le binranges[13] and Mi_PLE      gt binranges[12], N_ple)
if N ge 5 then begin
   oplot, z_full[w],    log_phi_full[w],    color=black+(clr_step*10), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*10), thick=linethick, linestyle=1
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[14] and Mi_bin_full gt binranges[13], N)
w_ple = where(Mi_PLE      le binranges[14] and Mi_PLE      gt binranges[13], N_ple)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*11), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*11), thick=linethick, linestyle=1
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[15] and Mi_bin_full gt binranges[14], N)
w_ple = where(Mi_PLE      le binranges[15] and Mi_PLE      gt binranges[14], N_ple)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w],       color=black+(clr_step*12), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*12), thick=linethick, linestyle=1
   print, N, mean(Mi_bin_full[w])
endif 

w     = where(Mi_bin_full le binranges[16] and Mi_bin_full gt binranges[15], N)
w_ple = where(Mi_PLE      le binranges[16] and Mi_PLE      gt binranges[15], N_ple)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*13), thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*13)+2, thick=linethick, linestyle=1
   print, N, mean(Mi_bin_full[w])
endif

w     = where(Mi_bin_full le binranges[17] and Mi_bin_full gt binranges[16], N)
w_ple = where(Mi_PLE      le binranges[17] and Mi_PLE      gt binranges[16], N_ple)
;; Here, N = 5 and clr_step*14 = 252 is 
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w],      color=black+(clr_step*14)+4, thick=linethick*solid
   oplot, z_PLE[w_ple], log_phi_PLE[w_ple], color=black+(clr_step*14)+4, thick=linethick, linestyle=1

   print, N, mean(Mi_bin_full[w])
endif

w = where(Mi_bin_full le binranges[18] and Mi_bin_full gt binranges[17], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*15), thick=thick
   print, N, mean(Mi_bin_full[w])
endif

w = where(Mi_bin_full le binranges[19] and Mi_bin_full gt binranges[18], N)
if N ge 5 then begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*16), thick=thick
   print, N, mean(Mi_bin_full[w])
endif

w = where(Mi_bin_full le binranges[20] and Mi_bin_full gt binranges[19], N)
if N ge 5 then begin 
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*17), thick=thick
   print, N, mean(Mi_bin_full[w])
endif

w = where(Mi_bin_full gt binranges[20], N)
if N ge 5 then  begin
   oplot, z_full[w], log_phi_full[w], color=black+(clr_step*18), thick=thick
   print, N, mean(Mi_bin_full[w])
endif




;;
;;  M O D E L S
;;

;;
;; P L E 
;;
;thick=thick*0.5
if choice_plot_PLE eq 'y' then begin
   w = where(mag_PLE gt binranges[1] and mag_PLE lt binranges[2], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black
   w = where(mag_PLE gt binranges[3] and mag_PLE lt binranges[4], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+32
   w = where(mag_PLE gt binranges[5] and mag_PLE lt binranges[6], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+64
   w = where(mag_PLE gt binranges[7] and mag_PLE lt binranges[8], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+96
   w = where(mag_PLE gt binranges[9] and mag_PLE lt binranges[10], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+128
   w = where(mag_PLE gt binranges[11] and mag_PLE lt binranges[12], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+160
   w = where(mag_PLE gt binranges[13] and mag_PLE lt binranges[14], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+196
   w = where(mag_PLE gt binranges[15] and mag_PLE lt binranges[16], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+224
   w = where(mag_PLE gt binranges[17] and mag_PLE lt binranges[18], N)
   oplot, Z_PLE[w], LOG_PHI_PLE[w], thick=thick, color=black+254
endif

thick=thick*2.0
if choice_plot_LDDE eq 'y' then begin
   w = where(mag_LDDE gt binranges[1] and mag_LDDE lt binranges[2], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black
   w = where(mag_LDDE gt binranges[3] and mag_LDDE lt binranges[4], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*1) 
   w = where(mag_LDDE gt binranges[5] and mag_LDDE lt binranges[6], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*2)
   w = where(mag_LDDE gt binranges[7] and mag_LDDE lt binranges[8], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*3)
   w = where(mag_LDDE gt binranges[9] and mag_LDDE lt binranges[10], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*4)
   w = where(mag_LDDE gt binranges[11] and mag_LDDE lt binranges[12], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*5)
   w = where(mag_LDDE gt binranges[13] and mag_LDDE lt binranges[14], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*6)
   w = where(mag_LDDE gt binranges[15] and mag_LDDE lt binranges[16], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*7)
   w = where(mag_LDDE gt binranges[17] and mag_LDDE lt binranges[18], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*8)
;; for boss21
   w = where(mag_LDDE gt binranges[19] and mag_LDDE lt binranges[20], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*9)
   w = where(mag_LDDE gt binranges[21] and mag_LDDE lt binranges[22], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*10)
   w = where(mag_LDDE gt binranges[23] and mag_LDDE lt binranges[24], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*11)
   w = where(mag_LDDE gt binranges[25] and mag_LDDE lt binranges[26], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+(clr_step*12)
endif

thick=thick*0.8
if choice_plot_LEDE eq 'y' then begin
   w = where(mag_LEDE gt binranges[1] and mag_LEDE lt binranges[2], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black
   w = where(mag_LEDE gt binranges[3] and mag_LEDE lt binranges[4], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+32 
   w = where(mag_LEDE gt binranges[5] and mag_LEDE lt binranges[6], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+64
   w = where(mag_LEDE gt binranges[7] and mag_LEDE lt binranges[8], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+96
   w = where(mag_LEDE gt binranges[9] and mag_LEDE lt binranges[10], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+128
   w = where(mag_LEDE gt binranges[11] and mag_LEDE lt binranges[12], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+160
   w = where(mag_LEDE gt binranges[13] and mag_LEDE lt binranges[14], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+196
   w = where(mag_LEDE gt binranges[15] and mag_LEDE lt binranges[16], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+224
   w = where(mag_LEDE gt binranges[17] and mag_LEDE lt binranges[18], N)
   oplot, Z_LEDE[w], LOG_PHI_LEDE[w], thick=thick, linestyle=2, color=black+254
endif

;;
;; LABELS
;;
charsize=charsize*0.9
charthick=charthick;*1.3
;xpos =  0.10
xpos =  1.20
ypos = -8.8
yoff =  0.20

xyouts, xpos, ypos+(yoff*0),  '-28.65',  charsize=charsize, charthick=charthick, color=black
xyouts, xpos, ypos+(yoff*1),  '-28.35',  charsize=charsize, charthick=charthick, color=black+(clr_step*1)
xyouts, xpos, ypos+(yoff*2),  '-28.05',  charsize=charsize, charthick=charthick, color=black+(clr_step*2)
xyouts, xpos, ypos+(yoff*3),  '-27.75',  charsize=charsize, charthick=charthick, color=black+(clr_step*3)
xyouts, xpos, ypos+(yoff*4),  '-27.45',  charsize=charsize, charthick=charthick, color=black+(clr_step*4)
xyouts, xpos, ypos+(yoff*5),  '-27.15',  charsize=charsize, charthick=charthick, color=black+(clr_step*5)
xyouts, xpos, ypos+(yoff*6),  '-26.85',  charsize=charsize, charthick=charthick, color=black+(clr_step*6)
xyouts, xpos, ypos+(yoff*7),  '-26.55',  charsize=charsize, charthick=charthick, color=black+(clr_step*7)
xyouts, xpos, ypos+(yoff*8),  '-26.25',  charsize=charsize, charthick=charthick, color=black+(clr_step*8)
xyouts, xpos, ypos+(yoff*9),  '-25.95',  charsize=charsize, charthick=charthick, color=black+(clr_step*9)
xyouts, xpos, ypos+(yoff*10), '-25.65',  charsize=charsize, charthick=charthick, color=black+(clr_step*10)
xyouts, xpos, ypos+(yoff*11), '-25.35',  charsize=charsize, charthick=charthick, color=black+(clr_step*11)
xyouts, xpos, ypos+(yoff*12), '-25.05',  charsize=charsize, charthick=charthick, color=black+(clr_step*12)
xyouts, xpos, ypos+(yoff*13), '-24.75',  charsize=charsize, charthick=charthick, color=black+(clr_step*13)
xyouts, xpos, ypos+(yoff*14), '-24.45',  charsize=charsize, charthick=charthick, color=black+(clr_step*14)
xyouts, xpos, ypos+(yoff*15), 'M!Ii!N[z=2]',  charsize=charsize, charthick=charthick, color=black


charthick=charthick;/1.3
charsize=charsize/0.9

;xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color


;;red, omega0=0.30, omegalambda=0.70, h100=0.70
zzz = [-0.01, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0]   
nticks=n_elements(zzz)
ttt = getage(zzz)   

;; zzz = -0.0100000     0.500000      1.00000      2.00000      3.00000      5.00000      10.0000
;; gives *where* on the top axis the labels should go. 
;; needs to be converted to age, and flipped around. 

zzz = ['0.0', '0.5', '1.0', '2', '3', '5', '7']
ttt_norm =  ttt / 13.608

xtickv  = alog10(ttt)
print, 'charthick', charthick

axis, xaxis=1,       xtickv=xtickv, $
      xtickname=zzz, $
      xticks=nticks-1, $
      xthick=xthick, charthick=charthick, charsize=charsize, $ 
;      XTICKFORMAT='(F5.2)',  $
      xtitle='!8z!3, redshift' 
;      xtitle='!8 Age of the Universe (Gyr)',charsize=2


device, /close
set_plot, 'X'     


end
