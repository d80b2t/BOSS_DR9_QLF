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

choice_plot_R06 = 'y'
;read, choice_plot_R06, PROMPT=' - Plot SDSS (R06) points?? y/n  '

choice_plot_boss ='y'
;read, choice_plot_boss, PROMPT=' - Plot BOSS (DR9) points?? y/n  '

choice_plot_boss21 = 'y' 
;read, choice_plot_boss, PROMPT=' - Plot boss21  points?? y/n  '

choice_boss_line  = 'n'
;read, choice_BOSS_line, PROMPT=' - Plot BOSS line?? y/n  '


choice_models    = 'y'

choice_plot_PLE  = 'n'
;read, choice_plot_PLE, PROMPT=' - Plot PLE models?? y/n  '

choice_plot_LDDE = 'y'
;read, choice_plot_LDDE, PROMPT=' - Plot LDDE models?? y/n  '

choice_plot_LEDE = 'n'
;read, choice_plot_LEDE, PROMPT=' - Plot LEDE models?? y/n  '


;;
;;  S D S S   D R 3
;;
readcol, '../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor

R06_delta_up    = alog10((10^log_PhiR06)+(sigma_Phi*1e-9))
R06_delta_up = (log_PhiR06 -R06_delta_up)

R06_delta_down = alog10((10^log_PhiR06)-(sigma_Phi*1e-9))
R06_delta_down = (log_PhiR06 - R06_delta_down)



;;
;;  B O S S   D R 9
;;
readcol, '../pro/qlf/My_QLF_iband_boss_temp.dat', $  ;; USE THIS FOR THE S82 plots!!!
         z_boss, Mi_z2_boss, blah_boss, log_Num_den_boss, raw_N_QSOs_boss, sigma_Phi_BOSS

log_Phi_BOSS =  log_Num_den_boss

boss_delta_up = alog10 ((10^log_Phi_BOSS+sigma_Phi_BOSS))
boss_delta_up = boss_delta_up - log_Phi_BOSS

boss_delta_down = alog10 ((10^log_Phi_BOSS-sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down - log_Phi_BOSS)


;;
;;  B O S S   2 1
;;

readcol, '../pro/qlf/My_QLF_iband_boss21_temp.dat', $  ;; USE THIS FOR THE boss21 chunk!!
         z_boss21, Mi_z2_boss21, blah_boss21, log_Num_den_boss21, raw_N_QSOs_boss21, sigma_Phi_BOSS21

log_Phi_BOSS21 =  log_Num_den_boss21

boss21_delta_up   = alog10 ((10^log_Phi_BOSS21 + sigma_Phi_BOSS21))
boss21_delta_up   = boss21_delta_up - log_Phi_BOSS21

boss21_delta_down = alog10 ((10^log_Phi_BOSS21 - sigma_Phi_BOSS21))
boss21_delta_down = abs(boss21_delta_down - log_Phi_BOSS21)




;;
;;  M o d e l   f i t s...
;;
if (choice_models eq 'y') then begin
   readcol, '../pro/models/Croom09b_PLE_temp.dat', $
            ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   LOG_PHI_PLE = alog10(PHI_PLE)
      
   readcol, '../pro/models/Croom09b_mLDDE_temp.dat', $
            ii_LDDE, z_LDDE, jj_LDDE, mag_LDDE, e_d_LDDE, zc_LDDE, Phi_LDDE
   LOG_PHI_LDDE = alog10(PHI_LDDE)
   
   readcol, '../pro/models/Croom09b_LEDE_temp.dat', $
            ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
   LOG_PHI_LEDE = alog10(PHI_LEDE)

endif



delta_bin = 0.30 
delta = delta_bin/2.
binranges=[-29.85, -29.55,-29.25, -28.95, -28.65, -28.35, -28.05, -27.75, -27.45, -27.15, -26.85, -26.55, -26.25, -25.95, -25.65, -25.35, -25.05, -24.75, -24.45, -24.15, -23.85, -23.55, -23.25, -22.95,  -22.65, -22.35]



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

charsize=2.8
charthick=7.2
thick=4.4
xthick=6.0
ythick=6.0
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; x-ranges
xmin = 0.0
xmax = 7.5

;; y-ranges
ymin = -10.0
ymax =  -5.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='downsizing_SDSS_BOSS_temp.eps', $
        xsize=12.0, ysize=12.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 1.4, /fill
plot, Z_R06, LOG_PHIR06, $
      psym=8, $
      position=[0.24, 0.20, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $
      xtitle=' redshift ', $
      ytitle='log!I10!N !7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
      color=black

plotsym, 0, 1.7, /fill

thick=thick*2.0
;;
;; SDSS DR3 from RIchards 06
;;
cstep=20

if choice_plot_R06 eq 'y' then begin
   
   w = where( Mi_z2 ge binranges[1] and Mi_z2 lt binranges[2], N)
   oplot, Z_R06[w], LOG_PHIR06[w], psym=8, color=black
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black, errcolor=black
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black, errcolor=black

   w = where( Mi_z2 ge binranges[3] and Mi_z2 lt binranges[4], N)
   oplot, Z_R06[w], LOG_PHIR06[w], psym=8, color=black+(cstep*1.)
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+(cstep*1.)32, errcolor=black+(cstep*1.)32
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+(cstep*1.)32, errcolor=black+(cstep*1.)32

   w = where( Mi_z2 ge binranges[5] and Mi_z2 lt binranges[6], N)
   oplot, Z_R06[w], LOG_PHIR06[w],  psym=8, color=black+(cstep*1.)64. 
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+(cstep*1.)64, errcolor=black+(cstep*1.)64
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+(cstep*1.)64, errcolor=black+(cstep*1.)64

   w = where( Mi_z2 ge binranges[7] and Mi_z2 lt binranges[8], N)
   oplot, Z_R06[w], LOG_PHIR06[w], psym=8, color=black+96.
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+96, errcolor=black+96
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+96, errcolor=black+96

   w = where( Mi_z2 ge binranges[9] and Mi_z2 lt binranges[10], N)
   oplot, Z_R06[w], LOG_PHIR06[w],  psym=8, color=black+128.
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+128, errcolor=black+128
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+128, errcolor=black+128

   w = where( Mi_z2 ge binranges[11] and Mi_z2 lt binranges[12], N)
   oplot, Z_R06[w], LOG_PHIR06[w], psym=8, color=black+160.
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+160, errcolor=black+160
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+160, errcolor=black+160

   w = where( Mi_z2 ge binranges[13] and Mi_z2 lt binranges[14], N)
   oplot, Z_R06[w], LOG_PHIR06[w],  psym=8, color=black+192
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+192, errcolor=black+192
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+192, errcolor=black+192

   w = where( Mi_z2 ge binranges[15] and Mi_z2 lt binranges[16], N)
   oplot, Z_R06[w], LOG_PHIR06[w], psym=8, color=black+224
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+224, errcolor=black+224
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+224, errcolor=black+224

   w = where( Mi_z2 ge binranges[17] and Mi_z2 lt binranges[18], N)
   oplot, Z_R06[w], LOG_PHIR06[w],  psym=8, color=black+254
   Z_R06_err = Z_R06[w] - Z_R06[w]
   oploterror, Z_R06[w], log_PHIR06[w], z_R06_err, R06_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+254, errcolor=black+254
   oploterror, Z_R06[w], log_PhiR06[w], z_R06_err, R06_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+254, errcolor=black+254

endif



;;
;;
;; BOSS DR9
;;
;;
if choice_plot_boss eq 'y' then begin
   
   w = where( Mi_z2_boss ge binranges[1] and Mi_z2_boss lt binranges[2] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w], psym=8, color=black
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black, errcolor=black
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black, errcolor=black
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black
   
   w = where( Mi_z2_boss ge binranges[3] and Mi_z2_boss lt binranges[4]and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w],  psym=8, color=black+32. 
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+32., errcolor=black+32.
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+32., errcolor=black+32.
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+32.

   w = where( Mi_z2_boss ge binranges[5] and Mi_z2_boss lt binranges[6] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w], psym=8, color=black+64.
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+64., errcolor=black+64.
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+64., errcolor=black+64.
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+64.

   w = where( Mi_z2_boss ge binranges[7] and Mi_z2_boss lt binranges[8] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w],  psym=8, color=black+96.
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+96., errcolor=black+96
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+96., errcolor=black+96
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+96.

   w = where( Mi_z2_boss ge binranges[9] and Mi_z2_boss lt binranges[10] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w], psym=8, color=black+128.
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+128., errcolor=black+128
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+128., errcolor=black+128
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+128.

   w = where( Mi_z2_boss ge binranges[11] and Mi_z2_boss lt binranges[12] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w],  psym=8, color=black+160
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+160, errcolor=black+160
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+160, errcolor=black+160
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+160.

   w = where( Mi_z2_boss ge binranges[13] and Mi_z2_boss lt binranges[14] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w], psym=8, color=black+192
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+192, errcolor=black+192
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+192, errcolor=black+192
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+192.

   w = where( Mi_z2_boss ge binranges[15] and Mi_z2_boss lt binranges[16] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w],  psym=8, color=black+224   
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+224, errcolor=black+224
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+224, errcolor=black+224
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+224.

   w = where( Mi_z2_boss ge binranges[17] and Mi_z2_boss lt binranges[18] and  Z_boss ge 2.10, N)
   oplot, Z_boss[w], LOG_PHI_boss[w],  psym=8, color=black+254  
   Z_boss_err = Z_boss[w] - Z_boss[w]
   oploterror, Z_boss[w], LOG_PHI_boss[w], z_boss_err, boss_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+254, errcolor=black+254
   oploterror, Z_boss[w], log_Phi_boss[w], z_boss_err, boss_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+254, errcolor=black+254
   if choice_boss_line eq 'y' then oplot, Z_boss[w], LOG_PHI_boss[w], thick=thick, color=black+254.

endif


;;
;;   boss21
;;
if choice_plot_boss21 eq 'y' then begin
   
;   w = where( Mi_z2_boss21 ge binranges[3] and Mi_z2_boss21 lt binranges[5] and  Z_boss21 ge 1.00 , N)
   w = where( Z_boss21 ge 1.00 , N)
   oplot, Z_boss21[w], LOG_PHI_boss21[w], psym=8, color=black+32
   Z_boss21_err = Z_boss21[w] - Z_boss21[w]
   oploterror, Z_boss21[w], LOG_PHI_boss21[w], z_boss21_err, boss21_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+32, errcolor=black+32
   oploterror, Z_boss21[w], log_Phi_boss21[w], z_boss21_err, boss21_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+32, errcolor=black+32

   w = where( Mi_z2_boss21 ge binranges[17] and Z_boss21 ge 1.00 , N)
   oplot, Z_boss21[w], LOG_PHI_boss21[w], psym=8, color=black+254
   Z_boss21_err = Z_boss21[w] - Z_boss21[w]
   oploterror, Z_boss21[w], LOG_PHI_boss21[w], z_boss21_err, boss21_delta_up[w], $
               /hibar, errthick=12, psym=8, color=black+254, errcolor=black+254
   oploterror, Z_boss21[w], log_Phi_boss21[w], z_boss21_err, boss21_delta_down[w], $
               /lobar, errthick=12, psym=8, color=black+254, errcolor=black+254
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
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+32 
   w = where(mag_LDDE gt binranges[5] and mag_LDDE lt binranges[6], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+64
   w = where(mag_LDDE gt binranges[7] and mag_LDDE lt binranges[8], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+96
   w = where(mag_LDDE gt binranges[9] and mag_LDDE lt binranges[10], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+128
   w = where(mag_LDDE gt binranges[11] and mag_LDDE lt binranges[12], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+160
   w = where(mag_LDDE gt binranges[13] and mag_LDDE lt binranges[14], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+196
   w = where(mag_LDDE gt binranges[15] and mag_LDDE lt binranges[16], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+224
   w = where(mag_LDDE gt binranges[17] and mag_LDDE lt binranges[18], N)
   oplot, Z_LDDE[w], LOG_PHI_LDDE[w], thick=thick, linestyle=1, color=black+254
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
charthick=charthick*1.4
xyouts, xmax*0.74, -9.5, '-29.50',  charsize=charsize, charthick=charthick, color=black
xyouts, xmax*0.74, -9.0, '-28.75',  charsize=charsize, charthick=charthick, color=black+32
xyouts, xmax*0.74, -8.5, '-28.25',  charsize=charsize, charthick=charthick, color=black+64
xyouts, xmax*0.74, -8.0, '-27.50',  charsize=charsize, charthick=charthick, color=black+96
xyouts, xmax*0.74, -7.5, '-27.00',  charsize=charsize, charthick=charthick, color=black+128
xyouts, xmax*0.74, -7.0, '-26.50',  charsize=charsize, charthick=charthick, color=black+160
xyouts, xmax*0.74, -6.5, '-25.75',  charsize=charsize, charthick=charthick, color=black+192
xyouts, xmax*0.74, -6.0, '-25.25',  charsize=charsize, charthick=charthick, color=black+224
xyouts, xmax*0.74, -5.5, '-24.50',  charsize=charsize, charthick=charthick, color=black+254


;xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

device, /close
set_plot, 'X'     


end
