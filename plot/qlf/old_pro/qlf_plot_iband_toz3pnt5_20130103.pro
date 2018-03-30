;+
; NAME:
;       qlf_plot
; 
; PURPOSE:
;        To plot Quasar Luminosity Functions
;
; CALLING SEQUENCE:
;       .run qlf_plot
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

red, omega0=0.30, omegalambda=0.70, h100=0.700

;; 
;;
;;  C h o i c e s...
;;
;;
print
print

choice_mydr3_points = 'n'
;read, choice_mydr3_points, PROMPT=' - Plot My DR3 points?? y/n  '

choice_BOSS_points   = 'y'
choice_S82_points    = 'n' 
choice_S82_errorbars = 'n'
choice_DR9_points    = 'y'
choice_DR9_errorbars = 'y'
choice_DR9_nowgt     = 'y'

;read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '
if choice_BOSS_points eq 'y' then begin
   choice_BOSS_line = 'y'
   ;read, choice_BOSS_line, PROMPT=' - Plot BOSS (data) line ?? y/n  '
;   read, choice_S82_points, PROMPT=' - Plot S82 points?? y/n  '
   if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
   read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
   if choice_DR9_points eq 'y' then begin
      ;read, choice_DR9_nowgt,     PROMPT=' - Plot DR9 points with NO WGT?  y/n  '
      ;read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
   endif

   choice_BOSS21_points = 'n'
;   read, choice_BOSS21_points, PROMPT=' - Plot BOSS21 points?? y/n  '

   choice_MMT_points = 'y'
  read, choice_MMT_points,    PROMPT=' - Plot boss21+MMT points?? y/n  '

   choice_2SLAQ_points = 'n'
;   read, choice_2SLAQ_points,    PROMPT=' - Plot 2SLAQ points?? y/n  '
endif   

print
print
choice_models = 'n'
;read, choice_models, PROMPT=' - Plot (PLE, LDDE and LEDE) model lines?? y/n  '

print
choice_models_HRH07  = 'n'
;read, choice_models_HRH07, PROMPT=' - Plot HRH07 models  ?? y/n  '
print
print

;; 
;;         D  A  T  A
;;
;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor
readcol, '../../pro/qlf/My_QLF_iband_dr3_final.dat', z_bin_dr3, Abs_mag_bin_dr3, blah_dr3, log_Num_den_dr3
print
print, '../../data/Richards06_Table06.dat READ-IN', N_elements(z_R06)
print

R06_delta_up    = alog10((10^log_PhiR06)+(sigma_Phi*1e-9))
R06_delta_up = (log_PhiR06 -R06_delta_up)
R06_delta_down = alog10((10^log_PhiR06)-(sigma_Phi*1e-9))
R06_delta_down = (log_PhiR06 - R06_delta_down)

R06_err = alog10(sigma_Phi*1e-9)

;;
;; There's also this:
;;
;; The sdssdr3fixyes.all from GTR
;readcol, '../../data/sdssdr3fixyes.all', $
 ;        z_dr3, Mi_dr3, num_dr3, numerr_dr3, N_Qcor_dr3, fill, N_Q, z_dr3_bin 
  ;       ;;0.490 -26.250 1.817604e-08 5.480553e-09 11.633511 1 11 0.59


;;
;;
;;   B O S S   D R 9
;;
;;
;; 
;;  Full DR9  (no correction)
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_temp.dat', $     
readcol, '../../pro/qlf/My_QLF_iband_boss_20120619.dat', $ ;; used in DR9 QLF paper!!
         z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss, raw_N_QSOs_boss, sigma_Phi_BOSS
print
print, '../pro/qlf/My_QLF_iband_boss_temp.dat  READ-IN', N_elements(z_bin_boss)
print

log_Phi_BOSS =  log_Num_den_boss
boss_delta_up = alog10 ((10^log_Phi_BOSS+sigma_Phi_BOSS))
boss_delta_up = boss_delta_up - log_Phi_BOSS
boss_delta_down = alog10 ((10^log_Phi_BOSS-sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down - log_Phi_BOSS)

;; 
;;  WITH correction/weighting (wgt)
;;
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $ 
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat', $ ;; used in the DR9 QLF paper!! 
        z_bin_boss_wgt, Abs_mag_bin_boss_wgt, blah_boss_wgt, log_Num_den_boss_wgt, raw_N_QSOs_boss_wgt, sigma_Phi_BOSS_wgt
print
print, '../pro/qlf/My_QLF_iband_boss_wgt_temp.dat  READ-IN', N_elements(z_bin_boss)
print

log_Phi_BOSS_wgt =  log_Num_den_boss_wgt

;; HORRIBLE LINE, just here FOR THE TIME BEING!!!
;sigma_Phi_BOSS_wgt = sigma_Phi_BOSS


boss_delta_up_wgt   = alog10((10^log_Phi_BOSS_wgt+sigma_Phi_BOSS_wgt))
boss_delta_up_wgt   = boss_delta_up_wgt - log_Phi_BOSS_wgt
boss_delta_down_wgt = alog10((10^log_Phi_BOSS_wgt-sigma_Phi_BOSS_wgt))
boss_delta_down_wgt = abs(boss_delta_down_wgt - log_Phi_BOSS_wgt)


;;
;; B O S S   S T R I P E   8 2
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_S82_stnd.dat', $ ;; USE THIS FOR THE S82 plots!!!
;readcol, '../../pro/qlf/My_QLF/My_QLF_iband_boss_S82_March2012.dat', $  ;; 
readcol, '../../pro/qlf/My_QLF_iband_boss_S82_McGkcor.dat', $  ;; 
        z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82
print
print, '../../pro/qlf/My_QLF_iband_boss_S82_March2012.dat  READ-IN', N_elements(z_bin_boss_S82)
print

log_Phi_BOSS_S82 =  log_Num_den_boss_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)


;;
;;  B O S S   2 1
;;
;readcol, '../../pro/qlf/My_QLF_iband_boss21_temp.dat', $  
;readcol, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_temp.dat', $  
;readcol, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_20120828.dat', $  ;; This used for the DR9 QLF paper!
readcol, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_McGkcor.dat', $  
         z_bin_boss21, Abs_mag_bin_boss21, blah_boss21, log_Num_den_boss21, raw_N_QSOs_boss21, sigma_Phi_BOSS21
print
print, '../pro/qlf/My_QLF_iband_boss21_wSelfn_20120305.dat  READ-IN', N_elements(z_bin_boss21)
print 

log_Phi_BOSS21 =  log_Num_den_boss21

boss21_delta_up   = alog10 ((10^log_Phi_BOSS21 + sigma_Phi_BOSS21))
boss21_delta_up   = boss21_delta_up - log_Phi_BOSS21
boss21_delta_down = alog10 ((10^log_Phi_BOSS21 - sigma_Phi_BOSS21))
boss21_delta_down = abs(boss21_delta_down - log_Phi_BOSS21)


;;
;;  B O S S   2 1   +   M M T
;;
readcol, '../../pro/qlf/Palanque-Delabrouille_2012.dat', $
         z_MMT, Abs_mag_MMT, nq_MMT, log_Phi_MMT, sigma_phi_MMT
;; M_i(z=2) = M_g(z=2) - 0.25
Abs_mag_MMT = Abs_mag_MMT - 0.25

MMT_delta_up   = sigma_phi_MMT
MMT_delta_down = sigma_phi_MMT

readcol, '../../data/Croom09b_2SLAQ_QLF.dat', $
         Mg_2SLAQ, z_bar, NQ_2SLAQ, log_phi_2SLAQ, delta_log_phi_lower, delta_log_phi_upper
Mi_2SLAQ = Mg_2SLAQ - 0.25

SLAQ_delta_up   = delta_log_phi_upper
SLAQ_delta_down = delta_log_phi_lower


;; Setting a sanity-check completeness limit. 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.

boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.80

red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]
dlums = DLUMINOSITY(red_bins) / 1e6
Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ; + kcor(fix(redphot/0.01))

iMag_limit_line = fltarr(11, 61)
limit_lines = (findgen(61)/10.)-10.0                            
for ii = 0L, 11-1 do iMag_limit_line[ii,*] = Abs_iMag_limit[ii]



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


charsize  = 4.8
charthick = 7.0
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.04
YTICKLEN  = 0.06

;; x-ranges
;x_min = -18.001
x_min = -22.001
x_max = -30.50

;; y-ranges
y_min = -9.20
y_max = -5.10

;; xy_outs 
x_xyouts = -26.00
y_xyouts =  -5.70

plot_sym_size_R06  = 2.0
plot_sym_size_BOSS = 2.8
plot_sym_size_MMT  = 1.4


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;  M o d e l   f i t s...
;;
;;
if (choice_models eq 'y') then begin
   readcol, '../../pro/models/Croom09b_PLE_temp.dat', $
            ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   
   readcol, '../../pro/models/Croom09b_mLDDE_temp.dat', $
            ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE
   
   readcol, '../../pro/models/Croom09b_LEDE_temp.dat', $
            ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;  Actually plotting stuff...
;;
;;
;;  Q L F  i-band  
;;   Richards+06 as a template...
;;
set_plot, 'ps'
!p.multi=[0,3,3]
device, filename='QLF_iband_R06_z2to3pnt5_temp.ps', $
        xsize=18.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;    z = 0.49
;;   0.30 < z < 0.68   for Richards06
;;   0.40 < z < 0.68   for Croom09
w  =  where(z_R06 gt 0.11 and z_R06 lt 0.68) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.70, 0.42, 0.98], $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xstyle=1, $
       ystyle=1, $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)'

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[0,*], limit_lines, color=160, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!80.30<z<0.68!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w], $
            /hibar, errthick=12, psym=8;, color=64
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], $
            /lobar, errthick=12, psym=8;, color=164

if (choice_mydr3_points eq 'y') then begin
   ww =  where(z_bin_dr3 gt 0.11 and z_bin_dr3 lt 0.68) 
   loadct, 6
   plotsym, 0, 1.5
   oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
endif


if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 0.30 and z_bin_boss21 lt 0.68 and RAW_N_QSOS_boss21 gt 0, N )
endif



if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 0.45 and z_ple lt   0.55, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 0.45 and z_mLDDE lt 0.55, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 0.45 and z_LEDE  lt 0.55, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

HRH07_color = 128 ;; green
if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z0.49.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif


loadct, 6
boss_color = 160 ;; the "dull turquiose, that I used a lot in plots"
boss_color =  60 ;; more or less red..
boss21_color = 60 ;; the "dull turquiose, that I used a lot in plots"
loadct, clr_table
MMT_color = blue 

if (choice_BOSS_points eq 'y') then begin
   loadct, 6
   if choice_DR9_points eq 'y' then begin
      www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 $
                   and Abs_mag_bin_boss lt -24.5 and RAW_N_QSOS_boss gt 0, N) 
      ;;NOTE TO NPR need to keep putting in this  RAW_N_QSOS_BOSS cut...
;     oplot, Abs_Mag_bin_boss[www],     log_Num_den_boss[www],     color=boss_color, thick=6
      oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   endif
   if (choice_S82_points eq 'y') then begin
      www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                   Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
      oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], color=60, thick=12
      xyouts, x_xyouts-1.0, y_xyouts-3.1, '!8z=2.4', $
              charsize=2.2, charthick=8., color=60.
   endif
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  z  =  0.87
;;
w  =  where(z_R06 gt 0.68 and z_R06 lt 1.06) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.42, 0.70, 0.70, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!80.68<z<1.06!3', charsize=2.2, charthick=6

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
   plotsym, 8, plot_sym_size_MMT, /fill

 w_MMT = where(z_MMT gt 0.68 and z_MMT lt 1.06 and NQ_MMT gt 0 and Abs_Mag_MMT lt -22.8, N)

   oplot, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], color=MMT_color, ps=8
   Mi_MMT_err = Abs_Mag_MMT[w_MMT] - Abs_Mag_MMT[w_MMT]
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_up[w_MMT], $
               /hibar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_down[w_MMT], $
               /lobar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
endif

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

if (choice_mydr3_points eq 'y') then begin
   ww =  where(z_bin_dr3 gt 0.68 and z_bin_dr3 lt 1.06) 
   loadct, 6
   plotsym, 0, 1.5
   ;oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
endif

if (choice_BOSS_points eq 'y') then begin
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[1,*], limit_lines, color=boss_color, linestyle=1, thick=12
   if (choice_DR9_points eq 'y') then begin
      www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                   Abs_mag_bin_boss lt -24.5 and RAW_N_QSOS_boss gt 0, N )
      loadct, 6
;     oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=boss_color, thick=6
      oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   endif

   if (choice_S82_points eq 'y') then begin
      www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                   Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
      oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=60, thick=12
   endif
endif

if (choice_BOSS21_points eq 'y') then begin 
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 0.68 and z_bin_boss21 lt 1.06 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -23.0, N )
   ;; Just really can't think I can justify this bin due to
   ;; e.g. host galaxy light...
;   oplot, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], color=boss21_color, ps=8
endif

if (choice_2SLAQ_points eq 'y') then begin
   loadct, 6
   plotsym, 8, plot_sym_size_R06/1.2, /fill
   w_2SLAQ = where(z_bar gt 0.68 and z_bar lt 1.06, N_2SLAQ)
   oplot, Mi_2SLAQ[w_2SLAQ], LOG_PHI_2SLAQ[w_2SLAQ], ps=8, color=132
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 0.80 and z_ple lt   0.90, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 0.80 and z_mLDDE lt 0.90, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 0.80 and z_LEDE  lt 0.90, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z0.87.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  z = 1.25
;;
w  =  where(z_R06 gt 1.06 and z_R06 lt 1.44) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.70, 0.70, 0.98, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!81.06<z<1.44!3', charsize=2.2, charthick=6

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
;   plotsym, 8, plot_sym_size_R06/1.2, /fill
   plotsym, 8, plot_sym_size_MMT, /fill

; w_MMT = where(z_MMT gt 1.06 and z_MMT lt 1.44 and NQ_MMT gt 0, N)
 w_MMT = where(z_MMT gt 1.06 and z_MMT lt 1.44 and NQ_MMT gt 0 and Abs_Mag_MMT lt -22.8, N)

   oplot, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], color=MMT_color, ps=8
   Mi_MMT_err = Abs_Mag_MMT[w_MMT] - Abs_Mag_MMT[w_MMT]
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_up[w_MMT], $
               /hibar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_down[w_MMT], $
               /lobar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
endif
loadct, 6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

if (choice_mydr3_points eq 'y') then begin
   ww =  where(z_bin_dr3 gt 1.06 and z_bin_dr3 lt 1.44) 
   loadct, 6
   plotsym, 0, 1.5
   oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
endif


loadct, 6
if (choice_BOSS_points eq 'y') then begin
   ;; Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[2,*], limit_lines, color=boss_color, linestyle=1, thick=12
   if (choice_DR9_points eq 'y') then begin
      www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                   Abs_mag_bin_boss lt -24.5 and RAW_N_QSOS_BOSS gt 0, N)  
;     oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=boss_color, thick=6
      oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   endif
   if (choice_S82_points eq 'y') then begin
      www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                   Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
      oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=60, thick=12
   endif
endif

if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 1.06 and z_bin_boss21 lt 1.44 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -22.8, N)
   oplot, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], color=boss21_color, ps=8
   Mi_boss_err = log_Num_den_boss21[w21] - log_Num_den_boss21[w21]
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_up[w21], /hibar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_down[w21], /lobar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
endif

loadct, 6
if (choice_2SLAQ_points eq 'y') then begin
   plotsym, 8, plot_sym_size_R06/1.2, /fill
   w_2SLAQ = where(z_bar gt 1.06 and z_bar lt 1.44, N_2SLAQ)
   oplot, Mi_2SLAQ[w_2SLAQ], LOG_PHI_2SLAQ[w_2SLAQ], ps=8, color=132
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 1.20 and z_ple lt   1.30, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 1.20 and z_mLDDE lt 1.30, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 1.20 and z_LEDE  lt 1.30, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z1.25.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   plot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           2 n d     R O W          
;; 
;;   1 . 4 4   <   z   <   1. 8 2
;;
w  =  where(z_R06 gt 1.44 and z_R06 lt 1.82) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.42, 0.42, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)'; , $
;       ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!81.44<z<1.82!3', charsize=2.2, charthick=6

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
;   plotsym, 8, plot_sym_size_R06/1.2, /fill
   plotsym, 8, plot_sym_size_MMT, /fill

   w_MMT = where(z_MMT gt 1.44 and z_MMT lt 1.68 and NQ_MMT gt 0 and $
                 Abs_Mag_MMT lt -22.8, N)

   oplot, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], color=MMT_color, ps=8
   Mi_MMT_err = Abs_Mag_MMT[w_MMT] - Abs_Mag_MMT[w_MMT]
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_up[w_MMT], $
               /hibar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_down[w_MMT], $
               /lobar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
endif
loadct, 6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

if (choice_mydr3_points eq 'y') then begin
   ww =  where(z_bin_dr3 gt 1.44 and z_bin_dr3 lt 1.82) 
   loadct, 6
   plotsym, 0, 1.5
   oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
endif


loadct, 6
plotsym, 0, 1.5, /fill
if (choice_BOSS_points eq 'y') then begin
   ;; Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[3,*], limit_lines, color=boss_color, linestyle=1, thick=12
   if (choice_DR9_points eq 'y') then begin
      ;;www =  where(z_bin_boss gt 1.44 and z_bin_boss lt 1.82) 
      ;;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=boss_color
      www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                   Abs_mag_bin_boss lt -24.5 and RAW_N_QSOS_BOSS gt 0, N)  
;      oplot, Abs_Mag_bin_boss[www],     log_Num_den_boss[www],     color=boss_color, thick=6
      oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   endif 
   if (choice_S82_points eq 'y') then begin
      www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                   Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
      oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], color=s82_color, thick=12
   endif
endif

if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 1.44 and z_bin_boss21 lt 1.86 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -23.8, N )
   oplot, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], color=boss21_color, ps=8
   Mi_boss_err = log_Num_den_boss21[w21] - log_Num_den_boss21[w21]
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_up[w21], /hibar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_down[w21], /lobar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
endif

loadct, 6
if (choice_2SLAQ_points eq 'y') then begin
   plotsym, 8, plot_sym_size_R06/1.2, /fill
   w_2SLAQ = where(z_bar gt 1.44 and z_bar lt 1.86, N_2SLAQ)
   oplot, Mi_2SLAQ[w_2SLAQ], LOG_PHI_2SLAQ[w_2SLAQ], ps=8, color=132
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 1.60 and z_ple lt  1.70, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 1.60 and z_mLDDE lt 1.70, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 1.60 and z_LEDE  lt 1.70, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z1.63.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif



;;;;;;;;;      z = 2.01         ;;;;;;;;;;;;;;;;;;;;;;;;;
w  =  where(z_R06 gt 1.82 and z_R06 lt 2.20) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.42, 0.42, 0.70, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!81.82<z<2.20!3', charsize=2.2, charthick=6

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
   ;plotsym, 8, plot_sym_size_R06/1.2, /fill
   plotsym, 8, plot_sym_size_MMT, /fill

   w_MMT = where(z_MMT gt 1.82 and z_MMT lt 2.20 and NQ_MMT gt 0 and $
                 Abs_Mag_MMT lt -22.8, N)

   oplot, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], color=MMT_color, ps=8
   Mi_MMT_err = Abs_Mag_MMT[w_MMT] - Abs_Mag_MMT[w_MMT]
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_up[w_MMT], $
               /hibar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_down[w_MMT], $
               /lobar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
endif
loadct, 6

if (choice_2SLAQ_points eq 'y') then begin
   plotsym, 8, plot_sym_size_R06/1.2, /fill
   w_2SLAQ = where(z_bar gt 1.82 and z_bar lt 2.20, N_2SLAQ)
   oplot, Mi_2SLAQ[w_2SLAQ], LOG_PHI_2SLAQ[w_2SLAQ], ps=8, color=132
endif

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

if (choice_mydr3_points eq 'y') then begin
   ww =  where(z_bin_dr3 gt 1.82 and z_bin_dr3 lt 2.20) 
   loadct, 6
   plotsym, 0, 1.5
   oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
endif


loadct, 6
plotsym, 0, 1.5, /fill
if (choice_BOSS_points eq 'y') then begin
   ;; Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[4,*], limit_lines, color=boss_color, linestyle=1, thick=12
   if (choice_DR9_points eq 'y') then begin
      ;; www =  where(z_bin_boss gt 1.82 and z_bin_boss lt 2.20) 
      ;; oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=boss_color
      www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                   Abs_mag_bin_boss lt -24.5 and RAW_N_QSOS_BOSS gt 0, N) 
;     oplot, Abs_Mag_bin_boss[www],     log_Num_den_boss[www],     color=boss_color, thick=6
      oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   endif
   if (choice_S82_points eq 'y') then begin
      www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                   Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
      oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=s82_color, thick=12
   endif
endif

if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   ;iMag_limit_line[4,*] =  -24.1705
   w21 =  where(z_bin_boss21 gt 1.82 and z_bin_boss21 lt 2.20 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -24.17, N) ; $
   oplot, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], color=boss21_color, ps=8
   Mi_boss_err = log_Num_den_boss21[w21] - log_Num_den_boss21[w21]
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_up[w21], /hibar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
   oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
               boss21_delta_down[w21], /lobar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE gt 1.95 and z_ple lt 2.05, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   
   w_mLDDE = where(z_mLDDE gt 1.95 and z_mLDDE lt 2.05, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   
   w_LEDE = where(z_LEDE gt 1.95 and z_LEDE lt 2.05, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif
if (choice_models eq 'y') then begin
   ;; Model lines...
   legend, pos=[-22.2,-7.5], ' ', box=0, thick=14, linestyle = 0, charsize=1.2
   xyouts, -25.0, -7.8, 'PLE', charsize=2.2, charthick=8.
   legend, pos=[-22.2,-8.0], ' ', box=0, thick=14, linestyle = 1, charsize=1.2
   xyouts, -25.0, -8.3, 'LDDE', charsize=2.2, charthick=8.
   legend, pos=[-22.2,-8.5], ' ', box=0, thick=14, linestyle = 2, charsize=1.2
   xyouts, -25.0, -8.8, 'LEDE', charsize=2.2, charthick=8.
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.01.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   z = 2.40     !!!  F i r s t   BOSS  BIN  !!!
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.70, 0.42, 0.98, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       ytickformat='(a1)'
;       xtickformat='(a1)' , $


xyouts, x_xyouts, y_xyouts, '!82.20<z<2.60!3', charsize=2.2, charthick=6

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
;   plotsym, 8, plot_sym_size_R06/1.2, /fill
   plotsym, 8, plot_sym_size_MMT, /fill

   w_MMT = where(z_MMT gt 2.20 and z_MMT lt 2.60 and NQ_MMT gt 0 and $
                 Abs_Mag_MMT lt -22.8, N)

;   oplot, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], color=MMT_color, ps=8
;   Mi_MMT_err = Abs_Mag_MMT[w_MMT] - Abs_Mag_MMT[w_MMT]
;   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_up[w_MMT], $
;               /hibar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
;   oploterror, Abs_Mag_MMT[w_MMT], log_Phi_MMT[w_MMT], Mi_MMT_err, MMT_delta_down[w_MMT], $
 ;              /lobar, errthick=12, psym=8, color=MMT_color, ERRcolor=MMT_color
endif
loadct, 6

if (choice_2SLAQ_points eq 'y') then begin
   plotsym, 8, plot_sym_size_R06/1.2, /fill
   w_2SLAQ = where(z_bar gt 2.20 and z_bar lt 2.60, N_2SLAQ)
   oplot, Mi_2SLAQ[w_2SLAQ], LOG_PHI_2SLAQ[w_2SLAQ], ps=8, color=132
endif

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

loadct, 6
if (choice_BOSS_points eq 'y') then begin
   ;;Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[5,*], limit_lines, color=boss_color, linestyle=1, thick=12
endif 

plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_DR9_points eq 'y') then begin
   
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt iMag_limit_line[5,*] and log_Phi_boss gt -500, N) 
   if choice_DR9_nowgt eq 'y' then begin
      plotsym, 0, plot_sym_size_BOSS
      oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=boss_color
   endif
   
   www =  where(z_bin_boss_wgt gt 2.20 and z_bin_boss_wgt lt 2.60 and $
                Abs_mag_bin_boss_wgt lt iMag_limit_line[5,*] and log_Phi_boss gt -500, N) 
   
   plotsym, 0, plot_sym_size_BOSS, /fill
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], psym=8, color=boss_color
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_up_wgt[www],  $
                  /hibar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_down_wgt[www], $
                  /lobar, errthick=12, psym=8, color=boss_color, errcolor=boss_color
   endif
endif

plot_sym_size_S82 = plot_sym_size_BOSS
plotsym, 0, plot_sym_size_S82, /fill
loadct, 6
s82_color = 160.
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  ge 2.20 and z_bin_boss_S82  lt 2.60 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=s82_color
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=s82_color, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=s82_color, ERRcolor=s82_color
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=s82_color, errcolor=s82_color
   endif
endif

if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 2.20 and z_bin_boss21 lt 2.60 and $ 
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -24.64, N) 
endif


if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE gt 2.35 and z_ple lt 2.45, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   
   w_mLDDE = where(z_mLDDE gt 2.35 and z_mLDDE lt 2.45, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   
   w_LEDE = where(z_LEDE gt 2.35 and z_LEDE lt 2.45, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.40.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; z = 2.80 
;;
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.14, 0.42, 0.42], $
       xstyle=1, ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, ythick=ythick, $
       charsize=charsize, charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN ; , $
;       xtickformat='(a1)';, ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!82.60<z<3.00!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] -Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8


plotsym, 0, 1.5
plotsym, 0, plot_sym_size_BOSS, /fill

loadct, 6
if (choice_BOSS_points eq 'y') then begin
   ;;Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[6,*], limit_lines, color=boss_color, linestyle=1, thick=12
endif 

plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 3.00 and $
                Abs_mag_bin_boss le iMag_limit_line[6,*] and log_Phi_boss gt -500, N) 

   if choice_DR9_nowgt eq 'y' then begin
      plotsym, 0, plot_sym_size_BOSS
      oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=boss_color
   endif
   
   www =  where(z_bin_boss_wgt ge 2.60 and z_bin_boss_wgt lt 3.00 and $
                Abs_mag_bin_boss_wgt le iMag_limit_line[6,*] and log_Phi_boss gt -500, N) 

   plotsym, 0, plot_sym_size_BOSS, /fill
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], psym=8, color=boss_color

   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_up_wgt[www],  $
                  /hibar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_down_wgt[www], $
                  /lobar, errthick=12, psym=8, color=boss_color, errcolor=boss_color
   endif
   ;; 2.20<z<2.60 line 
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
endif


plotsym, 0, plot_sym_size_S82, /fill
;s82_color = 60
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.60 and z_bin_boss_S82  lt 3.00 and $
                Abs_mag_bin_boss_S82  le -25.03 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=s82_color
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=s82_color, ERRcolor=s82_color
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=s82_color, errcolor=s82_color
   endif
   www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=s82_color, thick=12
endif


if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 2.60 and z_bin_boss21 lt 3.00 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -25.05, N )
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE gt 2.75 and z_ple lt 2.85, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   
   w_mLDDE = where(z_mLDDE gt 2.75 and z_mLDDE lt 2.85, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   
   w_LEDE = where(z_LEDE gt 2.75 and z_LEDE lt 2.85, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.80.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; z  = ~ 3.25 
;;
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.42, 0.14, 0.70, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
;       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

xyouts, x_xyouts, y_xyouts, '!83.00<z<3.50!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,      Mi_z2[w], log_PhiR06[w], psym=8
Mi_z2_err = Mi_z2[w] - Mi_z2[w]
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_up[w],   /hibar, errthick=12, psym=8
oploterror, Mi_z2[w], log_PhiR06[w], Mi_z2_err, R06_delta_down[w], /lobar, errthick=12, psym=8

loadct, 6  ;; red  = 60 
if (choice_BOSS_points eq 'y') then begin
   ;;Giving a *very rough* feel for where a BOSS mag-limit would be...
   oplot, iMag_limit_line[7,*], limit_lines, color=boss_color, linestyle=1, thick=12
endif 

plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_DR9_points eq 'y') then begin
   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
               Abs_mag_bin_boss lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 

   if choice_DR9_nowgt eq 'y' then begin
      plotsym, 0, plot_sym_size_BOSS
      oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=boss_color
     ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=boss_color, thick=6
   endif

   www = where(z_bin_boss_wgt gt 3.00 and z_bin_boss_wgt lt 3.50 and $
               Abs_mag_bin_boss_wgt lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 
   plotsym, 0, plot_sym_size_BOSS, /fill
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], psym=8, color=boss_color

   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_up_wgt[www],  $
                  /hibar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
      oploterror, Abs_Mag_bin_boss_wgt[www], log_Phi_boss_wgt[www], Mi_boss_err, boss_delta_down_wgt[www], $
                  /lobar, errthick=12, psym=8, color=boss_color, errcolor=boss_color
   endif

   ;; 2.20<z<2.60 line 
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
   oplot, Abs_Mag_bin_boss_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
endif

plotsym, 0, plot_sym_size_S82, /fill
;s82_color = 60
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.00 and z_bin_boss_S82  lt 3.50 and $
                Abs_mag_bin_boss_S82  lt iMag_limit_line[7,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=s82_color
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=s82_color, ERRcolor=s82_color
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=s82_color, errcolor=s82_color
   endif
   www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.60 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=s82_color, thick=12
endif


if (choice_BOSS21_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   w21 =  where(z_bin_boss21 gt 3.00 and z_bin_boss21 lt 3.50 and $
                RAW_N_QSOS_boss21 gt 0 and Abs_Mag_bin_boss21 lt -25.45, N )
   ;; Again, doesn't really add much...
   ;oplot, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], color=boss21_color, ps=8
   ;Mi_boss_err = log_Num_den_boss21[w21] - log_Num_den_boss21[w21]
   ;oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
   ;            boss21_delta_up[w21], /hibar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
   ;oploterror, Abs_Mag_bin_boss21[w21], log_Num_den_boss21[w21], Mi_boss_err, $
   ;            boss21_delta_down[w21], /lobar, errthick=12, psym=8, color=boss21_color, ERRcolor=boss21_color
endif

if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 3.20 and z_ple lt   3.30, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 3.20 and z_mLDDE lt 3.30, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 3.20 and z_LEDE  lt 3.30, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif
if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z3.25.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=HRH07_color
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  G L O B A L   A X E S 
;;
;; mind that the positions are tied to the ranges 
;; of the final panel...

print, 'charsize    ', charsize/1.2
print, 'charthick ', charthick

xyouts, -11.0, -8.0,  $
        'log!I10!N !7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
        ORIENTATION=90, charsize=charsize/1.2, $
       charthick=charthick

xyouts, -22.5, -10.40,  $
        'M!Ii!N[z=2]', $
       charthick=charthick, charsize=charsize/1.2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick       = 8.0
charthick_small = 4.6

charsize       = 2.8
charsize_small = 1.7

loadct, 0
plotsym, 8, 1.4, /fill
legend, 'SDSS DR3' , $
        position=[-31.00, -5.75], box=0, psym=8, color=0, $
        charsize=charsize, charthick=charthick
xyouts, -32.20, -6.53, 'Richards et al. 2006',  $
        color=0, charsize=charsize_small, charthick =charthick_small

if choice_DR9_points eq 'y' then begin
;if (choice_BOSS_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   legend, 'BOSS DR9' , $
           position=[-31.00, -6.60], box=0, psym=8, color=boss_color, $
           charsize=charsize, charthick=charthick

   plotsym, 0, 1.5
   legend, ' "   no selec. func.' , $
           position=[-31.00, -7.20], box=0, psym=8, color=boss_color, $
           charsize=charsize, charthick =charthick
;   xyouts, -32.30, -7.70, ' ', $
;           charsize=charsize_small, charthick =charthick_small, color=black

   legend, ' ' , $
           position=[-31.00, -8.85], box=0, linestyle=1, color=boss_color, $
           charsize=1.5, charthick =12., thick= 16.0
   xyouts, -34.30, -9.15,'i=21.8 limit', $
        color=boss_color, charsize=charsize, charthick =charthick
endif

if (choice_MMT_points eq 'y') then begin
   loadct, clr_table
   plotsym, 8, 1.2, /fill
   legend, 'boss21+MMT' , $
           position=[-31.00, -7.90], box=0, psym=8, color=MMT_color, $
           charsize=charsize, charthick =charthick
xyouts, -32.20, -8.70, 'Palanque-Delabrouille et al. 2012',  $
        color=0, charsize=charsize_small, charthick =charthick_small

endif

charsize=2.0
charthick=4.8

loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
