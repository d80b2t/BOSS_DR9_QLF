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
print, ' red, omega0=0.30, omegalambda=0.70, h100=0.700 '
print
print

;;
;;  S D S S   D R 3
;;
readcol, '../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent

;; My DR3...
;readcol, 'My_QLF_iband_20090717.dat', z_bin, Abs_mag_bin, blah, log_Num_den
;readcol, '../data/My_QLF_iband_20090826.dat',   z_bin, Abs_mag_bin, blah, log_Num_den
readcol, '../pro/qlf/My_QLF_iband_dr3_temp.dat', $
         z_bin_dr3, Abs_mag_bin_dr3, blah_dr3, log_Num_den_dr3

readcol, '../data/Richards06_Table04.dat', kcor_redshift, kcor, /silent
print, '  Richards06_kcor.dat READ-IN', n_elements(kcor)
print

;;
;;
;;  B O S S   D R 9
;;
;;
;; Just taking the "full" 2.2<z<2.6 points for that first bin/top left panel...
;;

;;  THE   S T R I P E   8 2   DATA  - 6,250 Quasars strong
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_S82.dat', $         ;; USE THIS FOR THE S82 plots!!!
readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_S82_wExtinc.dat', $  ;; 
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_temp.dat', $
         z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82

log_Phi_BOSS_S82 =  log_Num_den_boss_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)


;; THE `First'', 2.2<z<2.6 bin as a general guide/check
;; (would probably want this to be from the Stripe 82 data
;; until we're a bit happier with the Selec. Func. models...??!)
readcol, '../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $ 
         z_bin_boss_bin1, Abs_mag_bin_boss_bin1, blah_boss_bin1, log_Num_den_boss_bin1, raw_N_boss_QSOs_bin1

;; 
;; The  ("full")   D R 9   DATA 
;; 
readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_temp.dat', $
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_temp.dat', $
;
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_0pnt85comp.dat', $
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_0pnt85comp.dat', $
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_oldfiducial_grid_0pnt85.dat', $
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_oldfiducial_grid_newlinesvrev_0pnt85.dat', $
; readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_oldfiducial_grid_newlinesv5revx1d_0pnt85.dat', $
         z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss, raw_N_boss_QSOs,  sigma_Phi_BOSS

log_Phi_BOSS =  log_Num_den_boss

boss_delta_up   = alog10 ((10^log_Phi_BOSS +sigma_Phi_BOSS))
boss_delta_up   = boss_delta_up            - log_Phi_BOSS
boss_delta_down = alog10 ((10^log_Phi_BOSS - sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down      - log_Phi_BOSS)


;;
;; Setting a sanity-check completeness limit. 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.
;;
boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.80

;red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]
;
;red_bins = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]
red_bins = [2.40, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.90, 3.13, 3.35, 3.75, 4.00]

dlums = DLUMINOSITY(red_bins) / 1e6
Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ;+ kcor(fix(red_bins/0.01))

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


charsize  = 3.6
charthick = 4.8
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.04
YTICKLEN  = 0.06

;; For loadct, 6, 
;;  60 is red..
;; 160 is the teal colour...

loadct, 6
color_dr9 =  60
color_S82 = 160


;; x-ranges
;x_min = -18.001
x_min = -22.001
x_max = -30.50

;; y-ranges
y_min = -9.20
y_max = -5.10

;; xy_outs 
x_xyouts = -23.50
y_xyouts =  -5.70

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2
plot_sym_size_S82  = 1.6
plot_sym_size_DR9  = 1.6

print
choice_BOSS_points = 'y'
;read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '

choice_S82_points    = 'y'
choice_S82_errorbars = 'n'
choice_DR9_points    = 'n'
choice_DR9_errorbars = 'n'

if choice_BOSS_points eq 'y' then begin
   read, choice_S82_points, PROMPT=' - Plot S82 points?? y/n  '
   if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
   read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
   if choice_DR9_points eq 'y' then read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
endif


print
choice_models = 'n'
read, choice_models, PROMPT=' - Plot model lines?? y/n  '

print
print
print, ' choice_BOSS_points = ', choice_BOSS_points  , ' (default)'
print, ' choice_S82_points  = ', choice_S82_points
print, ' choice_DR9_points  = ', choice_DR9_points
print
print, ' choice_models      = ', choice_models , ' (default)'
print 
print

;;
;;
;;  M o d e l   f i t s...
;;
;;
choice_PLE_points    = 'n'
choice_LDDE_points   = 'y'
choice_LEDE_points   = 'n'
choice_models_HRH07  = 'n'
if (choice_models eq 'y') then begin
   read, choice_PLE_points,   PROMPT=' - Plot PLE  model??  y/n  '
   read, choice_LDDE_points,  PROMPT=' - Plot LDDE model??  y/n  '
   read, choice_LEDE_points,  PROMPT=' - Plot LEDE model??  y/n  '
   read, choice_models_HRH07, PROMPT=' - Plot HRH07     ??  y/n  '
   if choice_PLE_points eq 'y' then begin
      readcol, '../pro/models/Croom09b_PLE_temp.dat', $
               ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   endif
   
   if choice_LDDE_points eq 'y' then begin
      readcol, '../pro/models/Croom09b_mLDDE_temp.dat', $
               ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE
   endif
   
   if choice_LEDE_points eq 'y' then begin
      readcol, '../pro/models/Croom09b_LEDE_temp.dat', $
               ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
   endif
endif


;;
;;
;;  H R H 0 7   M O D E L S
;;
;;
;readcol, 'HRH07_Bband_z2.35.dat', $
 ;        band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07

;band_lum_sol = (10d^band_lum)/(3.9e33)
;Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_iband_R06_z2_narrowZ_temp.ps', $
        xsize=18.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_one
;;  ``First bin'', top, left,  
;;    giving a 1:1 comparison to  R06...
;; 
;;     2.20 < z < 2.60
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.70, 0.35, 0.98], $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xstyle=1, $
       ystyle=1, $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_BOSS_points eq 'y') then begin
   www =  where(z_bin_boss_bin1 gt 2.20 and z_bin_boss_bin1 lt 2.60 and Abs_mag_bin_boss_bin1 lt -24.5, N) 
   oplot, Abs_Mag_bin_boss_bin1[www], log_Num_den_boss_bin1[www], psym=8, color=color_dr9
endif

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      legend, pos=[-22.2,-7.5], ' ', box=0, thick=14, linestyle = 0, charsize=1.2
      xyouts, -25.8, -7.8, 'PLE', charsize=2.2, charthick=8.
   endif
   if choice_LDDE_points eq 'y' then begin
      legend, pos=[-22.2,-8.0], ' ', box=0, thick=14, linestyle = 1, charsize=1.2
      xyouts, -25.8, -8.3, 'LDDE', charsize=2.2, charthick=8.
   endif 
   if choice_LEDE_points eq 'y' then begin
      legend, pos=[-22.2,-8.5], ' ', box=0, thick=14, linestyle = 2, charsize=1.2
      xyouts, -25.8, -8.8, 'LEDE', charsize=2.2, charthick=8.
   endif
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[0,*], limit_lines, color=color_dr9, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.60!3', charsize=2.2, charthick=6



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_two
;;
;; 2.20 < z < 2.30
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.70, 0.56, 0.98], $
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

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[1,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.20<z<2.30!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.30 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_dr9, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.20 and z_ple lt 2.30, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.20 and z_mLDDE lt 2.30, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1 
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.20 and z_LEDE lt 2.30, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.25.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_three
;;
;;  2.30 < z < 2.40 
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.70, 0.77, 0.98], $
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

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[2,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.30<z<2.40!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8


loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.30 and z_bin_boss_S82  lt 2.40 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_dr9, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.30 and z_bin_boss lt 2.40 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

if (choice_models eq 'y') then begin
;;
;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.30 and z_ple lt 2.40, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.30 and z_mLDDE lt 2.40, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.30 and z_LEDE lt 2.40, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.35.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_four
;;
;;  2.40 < z < 2.50
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.77, 0.70, 0.98, 0.98], $
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

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[3,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.40<z<2.50!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;plotsym, 0, 1.5, /fill
plotsym, 0, plot_sym_size_BOSS, /fill

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.40 and z_bin_boss_S82  lt 2.50 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.40 and z_bin_boss lt 2.50 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif


if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.40 and z_ple lt 2.50, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.40 and z_mLDDE lt 2.50, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.40 and z_LEDE lt 2.50, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.45.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_five
;;
;;           2 n d     R O W          
;;
;;   2.50 < z 2.60    
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.42, 0.35, 0.70], $
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
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;plotsym, 0, 1.5, /fill
plotsym, 0, plot_sym_size_BOSS, /fill
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.50 and z_bin_boss_S82  lt 2.60 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.50 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif


;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[4,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.50<z<2.60!3', charsize=2.2, charthick=6

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.50 and z_ple lt 2.60, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.50 and z_mLDDE lt 2.60, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.50 and z_LEDE lt 2.60, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.55.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_six
;;
;;    2.60 < z < 2.70      
;;
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.42, 0.56, 0.70], $
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

plotsym, 0, plot_sym_size_BOSS, /fill

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[5,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.60<z<2.70!3', charsize=2.2, charthick=6

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.60 and z_bin_boss_S82  lt 2.70 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 2.70 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif


if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.60 and z_ple lt 2.70, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.60 and z_mLDDE lt 2.70, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.60 and z_LEDE lt 2.70, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.65.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_seven
;;
;;   2.70 < z < 2.80 
;;
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.42, 0.77, 0.70], $
       xstyle=1, ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, ythick=ythick, $
       charsize=charsize, charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 2.60 and z_bin_dr3 lt 3.00) 
loadct, 6
plotsym, 0, 1.5
;oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.70 and z_bin_boss_S82  lt 2.80 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.70 and z_bin_boss lt 2.80 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif


;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[6,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.70<z<2.80!3', charsize=2.2, charthick=6

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.70 and z_ple lt 2.80, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.70 and z_mLDDE lt 2.80, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.70 and z_LEDE lt 2.80, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.75.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   panel_eight
;;
;;  2.80 < z  < 3.00
;;
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.77, 0.42, 0.98, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6  ;; red  = 60 

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[7,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.80<z<3.00!3', charsize=2.2, charthick=6

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.80 and z_bin_boss_S82  lt 3.00 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.80 and z_bin_boss lt 3.00 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

if (choice_models eq 'y') then begin
   ;; Model lines...
   ;; (picking the z=2.90 model for the 2.80<z<3.00 bin...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 2.85 and z_ple lt 2.95, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 2.85 and z_mLDDE lt 2.95, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 2.85 and z_LEDE lt 2.95, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z2.90.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           3 r d     R O W          
;;
;;   3.00  <  z < 3.25
;;
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.14, 0.35, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[8,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.00<z<3.25!3', charsize=2.2, charthick=6


loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.00 and z_bin_boss_S82  lt 3.25 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.25 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 3.10 and z_ple lt 3.20, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 3.10 and z_mLDDE lt 3.20, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 3.10 and z_LEDE lt 3.20, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z3.15.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   panel_ten
;;
;;   3.25 < z< 3.50
;;
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.14, 0.56, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[9,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.25<z<3.50!3', charsize=2.2, charthick=6

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.25 and z_bin_boss_S82  lt 3.50 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.25 and z_bin_boss lt 3.50 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
   ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 3.30 and z_ple lt 3.40, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 3.30 and z_mLDDE lt 3.40, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 3.30 and z_LEDE lt 3.40, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z3.40.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_eleven
;;
;;      3.50 < z < 4.00
;;
w  =  where(z_R06 gt 3.50 and z_R06 lt 4.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.14, 0.77, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[10,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.50<z<4.00!3', charsize=2.2, charthick=6

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.50 and z_bin_boss_S82  lt 4.00 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Num_den_boss_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Num_den_boss_S82[www] - log_Num_den_boss_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   if N gt 0 then begin
      oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=color_dr9
      ;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www],         color=color_dr9, thick=12
      
      if (choice_DR9_errorbars eq 'y') then begin
         Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
         oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_up[www], $
                     /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
         oploterror, Abs_Mag_bin_boss[www], log_Phi_boss[www], Mi_boss_err, boss_delta_down[www], $
                     /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9
      endif
   endif
endif


if (choice_models eq 'y') then begin
   ;; Model lines...
   if choice_PLE_points eq 'y' then begin
      w_PLE   = where(z_PLE gt 3.70 and z_ple lt 3.80, N_PLE)
      oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0
   endif
   if choice_LDDE_points eq 'y' then begin
      w_mLDDE = where(z_mLDDE gt 3.70 and z_mLDDE lt 3.80, N_mLDDE)
      oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1
   endif
   if choice_LEDE_points eq 'y' then begin
      w_LEDE = where(z_LEDE gt 3.70 and z_LEDE lt 3.80, N_LEDE)
      oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
   endif
endif

if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/pdf/printed/Hopkins/HRH07_Bband_z3.75.dat', $
            band_lum_HRH07, a_HRH07, b_HRH07, bol_lum_HRH07, no_density_HRH07
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10, linestyle=0, color=128
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick=4.8
charsize=2.0
loadct, 0

plotsym, 8, 1.4, /fill
legend, 'SDSS QSOs' , $
        position=[-31.00, -6.00], box=0, psym=8, color=0, $
        charsize=charsize, charthick=charthick
xyouts, -31.20, -6.75,'(Richards+ 2006)',  color=0, charsize=charsize, charthick =charthick

loadct, 6
plotsym, 0, 1.5, /fill
;legend, 'R06 re-done by NPR' , $
;        position=[-31.00, -6.60], box=0, psym=8, color=60, $
;        charsize=charsize, charthick =charthick

if (choice_BOSS_points eq 'y') then begin
   if (choice_S82_points eq 'y') then begin
      loadct, 6
      plotsym, 0, 1.5, /fill
      legend, 'BOSS QSOs S82' , $
              position=[-31.00, -7.00], box=0, psym=8, color=color_S82, $
              charsize=charsize, charthick =charthick
   endif
   
   if (choice_DR9_points eq 'y') then begin
      loadct, 6
      plotsym, 0, 1.5, /fill
      legend, 'BOSS QSOs DR9' , $
              position=[-31.00, -7.50], box=0, psym=8, color=color_dr9, $
              charsize=charsize, charthick =charthick
   endif
   
   legend, ' ' , $
           position=[-31.00, -8.40], box=0, linestyle=1, color=color_dr9, $
           charsize=1.2, charthick =12., thick= 16.0
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
   xyouts, -34.00, -8.7,'i=21.8 limit', $
           color=color_dr9, charsize=charsize*1.4, charthick =charthick*1.4
endif


   
loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
