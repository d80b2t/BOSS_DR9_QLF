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

red, omega0=0.30, omegalambda=0.70, h100=0.700 

;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent

readcol, '../../data/Richards06_Table04.dat', kcor_redshift, kcor, /silent
print, '  Richards06_kcor.dat READ-IN', n_elements(kcor)
print


;;
;;
;;  B O S S   D R 9
;;
;;
;; Just taking the "full" 2.2<z<2.6 points for that first bin/top left panel...
;;

;;
;;  THE   S T R I P E   8 2   DATA  - 5,476 Quasars strong (z<3.50)
;;
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_20120807.dat', $  ;; Used for the DR9 QLF paper...
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_McGkcor.dat', $
readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_McGkcor_rv1.dat', $
         z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Phi_BOSS_S82, $
         raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)

;;
;; The ``first'', 2.2<z<2.6 bin as a general guide/check
;; (would probably want this to be from the Stripe 82 data
;; until we're a bit happier with the Selec. Func. models...??!)

;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $ 
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_fiducial_grid.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat', $ 
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_temp.dat', $ 
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor_rv1.dat', $ 
         z_bin_boss_bin1, Abs_mag_bin_boss_bin1, blah_boss_bin1, log_Num_den_boss_bin1, $
         raw_N_boss_QSOs_bin1

;; 
;; The  ("full")   D R 9   DATA 
;; 
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_temp.dat', $
;
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_20120619.dat', $
;readcol, ',../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_fiducial_grid.dat,' $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_expdust.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_McGkcor.dat', $
readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_McGkcor_rv1.dat', $
         z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss, $
         raw_N_boss_QSOs,  sigma_Phi_BOSS

log_Phi_BOSS =  log_Num_den_boss

boss_delta_up   = alog10 ((10^log_Phi_BOSS +sigma_Phi_BOSS))
boss_delta_up   = boss_delta_up            - log_Phi_BOSS
boss_delta_down = alog10 ((10^log_Phi_BOSS - sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down      - log_Phi_BOSS)

;; 
;; The  D R 9   DATA  but with NO selection function correction
;; 
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_NoWgt_1210.6389v1.dat', $
readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_NOwgt_McGkcor_rv1.dat', $
         z_bin_boss_nocor, Abs_mag_bin_boss_nocor, blah_boss, log_Num_den_boss_nocor, $
         raw_N_boss_QSOs_nocor,  sigma_Phi_BOSS_nocor

log_Phi_BOSS_nocor =  log_Num_den_boss_nocor

boss_delta_up_nocor   = alog10 ((10^log_Phi_BOSS_nocor + sigma_Phi_BOSS_nocor))
boss_delta_up_nocor   = boss_delta_up_nocor            - log_Phi_BOSS_nocor
boss_delta_down_nocor = alog10 ((10^log_Phi_BOSS_nocor - sigma_Phi_BOSS_nocor))
boss_delta_down_nocor = abs(boss_delta_down_nocor      - log_Phi_BOSS_nocor)




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
XTICKLEN  = 0.03
YTICKLEN  = 0.05

;; For loadct, 6, 
;;  60 is red..
;; 160 is the teal colour...

loadct, 6
color_dr9 =  60
color_S82 = 160


;; x-ranges
;x_min = -18.001
;x_min = -22.001
x_min = -24.3001
x_max = -30.250

;; y-ranges
y_min = -9.20
y_max = -4.80   ;; Was -5.10 for arXiv:1210.6389v1

;; xy_outs 
;x_xyouts = -23.50   ;; could also be e.g. -24.70, or something to put the label in the top-left corner more. 
x_xyouts = -25.50  
y_xyouts =  -5.40

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2
plot_sym_size_S82  = 1.6
plot_sym_size_DR9  = 1.6

print
choice_SDSS_points = 'n'
;read, choice_SDSS_points, PROMPT='- Plot SDSS points?? y/n  '

print
choice_BOSS_points = 'y'
read, choice_BOSS_points, PROMPT='- Plot BOSS points?? y/n  '

choice_S82_points       = 'n'
choice_S82_errorbars    = 'n'
choice_DR9_points       = 'n'
choice_DR9_errorbars    = 'n'
choice_DR9_points_nocor = 'n'
choice_DR9_line         = 'n'

if choice_BOSS_points eq 'y' then begin
   read, choice_S82_points, PROMPT='   - Plot S82 points?? y/n  '
;   if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
   if choice_S82_points eq 'y' then  choice_S82_errorbars = 'y'

   read, choice_DR9_points, PROMPT='   - Plot DR9 points?? y/n  '
;   if choice_DR9_points eq 'y' then read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
   if choice_DR9_points eq 'y' then choice_DR9_errorbars = 'y'
;  if choice_DR9_points eq 'y' then  read, choice_DR9_points_nocor, PROMPT=' - Plot the DR9 points W/O corrections?? y/n 

   read, choice_DR9_line, PROMPT='   - Plot DR9 line (2.2<z<2.3) ?? y/n  '
endif

dr9_line_thick = 8.  ;; used to be 22.

print
choice_models      = 'n'
read, choice_models, PROMPT=' - Plot model lines?? y/n  '

print
choice_first_panel = 'n'
;read, choice_first_panel, PROMPT=' - Plot first panel?? y/n  '

print
choice_last_panel = 'n'
;read, choice_first_panel, PROMPT=' - Plot first panel?? y/n  '


print
print
print, ' choice_BOSS_points = ', choice_BOSS_points  , ' (default)'
print, ' choice_S82_points  = ', choice_S82_points
print, ' choice_DR9_points  = ', choice_DR9_points
print, ' choice_DR9_line    = ', choice_DR9_line
print
print, ' choice_models      = ', choice_models , ' (default)'
print
print, ' choice_first_panel = ', choice_first_panel
print 
print

;;
;;
;;  M o d e l   f i t s...
;;
;;
choice_PLE_points    = 'n'
choice_LDDE_points   = 'n'
choice_LEDE_points   = 'n'
choice_models_HRH07  = 'n'
if (choice_models eq 'y') then begin
   read, choice_PLE_points,   PROMPT=' - Plot PLE  model??  y/n  '
;   read, choice_LDDE_points,  PROMPT=' - Plot LDDE model??  y/n  '
   read, choice_LEDE_points,  PROMPT=' - Plot LEDE model??  y/n  '
  ; read, choice_models_HRH07, PROMPT=' - Plot HRH07     ??  y/n  '

   if choice_PLE_points eq 'y' then begin
   ;   readcol, '../pro/models/Croom09b_PLE_temp.dat', $
;      readcol, '../../pro/models/Croom09b_PLE_toz5.dat', $
;      readcol, 'chi_sq_PLE_model_McG_2.2z3.5.dat', $
      readcol, 'chi_sq_PLE_model_McG_2.2z3.5_rv1.dat', $
               ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   endif
   
   if choice_LDDE_points eq 'y' then begin
      readcol, '../pro/models/Croom09b_mLDDE_temp.dat', $
               ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE
   endif
   
   if choice_LEDE_points eq 'y' then begin
      ;; Difference in format from Croom09b to Ross13 since
      ;; different equations and LEDE params between both...
;      readcol, '../pro/models/Croom09b_LEDE_temp.dat', $
 ;              ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
;      readcol, 'chi_sq_loglinear_LEDE_model_McG_2.2z3.5_temp.dat', $
      readcol, 'chi_sq_loglinear_LEDE_model_McG_2.2z3.5_rv1.dat', $
               ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, c1_LEDE, c2_LEDE,  alpha, Phi_LEDE
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
;!p.multi=[0,3,4]
!p.multi=[0,2,5]

device, filename='QLF_iband_R06_z2_narrowZ_temp.ps', $
;        xsize=18.0, ysize=14.0, $  ;; good for 3 by 4
        xsize=18.0, ysize=12.0, $  ;; good for 2 by 5
;        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_one
;;  ``First bin'', top, left,  
;;    giving a 1:1 comparison to  R06...
;; 
;;     2.20 < z < 2.60
;;
;;   THIS DOESN'T GET PLOTTED 
;;   IF choice_first_panel EQ N
;;   AND WANT THE 2x5 grid...
;;
if choice_first_panel eq 'y' then begin
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.70, 0.35, 0.98], $  ;; for 3 x 4 grid
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xstyle=1,      ystyle=1, $
       xthick=xthick, ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       /nodata, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', ytickformat='(a1)', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

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

endif

;; print, (.96-0.06)/5.
;; xmin_vals = [0.06, 0.24, 0.42, 0.60, 0.78, 0.96]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_two
;;
;; 2.20 < z < 2.30
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
;       position=[0.35, 0.70, 0.56, 0.98], $ ;; good for a 3x4 grid...
       position=[0.06, 0.56, 0.24, 0.96], $ ;; good for a 2x5 grid...
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       /nodata, $
       xthick=xthick,     ythick=ythick, $
       charsize=charsize, charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
       xtickformat='(a1)'

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[1,*], limit_lines, color=color_dr9, linestyle=1, thick=12
;xyouts, x_xyouts, y_xyouts, '!82.20<z<2.30!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8


;;
;;  S t r i p e     8 2     data 
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.20 and z_bin_boss_S82  lt 2.30 and $
                Abs_mag_bin_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_dr9, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9
;;
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.20 and z_bin_boss_nocor lt 2.30 and $
                Abs_mag_bin_boss_nocor lt -24.5 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif


;;
;;   DR9  l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
endif

;;
;;  M o d e l s
;;
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

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.30!3', charsize=2.2, charthick=6



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;    panel_three
;;
;;  2.30 < z < 2.40 
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
;       position=[0.56, 0.70, 0.77, 0.98], $
       position=[0.24, 0.56, 0.42, 0.96], $ ;; good for a 2x5 grid..
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       /nodata, $
       xstyle=1,          ystyle=1, $
       xthick=xthick,     ythick=ythick, $
       charsize=charsize, charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[2,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.30<z<2.40!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

;;
;;  S t r i p e    8 2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.30 and z_bin_boss_S82  lt 2.40 and $
                Abs_mag_bin_boss_S82  lt -24.5845 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_dr9, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.30 and z_bin_boss lt 2.40 and $
                Abs_mag_bin_boss lt -24.5845 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.30 and z_bin_boss_nocor lt 2.40 and $
                Abs_mag_bin_boss_nocor lt -24.5 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   DR9  2.2<z<2.3 l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
endif


;;
;;  M o d e l s... 
;;
if (choice_models eq 'y') then begin
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
;       position=[0.77, 0.70, 0.98, 0.98], $
       position=[0.42, 0.56, 0.60, 0.96], $ ;; good for a 2x5 grid...
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       /nodata, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[3,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.40<z<2.50!3', charsize=2.2, charthick=6

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;plotsym, 0, 1.5, /fill
plotsym, 0, plot_sym_size_BOSS, /fill

;;
;;   S t r i p e   8 2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.40 and z_bin_boss_S82  lt 2.50 and $
                Abs_mag_bin_boss_S82  lt -24.6945 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;   D R 9  
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.40 and z_bin_boss lt 2.50 and $
                Abs_mag_bin_boss lt -24.6945 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.40 and z_bin_boss_nocor lt 2.50 and $
                Abs_mag_bin_boss_nocor lt -24.5 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   DR9  l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.14, 0.42, 0.35, 0.70], $
       position=[0.60, 0.56, 0.78, 0.96], $ ;; good for a 2x5 grid...
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       /nodata, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
;       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'
       xtickformat='(a1)', ytickformat='(a1)'


plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;plotsym, 0, 1.5, /fill
plotsym, 0, plot_sym_size_BOSS, /fill
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.50 and z_bin_boss_S82  lt 2.60 and $
                Abs_mag_bin_boss_S82  lt -24.80 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9    
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.50 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.80 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.50 and z_bin_boss_nocor lt 2.60 and $
                Abs_mag_bin_boss_nocor lt -24.80 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   DR9  l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.35, 0.42, 0.56, 0.70], $   ;; for 3 x 4 panels...
       ;position=[0.35, 0.42, 0.56, 0.70], $   ;; for 3 x 3 panels...
       position=[0.78, 0.56, 0.96, 0.96], $ ;; good for a 2x5 grid...
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, ythick=ythick, $
       xstyle=1,      ystyle=1, $
       charsize=charsize, charthick=charthick, $
       /nodata, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
;       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8


plotsym, 0, plot_sym_size_BOSS, /fill

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[5,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.60<z<2.70!3', charsize=2.2, charthick=6

;;
;; Stripe 82...
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.60 and z_bin_boss_S82  lt 2.70 and $
                Abs_mag_bin_boss_S82  lt -24.9 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9  
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 2.70 and $
                Abs_mag_bin_boss lt -24.9 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.60 and z_bin_boss_nocor lt 2.70 and $
                Abs_mag_bin_boss_nocor lt -24.9 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   D R 9    l i n e ...
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.56, 0.42, 0.77, 0.70], $ ;; good for a 3x4 grid..
       position=[0.06, 0.16, 0.24, 0.56], $ ;; good for a 2x5 grid...
       xstyle=1, ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       /nodata, $
       xthick=xthick, ythick=ythick, $
       charsize=charsize, charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
plotsym, 0, 1.5
;oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

;;
;;  S t r i p e   8 2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.70 and z_bin_boss_S82  lt 2.80 and $
                Abs_mag_bin_boss_S82  lt -24.99 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.70 and z_bin_boss lt 2.80 and $
                Abs_mag_bin_boss lt -24.99 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.70 and z_bin_boss_nocor lt 2.80 and $
                Abs_mag_bin_boss_nocor lt -24.99 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   D R 9    l i n e ...
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.77, 0.42, 0.98, 0.70], $
       position=[0.24, 0.16, 0.42, 0.56], $ ;; good for a 2x5 grid...
       xrange=[x_min, x_max], yrange=[y_min, y_max], $
       xthick=xthick,         ythick=ythick, $
       xstyle=1,              ystyle=1, $
       /nodata, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6  ;; red  = 60 

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[7,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.80<z<3.00!3', charsize=2.2, charthick=6

;;
;;  S t r i p e   8  2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.80 and z_bin_boss_S82  lt 3.00 and $
                Abs_mag_bin_boss_S82  lt -25.1375 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9   
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.80 and z_bin_boss lt 3.00 and $
                Abs_mag_bin_boss lt -25.1375 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 2.80 and z_bin_boss_nocor lt 3.00 and $
                Abs_mag_bin_boss_nocor lt -25.1375 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   D R 9    l i n e ...
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.14, 0.14, 0.35, 0.42], $
       position=[0.42, 0.16, 0.60, 0.56], $ ;; good for a 2x5 grid...
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, ythick=ythick, $
       xstyle=1,      ystyle=1, $
       /nodata, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'


plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[8,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.00<z<3.25!3', charsize=2.2, charthick=6

;;
;;  S 8 2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.00 and z_bin_boss_S82  lt 3.25 and $
                Abs_mag_bin_boss_S82  lt -25.3368 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9    
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.25 and $
                Abs_mag_bin_boss lt -25.3368 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 3.00 and z_bin_boss_nocor lt 3.25 and $
                Abs_mag_bin_boss_nocor lt -25.3368 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif

;;
;;   DR9  l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
;       position=[0.35, 0.14, 0.56, 0.42], $
       position=[0.60, 0.16, 0.78, 0.56], $ ;; good for a 2x5 grid...
       xstyle=1, ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, ythick=ythick, $
       /nodata, $
       charsize=charsize, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'


plotsym, 8, plot_sym_size_R06, /fill
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[9,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.25<z<3.50!3', charsize=2.2, charthick=6

;;
;;  S 8 2
;;
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.25 and z_bin_boss_S82  lt 3.50 and $
                Abs_mag_bin_boss_S82  lt -25.5137 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
   ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

;;
;;  D R 9
;;
plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.25 and z_bin_boss lt 3.50 and $
                Abs_mag_bin_boss lt -25.5137 and RAW_N_BOSS_QSOS gt 0, N )
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

;;
;;  D R 9    N o   C o r r e c t i o n
;;
plotsym, 0, plot_sym_size_DR9
if (choice_DR9_points eq 'y' and choice_DR9_points_nocor eq 'y') then begin
   www =  where(z_bin_boss_nocor gt 3.25 and z_bin_boss_nocor lt 3.50 and $
                Abs_mag_bin_boss_nocor lt -25.5137 and RAW_N_BOSS_QSOS_nocor gt 0, N )
   oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www], psym=8, color=color_dr9, thick=12
   ;oplot, Abs_Mag_bin_boss_nocor[www], log_Num_den_boss_nocor[www],        color=color_dr9, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err_nocor = log_Num_den_boss_nocor[www] - log_Num_den_boss_nocor[www]
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_up_nocor[www], $
                  /hibar, errthick=1, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Abs_Mag_bin_boss_nocor[www], log_Phi_boss_nocor[www], Mi_boss_err_nocor, boss_delta_down_nocor[www], $
                  /lobar, errthick=1, psym=8, color=color_dr9, errcolor=color_dr9
   endif
endif



;;
;;   DR9  l i n e 
;;
if (choice_DR9_line eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], $
          linestyle=0, color=color_dr9, thick=dr9_line_thick
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
if choice_last_panel eq 'y' then begin
w  =  where(z_R06 gt 3.50 and z_R06 lt 4.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
;       position=[0.56, 0.14, 0.77, 0.42], $
       position=[0.78, 0.16, 0.96, 0.56], $ ;; good for a 2x5 grid...
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
if choice_SDSS_points eq 'y' then oplot,  Mi_z2[w], log_PhiR06[w], psym=8

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[10,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.50<z<4.00!3', charsize=2.2, charthick=6

loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.50 and z_bin_boss_S82  lt 4.00 and $
                Abs_mag_bin_boss_S82  lt -25.8060 and RAW_N_BOSS_QSOS_S82  gt 0, N )
   if N gt 1 then begin 
      oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www], psym=8, color=color_S82
                                ;oplot, Abs_Mag_bin_boss_S82[www], log_Phi_BOSS_S82[www],         color=color_S82, thick=12
      
      if (choice_S82_errorbars eq 'y') then begin
         Mi_boss_err_S82 = log_Phi_BOSS_S82[www] - log_Phi_BOSS_S82[www]
         oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_up_S82[www],  $
                     /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
         oploterror, Abs_Mag_bin_boss_S82[www], log_Phi_boss_S82[www], Mi_boss_err_S82, boss_delta_down_S82[www], $   
                     /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
      endif
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00 and $
                Abs_mag_bin_boss lt -25.8060 and RAW_N_BOSS_QSOS gt 0, N )
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

;; if choice_last_panel eq 'y' then begin
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick=4.8
charsize=2.0
factor=1.2
loadct, 0

if choice_SDSS_points eq 'y' then begin
   plotsym, 8, 1.4, /fill
   legend, 'SDSS QSOs' , $
           position=[-31.00, -6.00], box=0, psym=8, color=0, $
           charsize=charsize*factor, charthick=charthick*(factor^2)
   xyouts, -31.20, -6.95,'(Richards+ 2006)',  color=0, $
           charsize=charsize*factor, charthick =charthick*(factor^2)
endif

loadct, 6
if (choice_BOSS_points eq 'y') then begin
   if (choice_DR9_points eq 'y') then begin
      loadct, 6
      plotsym, 0, 1.6, /fill
      legend, 'BOSS QSOs DR9' , $
              ;; Good when SDSS R06 is also plotted...
              ;position=[-31.00, -7.20], box=0, psym=8, color=color_dr9, $
              position=[-30.20, -5.60], box=0, psym=8, color=color_dr9, $
              charsize=charsize*factor, charthick =charthick*(factor^2)
      xyouts, -30.80, -6.45,'(color-selected)',  color=color_dr9, $
              charsize=charsize*factor, charthick =charthick*(factor^2)
   endif

   if (choice_S82_points eq 'y') then begin
      loadct, 6
      plotsym, 0, 1.6, /fill
      legend, 'BOSS QSOs S82' , $
              ;; Good when SDSS R06 is also plotted...
              ;position=[-31.00, -7.75], box=0, psym=8, color=color_S82, $
              position=[-30.20, -6.80], box=0, psym=8, color=color_S82, $
              charsize=charsize*factor, charthick =charthick*(factor^2)
          xyouts, -30.80, -7.65,'(variability-selected)',  color=color_S82, $
              charsize=charsize*factor, charthick =charthick*(factor^2)
   endif
   
;   xyouts, -31.50, -8.40,' fiducial_linetweak', $
;;   xyouts, -31.50, -8.40,' fiducial_grid', $
;   xyouts, -31.50, -8.40,' McGreer k-corr', $
;           xyouts, -31.50, -8.40,'expdust_grid', $
 ;         color=color_dr9, charsize=charsize*1.4, charthick =charthick*1.4

   legend, ' ' , $
           position=[-30.20, -8.10], box=0, linestyle=0, color=color_dr9, $
           charsize=1.1, charthick =12., thick= 16.0
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
   xyouts, -32.10, -8.30,'2.2<z<2.3 data', $
           color=color_dr9, charsize=charsize*1.2, charthick =charthick*1.4

 ;  legend, ' ' , $
  ;         position=[-31.00, -8.70], box=0, linestyle=1, color=color_dr9, $
   ;        charsize=1.1, charthick =12., thick= 16.0
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
;   xyouts, -34.00, -9.0,'i=21.8 limit', $
    ;       color=color_dr9, charsize=charsize*1.4, charthick =charthick*1.4

   legend, ' ' , $
           position=[-30.20, -8.70], box=0, linestyle=0, color=black, $
           charsize=1.1, charthick =12., thick= 16.0
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
   xyouts, -32.10, -9.0,'PLE', $
           color=black, charsize=charsize*1.4, charthick =charthick*1.4
   legend, ' ' , $
           position=[-33.50, -8.70], box=0, linestyle=1, color=black, $
           charsize=1.1, charthick =12., thick= 16.0
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
   xyouts, -35.40, -9.0,'LEDE', $
           color=black, charsize=charsize*1.4, charthick =charthick*1.4

endif


   
loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
