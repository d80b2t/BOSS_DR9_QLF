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

print
choice_BOSS_points = 'y'
;read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '

choice_S82_points    = 'y'
choice_S82_errorbars = 'y'
choice_DR9_points    = 'y'
choice_DR9_errorbars = 'y'

if choice_BOSS_points eq 'y' then begin
;   read, choice_S82_points, PROMPT=' - Plot S82 points?? y/n  '
 ;  if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
  ; read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
  ; if choice_DR9_points eq 'y' then read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
endif

print
print
print, ' choice_BOSS_points = ', choice_BOSS_points  , ' (default)'
print, ' choice_S82_points  = ', choice_S82_points
print, ' choice_DR9_points  = ', choice_DR9_points
print
print


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
;;  THE   S T R I P E   8 2   DATA  - 6,250 Quasars strong
;;
;;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_S82.dat', $         ;; 
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_wExtinc_20120622.dat', $  ;; Used this for the DR9 QLF paper...
readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_20120807.dat', $  ;; Us
         z_bin_boss_s82, Mi_boss_s82, blah_boss_s82, log_Phi_BOSS_S82, $
         raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)

;;
;; The ``first'', 2.2<z<2.6 bin as a general guide/check
;; (would probably want this to be from the Stripe 82 data
;; until we're a bit happier with the Selec. Func. models...??!)

readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $ 
         z_bin_boss_bin1, Mi_boss_bin1, blah_boss_bin1, log_Phi_BOSS_bin1, raw_N_boss_QSOs_bin1

readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_fiducial_grid.dat', $
         z_bin_boss_bin2, Mi_boss_bin2, blah_boss_bin2, log_Phi_BOSS_bin2, raw_N_boss_QSOs_bin2

readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_expdust.dat', $ 
         z_bin_boss_bin3, Mi_boss_bin3, blah_boss_bin3, log_Phi_BOSS_bin3, raw_N_boss_QSOs_bin3


;; 
;; The  ("full")   D R 9   DATA 
;; 
;readcol, '../pro/qlf/My_QLF_iband_boss_narrowZ_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_temp.dat', $
;
readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_20120619.dat', $
         z_bin_boss, Mi_boss, blah_boss, log_Phi_BOSS, raw_N_boss_QSOs,  sigma_Phi_BOSS

boss_delta_up   = alog10 ((10^log_Phi_BOSS +sigma_Phi_BOSS))
boss_delta_up   = boss_delta_up            - log_Phi_BOSS
boss_delta_down = alog10 ((10^log_Phi_BOSS - sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down      - log_Phi_BOSS)

readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_fiducial_grid.dat', $
         z_bin_boss, Mi_boss, blah_boss, log_Phi_BOSS_fid, raw_N_boss_QSOs, sigma_Phi_BOSS_fid
boss_delta_up_fid   = alog10 ((10^log_Phi_BOSS_fid + sigma_Phi_BOSS_fid))
boss_delta_up_fid   = boss_delta_up_fid            - log_Phi_BOSS_fid
boss_delta_down_fid = alog10 ((10^log_Phi_BOSS_fid - sigma_Phi_BOSS_fid))
boss_delta_down_fid = abs(boss_delta_down_fid      - log_Phi_BOSS_fid)

readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_wgt_expdust.dat', $
         z_bin_boss, Mi_boss, blah_boss, log_Phi_BOSS_exp, raw_N_boss_QSOs, sigma_Phi_BOSS_exp
boss_delta_up_exp   = alog10 ((10^log_Phi_BOSS_exp + sigma_Phi_BOSS_exp))
boss_delta_up_exp   = boss_delta_up_exp            - log_Phi_BOSS_exp
boss_delta_down_exp = alog10 ((10^log_Phi_BOSS_exp - sigma_Phi_BOSS_exp))
boss_delta_down_exp = abs(boss_delta_down_exp      - log_Phi_BOSS_exp)


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
;; limit_lines = (findgen(61)/10.)-10.0                            
limit_lines = (findgen(61)/20.) -1.5
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

color_dr9_fid =  135
color_dr9_exp =  220

;; x-ranges
;x_min = -18.001
x_min = -24.001
x_max = -30.50

;; y-ranges
;y_min = -9.20
;y_max = -5.10
y_min = -.9999
y_max =  1.00

;; xy_outs 
x_xyouts = -26.20   ;; could also be e.g. -24.70, or something to put the label in the top-left corner more. 
y_xyouts = 0.70

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2
plot_sym_size_S82  = 1.6
plot_sym_size_DR9  = 1.6

plot_sym_fid = 4 
plot_sym_exp = 8 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_iband_R06_z2_narrowZ_selfunc_divmodel_temp.ps', $
        xsize=18.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
       charsize=charsize*1.4, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
;       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'
       ytitle='log!I10!N (!7U!3!Iobs!N!3/!7U!3!Imodel!N!3)'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

divmodel = best_fit_PLE(2.25, Mi_z2[w])
oplot, Mi_z2[w], alog10((10^log_PhiR06[w]) /divmodel), psym=8


loadct, 6
plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_BOSS_points eq 'y') then begin
   www =  where(z_bin_boss_bin1 gt 2.20 and z_bin_boss_bin1 lt 2.60 and Mi_boss_bin1 lt -24.5, N) 
;   oplot, Mi_boss_bin1[www], log_Phi_BOSS_bin1[www], psym=8, color=color_dr9

   plotsym, 0, plot_sym_size_BOSS, /fill
   divmodel = best_fit_PLE(2.25, Mi_BOSS_bin1[www])
   oplot, Mi_boss_bin1[www], alog10((10^log_Phi_BOSS_bin1[www]) /divmodel), psym=8, color=color_dr9

   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   www =  where(z_bin_boss_bin2 gt 2.20 and z_bin_boss_bin2 lt 2.60 and Mi_boss_bin2 lt -24.5, N) 
   oplot, Mi_boss_bin2[www], alog10((10^log_Phi_BOSS_bin2[www]) /divmodel), psym=8, color=color_dr9_fid
   
   plotsym, plot_sym_exp, plot_sym_size_DR9, /fill
   www =  where(z_bin_boss_bin3 gt 2.20 and z_bin_boss_bin3 lt 2.60 and Mi_boss_bin3 lt -24.5, N) 
   oplot, Mi_boss_bin3[www], alog10((10^log_Phi_BOSS_bin3[www]) /divmodel), psym=8, color=color_dr9_exp
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
;       yrange=[-1.0,1.0], $
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
                Mi_boss_S82  lt iMag_limit_line[1,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.25, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.30 and $
                Mi_boss lt iMag_limit_line[1,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.25, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp

   endif
endif

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.30!3', charsize=2.2, charthick=6



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
                Mi_boss_S82  lt iMag_limit_line[2,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.35, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.30 and z_bin_boss lt 2.40 and $
                Mi_boss lt iMag_limit_line[2,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.35, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp
   endif
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
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.40 and z_bin_boss_S82  lt 2.50 and $
                Mi_boss_S82  lt iMag_limit_line[3,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.45, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.40 and z_bin_boss lt 2.50 and $
                Mi_boss lt iMag_limit_line[3,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.45, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp
   endif
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
       charsize=charsize*1.4, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat='(a1)', $
;       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'
       ytitle='log!I10!N (!7U!3!Iobs!N!3/!7U!3!Imodel!N!3)'


plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[4,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.50<z<2.60!3', charsize=2.2, charthick=6
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.50 and z_bin_boss_S82  lt 2.60 and $
                Mi_boss_S82  lt iMag_limit_line[4,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.55, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.50 and z_bin_boss lt 2.60 and $
                Mi_boss lt iMag_limit_line[4,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.55, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp

   endif
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
                Mi_boss_S82  lt iMag_limit_line[5,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.65, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 2.70 and $
                Mi_boss lt iMag_limit_line[5,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.65, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp

   endif
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

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[6,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!82.70<z<2.80!3', charsize=2.2, charthick=6


loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 2.70 and z_bin_boss_S82  lt 2.80 and $
                Mi_boss_S82  lt iMag_limit_line[6,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.75, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.70 and z_bin_boss lt 2.80 and $
                Mi_boss lt iMag_limit_line[6,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.75, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp

   endif
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
                Mi_boss_S82  lt iMag_limit_line[7,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.90, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.80 and z_bin_boss lt 3.00 and $
                Mi_boss lt iMag_limit_line[7,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.90, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp

   endif
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
       charsize=charsize*1.4, $
       charthick=charthick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtitle='!3M!Ii!N[z=2]', $
       ytitle='log!I10!N (!7U!3!Iobs!N!3/!7U!3!Imodel!N!3)'
;       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[8,*], limit_lines, color=color_dr9, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.00<z<3.25!3', charsize=2.2, charthick=6
loadct, 6
plotsym, 0, plot_sym_size_S82, /fill
if (choice_S82_points eq 'y') then begin
   www =  where(z_bin_boss_S82  gt 3.00 and z_bin_boss_S82  lt 3.25 and $
                Mi_boss_S82  lt iMag_limit_line[8,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(2.90, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.25 and $
                Mi_boss lt iMag_limit_line[8,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(2.90, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp
   endif
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
                Mi_boss_S82  lt iMag_limit_line[9,0] and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(3.35, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif

plotsym, 0, plot_sym_size_DR9, /fill
if (choice_DR9_points eq 'y') then begin
   www =  where(z_bin_boss gt 3.25 and z_bin_boss lt 3.50 and $
                Mi_boss lt iMag_limit_line[9,0] and RAW_N_BOSS_QSOS gt 0, N )
   
   divmodel = best_fit_PLE(3.35, Mi_boss[www])
   
   oplot, Mi_boss[www], alog10((10^log_Phi_boss[www])/divmodel), psym=8, color=color_dr9, thick=12
   plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_fid[www])/divmodel), psym=8, color=color_dr9_fid, thick=12
   plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
   oplot, Mi_boss[www], alog10((10^log_Phi_boss_exp[www])/divmodel), psym=8, color=color_dr9_exp, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Phi_boss[www] - log_Phi_boss[www]
      ratio     =  alog10((10^log_Phi_boss[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_up[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss[www] + boss_delta_down[www])) / divmodel)

      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9, ERRcolor=color_dr9
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9, errcolor=color_dr9

      ratio     =  alog10((10^log_Phi_boss_fid[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_fid[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_fid, plot_sym_size_DR9, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_fid, ERRcolor=color_dr9_fid
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_fid, errcolor=color_dr9_fid

      ratio     =  alog10((10^log_Phi_boss_exp[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_up_fid[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_exp[www] + boss_delta_down_fid[www])) /divmodel)
      plotsym, plot_sym_exp, plot_sym_size_DR9;, /fill
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_up, $
                  /hibar, errthick=12, psym=8, color=color_dr9_exp, ERRcolor=color_dr9_exp
      oploterror, Mi_boss[www], ratio, Mi_boss_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=color_dr9_exp, errcolor=color_dr9_exp
   endif
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
                Mi_boss_S82  lt -24.5 and RAW_N_BOSS_QSOS_S82  gt 0, N )

   divmodel = best_fit_PLE(3.75, Mi_BOSS_S82[www])
   oplot, Mi_boss_S82[www], alog10((10^log_Phi_BOSS_S82[www]) /divmodel), psym=8, color=color_S82
   
   if (choice_S82_errorbars eq 'y') then begin
      Mi_boss_err_S82 = log_Phi_boss_S82[www] - log_Phi_boss_S82[www]
      
      ratio     =  alog10((10^log_Phi_boss_S82[www]) /divmodel)
      diff_up   =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_up_S82[www])) / divmodel)
      diff_down =  ratio - alog10((10^(log_Phi_boss_S82[www] + boss_delta_down_S82[www])) / divmodel)
      
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_up,  $
                  /hibar, errthick=12, psym=8, color=color_S82, ERRcolor=color_S82
      oploterror, Mi_boss_S82[www], ratio, Mi_boss_err_S82, diff_down, $   
                  /lobar, errthick=12, psym=8, color=color_S82, errcolor=color_S82
   endif
endif





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick=4.8*1.4
charsize=2.0*1.2
loadct, 0

plotsym, 8, 1.4, /fill
;legend, 'SDSS QSOs' , $
 ;       position=[-31.00, -6.00], box=0, psym=8, color=0, $
  ;      charsize=charsize, charthick=charthick
;xyouts, -31.20, -0.75,'(Richards+ 2006)',  color=0, charsize=charsize, charthick =charthick


loadct, 6
plotsym, 0, 1.5, /fill
legend, 'BOSS QSOs S82' , $
        position=[-31.00, 0.50], box=0, psym=8, color=color_S82, $
        charsize=charsize, charthick =charthick

loadct, 6
plotsym, 0, 1.5, /fill
legend, 'fiducial_tweak' , $
        position=[-31.00, 0.20], box=0, psym=8, color=color_dr9, $
        charsize=charsize, charthick =charthick

plotsym, plot_sym_fid, 1.5, /fill
legend, 'fiducial' , $
        position=[-31.00, -0.10], box=0, psym=8, color=color_dr9_fid, $
           charsize=charsize, charthick =charthick

plotsym, plot_sym_exp, 1.5, /fill
legend, 'exp_dust' , $
        position=[-31.00, -0.40], box=0, psym=8, color=color_dr9_exp, $
           charsize=charsize, charthick =charthick

legend, ' ' , $
        position=[-30.70, -0.90], box=0, linestyle=1, color=color_dr9, $
        charsize=1.1, charthick =12., thick= 16.0
xyouts, -32.50, -1.0,'i=21.8 limit', $
        color=color_dr9, charsize=charsize*1.2, charthick =charthick*1.2


   
loadct, 0
!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
