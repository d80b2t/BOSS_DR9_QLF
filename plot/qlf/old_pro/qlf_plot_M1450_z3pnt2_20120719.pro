;+
; NAME:
;    qlf_plot
; 
; PURPOSE:
;    To plot Quasar Luminosity Functions
;    in M_1450 at z~3.2  
;
;    So, it seems that one key thing to do would be to plot, 
;    and compare, M_1450 magnitudes and results. 
;    The key equation seems simply to be that of R06, Eq. (3).  
;
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
choice_DR9_points    = 'n'
choice_DR9_errorbars = 'n'
choice_DR9_nowgt     = 'n'

;read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '
if choice_BOSS_points eq 'y' then begin
   choice_BOSS_line = 'y'
   ;read, choice_BOSS_line, PROMPT=' - Plot BOSS (data) line ?? y/n  '
   ;read, choice_S82_points, PROMPT=' - Plot S82 points?? y/n  '
   if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '

   read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
   if choice_DR9_points eq 'y' then begin
;      read, choice_DR9_nowgt,     PROMPT=' - Plot DR9 points with NO WGT?  y/n  '
      read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
   endif
endif   


print
print
choice_models = 'y'
read, choice_models, PROMPT=' - Plot (PLE, LDDE and LEDE) model lines?? y/n  '

print
choice_models_HRH07  = 'n'
read, choice_models_HRH07, PROMPT=' - Plot HRH07 models  ?? y/n  '
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
print
print, '../../data/Richards06_Table06.dat READ-IN', N_elements(z_R06)
print

;; Key Line, Fix for M_1450..., 
;; e.g. Equation (3) or Richards06
M_SDSS = Mi_z2 +1.486

R06_delta_up    = alog10((10^log_PhiR06)+(sigma_Phi*1e-9))
R06_delta_up = (log_PhiR06 -R06_delta_up)
R06_delta_down = alog10((10^log_PhiR06)-(sigma_Phi*1e-9))
R06_delta_down = (log_PhiR06 - R06_delta_down)

R06_err = alog10(sigma_Phi*1e-9)

;;
;;
;;   B O S S   D R 9
;;
;;
;; 
;;  Full DR9  (no correction)
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_temp.dat', $     
readcol, '../../pro/qlf/My_QLF_iband_boss_20120619.dat', $ ;; used in the DR9 QLF paper!!
         z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss, raw_N_QSOs_boss, sigma_Phi_BOSS
print
print, '../../pro/qlf/My_QLF_iband_boss_temp.dat  READ-IN', N_elements(z_bin_boss)
print

M_BOSS = Abs_mag_bin_boss +1.486

log_Phi_BOSS =  log_Num_den_boss
boss_delta_up = alog10 ((10^log_Phi_BOSS+sigma_Phi_BOSS))
boss_delta_up = boss_delta_up - log_Phi_BOSS
boss_delta_down = alog10 ((10^log_Phi_BOSS-sigma_Phi_BOSS))
boss_delta_down = abs(boss_delta_down - log_Phi_BOSS)

;; 
;;  WITH correction/weighting (wgt)
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_wgt_temp.dat', $
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $   ;; used in the DR9 QLF paper!!
         z_bin_boss_wgt, Abs_mag_bin_boss_wgt, blah_boss_wgt, log_Num_den_boss_wgt, raw_N_QSOs_boss_wgt, sigma_Phi_BOSS_wgt
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_temp.dat  READ-IN', N_elements(z_bin_boss)
print

M_BOSS_wgt = Abs_mag_bin_boss +1.486

log_Phi_BOSS_wgt =  log_Num_den_boss_wgt

boss_delta_up_wgt   = alog10((10^log_Phi_BOSS_wgt+sigma_Phi_BOSS_wgt))
boss_delta_up_wgt   = boss_delta_up_wgt - log_Phi_BOSS_wgt
boss_delta_down_wgt = alog10((10^log_Phi_BOSS_wgt-sigma_Phi_BOSS_wgt))
boss_delta_down_wgt = abs(boss_delta_down_wgt - log_Phi_BOSS_wgt)


;;
;; B O S S   S T R I P E   8 2
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_S82_stnd.dat', $ ;; USE THIS FOR THE S82 plots!!!
readcol, '../../pro/qlf/My_QLF_iband_boss_S82_March2012.dat', $  ;; 
        z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82
print
print, '../pro/qlf/My_QLF_iband_boss_S82_March2012.dat  READ-IN', N_elements(z_bin_boss_S82)
print

M_S82 = Abs_mag_bin_boss_S82 +1.486

log_Phi_BOSS_S82 =  log_Num_den_boss_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)

;;
;;  S W I R E    from  Siana et al. (2008)   results..
;;
readcol, 'SWIRE_Siana2008_z3pnt2.dat', M_SWIRE, phi_SWIRE

log_phi_SWIRE = alog10(phi_SWIRE)

;;
;;  V V D S      from    Bongiorno et al. (2007)   
;;
;readcol, 'VVDS_Bongiorno2007.dat', M_VVDS, phi_VVDS
;log_phi_VVDS = alog10(phi_VVDS)


;;
;;  C O S M O S   from      Masters et al. 1207.2154v1)
;;
readcol, 'COSMOS_Masters12_z3pnt2.dat', M_COSMOS, phi_COSMOS, sigma_phi_cosmos

log_phi_COSMOS = alog10(phi_COSMOS)

cosmos_delta_up   = alog10 ((10^log_Phi_cosmos + sigma_Phi_cosmos))
cosmos_delta_up   =     cosmos_delta_up        - log_Phi_cosmos
cosmos_delta_down = alog10 ((10^log_Phi_cosmos - sigma_Phi_cosmos))
cosmos_delta_down = abs(cosmos_delta_down      - log_Phi_cosmos)



;; Setting a sanity-check completeness limit. 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.

boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.80

red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]
;;dlums = DLUMINOSITY(red_bins) / 1e6
;;Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ; + kcor(fix(redphot/0.01))

;; For BOSS in M_1450 at z=1.25, 2.40, 3.24, 4.25]
M_limit=[-21.413717, -23.154132, -23.948839,-24.642529] 

;;iMag_limit_line = fltarr(11, 61)
;;limit_lines = (findgen(61)/10.)-10.0                            
;;for ii = 0L, 11-1 do iMag_limit_line[ii,*] = Abs_iMag_limit[ii] + 1.486



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


charsize  = 2.6
charthick = 6.8
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.04
YTICKLEN  = 0.06

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



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;  Actually plotting stuff...
;;
;;
;;  Q L F  i-band  
;;   Richards+06 as a template...
;;

;; x-ranges
;x_min = -18.001
x_min = -18.250
x_max = -30.50

;; y-ranges
y_min = -9.20
y_max = -5.00

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

set_plot, 'ps'
!p.multi=0
device, filename='QLF_M1450_z3pnt2_temp.ps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

boss_color=60
;;S82_color=160  ;; set below...

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; z  = ~ 3.25 
;;
plot,  M_SDSS, log_PhiR06, $
       position=[0.20, 0.20, 0.96, 0.96], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
       ytickformat='(a1)', $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       /nodata, $
       xtitle='!3M!I1450!N' ;, $
;       ytitle='log!I10!N !7U!3(M!I1450!N) [Mpc!E-3!N mag!E-1!N]'

;;
;;   S  D  S  S    R 0 6
;;
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
plotsym, 8, plot_sym_size_R06/1.2, /fill
;oplot,  M_SDSS[w], log_PhiR06[w], psym=8
M_SDSS_err = M_SDSS[w] -M_SDSS[w]
oploterror, M_SDSS[w], log_PhiR06[w], M_SDSS_err, R06_delta_up[w],   $
            /hibar, errthick=12, psym=8, color=black, ERRcolor=black
oploterror, M_SDSS[w], log_PhiR06[w], M_SDSS_err, R06_delta_down[w], $
            /lobar, errthick=12, psym=8, color=black, ERRcolor=black



loadct, 6  ;; red  = 60 

plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_DR9_points eq 'y') then begin
;   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
;               M_BOSS lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 
   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
               log_Phi_boss gt -500 and M_BOSS_wgt le -23.50, N) 
   
   oplot, M_BOSS_wgt[www], log_Num_den_boss_wgt[www], psym=8, color=boss_color
   ;oplot, M_BOSS_wgt[www], log_Num_den_boss_wgt[www], color=boss_color, thick=12
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = log_Num_den_boss[www] - log_Num_den_boss[www]
      oploterror, M_BOSS_wgt[www], log_Phi_boss_wgt[www], M_SDSS_err, boss_delta_up_wgt[www],  $
                  /hibar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
      oploterror, M_BOSS_wgt[www], log_Phi_boss_wgt[www], M_SDSS_err, boss_delta_down_wgt[www], $
                  /lobar, errthick=12, psym=8, color=boss_color, errcolor=boss_color
   endif
endif


;;
;;  S  W  I  R  E  
;; 
;;   from Siana et al. 2008, trying for green SWIRE squares...
PLOTSYM, 8, 1.8, /fill
oplot, M_swire, log_phi_swire, psym=8,  color = green

;;
;; VVDS points (from Bongiorno et al. 2007)
;;    (trying for yellow VVDS triangles...)
PLOTSYM, 4, 1.4, /fill
;oplot, M_vvds, log_phi_vvds, psym=8,  color = yellow

;;
;; COSMOS points (from Masters et al. 1207.2154v1)
;;
cosmos_sym = 4
loadct, clr_table
PLOTSYM, cosmos_sym, 1.4, /fill
oplot, M_cosmos, log_phi_cosmos, psym=4,  thick=28, color = purple, SYMSIZE=3.


;; xy_outs 
x_xyouts = -25.25
y_xyouts =  -5.70
;xyouts, x_xyouts, y_xyouts, '!83.00<z<3.50!3', charsize=2.2, charthick=6

;;
;;  L A B E L S 
;;
x_xyouts = -20.25
y_xyouts =  -9.10
offset = 0.34
;xyouts, x_xyouts, y_xyouts+(4*offset), '!8z~0.5', charsize=2.2, charthick=6, color=black
;xyouts, x_xyouts, y_xyouts+(3*offset), '!8z~1.25',  charsize=2.4, charthick=6, color=purple
;xyouts, x_xyouts, y_xyouts+(offset), '!8z~2.4!3', charsize=3.2, charthick=8, color=purple
xyouts, x_xyouts, y_xyouts+offset,     '!8z~3.2!3', charsize=2.4, charthick=6, color=blue
;xyouts, x_xyouts, y_xyouts,            '!8z~4!3',   charsize=2.4, charthick=6, color=green



x_xyouts = -19.50
y_xyouts = -6.5
offset = 0.40

;loadct, 6
plotsym, 0, plot_sym_size_BOSS, /fill
legend, '', position=[x_xyouts, y_xyouts], box=0, psym=8, color=red, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25, 'BOSS', $
        charsize=charsize, charthick =charthick*1.2, color=red

plotsym, 8, plot_sym_size_R06, /fill
legend, '', position=[x_xyouts, y_xyouts-offset], box=0, psym=8, color=black, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25-offset, 'SDSS', $
        charsize=charsize, charthick =charthick*1.2, color=black

;plotsym, sym_2slaq, plot_sym_size_R06, /fill
;legend, '' , position=[x_xyouts, y_xyouts-(2*offset)], box=0, psym=8, color=blue, $
;        charsize=charsize/1.4, charthick =charthick*1.2
;xyouts, x_xyouts-1.2, y_xyouts-0.25-(2*offset), '2SLAQ', $
;        charsize=charsize, charthick =charthick*1.2, color=blue

plotsym, 4, 1.6, /fill
legend, '' , position=[x_xyouts, y_xyouts+0.073-(2*offset)], box=0, psym=8, color=purple, $
        charsize=charsize/1.4, charthick =charthick*20.
plotsym, 5, 1.6, /fill
legend, '' , position=[x_xyouts, y_xyouts-0.073-(2*offset)], box=0, psym=8, color=purple, $
        charsize=charsize/1.4, charthick =charthick*20.
xyouts, x_xyouts-1.2, y_xyouts-0.25-(2*offset), 'COSMOS', $
        charsize=charsize, charthick =charthick*1.2, color=purple

loadct, 6  ;; red  = 60 
plotsym, 8, 1.6, /fill
legend, '', position=[x_xyouts, y_xyouts-(3*offset)], box=0, psym=8, color=green, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25-(3*offset), 'SWIRE', $
        charsize=charsize, charthick =charthick*1.2, color=green

loadct, clr_table

if (choice_models_HRH07 eq 'y') then begin
   legend, '' , position=[x_xyouts+0.5, y_xyouts-0.05-(4*offset)], box=0, line=0, $
           color=TURQUIOSE, thick =thick*2.
   xyouts, x_xyouts-1.2, y_xyouts-0.25-(4*offset), 'HRH07', $
           charsize=charsize, charthick =charthick*1.2, color=TURQUIOSE
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  H R H 0 7    m o d e l s
;;
if (choice_models_HRH07 eq 'y') then begin
   readcol, '/cos_pc19a_npr/BOSS/QLF/pro/models/HRH07/HRH07_1450A_z3.2.dat', $
            band_lum_HRH07, ABMag_HRH07, Monoflux_HRH07, bol_lum_HRH07, phi_HRH07_z3pnt2
   
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
;   oplot, Bmag_HRH07, alog10(no_density_HRH07), thick=10,
;   linestyle=0, color=HRH07_color

   per_mag = alog10(2.5)
   oplot, ABMag_HRH07, alog10(phi_HRH07_z3pnt2)-per_mag, thick=8, linestyle=0, color=turquiose ;;red
endif



if choice_models eq 'y' then begin
   w_PLE   = where(z_PLE   gt 3.20 and z_ple lt   3.30, N_PLE)
   oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=8, linestyle = 0

   w_mLDDE = where(z_mLDDE gt 3.20 and z_mLDDE lt 3.30, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=8, linestyle = 1

   w_LEDE = where(z_LEDE   gt 3.20 and z_LEDE  lt 3.30, N_LEDE)
   oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=8, linestyle = 2
endif



;;
;; SINGLE REDSHIFT BIN => NO REDSHIFT DEPENDENCE!! 
;; 
;;
;; Masters et al (1207.2154v1) best-fit (PLE) params
;; e.g. Table 2. 
;;
alpha        = -2.98
beta         = -1.73
Mstar_1450   = -25.54
Phi_star     = 2.65000e-7
log_Phi_star = alog10(2.65e-7)

Mstar_1450_faint  = -25.54 +0.68
Mstar_1450_bright = -25.54 -0.68


;; No. of mag bins
mag_bins = 60

;; Set up phi_fit...
Phi_fit        = fltarr(mag_bins)
Phi_fit_faint  = fltarr(mag_bins)
Phi_fit_bright = fltarr(mag_bins)

mag_fit = fltarr(mag_bins)

for jj=0L, mag_bins-1 do begin
   mag_fit[jj] = -32.0 + (jj*0.25)
   dmag        = mag_fit[jj] - Mstar_1450
   Phi_fit[jj] = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

   dmag        = mag_fit[jj] - Mstar_1450_faint
   Phi_fit_faint[jj]  = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
   dmag        = mag_fit[jj] - Mstar_1450_bright
   Phi_fit_bright[jj] = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
endfor

log_phi_fit    = alog10(phi_fit)
log_phi_faint  = alog10(phi_fit_faint)
log_phi_bright = alog10(phi_fit_bright)

loadct, 13
oplot, mag_fit, log_Phi_fit,    thick=thick*2.2, color=dark_blue
;oplot, mag_fit, log_Phi_faint,  thick=thick*1.6, color=dark_blue
;oplot, mag_fit, log_Phi_bright, thick=thick*1.6, color=dark_blue






   

loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
