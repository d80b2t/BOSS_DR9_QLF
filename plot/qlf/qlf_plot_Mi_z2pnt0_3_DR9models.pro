;+
; NAME:
;    qlf_plot_Mi_z2pnt0
; 
; PURPOSE:
;    To plot Quasar Luminosity Functions
;    in M_1450 at z~2.0  
;
;    So, it seems that one key thing to do would be to plot, 
;    and compare, M_i(z=2) magnitudes and results. 
;    The key equation seems simply to be that of R06, Eq. (3).  
;
; CALLING SEQUENCE:
;       .run qlf_plot_Mi_z2pnt0
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

print
print
Mi2_offset = (-0.25)   ;; Conroy&White12 Eq. 6, by default
print, ' M_i(z=2)        = M_g(z=2) -  0.25  set as defaut....'
print, ' M_i(z=2) + 0.25 = M_g(z=2)          set as defaut....'
;read,  Mi2_offset, PROMPT=' - Value of  M_i(z=2) to M_g(z=2)  offset??  '
print


;; 
;;
;;  C h o i c e s...
;;
;;
print
print

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


;; 
;;         D  A  T  A
;;
;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_SDSS, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor
print
print, '../../data/Richards06_Table06.dat READ-IN', N_elements(z_R06)
print

R06_delta_up    = alog10((10^log_PhiR06)+(sigma_Phi*1e-9))
R06_delta_up = (log_PhiR06 -R06_delta_up)
R06_delta_down = alog10((10^log_PhiR06)-(sigma_Phi*1e-9))
R06_delta_down = (log_PhiR06 - R06_delta_down)

R06_err = alog10(sigma_Phi*1e-9)


;;
;;   2 S L A Q    Q S O
;;
readcol, '../../data/Croom09b_2SLAQ_QLF.dat', $
         Mg_2SLAQ, z_2SLAQ, NQ, log_phi_2SLAQ, log_phi_2SLAQ_down, log_phi_2SLAQ_up
;; Mi2_offset = (-0.25)
Mi_2SLAQ = Mg_2SLAQ + Mi2_offset 


;;
;;
;;   B O S S   D R 9
;;
;;
;;  WITH correction/weighting (wgt)
;;
;readcol, '../pro/qlf/My_QLF_iband_boss_wgt_temp.dat', $
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_20120619.dat', $ 
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat', $   ;; used in the DR9 QLF paper!!
         z_bin_boss_wgt, Mi_BOSS_wgt, blah_boss_wgt, log_Phi_BOSS_wgt, raw_N_QSOs_boss_wgt, sigma_Phi_BOSS_wgt
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_temp.dat  READ-IN', N_elements(z_bin_boss_wgt)
print

boss_delta_up_wgt   = alog10((10^log_Phi_BOSS_wgt+sigma_Phi_BOSS_wgt))
boss_delta_up_wgt   = boss_delta_up_wgt - log_Phi_BOSS_wgt
boss_delta_down_wgt = alog10((10^log_Phi_BOSS_wgt-sigma_Phi_BOSS_wgt))
boss_delta_down_wgt = abs(boss_delta_down_wgt - log_Phi_BOSS_wgt)


;;
;;  ``VdB model..."
;;
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_fiducial_grid.dat', $
         z_bin_boss_VdB, Mi_BOSS_VdB, blah_boss_VdB, log_Phi_BOSS_VdB, raw_N_QSOs_boss_VdB, sigma_Phi_BOSS_VdB
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_fiducial_grid.dat  READ-IN', N_elements(z_bin_boss_VdB)
print

boss_delta_up_VdB   = alog10((10^log_Phi_BOSS_VdB+sigma_Phi_BOSS_VdB))
boss_delta_up_VdB   = boss_delta_up_VdB - log_Phi_BOSS_VdB
boss_delta_down_VdB = alog10((10^log_Phi_BOSS_VdB-sigma_Phi_BOSS_VdB))
boss_delta_down_VdB = abs(boss_delta_down_VdB - log_Phi_BOSS_VdB)


;;
;;  ``Exp dust model..."
;;
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_expdust.dat', $
         z_bin_boss_Exp, Mi_BOSS_Exp, blah_boss_Exp, log_Phi_BOSS_Exp, raw_N_QSOs_boss_Exp, sigma_Phi_BOSS_Exp
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_expdust.dat  READ-IN', N_elements(z_bin_boss_Exp)
print

boss_delta_up_VdB   = alog10((10^log_Phi_BOSS_VdB+sigma_Phi_BOSS_VdB))
boss_delta_up_VdB   = boss_delta_up_VdB - log_Phi_BOSS_VdB
boss_delta_down_VdB = alog10((10^log_Phi_BOSS_VdB-sigma_Phi_BOSS_VdB))
boss_delta_down_VdB = abs(boss_delta_down_VdB - log_Phi_BOSS_VdB)




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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
x_max = -32.50  ;;  -30.50

;; y-ranges
y_min = -10.00      ;; -9.20 
y_max = -2.70      ;; -4.70   Used to be -5.00 for the z=2.6 plot so take care here!!

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

set_plot, 'ps'
!p.multi=0
device, filename='QLF_Mi_z2pnt0_3_DR9models_temp.ps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

boss_color=60
;;S82_color=160  ;; set below...

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; z  = ~ 2.00 
;;
plot,  Mi_SDSS, log_PhiR06, $
       position=[0.20, 0.20, 0.96, 0.96], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
;       ytickformat='(a1)', $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       /nodata, $
       xtitle='!3M!Ii!N(z=2)', $
       ytitle='log!I10!N !7U!3(M!Ii!N(z=2)) [Mpc!E-3!N mag!E-1!N]'

per_mag = alog10(2.5)


;;
;;   S  D  S  S    R 0 6
;;
w  =  where(z_R06 gt 1.80 and z_R06 lt 2.20) 
plotsym, 8, plot_sym_size_R06/1.2, /fill
;oplot,  M_SDSS[w], log_PhiR06[w], psym=8
Mi_SDSS_err = Mi_SDSS[w] - Mi_SDSS[w]
oploterror, Mi_SDSS[w], log_PhiR06[w], Mi_SDSS_err, R06_delta_up[w],   $
            /hibar, errthick=12, psym=8, color=black, ERRcolor=black
oploterror, Mi_SDSS[w], log_PhiR06[w], Mi_SDSS_err, R06_delta_down[w], $
            /lobar, errthick=12, psym=8, color=black, ERRcolor=black


;;
;;   2 S L A Q    Q S Os
;;
;;        Mg_2SLAQ, z_2SLAQ, NQ, log_phi_2SLAQ, log_phi_2SLAQ_down, log_phi_2SLAQ_up
w = where(z_2SLAQ ge 1.80 and z_2SLAQ lt 2.2, N)

loadct, clr_table
if choice_models eq 'y' or choice_fits eq 'y' then loadct, 0
color_2slaq = light_blue+32
sym_2slaq = 4
plotsym, sym_2slaq, plot_sym_size_R06*1.35, /fill
oplot,       Mi_2SLAQ[w], log_phi_2SLAQ[w], psym=8, color=black

plotsym, sym_2slaq, plot_sym_size_R06, /fill
Mi_2SLAQ_err = Mi_2SLAQ[w] - Mi_2SLAQ[w]
oploterror, Mi_2SLAQ[w], log_phi_2SLAQ[w], Mi_2SLAQ_err,  log_phi_2SLAQ_up[w], $
            /hibar, errthick=12, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, Mi_2SLAQ[w], log_phi_2SLAQ[w], Mi_2SLAQ_err, log_phi_2SLAQ_down[w], $
            /lobar, errthick=12, psym=8, color=color_2slaq, errcolor=color_2slaq





;;
;;  L A B E L S 
;;
if (choice_models eq 'n' and choice_fits eq 'n') then begin

x_xyouts = -19.30
y_xyouts = -6.5
offset = 0.40

;loadct, 6
plotsym, 0, plot_sym_size_BOSS/1.2, /fill
legend, '', position=[x_xyouts, y_xyouts], box=0, psym=8, color=red, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25, 'boss21', $
        charsize=charsize, charthick =charthick*1.2, color=red

plotsym, 8, plot_sym_size_R06/1.2, /fill
legend, '', position=[x_xyouts, y_xyouts-(1*offset)], box=0, psym=8, color=black, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25-(1*offset), 'SDSS', $
        charsize=charsize, charthick =charthick*1.2, color=black

plotsym, sym_2slaq, plot_sym_size_R06, /fill
legend, '' , position=[x_xyouts, y_xyouts-(2*offset)], box=0, psym=8, color=color_2slaq, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-1.2, y_xyouts-0.25-(2*offset), '2SLAQ', $
        charsize=charsize, charthick =charthick*1.2, color=color_2slaq

endif



;;
;;    z  ~  2 . 0
;;
x_xyouts = -28.25   ;  -26.75
y_xyouts = -5.0     ;  -5.50
offset = 0.00
xyouts, x_xyouts, y_xyouts+offset, '!8z~2.0!3', charsize=2.8, charthick=8, color=black



!p.multi=0
device, /close
set_plot, 'X'


close, /all

end
