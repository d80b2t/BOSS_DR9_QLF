;+
; NAME:
;    qlf_plot_Mi_z2pnt4
; 
; PURPOSE:
;    To plot Quasar Luminosity Functions
;    in M_i(z=~2.4)
;
; CALLING SEQUENCE:
;       .run qlf_plot_Mi_z2pnt4
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
;;  C H O I C E S ... 
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
   ;if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
   read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
   if choice_DR9_points eq 'y' then begin
     ;read, choice_DR9_nowgt,     PROMPT=' - Plot DR9 points with NO WGT?  y/n  '
      read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
   endif
endif   


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
;;  N. B. ::  ALL   WITH correction/weighting (wgt)
;;
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat', $   ;; used in the DR9 QLF paper!!
         z_bin_boss_wgt, Abs_mag_bin_boss, blah_boss_wgt, log_Phi_BOSS_wgt, raw_N_QSOs_boss_wgt, sigma_Phi_BOSS_wgt
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat  READ-IN', N_elements(z_bin_boss_wgt)
print

Mi_BOSS_wgt = Abs_mag_bin_boss 

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

boss_delta_up_Exp   = alog10((10^log_Phi_BOSS_Exp+sigma_Phi_BOSS_Exp))
boss_delta_up_Exp   = boss_delta_up_Exp - log_Phi_BOSS_Exp
boss_delta_down_Exp = alog10((10^log_Phi_BOSS_Exp-sigma_Phi_BOSS_Exp))
boss_delta_down_Exp = abs(boss_delta_down_Exp - log_Phi_BOSS_Exp)




;; 
;; Setting a sanity-check completeness limit. 
;; 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.
;;
boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.80

red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]

;; For BOSS in M_1450 at z=1.25, 2.40, 3.24, 4.25]
M_limit=[-21.413717, -23.154132, -23.948839,-24.642529] 

Mg_limit = M_limit 


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
x_max =  -32.50    ; -30.50

;; y-ranges
y_min = -10.00    ;; -9.20 
y_max = -2.70     ;; -4.70      ;; Used to be -5.00, so take care here!!

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

set_plot, 'ps'
!p.multi=0
device, filename='QLF_Mi_z2pnt4_3_DR9models_temp.ps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

boss_color=60
;;S82_color=160  ;; set below...
color_dr9_VdB =  green
color_dr9_exp =  blue


plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2
plot_sym_size_S82  = 1.6
plot_sym_size_DR9  = 1.6

      ps_fid = 3
plot_sym_VdB = 4 
plot_sym_exp = 8 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  z  ~  2.4
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
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       /nodata, $
       xtitle='!3M!Ii!N(z=2)', $
       ytitle='log!I10!N !7U!3(M!Ii!N(z=2)) [Mpc!E-3!N mag!E-1!N]'

per_mag = alog10(2.5)

loadct, clr_table
;;
;;   S  D  S  S    R 0 6
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plotsym, 8, plot_sym_size_R06/1.2, /fill
;oplot,  Mi_SDSS[w], log_PhiR06[w], psym=8
Mi_SDSS_err = Mi_SDSS[w] - Mi_SDSS[w]
oploterror, Mi_SDSS[w], log_PhiR06[w], Mi_SDSS_err, R06_delta_up[w],   $
            /hibar, errthick=12, psym=8, color=black, ERRcolor=black
oploterror, Mi_SDSS[w], log_PhiR06[w], Mi_SDSS_err, R06_delta_down[w], $
            /lobar, errthick=12, psym=8, color=black, ERRcolor=black


;;
;;    B O S S   D R 9
;;
plotsym, 0, plot_sym_size_BOSS, /fill
boss_color = red


;;
;; ''fiducial'' model (technically, the  "fiducial_tweak" model...
;; 
if (choice_DR9_points eq 'y') then begin
;   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
;               M_BOSS lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 
   www = where(z_bin_boss_wgt gt 2.20 and z_bin_boss_wgt lt 2.60 and $
               log_Phi_boss_wgt gt -50 and Mi_BOSS_wgt le -24.2, N) 
   
   oplot, Mi_BOSS_wgt[www], log_Phi_boss_wgt[www], psym=8,   color=boss_color
   ;oplot, M_gBOSS_wgt[www], log_Phi_boss_wgt[www], thick=12, color=boss_color
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = Mi_BOSS_wgt[www] - Mi_BOSS_wgt[www]
      oploterror, Mi_BOSS_wgt[www], log_Phi_boss_wgt[www], $
                  Mi_BOSS_err, boss_delta_up_wgt[www], $
                  /hibar, errthick=12, psym=8, color=boss_color, $
                  ERRcolor=boss_color
      oploterror, Mi_BOSS_wgt[www], log_Phi_boss_wgt[www], $
                  Mi_BOSS_err, boss_delta_down_wgt[www], $
                  /lobar, errthick=12, psym=8, color=boss_color, $
                  errcolor=boss_color
   endif
endif

;;
;; the older VdB model 
;; 
if (choice_DR9_points eq 'y') then begin
;   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
;               M_BOSS lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 
   www = where(z_bin_boss_VdB gt 2.20 and z_bin_boss_VdB lt 2.60 and $
               log_Phi_boss_VdB gt -50 and Mi_BOSS_VdB le -24.2, N) 
   
   oplot, Mi_BOSS_VdB[www], log_Phi_boss_VdB[www], psym=8,   color=green
   ;oplot, M_gBOSS_VdB[www], log_Phi_boss_VdB[www], thick=12, color=green
   plotsym, plot_sym_VdB, plot_sym_size_DR9, /fill

   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = Mi_BOSS_VdB[www] - Mi_BOSS_VdB[www]
      oploterror, Mi_BOSS_VdB[www], log_Phi_boss_VdB[www], $
                  Mi_BOSS_err, boss_delta_up_VdB[www], $
                  /hibar, errthick=12, psym=8, color=green, $
                  ERRcolor=green
      oploterror, Mi_BOSS_VdB[www], log_Phi_boss_VdB[www], $
                  Mi_BOSS_err, boss_delta_down_VdB[www], $
                  /lobar, errthick=12, psym=8, color=green, $
                  errcolor=green
   endif
endif

;;
;; the EXP-DUST model
;; 
if (choice_DR9_points eq 'y') then begin
;   www = where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and $
;               M_BOSS lt iMag_limit_line[7,0] and log_Phi_boss gt -500, N) 
   www = where(z_bin_boss_Exp gt 2.20 and z_bin_boss_Exp lt 2.60 and $
               log_Phi_boss_Exp gt -50 and Mi_BOSS_Exp le -24.2, N) 

   plotsym, plot_sym_exp, plot_sym_size_DR9, /fill
   oplot, Mi_BOSS_Exp[www], log_Phi_boss_Exp[www], psym=8,   color=blue
   ;oplot, M_gBOSS_Exp[www], log_Phi_boss_Exp[www], thick=12, color=blue
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = Mi_BOSS_Exp[www] - Mi_BOSS_Exp[www]
      oploterror, Mi_BOSS_Exp[www], log_Phi_boss_Exp[www], $
                  Mi_BOSS_err, boss_delta_up_Exp[www], $
                  /hibar, errthick=12, psym=8, color=blue, $
                  ERRcolor=blue
      oploterror, Mi_BOSS_Exp[www], log_Phi_boss_Exp[www], $
                  Mi_BOSS_err, boss_delta_down_Exp[www], $
                  /lobar, errthick=12, psym=8, color=blue, $
                  errcolor=blue
   endif
endif


;;
;;   2 S L A Q    Q S Os
;;
;;        Mg_2SLAQ, z_2SLAQ, NQ, log_phi_2SLAQ, log_phi_2SLAQ_down,
;;        log_phi_2SLAQ_up
w = where(z_2SLAQ ge 1.80 and z_2SLAQ lt 2.2, N)
;w = where(z_2SLAQ ge 2.2 and z_2SLAQ lt 2.6, N)
loadct, clr_table

color_2slaq = light_blue+32
sym_2slaq = 4

plotsym, sym_2slaq, plot_sym_size_R06*1.35, /fill
oplot,       Mi_2SLAQ[w], log_phi_2SLAQ[w], psym=8, color=color_2slaq

plotsym, sym_2slaq, plot_sym_size_R06*1.35, /fill
Mi_2SLAQ_err = Mi_2SLAQ[w] - Mi_2SLAQ[w]
oploterror, Mi_2SLAQ[w], log_phi_2SLAQ[w], $
            Mi_2SLAQ_err, log_phi_2SLAQ_up[w], $
           /hibar, errthick=12, psym=8, $
           color=color_2slaq, ERRcolor=color_2slaq
oploterror, Mi_2SLAQ[w], log_phi_2SLAQ[w], $
          Mi_2SLAQ_err, log_phi_2SLAQ_down[w], $
           /lobar, errthick=12, psym=8, $
         color=color_2slaq, errcolor=color_2slaq

charthick_lab = charthick
charsize_lab  = charsize/1.2

x_xyouts = -19.50    ;; to match the z~2.4 plot...
y_xyouts = -6.50     ;; -6.20 if plotting 6 models;  -6.40 otherwise...
x_off    =  0.50
y_off    =  0.35     ;; 0.25
offset   =  0.60

plotsym, sym_2slaq, 1.3, /fill
legend, '2SLAQ ' , $
        position=[x_xyouts+offset, y_xyouts-0.05+(5*offset)], $
        box=0, psym=8, color=color_2slaq, $
        charsiz=charsize_lab, charthick =charthick_lab

charsize_lab  = charsize/1.4
xyouts, x_xyouts-(0.8*offset), y_xyouts-0.05+(3.6*offset), $
        '!8(1.82<z<2.2)!6', $
        color=black, charsize=charsize_lab/1.1, charthick =charthick_lab/1.2

;legend, '!8(1.82<z<2.2)!6' , $
;        position=[x_xyouts+offset, y_xyouts-0.05+(4.5*offset)], $
;        box=0, psym=8, color=color_2slaq, $
 ;       charsize=charsize_lab, charthick =charthick_lab




;;
;;   z  ~  2 . 4 
;;
x_xyouts = -28.25 ;; -26.75
y_xyouts = -5.0   ;;  -5.50
offset = 0.00
;xyouts, x_xyouts, y_xyouts+(offset), '!8z~2.4!3', charsize=2.8, charthick=8, color=black
xyouts, x_xyouts+1.2, y_xyouts+(offset), '!82.2<z<2.6!3', charsize=2.6, charthick=8/1.2, color=black

charthick_lab = charthick
charsize_lab  = charsize

x_xyouts = -19.30    ;; to match the z~2.4 plot...
y_xyouts = -6.50     ;; -6.20 if plotting 6 models;  -6.40 otherwise...
x_off    =  0.50
y_off    =  0.35     ;; 0.25
offset   =  0.60

loadct, clr_table

plotsym, 0, 1.3, /fill
legend, 'fiducial' , $
        ;position=[xpos, ypos+(3*y_off)], $
        position=[x_xyouts+offset, y_xyouts-0.05-(1*offset)], $
        box=0, psym=8, color=boss_color, $
        charsiz=charsize_lab, charthick =charthick_lab

plotsym, plot_sym_exp, 1.3, /fill
legend, 'exp dust' , $
        position=[x_xyouts+offset, y_xyouts-0.05-(2*offset)], $
        box=0, psym=8, color=light_blue, $
        charsize=charsize_lab, charthick =charthick_lab

plotsym,  plot_sym_VdB, 1.3, /fill
legend, 'VdB lines' , $
        position=[x_xyouts+offset, y_xyouts-0.05-(3*offset)], $
        box=0, psym=8, color=color_dr9_VdB, $
        charsize=charsize_lab, charthick =charthick_lab



loadct, 0
!p.multi=0

device, /close
set_plot, 'X'

close, /all

end
