;+
; NAME:
;       qlf_plot_gband
; 
; PURPOSE:
;        To plot Quasar Luminosity Functions in the g-band
;
; CALLING SEQUENCE:
;       .run qlf_plot_gband
;
; INPUTS:
;        Croom09b_2SLAQ_QLF.dat  
;             - The QLF from 2SLAQ and SDSS, from Croom et al. (2009b)
;
; OUTPUTS:
;       .ps file
;
; COMMENTS:
;       /usr/common/rsi/lib/general/LibAstro/ 
;         
;-

red, omega0=0.30, omegalambda=0.70, h100=0.700


readcol, '../data/Croom09b_2SLAQ_QLF.dat', $ 
         Mg_2slaq, z_2slaq, NQ_2slaq, log_phi_2SLAQ, delta_log_phi_upper, delta_log_phi_lower


readcol, '../pro/My_QLF_gband_boss_temp.dat',  z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss

;; Setting a sanity-check completeness limit. 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.

;boss_glimit = 21.85
boss_glimit = 22.0
;boss_rlimit = 21.85
;boss_ilimit = 21.79

red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]
dlums = DLUMINOSITY(red_bins) / 1e6
Abs_gMag_limit = boss_glimit - (5 * alog10(dlums))  - 25.00 ; + kcor(fix(redphot/0.01))

gMag_limit_line = fltarr(11, 61)
limit_lines = (findgen(61)/10.)-10.0                            
for ii = 0L, 11-1 do gMag_limit_line[ii,*] = Abs_gMag_limit[ii]



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

;; x-ranges
;x_min = -18.001
x_min = -20.001
x_max = -30.50

;; y-ranges
y_min = -9.20
y_max = -5.10

;; xy_outs 
x_xyouts = -21.50
y_xyouts =  -5.70

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

choice_BOSS_points = 'y'
read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '

choice_BOSS_line = 'y'
;read, choice_BOSS_line, PROMPT=' - Plot BOSS line ?? y/n  '


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  RE-DOING   
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_gband_C09_z2_temp.ps', $
        xsize=14.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

;;    z = 0.49
;;   0.30 < z < 0.68   for Richards06
;;   0.40 < z < 0.68   for Croom09
w  =  where(z_2slaq gt 0.40 and z_2slaq lt 0.68) 
ww  =  where(z_2slaq gt 1.82 and z_2slaq lt 2.20) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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

loadct, 6
plotsym, 0, 1.4, /fill
oplot, Mg_2slaq[w],  log_phi_2SLAQ[w],  psym=8, color=60.
oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 0.11 and z_bin_boss lt 0.68) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[0,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!80.40<z<0.68!3', charsize=2.2, charthick=6

;;
;; z = 0.87
;;
w  =  where(z_2slaq gt 0.68 and z_2slaq lt 1.06) 

plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Mg_2slaq[w],  log_phi_2SLAQ[w],  psym=8, color=60
oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 0.68 and z_bin_boss lt 1.06) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[1,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!80.68<z<1.06!3', charsize=2.2, charthick=6

;;
;;  z = 1.25
w  =  where(z_2slaq gt 1.06 and z_2slaq lt 1.44) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
oplot,  Mg_2slaq[w],  log_phi_2SLAQ[w],  psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 1.06 and z_bin_boss lt 1.44) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[2,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.06<z<1.44!3', charsize=2.2, charthick=6

;;
;; z = 1.63
;;
w  =  where(z_2slaq gt 1.44 and z_2slaq lt 1.82)
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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

loadct, 6
plotsym, 0, 1.5, /fill
oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 1.44 and z_bin_boss lt 1.82) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[3,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.44<z<1.82!3', charsize=2.2, charthick=6


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           2 n d     R O W          
;;
;;;;;;;;;      z = 2.01         ;;;;;;;;;;;;;;;;;;;;;;;;;
w  =  where(z_2slaq gt 1.82 and z_2slaq lt 2.20) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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

loadct, 6
plotsym, 0, 1.5, /fill
oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 1.82 and z_bin_boss lt 2.20) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[4,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.82<z<2.20!3', charsize=2.2, charthick=6

;;
;; z = 2.40 
;;
w  =  where(z_2slaq gt 2.20 and z_2slaq lt 2.60) 
;ww  =  where(z_2slaq gt 1.82 and z_2slaq lt 2.20) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[5,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.60!3', charsize=2.2, charthick=6


;; z = 2.80 
ww  =  where(z_2slaq gt 1.82 and z_2slaq lt 2.20) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
       position=[0.56, 0.42, 0.77, 0.70], $
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
loadct, 6
plotsym, 0, 1.5, /fill
;; no 2SLAQ QSOs > 2.60!!
;oplot, Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60  
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.
www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 3.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[6,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!82.60<z<3.00!3', charsize=2.2, charthick=6


;; z  = ~ 3.25 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
;oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.
www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[7,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!83.00<z<3.50!3', charsize=2.2, charthick=6



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           3 r d     R O W          
;;
;;           z = 3.75
;;
;w  =  where(z_2slaq gt 2.20 and z_2slaq lt 2.60) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
;oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.
www =  where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[8,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!83.50<z<4.00!3', charsize=2.2, charthick=6

;;
;;
;;           z = ~ 4.25
;w  =  where(z_2slaq gt 2.20 and z_2slaq lt 2.60) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
;oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.
www =  where(z_bin_boss gt 4.00 and z_bin_boss lt 4.50) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[9,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!84.00<z<4.50!3', charsize=2.2, charthick=6


;;
;;           z = ~ 4.75
;;
;w  =  where(z_2slaq gt 2.20 and z_2slaq lt 2.60) 
plot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], $
       /nodata, $
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
loadct, 6
plotsym, 0, 1.5, /fill
;oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
oplot,  Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.
www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, gMag_limit_line[10,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!84.50<z<5.00!3', charsize=2.2, charthick=6


charsize  = 1.6
charthick = 6.2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Labels in the bottom right hand corner
;;
loadct, 6
plotsym, 8, 1.4, /fill
legend, 'Croom et al. (2009b)' , $
        position=[-30.80, -6.20], box=0, psym=8, color=60, $
        charsize=charsize, charthick=charthick
xyouts, -34.00, -6.80,' 1.82<z<2.20 line ', charsize=1.2, charthick =6.2
legend, ' ' , $
        position=[-30.80, -6.60], box=0, linestyle=0, color=60, $
        charsize=1.1, charthick =6.2, thick=4.0

;loadct, 6
;plotsym, 0, 1.5, /fill
;legend, 'R06 by NPR (hmmm...)' , $
 ;       position=[-31.00, -6.60], box=0, psym=8, color=60, $
  ;      charsize=1.6, charthick =4.2

loadct, 6
plotsym, 0, 1.5, /fill
legend, 'BOSS QSOs' , $
        position=[-30.80, -7.20], box=0, psym=8, color=160, $
        charsize=charsize, charthick=charthick

legend, '(Vmax) ' , $
        position=[-31.10, -7.60], box=0, color=160, $
        charsize=charsize, charthick=charthick


xyouts, -34.00, -8.50,' g=22.0 limit',   color=160, charsize=2.0, charthick =6.2
legend, ' ' , $
        position=[-30.80, -8.30], box=0, linestyle=1, color=160, $
        charsize=1.1, charthick =6.2, thick=6.0


loadct, 0




device, /close
set_plot, 'X'
close, /all

end
