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

red, omega0=0.30, omegalambda=0.70, h100=0.700


readcol, '../data/Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor

;readcol, 'My_QLF_iband_20090717.dat', z_bin, Abs_mag_bin, blah, log_Num_den
;readcol, '../data/My_QLF_iband_20090826.dat',   z_bin, Abs_mag_bin, blah, log_Num_den
readcol, '../pro/My_QLF_iband_dr3_temp.dat',  z_bin_dr3, Abs_mag_bin_dr3, blah_dr3, log_Num_den_dr3

readcol, '../pro/My_QLF_iband_boss_temp.dat',  z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss


;; Setting a sanity-check completeness limit. 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.

boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.79

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


charsize  = 3.6
charthick = 4.8
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



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  RE-DOING   R I C H A R D S   '06
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_R06_z2_temp.ps', $
        xsize=14.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


;; z = 0.49
w  =  where(z_R06 gt 0.11 and z_R06 lt 0.68) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 0.11 and z_bin_dr3 lt 0.68) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 0.11 and z_bin_boss lt 0.68) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[0,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -26.50, -5.70, '!8z=0.49!3', charsize=2.2, charthick=6


;; z = 0.87
w  =  where(z_R06 gt 0.68 and z_R06 lt 1.06) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 0.68 and z_bin_dr3 lt 1.06) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 0.68 and z_bin_boss lt 1.06) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[1,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -26.50, -5.70, '!8z=0.87!3', charsize=2.2, charthick=6


;; z = 1.25
w  =  where(z_R06 gt 1.06 and z_R06 lt 1.44) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 1.06 and z_bin_dr3 lt 1.44) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 1.06 and z_bin_boss lt 1.44) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[2,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -26.50, -5.70, '!8z=1.25!3', charsize=2.2, charthick=6


;; z = 1.63
w  =  where(z_R06 gt 1.44 and z_R06 lt 1.82) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 1.44 and z_bin_dr3 lt 1.82) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 1.44 and z_bin_boss lt 1.82) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[3,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -23.50, -5.70, '!8z=1.63!3', charsize=2.2, charthick=6


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           2 n d     R O W          
;;
;;;;;;;;;      z = 2.01         ;;;;;;;;;;;;;;;;;;;;;;;;;
w  =  where(z_R06 gt 1.82 and z_R06 lt 2.20) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 1.82 and z_bin_dr3 lt 2.20) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 1.82 and z_bin_boss lt 2.20) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[4,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -23.50, -5.70, '!8z=2.01', charsize=2.2, charthick=6


;; z = 2.40 
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 2.20 and z_bin_dr3 lt 2.60) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[5,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -23.50, -5.70, '!8z=2.40!3', charsize=2.2, charthick=6


;; z = 2.80 
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 2.60 and z_bin_dr3 lt 3.00) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 3.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[6,*], limit_lines, color=160, linestyle=1, thick=4

;xyouts, -23.50, -8.40, '!8z!9A!82.80!3', charsize=2.2, charthick=6
xyouts, -23.50, -5.70, '!8z=2.80!3', charsize=2.2, charthick=6


;; z  = ~ 3.25 
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 3.00 and z_bin_dr3 lt 3.50) 
loadct, 6
red  = 60
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[7,*], limit_lines, color=160, linestyle=1, thick=4


xyouts, -23.50, -5.70, '!8z=3.25!3', charsize=2.2, charthick=6


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;           3 r d     R O W          
;;
;;           z = 3.75
;;
w  =  where(z_R06 gt 3.50 and z_R06 lt 4.00) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 3.50 and z_bin_dr3 lt 4.00) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[8,*], limit_lines, color=160, linestyle=1, thick=4

www =  where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

xyouts, -24.00, -5.70, '!8z=3.75!3', charsize=2.2, charthick=6


;;           z = ~ 4.25
w  =  where(z_R06 gt 4.00 and z_R06 lt 4.50) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 4.00 and z_bin_dr3 lt 4.50) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 4.00 and z_bin_boss lt 4.50) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[9,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -24.00, -5.70, '!8z=4.25!3', charsize=2.2, charthick=6


;;
;;           z = ~ 4.75
;;
w  =  where(z_R06 gt 4.50 and z_R06 lt 5.00) 
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

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin_dr3 gt 4.50 and z_bin_dr3 lt 5.00) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin_dr3[ww], log_Num_den_dr3[ww], psym=8, color=60

www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00) 
oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[10,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, -24.00, -5.70, '!8z=4.75!3', charsize=2.2, charthick=6

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Labels in the bottom right hand corner
;;
loadct, 0
plotsym, 8, 1.4, /fill
legend, 'Richards et al. (2006)' , $
        position=[-31.00, -5.90], box=0, psym=8, color=0, $
        charsize=1.6, charthick =4.2

loadct, 6
plotsym, 0, 1.5, /fill
legend, 'R06 by NPR (hmmm...)' , $
        position=[-31.00, -6.60], box=0, psym=8, color=60, $
        charsize=1.6, charthick =4.2

loadct, 6
plotsym, 0, 1.5, /fill
legend, 'BOSS QSOs' , $
        position=[-31.00, -7.50], box=0, psym=8, color=160, $
        charsize=1.6, charthick =4.2

xyouts, -31.00, -7.90,'g=22 limit',   color=160, charsize=1.0, charthick =4.2
legend, ' ' , $
        position=[-31.00, -7.90], box=0, linestyle=1, color=160, $
        charsize=1.0, charthick =4.2, thick=3.0


loadct, 0




device, /close
set_plot, 'X'
close, /all

end
