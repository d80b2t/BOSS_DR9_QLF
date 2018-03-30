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
readcol, '../pro/qlf/My_QLF_iband_dr3_temp.dat', z_bin_dr3, Abs_mag_bin_dr3, blah_dr3, log_Num_den_dr3

;readcol, '../pro/qlf/My_QLF_iband_boss_S82_stnd.dat', $
readcol, '../pro/qlf/My_QLF_iband_boss_wgt_temp.dat', $
;readcol, '../pro/qlf/My_QLF_iband_boss_temp.dat', $
         z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss, raw_N_boss_QSOs


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

;; xy_outs 
x_xyouts = -23.50
y_xyouts =  -5.70

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

choice_BOSS_points = 'y'
read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '

choice_BOSS_line = 'y'
;read, choice_BOSS_line, PROMPT=' - Plot BOSS line ?? y/n  '


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  RE-DOING   R I C H A R D S   '06
;;
set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_iband_R06_z2_temp.ps', $
        xsize=14.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


;;    z = 0.49
;;   0.30 < z < 0.68   for Richards06
;;   0.40 < z < 0.68   for Croom09
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
;; Want to oplot the BOSS 2.2<z<2.6 line
;; (and not the BOSS 0.3<z<0.7 line!
if (choice_BOSS_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 $
                and Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N) 
;;NOTE TO NPR need to keep putting in this  RAW_N_BOSS_QSOS cut...
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[0,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!80.30<z<0.68!3', charsize=2.2, charthick=6


;;
;; z = 0.87
;;
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
if (choice_BOSS_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N )
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[1,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!80.68<z<1.06!3', charsize=2.2, charthick=6

;;
;;  z = 1.25
;;
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
if (choice_BOSS_points eq 'y') then begin
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N)  
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[2,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.06<z<1.44!3', charsize=2.2, charthick=6

;;
;; z = 1.63
;;
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
plotsym, 0, 1.5, /fill
if (choice_BOSS_points eq 'y') then begin
   ;;www =  where(z_bin_boss gt 1.44 and z_bin_boss lt 1.82) 
   ;;oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N)  
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[3,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.44<z<1.82!3', charsize=2.2, charthick=6



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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6
plotsym, 0, 1.5, /fill
if (choice_BOSS_points eq 'y') then begin
   ;; www =  where(z_bin_boss gt 1.82 and z_bin_boss lt 2.20) 
   ;; oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and $
                Abs_mag_bin_boss lt -24.5 and RAW_N_BOSS_QSOS gt 0, N) 
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[4,*], limit_lines, color=160, linestyle=1, thick=4

xyouts, x_xyouts, y_xyouts, '!81.82<z<2.20!3', charsize=2.2, charthick=6

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   z = 2.40     !!!  F i r s t   BOSS  BIN  !!!
;;
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

plotsym, 0, plot_sym_size_BOSS, /fill
www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[5,*], limit_lines, color=160, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!82.20<z<2.60!3', charsize=2.2, charthick=6


;;
;; z = 2.80 
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

plotsym, 0, 1.5
plotsym, 0, plot_sym_size_BOSS, /fill
www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 3.00 and Abs_mag_bin_boss lt iMag_limit_line[6,0]) 
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   ;; 2.20<z<2.60 line...
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60  and Abs_mag_bin_boss lt -24.5, N)  
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif


;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[6,*], limit_lines, color=160, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!82.60<z<3.00!3', charsize=2.2, charthick=6


;;
;; z  = ~ 3.25 
;;
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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6  ;; red  = 60 

plotsym, 0, plot_sym_size_BOSS, /fill
www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and Abs_mag_bin_boss lt iMag_limit_line[7,0]) 
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   
   ;; 2.20<z<2.60 line...
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60  and Abs_mag_bin_boss lt -24.5, N)  
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[7,*], limit_lines, color=160, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!83.00<z<3.50!3', charsize=2.2, charthick=6



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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6

;; PLOTTING BOSS POINTS, LINES
plotsym, 0, 1.5
plotsym, 0, plot_sym_size_BOSS, /fill
www = where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00 and Abs_mag_bin_boss lt iMag_limit_line[8,0])  

if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   ;; 2.20<z<2.60 line..
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N)  
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[8,*], limit_lines, color=160, linestyle=1, thick=12
xyouts, x_xyouts, y_xyouts, '!83.50<z<4.00!3', charsize=2.2, charthick=6


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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

loadct, 6

plotsym, 0, 1.5
plotsym, 0, plot_sym_size_BOSS, /fill 
www =  where(z_bin_boss gt 4.00 and z_bin_boss lt 4.50 and Abs_mag_bin_boss lt iMag_limit_line[9,0], N) 
if (choice_BOSS_points eq 'y') then begin
   if N gt 0 then oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
   ;; 2.20<z<2.60 line..
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif

;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[9,*], limit_lines, color=160, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!84.00<z<4.50!3', charsize=2.2, charthick=6


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

plotsym, 8, plot_sym_size_R06, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

plotsym, 0, plot_sym_size_BOSS, /fill
www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00 and Abs_mag_bin_boss lt iMag_limit_line[10,0], N)
;www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00, N)
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8,color=160
   ;; 2.20<z<2.60..
   ;; 2.20<z<2.60 line..
   www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], color=160, thick=12
endif
;; Giving a *very rough* feel for where a BOSS mag-limit would be...
oplot, iMag_limit_line[10,*], limit_lines, color=160, linestyle=1, thick=12

xyouts, x_xyouts, y_xyouts, '!84.50<z<5.00!3', charsize=2.2, charthick=6



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  L A B E L S   in the bottom right hand corner
;;
charthick=4.8
charsize=2.0
loadct, 0
plotsym, 8, 1.4, /fill
legend, 'SDSS QSOs' , $
        position=[-31.00, -6.10], box=0, psym=8, color=0, $
        charsize=charsize, charthick=charthick

xyouts, -31.20, -7.00,'(Richards+ 2006)',  color=0, charsize=charsize, charthick =charthick


loadct, 6
plotsym, 0, 1.5, /fill
;legend, 'R06 re-done by NPR' , $
;        position=[-31.00, -6.60], box=0, psym=8, color=60, $
;        charsize=charsize, charthick =charthick

if (choice_BOSS_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   legend, 'BOSS QSOs' , $
;        position=[-31.00, -7.50], box=0, psym=8, color=160, $
           position=[-31.00, -7.25], box=0, psym=8, color=160, $
           charsize=charsize, charthick =charthick
   
   legend, ' ' , $
           position=[-31.00, -8.10], box=0, linestyle=1, color=160, $
           charsize=1.5, charthick =12., thick= 16.0
endif
   
;xyouts, -31.20, -9.0,'i!9A!321.8 limit', $
xyouts, -31.20, -9.0,'i=21.8 limit', $
        color=160, charsize=charsize*1.4, charthick =charthick*1.4

loadct, 0




device, /close
set_plot, 'X'
close, /all

end
