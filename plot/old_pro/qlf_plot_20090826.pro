;;+
;; NAME:
;;       qlf_plot
;; 
;; PURPOSE:
;;        To plot Quasar Luminosity Functions
;;
;; CALLING SEQUENCE:
;;       .run qlf_plot
;;
;; INPUTS:
;;        Richards06_Table06.dat   - QLF from Richards06, z=2.01 bin only
;;
;; OUTPUTS:
;;       .ps file
;;
;; COMMENTS:
;;       /usr/common/rsi/lib/general/LibAstro/ 
;;         
;;-

readcol, 'Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor

;readcol, 'My_QLF_iband_20090717.dat', z_bin, Abs_mag_bin, blah, log_Num_den
readcol, 'My_QLF_iband_temp.dat',  z_bin, Abs_mag_bin, blah, log_Num_den
      
x_min = -18.001
x_max = -30.50
y_min = -8.90
y_max = -4.60

char_size  = 3.6
char_thick = 5.0   

set_plot, 'ps'
!p.multi=[0,3,4]
device, filename='QLF_R06_z2_temp.ps', $
        xsize=14, ysize=12, $
        /inches, /color


;; z = 0.49
w  =  where(z_R06 gt 0.11 and z_R06 lt 0.68) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.70, 0.35, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 0.11 and z_bin lt 0.68) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -26.50, -6.40, '!8z=0.49!3', charsize=2.2, charthick=6


;; z = 0.87
w  =  where(z_R06 gt 0.68 and z_R06 lt 1.06) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.70, 0.56, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 0.68 and z_bin lt 1.06) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -26.50, -6.40, '!8z=0.87!3', charsize=2.2, charthick=6


;; z = 1.25
w  =  where(z_R06 gt 1.06 and z_R06 lt 1.44) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.70, 0.77, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $ 
       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 1.06 and z_bin lt 1.44) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z=1.25!3', charsize=2.2, charthick=6


;; z = 1.63
w  =  where(z_R06 gt 1.44 and z_R06 lt 1.82) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.77, 0.70, 0.98, 0.98], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 1.44 and z_bin lt 1.82) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z=1.63!3', charsize=2.2, charthick=6


;;;;;;;;;       2nd ROW         ;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; 
;;;;;;;;;      z = 2.01         ;;;;;;;;;;;;;;;;;;;;;;;;;
w  =  where(z_R06 gt 1.82 and z_R06 lt 2.20) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.42, 0.35, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 1.82 and z_bin lt 2.2) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z=2.01', charsize=2.2, charthick=6


;; z = 2.40 
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.42, 0.56, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 2.20 and z_bin lt 2.60) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z=2.40!3', charsize=2.2, charthick=6


;; z = 2.80 
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.42, 0.77, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtickformat='(a1)', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 2.60 and z_bin lt 3.00) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z!9A!82.80!3', charsize=2.2, charthick=6



;; z  = ~ 3.25 
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.77, 0.42, 0.98, 0.70], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 3.00 and z_bin lt 3.50) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -23.50, -8.40, '!8z!9A!83.25!3', charsize=2.2, charthick=6



;;;;;;;;;       3rd ROW         ;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; 
;;           z = ~ 3.75
w  =  where(z_R06 gt 3.50 and z_R06 lt 4.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.14, 0.14, 0.35, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtitle='!3M!Ii!N[z=2]', $
       ytitle='!7U!3(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 3.50 and z_bin lt 3.90) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

ww =  where(z_bin gt 3.90 and z_bin lt 4.20) 
loadct, 6
plotsym, 0, 1.5
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -24.00, -6.40, '!8z!9A!83.75!3', charsize=2.2, charthick=6


;;           z = ~ 4.25
w  =  where(z_R06 gt 4.00 and z_R06 lt 4.50) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.35, 0.14, 0.56, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'

plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 4.00 and z_bin lt 4.50) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -24.00, -6.40, '!8z!9A!84.25!3', charsize=2.2, charthick=6


;;           z = ~ 4.75
w  =  where(z_R06 gt 4.50 and z_R06 lt 5.00) 
plot,  Mi_z2[w], log_PhiR06[w], $
       position=[0.56, 0.14, 0.77, 0.42], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=4, $
       ythick=4, $
       charsize=char_size, $
       charthick=char_thick, $
       xtitle='!3M!Ii!N[z=2]', $
       ytickformat='(a1)'


plotsym, 8, 1.4, /fill
oplot,  Mi_z2[w], log_PhiR06[w], psym=8

ww =  where(z_bin gt 4.50 and z_bin lt 5.00) 
loadct, 6
plotsym, 0, 1.5, /fill
oplot, Abs_Mag_bin[ww], log_Num_den[ww], psym=8, color=60

xyouts, -24.00, -6.40, '!8z!9A!84.75!3', charsize=2.2, charthick=6


;; Labels in the bottom right hand corner
loadct, 0
plotsym, 8, 1.4, /fill
legend, 'Richards et al. (2006)' , $
        position=[-31.00, -6.50], box=0, psym=8, color=0, $
        charsize=1.6, charthick =4.2

loadct, 6
plotsym, 0, 1.5, /fill
legend, 'Stripe 82 (this paper)' , $
        position=[-31.00, -7.00], box=0, psym=8, color=60, $
        charsize=1.6, charthick =4.2
loadct, 0









device, /close
set_plot, 'X'
close, /all

end
