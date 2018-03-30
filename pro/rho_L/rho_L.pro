;+
; NAME:
;   program_name
;
; PURPOSE:
;   Purpose here. 
;
; CALLING SEQUENCE:
;    program_name, [ option= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine ...
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   readspec
;   splog
;
; NOTES:
;
; REVISION HISTORY:
;   11-Jan-2011  v0.0.1     NPR
;-


readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_R06, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor, /silent

;; -26.0 M* from Table 7 of Richards et al. 
rho_L_R06 =  (2.51189)^((-26.0) - Mi_R06 )  * (10^(log_PhiR06))

;;
;; BOSS DR9 
;;
readcol, '../../pro/qlf/My_QLF_iband_boss_20120619.dat', $ ;; used in DR9 QLF paper!!
         z_DR9, Mi_DR9, blah_boss, log_PhiDR9, raw_N_QSOs_boss, sigma_Ph

rho_L_DR9 =  (2.51189)^((-26.0) - Mi_DR9 )  * (10^(log_PhiDR9))

;;
;; For Richards 06
;;
red_bins_r06 = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]

rho_L_R06_zbinned = fltarr(n_elements(red_bins_r06))

w  =  where(z_R06 gt 0.11 and z_R06 lt 0.68) 
rho_L_R06_zbinned[0] = total(rho_L_R06[w])
w  =  where(z_R06 gt 0.68 and z_R06 lt 1.06) 
rho_L_R06_zbinned[1] = total(rho_L_R06[w])
w  =  where(z_R06 gt 1.06 and z_R06 lt 1.44) 
rho_L_R06_zbinned[2] = total(rho_L_R06[w])
w  =  where(z_R06 gt 1.44 and z_R06 lt 1.82) 
rho_L_R06_zbinned[3] = total(rho_L_R06[w])
w  =  where(z_R06 gt 1.82 and z_R06 lt 2.20) 
rho_L_R06_zbinned[4] = total(rho_L_R06[w])
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
rho_L_R06_zbinned[5] = total(rho_L_R06[w])
w  =  where(z_R06 gt 2.60 and z_R06 lt 3.00) 
rho_L_R06_zbinned[6] = total(rho_L_R06[w])
w  =  where(z_R06 gt 3.00 and z_R06 lt 3.50) 
rho_L_R06_zbinned[7] = total(rho_L_R06[w])
w  =  where(z_R06 gt 3.50 and z_R06 lt 4.00) 
rho_L_R06_zbinned[8] = total(rho_L_R06[w])
w  =  where(z_R06 gt 4.00 and z_R06 lt 4.50) 
rho_L_R06_zbinned[9] = total(rho_L_R06[w])
w  =  where(z_R06 gt 4.50 and z_R06 lt 5.00) 
rho_L_R06_zbinned[10] = total(rho_L_R06[w])


;;
;; DR9 
;;
red_bins_dr9 = [2.40, 2.80, 3.25]
rho_L_dr9_zbinned = fltarr(n_elements(red_bins_dr9))

w  =  where(z_dr9 ge 2.20 and z_dr9 lt 2.60) 
rho_L_dr9_zbinned[0] = total(rho_L_dr9[w])
w  =  where(z_dr9 ge 2.60 and z_R06 lt 3.00) 
rho_L_dr9_zbinned[1] = total(rho_L_dr9[w])
w  =  where(z_dr9 ge 3.00 and z_R06 lt 3.50) 
rho_L_dr9_zbinned[2] = total(rho_L_dr9[w])



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

charsize=3.0
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; x-ranges
xmin = 0.0
xmax = 5.2

;; y-ranges
ymin = -7.0
ymax = -5.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='rho_L_temp.eps', $
        xsize=12.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 0.4, /fill
plot, red_bins_r06, $
      alog10(rho_L_R06_zbinned), $
;      psym=8, $
      position=[0.20, 0.20, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
;      /nodata, $
      xtitle='!8z!6, redshift', $
      ytitle='!6!9? !6 log!D10!N(!7q!D!6L!N)', $
;      ytitle='!6!9? !6 log!D10!N(luminosity density)', $
      color=black

plotsym, 0, 2, /fill
oplot, red_bins_dr9, alog10(rho_L_dr9_zbinned), $
       ps=8, color=red


;xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

;legend, pos=[-22.2, -7.5], ' ',   box=0, thick=14, linestyle = 0, charsize=1.2
;xyouts,      -25.0, -7.8,  'PLE', charsize=2.2, charthick=8.


device, /close
set_plot, 'X'     



end
