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


readcol, 'QLF_from_McGreer13.dat', $
         mag_z5, phi_z5, log_Phi_z5

;;
;; Colour Table
;; http://ham.space.umn.edu/johnd/ct/ct-names.html
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

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03

;; positions...
xpos_min = 0.20
xpos_max = 0.98
ypos_min = 0.20
ypos_max = 0.98

;; x-ranges
xmin = -23.00
xmax = -29.00

;; y-ranges
ymin = 
ymax = 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 0.4, /fill
plot, mag_z5, phi_z5, $
      psym=8, $
      position=[xpos_min, ypos_min, xpos_max, ypos_max], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
;      /nodata, $
      xtitle='!8z, !6redshift', $
      ytitle='!6N(!8z!6), !8z!6=0.05 bins', $
      color=black

xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

legend, pos=[-22.2, -7.5], ' ',   box=0, thick=14, linestyle = 0, charsize=1.2
xyouts,      -25.0, -7.8,  'PLE', charsize=2.2, charthick=8.


device, /close
set_plot, 'X'     

end
