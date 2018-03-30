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

print
print
print, '------------------------------------------------------'
print, '------------------------------------------------------'
print, '!!!! You have to run the "red" routine before this !!!!'
print, '          '
print, '   red, omega0 = 0.228, omegalambda = 0.772, h100=0.702  ' ;; (Komatsu et al. 2011) 
;print, '  red, omega0 = 0.300, omegalambda = 0.700, h100=0.70   '
print, '  '
print, '------------------------------------------------------'
print, '------------------------------------------------------'
print
;red, omega0=0.30, omegalambda=0.70, h100=0.700


readcol, '../data/Croom09b_2SLAQ_QLF.dat', $ 
         Mg_2slaq, z_2slaq, NQ_2slaq, log_phi_2SLAQ, delta_log_phi_upper, delta_log_phi_lower


readcol, '../pro/qlf/My_QLF_gband_boss_temp.dat',  z_bin_boss, Abs_mag_bin_boss, blah_boss, log_Num_den_boss

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


;;  Below I list a piece of Python code which gives the Croton09 fit to the QLF as a function of bJ magnitude and redshift.  This is a modified version of the Croom++04 fit.  I turned the densities into number per comoving (Mpc/h)^3 per magnitude from their per Mpc^3 per magnitude by dividing by hub**3 in the last line ... this of course assumes you've defined "hub" to be something useful like 0.72!  I was wondering whether we could overplot that line on your figure?  My by-eye comparisons show that Darren's modification does a pretty good job!

;;   lumfn_C09(mag, zz):
;;   The QSO luminosity function at redshift zz: dn/dM with magnitudes
;;   in the bJ system.
;;   Parameters for the LF fit from Croom++04, as modified by Croton09.
;;   This agrees with the C04 for z<3.
;;   """
;; 
;;

;; 
;; Copying straight from Darren's Code...
;;
numz = 100
zz = findgen(numz)/10    ;;  Redshift bins..

mvir = 10.0^((findgen(60)/10)+10.0)  ;; Virial Mass, 1e10 -> 1e15, in "alog10" steps ;-) 
numm = n_elements(mvir)

;Abs_Mag_model = fltarr(numm)
Abs_Mag_model = (findgen(numm)/4.)-32.

mass  = fltarr(numz, numm)
;mag   = fltarr(numz, numm)
Phi   = fltarr(numz, numm)
Phi_h = fltarr(numz, numm)

px = fltarr(numz)
py = fltarr(numz)

G   =  4.3e-9   ;; big G in units of (km s-1)^2 Mpc Msol^-1
H0  = 100.0     ;; H-nought
hub = 0.70      ;; Something useful!   

;print, alog10(1.67e-6)
;     -5.77728
;; MW's modified numbers
PhiStar   =   1.67e-6  ;;
Mstar0_bJ = -21.61     ;; b_J magnitude value
Mstar0_g  = -22.24     ;; converted to g-band (Croom09b

beta  =  -1.09
;; z< 3
;   k1,k2,a = 1.39,-0.29,-3.31
;;z>3
;   k1,k2,a = 1.22,-0.23,-3.31+0.5*(zz-3.0)

;; Croom et al. 2009 QLF
;; LEDE Model, Table 4
;PhiStar = 1.62e-6
;alpL = -3.48 & betL = -1.4
;Mstar0 = -22.24
;k1 = 1.23 & k2 = -0.206

for ii=0L, numz-1 do begin
   Hz = H0 * sqrt(0.3*(1.0+zz[ii])^3.0 + 0.7)
   ;Hz = H0 * sqrt(0.25*(1.0+zz[i])^3.0 + 0.75)
   
   if zz[ii] lt 3.0 then begin
      ;; Actual numbers from Croom09b,Table 4
      ;k1    =  1.23
      ;k2    = -0.206
      ;alpha = -3.48
      ;; Modified no.s from Croton/Martin
      k1    =  1.39
      k2    = -0.29
      alpha = -3.31
   endif else begin
      ;; Croton modified values.... 
      ;; (no values for Croom09b above z=2.6...)
      k1    =  1.22
      k2    = -0.23
      alpha = -3.31 + (0.5*(zz[ii]-3.0))
   endelse
   
   ;; For both bJ and g-band...
   ;Mstar_bJ = Mstar0_bJ - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   
   for jj=0L, numm-1 do begin
      mag = -32.0 + (jj*0.25)
;      Abs_Mag_model[jj] = mag
      dmag   = mag - Mstar_g 
;      Phi[ii, *] = PhiStar / ( (10.0^(0.4*(alpL+1)*(magQ-Mstar))) + (10.0^(0.4*(betL+1)*(magQ-Mstar))) )
      Phi[ii, jj] = PhiStar / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

      print, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
   endfor
   
;   phi   = phi0/(10.0**(0.4*(a+1)*dm)+10.0**(0.4*(b+1)*dm));
;   phi  /= hub**3  ;; Convert from Mpc to Mpc/h volumes.
;   print, ii, zz[ii], 
endfor
Phi_h = Phi / (hub^3)
print

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

choice_2SLAQ_line = 'y'
read, choice_2SLAQ_line, PROMPT=' - Plot 2SLAQ line?? y/n  '

choice_BOSS_points = 'y'
read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '

choice_BOSS_line = 'y'
;read, choice_BOSS_line, PROMPT=' - Plot BOSS line ?? y/n  '

choice_Croton09_model = 'y'
read, choice_Croton09_model, PROMPT=' - Overplot the Croton (2009) model?  y/n  '
Croton_color = 220
Croton_thick = 8.0
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
       ytitle='!7U!3(M!Ig!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

loadct, 6
plotsym, 0, 1.4, /fill
oplot, Mg_2slaq[w],  log_phi_2SLAQ[w],  psym=8, color=60.

if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12
;oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(5,*)), $
   color=Croton_color, thick=Croton_thick





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

if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(8,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(13,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12


if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(16,*)), $
   color=Croton_color, thick=Croton_thick

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
       ytitle='!7U!3(M!Ig!N[z=2]) [Mpc!E-3!N mag!E-1!N]'

loadct, 6
plotsym, 0, 1.5, /fill
oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(20,*)), $
      color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12
;oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

;www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt -24.5, N) 
www =  where(z_bin_boss gt 2.20 and z_bin_boss lt 2.60 and Abs_mag_bin_boss lt  gMag_limit_line[5,0], N) 
if (choice_BOSS_points eq 'y') then begin
   plotsym, 0, 1.5, /fill
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(24,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

www =  where(z_bin_boss gt 2.60 and z_bin_boss lt 3.00 and Abs_mag_bin_boss lt gMag_limit_line[6,0] ) 
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(28,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

www =  where(z_bin_boss gt 3.00 and z_bin_boss lt 3.50 and Abs_mag_bin_boss lt gMag_limit_line[7,0] )  
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(33,*)), $
   color=Croton_color, thick=Croton_thick

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
       ytitle='!7U!3(M!Ig!N[z=2]) [Mpc!E-3!N mag!E-1!N]'
loadct, 6
plotsym, 0, 1.5, /fill
;oplot,  Mg_2slaq[w], log_phi_2SLAQ[w], psym=8, color=60
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12

www =  where(z_bin_boss gt 3.50 and z_bin_boss lt 4.00 and Abs_mag_bin_boss lt gMag_limit_line[8,0]) 
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(38,*)), $
      color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12
;oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 4.00 and z_bin_boss lt 4.50 and Abs_mag_bin_boss lt gMag_limit_line[9,0])  
if (choice_BOSS_points eq 'y') then begin
   oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(43,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60, thick=12
;oplot, Mg_2slaq[ww], log_phi_2SLAQ[ww], color=60.

www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00 and Abs_mag_bin_boss lt gMag_limit_line[10,0], N) 
if (choice_BOSS_points eq 'y') then begin
;www =  where(z_bin_boss gt 4.50 and z_bin_boss lt 5.00) 
   if N gt 0 then oplot, Abs_Mag_bin_boss[www], log_Num_den_boss[www], psym=8, color=160
endif 

if choice_Croton09_model eq 'y' then oplot, Abs_Mag_model, alog10(phi(48,*)), $
   color=Croton_color, thick=Croton_thick

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
if choice_2SLAQ_line eq 'y' then begin
   xyouts, -34.00, -6.80,' 1.82<z<2.20 line ', charsize=1.2, charthick =6.2
   legend, ' ' , $
           position=[-30.80, -6.60], box=0, linestyle=0, color=60, $
           charsize=1.1, charthick =6.2, thick=4.0
endif

if (choice_BOSS_points eq 'y') then begin
   loadct, 6
   plotsym, 0, 1.5, /fill
   legend, 'BOSS QSOs' , $
           position=[-30.80, -7.20], box=0, psym=8, color=160, $
           charsize=charsize, charthick=charthick
   
;   legend, '(Vmax) ' , $
 ;          position=[-31.10, -7.60], box=0, color=160, $
;           charsize=charsize, charthick=charthick
endif  

if CHOICE_CROTON09_MODEL eq 'y' then begin
   xyouts, -34.00, -7.90,' Croton (2009) model ', charsize=1.2, charthick =6.2
   legend, ' ' , $
           position=[-30.80, -7.70], box=0, linestyle=0, color=Croton_color, $
           charsize=1.1, charthick =Croton_thick, thick=4.0*2.
endif

xyouts, -34.00, -8.50,' g=22.0 limit',   color=160, charsize=2.0, charthick =6.2
legend, ' ' , $
        position=[-30.80, -8.30], box=0, linestyle=1, color=160, $
        charsize=1.1, charthick =6.2, thick=6.0




loadct, 0

!p.multi=0


device, /close
set_plot, 'X'
close, /all

end
