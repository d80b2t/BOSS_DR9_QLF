;+
; NAME:
;   number_counts_plot.pro
;
; PURPOSE:
;   Plots i-band QSO number counts
;
; INPUTS:
;
; OUTPUTS:
;
; COMMENTS:
;   This routine ...
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   15-May-2012  v0.9.1     NPR
;-

print
print

choice_plot_SDSS = 'n'
;read,  choice_plot_SDSS, prompt=' Plot SDSS DR3 R06 number counts?? y/n '
print
print

;;
;;
;; Richards et al. (2006)
;;
;; Quasar Number Counts (0.3 < z < 2.2)
;;
readcol, '../../pro/number_counts/Richards06_Table2_0pnt3z2pnt2.dat', $
         SDSS_mag, $
         SDSS_N_g, SDSS_N_gerr, SDSS_N_lt_g, SDSS_N_lt_gerr, N_Q_gmag, $
         SDSS_N_i, SDSS_N_ierr, SDSS_N_lt_i, SDSS_N_lt_ierr, N_Q_imag 
;; print, mean(N_Q_imag/ SDSS_N_i)
;;      1554.20

;; Don't like 0.00 errorbars, even if this is what's quoted in R06!!
w = where(SDSS_N_ierr eq 0.000, N)
;; as good a estimate as anything!
SDSS_N_ierr[w] = (sqrt(N_Q_imag[w])/1550.)

w = where(SDSS_N_lt_ierr eq 0.0, N)
SDSS_N_lt_ierr[w] = (sqrt(N_Q_imag[w])/1550.)

;;
;; Richards et al. (2006)
;; Quasar Number Counts (3 < z < 5)
;;
readcol, '../../pro/number_counts/Richards06_Table3_3z5.dat', $
         SDSS_mag_hiz, $
         SDSS_N_i_hiz, SDSS_N_ierr_hiz, SDSS_N_lt_i_hiz, SDSS_N_lt_ierr_hiz, $
         N_Q_imag_hiz, N_Q_corr_imag_hiz
;; print, mean(N_Q_corr_imag_hiz/  SDSS_N_i_hiz)
;;      1597.36
w = where(SDSS_N_ierr_hiz eq 0.000, N)
;; as good a estimate as anything!
SDSS_N_ierr_hiz[w] = (sqrt(N_Q_corr_imag_hiz[w])/1550.)

w=where(SDSS_N_lt_ierr_hiz eq 0.0, N)
SDSS_N_lt_ierr_hiz[w] = (sqrt(N_Q_corr_imag_hiz[w])/1550.)

;;
;; Mi(z = 2) = Mg(z = 2) −0.25	(6) 
;;           = M1450     − 0.29	(7) 
;;           = MbJ       −0.71
;;
;; If you wanted to cf.. with the e.g. 2QZ or 2SLAQ...



;;
;;  D A T A    F O R   B O S S   D R 9  
;;
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_DR9_XDQSO_iband_numcounts_20120518.dat', $
         BOSS_imag, BOSS_N_DR9, BOSS_DR9_nm, BOSS_N_DR9_sum, BOSS_DR9_nm_sum
   
area_DR9=2236.
BOSS_DR9_nm_err     = sqrt(BOSS_N_DR9)/area_DR9
BOSS_DR9_nm_sum_err = sqrt(BOSS_N_DR9_sum)/area_DR9

;; 
;;  C O M P L E T E N E S S 
;;
readcol, '../../completeness/fiducial_linetweak_grid_i_sel.txt', $
        i_start, i_end, z_start, z_end, sel

N_sel_bins = (max(i_start)-min(i_start))/0.1
BOSS_DR9_nm_corr         = fltarr(N_elements(BOSS_DR9_nm))
BOSS_DR9_nm_sum_corr     = fltarr(N_elements(BOSS_DR9_nm))
BOSS_DR9_nm_sum_err_corr = fltarr(N_elements(BOSS_DR9_nm))

openw, 10, 'nm_cuml_BOSS_R06_imag_temp.dat'
printf, 10, '# BOSS_imag,  BOSS_N_DR9, BOSS_DR9_nm_sum, BOSS_DR9_nm_sum_corr, BOSS_DR9_nm_sum_err_corr'

for ii=0ll, N_elements(BOSS_DR9_nm)-2 do begin
;   print, ii, BOSS_imag[ii], BOSS_imag[ii] , BOSS_DR9_nm[ii], N, corr_factor, BOSS_DR9_nm_corr[ii]
   BOSS_DR9_nm_corr[ii] = BOSS_DR9_nm[ii]
;   w = where(i_start le BOSS_imag[ii] and i_start le BOSS_imag[ii+1] and , N)
   w = where(i_start le BOSS_imag[ii] and i_end ge BOSS_imag[ii] and z_start ge 2.20 and z_end le 3.50, N)
   if (N gt 0) then begin
      if (max(sel[w]) gt 0) then begin
         corr_factor = 1./mean(sel[w])
         BOSS_DR9_nm_corr[ii]           = BOSS_DR9_nm[ii]*corr_factor
         BOSS_DR9_nm_sum_corr[ii]       = total(BOSS_DR9_nm_corr)
;         BOSS_DR9_nm_sum_err_corr[ii] =  sqrt(BOSS_DR9_sum[ii]*corr_factor)/area_DR9
         BOSS_DR9_nm_sum_err_corr[ii] =  BOSS_DR9_nm_sum_err[ii]*corr_factor

         print, ii, BOSS_imag[ii], BOSS_imag[ii+1], BOSS_N_DR9[ii], BOSS_DR9_nm[ii], corr_factor, BOSS_DR9_nm_corr[ii], BOSS_DR9_nm_sum[ii], BOSS_DR9_nm_sum_corr[ii]
         printf, 10, 0.5*( BOSS_imag[ii]+ BOSS_imag[ii+1]),  BOSS_N_DR9[ii],  BOSS_DR9_nm_sum[ii], BOSS_DR9_nm_sum_corr[ii],  BOSS_DR9_nm_sum_err_corr[ii]
      endif
   endif
endfor
close, 10


;;
;;  D A T A    F O R    S T R I P E  8 2
;;
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_DR9_Stripe82_iband_numcounts_20120518.dat', $
         BOSS_imag, BOSS_N_S82, BOSS_S82_nm, BOSS_N_S82_sum, BOSS_N_S82_nm_sum
area_S82=220.
BOSS_S82_nm_err     = sqrt(BOSS_N_S82)/area_S82
BOSS_S82_nm_sum_err = sqrt(BOSS_N_S82_sum)/area_S82

;;
;;  D A T A    F O R    ``b o s s  2 1 ''
;;
;readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_boss21_iband_numcounts.dat', $
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_DR9_boss21_iband_numcounts_20120518.dat', $
         BOSS21_imag, BOSS_N_boss21, BOSS_boss21_nm,  BOSS_N_boss21_sum, BOSS_boss21_nm_sum
area_boss21 =14.5
BOSS_boss21_nm_err  = sqrt(BOSS_N_boss21)/area_boss21

BOSS_boss21_nm_sum_err = sqrt(BOSS_N_boss21_sum)/area_boss21

openw, 11, 'nm_cuml_boss21_imag_temp.dat'
printf, 11, '# boss21_imag, boss21_N, BOSS_boss21_nm_sum, BOSS_boss21_nm_sum_err'
for ii=0ll, n_elements( BOSS21_imag) -1 do begin
   printf, 11, BOSS21_imag[ii], BOSS_N_boss21[ii], BOSS_boss21_nm_sum[ii], BOSS_boss21_nm_sum_err[ii]
endfor
close, 11


;;
;;  M A G   B I N   S I Z E ,    will vaery between  
;;
mag_bin = 0.10


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

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.05
YTICKLEN  = 0.05



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   i-band  N(M)'s
;;
set_plot, 'ps'
device, filename='nm_BOSS_R06_imag_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

;; x-ranges
xmin = 14.21   ;; 15.01 is nice...
xmax = 24.99   ;; 23.49 is okay

;; y-range
ymin = 0.00002   ;; 0.0002 is okay...
ymax = 99.9     ;; 10.0 is okay w/o boss21...

plotsym, 0,1, /fill
plot, SDSS_mag,  SDSS_N_i, $
      position=[0.22, 0.14, 0.96, 0.96], $
      ps=2, $
      /ylog, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], ystyle=1, xstyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick, thick=thick,$ 
      YTICKNAME = ['10!E-4!N','10!E-3!N', '0.01', '0.1', '1', '10', '100'], $
      yminor=10., $
      color=black, $
      xtitle='i-band PSF magnitude', $
      ytitle='N(<i) (deg!E-2!N)'
;      ytitle='N(i) (deg!E-2!N [0.25 mag]!E-1!N)'

;;
;; errorbar styles
;; http://idlastro.gsfc.nasa.gov/ftp/pro/plot/oploterror.pro
;;
sdss_color = black
boss_color = red
boss21_color = blue
errthick = 1.0

SDSS_mag_xerr = SDSS_mag - SDSS_mag
oploterror, SDSS_mag, SDSS_N_i, SDSS_mag_xerr, SDSS_N_ierr, $ 
            /hibar, errthick=errthick, psym=2, color=sdss_color, ERRcolor=sdss_color
oploterror, SDSS_mag, SDSS_N_i, SDSS_mag_xerr, SDSS_N_ierr, $ 
            /lobar, errthick=errthick, psym=2, color=sdss_color, ERRcolor=sdss_color

;;
;; SDSS  3 < z < 5
;;
oplot, SDSS_mag_hiz,  SDSS_N_i_hiz, ps=8, color=black
oploterror, SDSS_mag_hiz, SDSS_N_i_hiz, SDSS_mag_xerr, SDSS_N_ierr_hiz, $
            /hibar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color
oploterror, SDSS_mag_hiz, SDSS_N_i_hiz, SDSS_mag_xerr, SDSS_N_ierr_hiz, $
            /lobar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color

s82_color = blue
errthick  = 6.
;;
;;  B O S S    S T R I P E   8 2
;;
w = where(BOSS_imag lt 22.85, N)
;oplot,  BOSS_imag[w], BOSS_S82_nm[w]*(1./mag_bin), ps=1, color=red
BOSS_mag_err = BOSS_imag - BOSS_imag
oploterror, BOSS_imag[w],    BOSS_S82_nm[w]     * (1./mag_bin), $
            BOSS_mag_err[w], BOSS_S82_nm_err[w] * (1./mag_bin), $
            /hibar, errthick=errthick, psym=1, color=s82_color, ERRcolor=s82_color
oploterror, BOSS_imag[w],    BOSS_S82_nm[w]      *(1./mag_bin), $
            BOSS_mag_err[w], BOSS_S82_nm_err[w] * (1./mag_bin), $
            /lobar, errthick=errthick, psym=1, color=s82_color, ERRcolor=s82_color

;;
;;  B O S S    D R 9 
;;
plotsym, 0, 1.4
oplot,  BOSS_imag, BOSS_DR9_nm*(1./mag_bin),      ps=8, color=red

;; with sel. fn. correction applied...
plotsym, 0, 1.4, /fill
w = where(BOSS_imag lt max(i_end), N)
oplot,  BOSS_imag[w], BOSS_DR9_nm_corr[w]*(1./mag_bin), ps=8, color=red

;; not much point really, since they are ~same size as symbols!
BOSS_mag_err = BOSS_imag - BOSS_imag
w = where(BOSS_N_DR9 gt 0, N)
;oploterror, BOSS_imag[w],    BOSS_DR9_nm[w]    * (1./mag_bin), $
;            BOSS_mag_err[w], BOSS_DR9_nm_err[w]* (1./mag_bin), $
;            /hibar, errthick=errthick, psym=2, color=boss_color, ERRcolor=boss_color
;oploterror, BOSS_imag[w],    BOSS_DR9_nm[w]    * (1./mag_bin), $
;            BOSS_mag_err[w], BOSS_DR9_nm_err[w]* (1./mag_bin), $
;            /lobar, errthick=errthick, psym=2, color=boss_color, ERRcolor=boss_color



;;
;;   `` B O S S   2 1''
;;
oplot,  BOSS21_imag, BOSS_boss21_nm*(1./mag_bin),  ps=8,  color=blue
oploterror, BOSS_imag,    BOSS_boss21_nm*(1./mag_bin), $
            BOSS_mag_err, BOSS_boss21_nm_err, $
            /hibar, errthick=errthick, psym=8, color=boss21_color, ERRcolor=boss21_color
oploterror, BOSS_imag,    BOSS_boss21_nm*(1./mag_bin), $
            BOSS_mag_err, BOSS_boss21_nm_err, $
            /lobar, errthick=errthick, psym=8, color=boss21_color, ERRcolor=boss21_color


;;
;; POWER-LAWs...
;;
amplitude = 18.5
pl_slope  = 0.94
slope_x = findgen(5)/1.+16.
slope_y = (10^(slope_x-amplitude))^pl_slope
oplot, slope_x, slope_y, thick=6, linestyle=2

;;
;; Smoothed power-law model (a al Myers et al. (2003) and every other
;; QSO paper before and after that... ;-) 
;;

N_zero   = 10.0
alpha_d  = 0.98
beta_d   = 0.15
m_nought = 19.10
m_running = (findgen(100)/20.)+16.5 

dmag  = m_running - m_nought

double_power_law = (10^((-1.)*alpha_d*dmag) ) +  (10^((-1.)*beta_d*dmag)  )
dN_dm = N_zero / double_power_law

oplot, m_running, dN_dm, thick=6, linestyle=2, color=purple


;; 
;;  L a b e l s....
;;

;; PL label first
xyouts, 15.1, 5., 'N(i)!9?!3i!U0.94!N', charthick=4.8, charsize=2.6


up = 0.
   top     = 10^(-2.5+up)
middle     = 10^(-2.9+up)
bottom     = 10^(-3.3+up)
below_btm  = 10^(-3.7+up)
below_btmm = 10^(-4.1+up)

up = 0.20
   top_leg     = 10^(-2.5+up)
middle_leg     = 10^(-2.9+up)
bottom_leg     = 10^(-3.3+up)
below_btm_leg  = 10^(-3.7+up)
below_btmm_leg = 10^(-4.1+up)

plotsym, 0, 1.4, /fill
;legend, [' '], linestyle=lines, color=32, pos=[0.1,9.4], box=0, thick=8, charthick=4.2 ;, charsize=1.8
legend, [' '], pos=[18.7, top_leg],         box=0, psym=2, color=black
legend, [' '], pos=[18.7, middle_leg],      box=0, psym=8, color=black

plotsym, 0, 1.4
legend, [' '], pos=[18.7, bottom_leg],      box=0, psym=8, color=red
plotsym, 0, 1.4, /fill
legend, [' '], pos=[18.7, below_btm_leg],  box=0, psym=8, color=red
legend, [' '], pos=[18.7, below_btmm_leg], box=0, psym=1, color=S82_color
;legend, [' '], pos=[19.0, below_btmm_leg], box=0, psym=8, color=blue

charsize  = charsize/1.2
charthick = charthick*1.2
xyouts, 19.3, top,       'SDSS 0.3<z<2.2', charsize=charsize, charthick=charthick, color=black
xyouts, 19.3, middle,    'SDSS   3<z<5  ', charsize=charsize, charthick=charthick, color=black
xyouts, 19.3, bottom,    'BOSS  2.2<z<3.5', charsize=charsize, charthick=charthick, color=red
xyouts, 19.3, below_btm, ' "  w/ correction', charsize=charsize, charthick=charthick, color=red
xyouts, 19.3, below_btmm, 'Stp82 2.2<z<3.5', charsize=charsize, charthick=charthick, color=S82_color
;xyouts, 19.6, below_btmm, 'BOSS 1.0<z<2.2', charsize=charsize, charthick=charthick, color=blue
;xyouts, 18., below_btm, 'boss21 0.3<z<2.2', charsize=charsize, charthick=charthick, color=light_blue+8.
charsize  = charsize/1.2
charthick=charthick/1.2

device, /close
set_plot, 'X'     






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  
;;    i - b a n d     N ( < M ) 's
;;
;;

set_plot, 'ps'
if choice_plot_SDSS eq 'y' then begin
   device, filename='nm_cuml_BOSS_R06_imag_temp.eps', $
           xsize=8.0, ysize=8.0, $
           xoffset=0.2, yoffset=0.2, $
           /inches, /color, /encapsulated
endif

if choice_plot_SDSS eq 'n' then begin
   device, filename='nm_cuml_BOSS_imag_temp.eps', $
           xsize=8.0, ysize=8.0, $
           xoffset=0.2, yoffset=0.2, $
           /inches, /color, /encapsulated
endif

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.04
YTICKLEN  = 0.04

;; x-ranges
if choice_plot_SDSS eq 'y' then xmin = 14.21   ;; 15.01 is nice...
if choice_plot_SDSS eq 'n' then xmin = 16.21 
 xmax = 24.49        ;; 24.99   ;; 23.49 is okay

;; y-range
if choice_plot_SDSS eq 'y' then ymin = 0.00009   ;; 0.0002 is okay...
if choice_plot_SDSS eq 'n' then ymin = 0.0020   ;; 0.0002 is okay...
ymax = 150.0          ;; 599.99     ;; 10.0 is okay w/o boss21...

;readcol, '../../pro/number_counts/Richards06_Table2_0pnt3z2pnt2.dat', $
 ;        SDSS_mag, $
  ;       SDSS_N_g, SDSS_N_gerr, SDSS_N_lt_g, SDSS_N_lt_gerr, N_Q_gmag, $
   ;      SDSS_N_i, SDSS_N_ierr, SDSS_N_lt_i, SDSS_N_lt_ierr, N_Q_imag 

plot_sym_size = 1.6

if choice_plot_SDSS eq 'y' then begin
plotsym, 0,plot_sym_size, /fill
plot, SDSS_mag,  SDSS_N_lt_i, $
      position=[0.22, 0.14, 0.96, 0.98], $
      ps=2, $
      /ylog, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], ystyle=1, xstyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick, thick=thick,$ 
;      YTICKNAME = ['10!E-4!N','0.001', '0.01', '0.1', '1', '10', '100'], $
      YTICKNAME = ['10!E-4!N','10!E-3!N', '0.01', '0.1', '1', '10', '100'], $
      yminor=10., $
      /nodata, $
      color=black, $
      xtitle=' i ', $
      ytitle='N(<i) (deg!E-2!N)'
endif

if choice_plot_SDSS eq 'n' then begin
plotsym, 0,plot_sym_size, /fill
plot, SDSS_mag,  SDSS_N_lt_i, $
      position=[0.18, 0.14, 0.96, 0.98], $
      ps=2, $
      /ylog, $
      xrange=[xmin, xmax], yrange=[ymin, ymax], ystyle=1, xstyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick, thick=thick,$ 
;      YTICKNAME = ['10!E-4!N','0.001', '0.01', '0.1', '1', '10', '100'], $
      YTICKNAME = ['0.01', '0.1', '1', '10', '100'], $
      yminor=10., $
      /nodata, $
      color=black, $
;      ytitle='N(<i) (deg!E-2!N)', $
      xtitle='i-band PSF magnitude'
endif


;;
;; errorbar styles
;;
sdss_color = black
boss_color = red
boss21_color = blue
errthick = 4.0


if choice_plot_SDSS eq 'y' then begin
;;CALLING SEQUENCE:
;      oploterror, [ x,]  y, [xerr], yerr,   
;            [ /NOHAT, HATLENGTH= , ERRTHICK =, ERRSTYLE=, ERRCOLOR =, 
;              /LOBAR, /HIBAR, NSKIP = , NSUM = , /ADDCMD, ... OPLOT keywords ]
; http://idlastro.gsfc.nasa.gov/ftp/pro/plot/oploterror.pro

;SDSS_mag_xerr = SDSS_mag - SDSS_mag
plotsym, 3, plot_sym_size, /fill
oploterror, SDSS_mag, SDSS_N_lt_i, SDSS_mag_xerr, SDSS_N_lt_ierr, $
            /hibar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color
oploterror, SDSS_mag, SDSS_N_lt_i, SDSS_mag_xerr, SDSS_N_lt_ierr, $
            /lobar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color

;;
;; SDSS  3 < z < 5
;;
plotsym, 0,plot_sym_size, /fill
oplot,      SDSS_mag_hiz, SDSS_N_lt_i_hiz, ps=8, color=black
oploterror, SDSS_mag_hiz, SDSS_N_lt_i_hiz, SDSS_mag_xerr, SDSS_N_lt_ierr_hiz, $
            /hibar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color
oploterror, SDSS_mag_hiz, SDSS_N_lt_i_hiz, SDSS_mag_xerr, SDSS_N_lt_ierr_hiz, $
            /lobar, errthick=errthick, psym=8, color=sdss_color, ERRcolor=sdss_color

endif


;;
;;  B O S S    D R 9 
;;
plotsym, 0, plot_sym_size
w = where(BOSS_imag le 22.0, N)
oplot,  BOSS_imag[w], BOSS_DR9_nm_sum[w], ps=8, color=red
plotsym, 0, plot_sym_size, /fill
oplot,  BOSS_imag[w], BOSS_DR9_nm_sum_corr[w], ps=8, color=red


BOSS_mag_err = BOSS_imag - BOSS_imag
w = where(BOSS_N_DR9 gt 0, N)
;oploterror, BOSS_imag[w],    BOSS_DR9_nm_sum_corr[w], $
 ;           BOSS_mag_err[w], BOSS_N_DR9_nm_sum_err_corr[w], $
  ;          /hibar, errthick=errthick, psym=2, color=boss_color, ERRcolor=boss_color
;oploterror, BOSS_imag[w],    BOSS_N_DR9_nm_sum_err, $
 ;           BOSS_mag_err[w], BOSS_N_DR9_nm_sum_err_corr[w], $
  ;          /lobar, errthick=errthick, psym=2, color=boss_color, ERRcolor=boss_color


;;
;;  B O S S    S T R I P E   8 2
;;
w = where(BOSS_imag le 22.0, N)
oplot,      BOSS_imag[w], BOSS_N_S82_nm_sum[w], ps=1, color=red, symsize=2.0, thick=6.0

;w = where(BOSS_imag le 19.1, N)
;plotsym, 0, plot_sym_size
BOSS_imag_err = BOSS_imag[w] - BOSS_imag[w]
oploterror, BOSS_imag[w], BOSS_N_S82_nm_sum[w], $
            BOSS_mag_err, BOSS_S82_nm_sum_err[w], $
            /hibar, errthick=errthick, psym=1, color=boss_color, ERRcolor=boss_color
oploterror, BOSS_imag[w], BOSS_N_S82_nm_sum[w], $
            BOSS_mag_err, BOSS_S82_nm_sum_err[w], $
            /lobar, errthick=errthick, psym=1, color=boss_color, ERRcolor=boss_color


;;
;;   `` B O S S   2 1''
;;
plotsym, 0, plot_sym_size, /fill
w = where(BOSS_imag le 23.0, N)
;oplot,  BOSS_imag[w], BOSS_BOSS21_NM_SUM[w],  ps=8,  color=blue

w = where(BOSS_boss21_nm_sum_err gt 0 and BOSS_imag le 23.0, N)
oploterror, BOSS_imag[w],    BOSS_BOSS21_NM_SUM[w], $
            BOSS_mag_err[w], BOSS_boss21_nm_sum_err[w], $
            /hibar, errthick=errthick, psym=8, color=boss21_color, ERRcolor=boss21_color
oploterror, BOSS_imag[w],    BOSS_boss21_nm_SUM[w], $
            BOSS_mag_err[w], BOSS_boss21_nm_SUM_Err[w], $
            /lobar, errthick=errthick, psym=8, color=boss21_color, ERRcolor=boss21_color


;;
;; POWER-LAW FITS
;;
;;  S I N G L E   POWER-LAW....
amplitude = 19.0
pl_slope  = 0.80

slope_x = findgen(10)/1.+14.
slope_x = findgen(5)/1.+17.
;slope_y = 10^(slope_x-18)
slope_y = (10^(slope_x-amplitude))^pl_slope
;oplot, slope_x, slope_y, thick=6, linestyle=2


;; D O U B L E   POWER-LAW....
N_zero   = 10.0
alpha_d  = 0.98
beta_d   = 0.15
m_nought = 19.10
m_running = (findgen(100)/10.)+15.0 

dmag  = m_running - m_nought
double_power_law = (10^((-1.)*alpha_d*dmag) ) +  (10^((-1.)*beta_d*dmag)  )
dN_dm = N_zero / double_power_law
;oplot, m_running, dN_dm, thick=8, linestyle=2, color=purple

;;
;; ACTUALLY FIT 
;; to the Corrected DR9 data...
 N_zero   =       2.63027
 alpha_d  =       1.50000
 beta_d   =      0.400000
 m_nought =       19.0000
dmag  = m_running - m_nought

double_power_law = (10^((-1.)*alpha_d*dmag) ) +  (10^((-1.)*beta_d*dmag)  )
dN_dm = N_zero / double_power_law
oplot, m_running, dN_dm, thick=8, linestyle=0, color=red


;; Fit to the BOSS21 data...
N_zero   =       43.6516
 alpha_d  =       0.10000
 beta_d   =      0.800000
 m_nought =       20.4000

dmag  = m_running - m_nought

double_power_law = (10^((-1.)*alpha_d*dmag) ) +  (10^((-1.)*beta_d*dmag)  )
dN_dm = N_zero / double_power_law
oplot, m_running, dN_dm, thick=8, linestyle=0, color=blue



;;
;;  G L O B A L    L A B E L S
;;
x_xyouts = 15.0
y_xyouts = 0.1
xyouts, x_xyouts, y_xyouts,'N(<i) (deg!E-2!N)', $
        charthick=charthick, charsize=charsize, ORIENTATION=90



;;
;;  L A B E L S
;;

;; PL label first
;xyouts, 15.1, 50., 'N(<i)!9?!3i!U0.80!N', charthick=4.2, charsize=2.6

x_pos = 18.4
x_off =  0.5

y_pos = -4.1
y_off = 0.35

up = 1.35
   top     = 10^(y_pos+(4*y_off)+up)
middle     = 10^(y_pos+(3*y_off)+up)
bottom     = 10^(y_pos+(2*y_off)+up)
below_btm  = 10^(y_pos+(  y_off)+up)
below_btmm = 10^(y_pos+up)

up = 1.50
   top_legd     = 10^(y_pos+(4*y_off)+up)
middle_legd     = 10^(y_pos+(3*y_off)+up)
bottom_legd     = 10^(y_pos+(2*y_off)+up)
below_btm_legd  = 10^(y_pos+( y_off)+up)
below_btmm_legd = 10^(y_pos+up)


;legend, [' '], linestyle=lines, color=32, pos=[0.1,9.4], box=0,
;thick=8, charthick=4.2 ;, charsize=1.8
if choice_plot_SDSS eq 'y' then begin
plotsym, 3, plot_sym_size, /fill
legend, [' '], pos=[x_pos, top_legd],         box=0, psym=8, color=black
plotsym, 0, plot_sym_size, /fill
legend, [' '], pos=[x_pos, middle_legd],      box=0, psym=8, color=black
endif

plotsym, 0, plot_sym_size
legend, [' '], pos=[x_pos, top_legd],      box=0, psym=8, color=red
plotsym, 0, plot_sym_size, /fill
legend, [' '], pos=[x_pos, middle_legd],      box=0, psym=8, color=red
legend, [' '], pos=[x_pos, bottom_legd],  box=0, psym=1, color=red, symsize=2.0, thick=6

plotsym, 0, plot_sym_size, /fill
legend, [' '], pos=[x_pos, below_btm_legd], box=0, psym=8, color=blue


charsize  = charsize/1.2
charthick = charthick*1.2
if choice_plot_SDSS eq 'y' then begin
xyouts, x_pos+x_off, top,       'SDSS 0.3<z<2.2', charsize=charsize, charthick=charthick, color=black
xyouts, x_pos+x_off, middle,    'SDSS   3<z<5  ', charsize=charsize, charthick=charthick, color=black
endif 
xyouts, x_pos+x_off,     top, 'BOSS', charsize=charsize, charthick=charthick, color=red
xyouts, x_pos+(4.5*x_off), top, '2.2<z<3.5', charsize=charsize, charthick=charthick, color=red

xyouts, x_pos+x_off,     middle, 'BOSS', charsize=charsize, charthick=charthick, color=red
xyouts, x_pos+(4.5*x_off), middle, '2.2<z<3.5 w/cor', charsize=charsize, charthick=charthick, color=red

xyouts, x_pos+x_off,     bottom, 'Stp82', charsize=charsize, charthick=charthick, color=red
xyouts, x_pos+(4.5*x_off), bottom, '2.2<z<3.5', charsize=charsize, charthick=charthick, color=red
;xyouts, x_pos+x_off, below_btm, 'boss21+MMT 1.0<z<2.2', charsize=charsize, charthick=charthick, color=blue
xyouts, x_pos+x_off,     below_btm, 'boss21',    charsize=charsize, charthick=charthick, color=blue
xyouts, x_pos+(4.5*x_off), below_btm, '1.0<z<2.2', charsize=charsize, charthick=charthick, color=blue
;xyouts, 18., below_btm, 'boss21 0.3<z<2.2', charsize=charsize, charthick=charthick, color=light_blue+8.
charsize  = charsize/1.2
charthick=charthick/1.2

device, /close
set_plot, 'X'     





end
