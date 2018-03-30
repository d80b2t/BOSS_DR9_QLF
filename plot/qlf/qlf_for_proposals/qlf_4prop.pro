;+
; NAME:
;       qlf_4prop.pro
; 
; PURPOSE:
;        To plot Quasar Luminosity Functions, as several M_i lines on 
;        the redshift-number density plane, a la e.g. Boyle et
;        al. 2000.
;
; CALLING SEQUENCE:
;       .run downsizing
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


readcol, '../../../pro/qlf/QLF_iband_sdss_boss_boss21_wSelfn_formatted_McGkcor.dat', $
;readcol, '../../../pro/qlf/QLF_iband_sdss_S82_narrow_boss21_wSelfn_formatted_McGkcor.dat', $
         z_full,  fo,  Mi_mean_full, foo, Mi_bin_full, fooo, $ 
         NQ_full, ba,  log_Phi_full, bar, error_full, $
         format='(d,a, d,a, d,a,  d,a, d,a, d)', /silent

red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]
red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50]; 4.00];, 4.50, 5.00]

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

 clr =[red, blue, green, yellow+8]

expand = 1.0

charsize  =  3.2*expand
charthick =  7.2*expand
xthick    =  6.0
ythick    =  6.0
XTICKLEN  = 0.03
YTICKLEN  = 0.03
xstyle = 1
ystyle = 1

;; x-ranges
;x_min = -18.001
xmin = -22.001
xmax = -30.50

;; y-ranges
ymin = -9.20
ymax = -5.10

;; xy_outs 
x_xyouts = -25.25
y_xyouts =  -5.70


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;
;;
;;
set_plot, 'ps'
device, filename='qlf_4props_temp.eps', $
        xsize=12.0, ysize=12.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

print, 'charthick', charthick

plotsym, 0, 1.4, /fill
plot, z_full, log_phi_full, $
      psym=8, $
      position=[0.20, 0.20, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=xstyle, ystyle=ystyle, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      /nodata, $
      xtitle='!6M!Ii!N[z=2]', $
      ytitle='!6log!I10!N !7U!6(M!Ii!N[z=2]) [Mpc!E-3!N mag!E-1!N]', $
      color=black

clr_step  = 18
linethick = 30

plotsym, 0, 1.4, /fill
sym_no = [0,4,8]
for ii=0ll, n_elements(red_bins)-2 do begin
   
   res = ii mod 2
   res_sym = ii mod 3
 
   if res eq 0 then  plotsym,  sym_no[res_sym], 2.8, /fill
   if res eq 1 then  plotsym,  sym_no[res_sym], 2.8, /fill
   

   w = where(z_full ge red_bins[ii] and z_full le red_bins[ii+1], N)

   res_clr = ii mod 4
   ;oplot, Mi_mean_full[w], log_Phi_full[w],   color=clr[res_clr], thick=6.8
   ;oplot, Mi_mean_full[w], log_Phi_full[w], ps=8, color=clr[res_clr]

   oplot, Mi_mean_full[w], log_Phi_full[w],       color=(ii*36), thick=7.4
   oplot, Mi_mean_full[w], log_Phi_full[w], ps=8, color=(ii*36.)
   
   print, ii, res, res_sym, ii*37

   red_bins_string = string(red_bins, format='(d5.2)')

;   legend, pos=[-22.4, -7.40+(ii*(-0.2))], ' ',   box=0, thick=14, psym = 8, color=clr[res_clr]
; 
;   xyouts, -22.6, -7.5+(-0.2*ii), red_bins_string[ii], charsize=2.82, charthick=8.,color=clr[res_clr]
;   xyouts, -23.62, -7.5+(-0.2*ii), '<z<', charsize=2.82, charthick=8.,color=clr[res_clr]
;   xyouts, -24.2, -7.5+(-0.2*ii), red_bins_string[ii+1], charsize=2.82, charthick=8.,color=clr[res_clr]

   legend, pos=[-22.4, -7.40+(ii*(-0.2))], ' ',   box=0, thick=14, psym = 8, color=(ii*36)
 
   xyouts, -22.6, -7.5+(-0.2*ii), red_bins_string[ii], charsize=2.82, charthick=8.,color=(ii*36)
   xyouts, -23.62, -7.5+(-0.2*ii), '!8<z<!6', charsize=2.82, charthick=8.,color=(ii*36)
   xyouts, -24.2, -7.5+(-0.2*ii), red_bins_string[ii+1], charsize=2.82, charthick=8.,color=(ii*36)

endfor

charthick = 7.0
xyouts, -24.0, -5.4,  '!6Quasar Luminosity Function', charsize=3.0, charthick=charthick,color=black
xyouts, -27.0, -5.65, '!6from SDSS DR3', charsize=3.0, charthick=charthick,color=black
xyouts, -27.0, -5.9,  '!6and  BOSS DR9', charsize=3.0, charthick=charthick,color=black
   

device, /close
set_plot, 'X'     






end
