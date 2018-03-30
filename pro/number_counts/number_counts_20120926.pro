;+
; NAME:
;   number_counts.pro
;
; PURPOSE:
;   This is (hopefully!) just a fairly simple bit of code, that's
;   gonnae work out some QSO number counts for from SDSS....
;
; CALLING SEQUENCE:
;    .run number_counts
;
; INPUTS:
;   SDSS_DR8_field_areas.dat    
;       ;; the RUN number and areas of all the SDSS where 
;       ;; "RUN6" is the six-digit imaging run num- ber
;
;   xdcore_targets.sweeps.fits
;       ;; Kindly provided by ADM.  
;       ;; The Bovy11-XDQSO catalog, cut at P(QSO)>0.400.
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; USEFUL URLs:
;     http://data.sdss3.org/sas/dr8/groups/boss/photoObj/xdqso/xdcore/
;     http://data. sdss3.org/datamodel/files/BOSS_PHOTOOBJ/xdqso/xdcore/xdcore_RUN6.html
;   
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  
; NOTES:
;
; REVISION HISTORY:
;   18-May-2011  v0.0.1     NPR
;-

;;
;; CHOICES..
;;
choice_plot_powerlaws = 'y'
read, choice_plot_powerlaws, PROMPT=' - Plot (v. simple!!) power law fits??  y/n  '

;data = hogg_mrdfits('../data/xdcore_pqsogt0.4_targets.sweeps.fits',nchunk=100,1)


;boss = mrdfits('../../bossqsomask/trunk/data/spall/mini-spAll-v5_4_31.fits',1)
boss = mrdfits('../../bossqsomask/trunk/data/spall/mini-spAll-v5_4_31_wPSFmags.fits',1)

;; print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19
;;     1047552
;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41
;;          3298535930880
;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41+2LL^42+2LL^43
;;         16492675464192
;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41+2LL^42+2LL^43+2LL^44
;;         34084861508608LL
;target_flag = 3298535930880LL ;; 
;target_flag = 16492675464192LL ;; 
target_flag = 34084861508608LL

readcol, '../../Stripe82/Stripe_82_QSOs_unique.dat', $
         ra_S82_full, dec_S82_full, z_S82_full, gmag_S82_full
w = where( ((ra_S82_full le 45.0) or (ra_S82_full gt 317.0)) and $
           dec_S82_full ge -1.25 and dec_S82_full le 1.25, N) 
 ra_S82  =  ra_S82_full[w] 
dec_S82  = dec_S82_full[w] 
 z_S82   = z_S82_full[w]
gmag_S82 = gmag_S82_full[w]


readcol, 'SDSS_DR8_field_areas.dat', field, area

print
print
print, ' T O T A L   A R E A   F R O M   D R 8  ', total(area)
print
print

;;
;; Think this is the relatively simple, but key, loop...
;;
fields_total = fltarr(N_elements(field))
fields_lowz = fltarr(N_elements(field))
fields_midz = fltarr(N_elements(field))
fields_hiz = fltarr(N_elements(field))

faint_mag_limit  = 23.50
bright_mag_limit = 17.00
mag_binsize = 0.25
range   = (faint_mag_limit - bright_mag_limit)
no_bins = (faint_mag_limit - bright_mag_limit)/mag_binsize

mid_mag                 = fltarr(no_bins)
mid_mag_limit           = fltarr(no_bins)
no_count_bin_boss       = fltarr(no_bins)
no_count_bin_S82        = fltarr(no_bins)
no_count_bin_S82_err    = fltarr(no_bins)

no_count_bin_boss_noS82 = fltarr(no_bins)
;; Cumulative counts
no_count_bin_S82_limit  = fltarr(no_bins)

;for ii=0L, n_elements(field)-1 do begin
print, 'jj, bright_mag, faint_mag, mid_mag[jj], no_count_bin_S82[jj],  no_count_bin_S82_limit[jj]'

for jj=0L, no_bins-1 do begin
   ;; Total Number Counts in 17.75 􏰑 i < 22.45 ;; in 0.25 mag bins...
   bright_mag  = bright_mag_limit+(jj*0.25)
   faint_mag   = bright_mag_limit+((jj+1)*0.25)
   mid_mag[jj] = (bright_mag+faint_mag)/2.
   
   specz        = where(boss.z ge 2.20 and boss.z lt 3.50 and boss.zwarning eq 0,  N_midz)
   specz_magbin = where(boss.z ge 2.20 and boss.z lt 3.50 and boss.zwarning eq 0 and $
                        boss.PSFMAG[1] ge bright_mag and boss.PSFMAG[1] lt faint_mag, N_counts)
   no_count_bin_boss[jj] = N_counts 
   
   specz_magbin = where(z_S82 ge 2.20 and z_S82 lt 3.50 and $
                        gmag_S82 ge bright_mag and gmag_S82 lt faint_mag, $
                        N_counts_S82)
   specz_magbin = where(z_S82 ge 2.20 and z_S82 lt 3.50 and $
                        gmag_S82 lt faint_mag,    N_counts_S82_limit)
   mid_mag_limit[jj] = faint_mag-0.25

   no_count_bin_S82[jj]     = N_counts_S82
   no_count_bin_S82_err[jj] = sqrt(N_counts_S82)

   no_count_bin_S82_limit[jj] = N_counts_S82_limit

   print, jj, bright_mag, faint_mag, mid_mag[jj], no_count_bin_S82[jj],  no_count_bin_S82_limit[jj]
endfor
print
print

print, ' (Kinda rough) area for BOSS DR9 ',  0.715787*((180/!dpi)^2)
print
print
area_DR9      = 2349.7901  ;;
area_Stripe82 =  220.0

;; Errors
no_count_bin_S82_err_hi = no_count_bin_S82+no_count_bin_S82_err
no_count_bin_S82_err_lo = no_count_bin_S82-no_count_bin_S82_err

charsize=2.8
charthick=5.0
thick=4.2
xthick=4
ythick=4
XTICKLEN  = 0.025
YTICKLEN  = 0.025


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


;; Total Number Counts in 17.75 􏰑 i < 22.45
;; x-ranges
xmin = 17.00      ; i-band PSF mag
xmax = 24.00

;; y-ranges
ymin =  0.11    ;; To be close to 2SLAQ QSO, Croom09b
;ymax = 75.0       ;; Linear axis
ymax = 179.0       ;; Linear axis

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='BOSS_number_counts_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 0.4, /fill
plot, mid_mag, ((no_count_bin_S82*(1./mag_binsize))/area_Stripe82), $
;      psym=8, $
      position=[0.23, 0.20, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      /ylog, $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      ytickformat='(i3)', $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick*2.,$ 
;      /nodata, $
      xtitle='!8g!3-PSF mag', $
      ytitle='No. of objects /deg!E-2!3'
;      ytitle=' log10(No. objects) /deg!E-2!3'

oplot, mid_mag, ((no_count_bin_S82*(1./mag_binsize))/area_Stripe82),  psym=8
;; http://idlastro.gsfc.nasa.gov/ftp/pro/plot/oploterror.pro
;oploterror, mid_mag, ((no_count_bin_S82*(1./mag_binsize))/area_Stripe82), ((no_count_bin_S82_err*(1./mag_binsize))/area_Stripe82), psym=8, thick=thick*2.

oplot, mid_mag_limit, ((no_count_bin_S82_limit)/area_Stripe82), thick=thick*2., color=red
oplot, mid_mag_limit, ((no_count_bin_S82_limit)/area_Stripe82), thick=thick*2., color=red, psym=8

;; Power-law slopes...
pl1 = (findgen(15)*0.25)+17.0
pl2 = (findgen(10)*0.25)+20.0

if choice_plot_powerlaws eq 'y' then begin
   ;oplot, pl1, (pl1^(0.8)), thick=thick*2., color=blue
   ;oplot, pl1, (pl1^(0.3)), thick=thick*2., color=blue
   oplot, pl1, (exp((pl1*(1.3))-24.2)),  thick=thick*2., color=blue
   ;oplot, pl1, (exp((pl1*(0.8))-14.2)),  thick=thick*2., color=light_blue
   oplot, pl2, (exp((pl2*(0.3))-3.8)),  thick=thick*2., color=light_blue+24

   xyouts, 17.2, 8.0, 'log(N)!9?!31.3!8m!3 ',  charsize=charsize/1.4, charthick=charthick, color=blue
   xyouts, 17.2, 4.0, 'log(N)!9?!30.3!8m!3 ',  charsize=charsize/1.4, charthick=charthick, color=light_blue+24
endif


;; Labels 
;; FOR LINEAR (Y)-AXIS
;xyouts, xmax*0.75, ymax*0.80, '2.20<!8z!3<3.50 Quasars',  charsize=charsize, charthick=charthick
;xyouts, xmax*0.75, ymax*0.70, 'BOSS Stripe82',  charsize=charsize, charthick=charthick

;legend, '', position=[xmax*0.74, ymax*0.60], box=0, linestyle=0, color=red, thick=thick, charsize=1.3
;legend, 'Cumulative', position=[xmax*0.78, ymax*0.64], box=0, textcolors=red, charsize=charsize, charthick=charthick


xyouts, xmax*0.75, 95., '2.20<!8z!3<3.50 Quasars',  charsize=charsize, charthick=charthick*1.4
xyouts, xmax*0.75, 54., 'BOSS Stripe82',  charsize=charsize, charthick=charthick*1.4

legend, '', position=[xmax*0.74, 45.], box=0, linestyle=0, color=red, thick=thick, charsize=1.2
legend, 'cumulative', position=[xmax*0.78, 55.], box=0, textcolors=red, $
         charsize=charsize/1.1, charthick=charthick

;xyouts, xmax*0.40, ymax*0.60, 'All BOSS',  charsize=charsize, charthick=charthick, color=color

device, /close
set_plot, 'X'     


end
