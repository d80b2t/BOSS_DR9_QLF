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

data = hogg_mrdfits('../data/xdcore_pqsogt0.4_targets.sweeps.fits',nchunk=100,1)

psfmag_u = 22.5-2.5*alog10(data.psfflux[0] > 0.001)                                  
psfmag_g = 22.5-2.5*alog10(data.psfflux[1] > 0.001)                                  
psfmag_r = 22.5-2.5*alog10(data.psfflux[2] > 0.001)                                  
psfmag_i = 22.5-2.5*alog10(data.psfflux[3] > 0.001)                                  
psfmag_z = 22.5-2.5*alog10(data.psfflux[4] > 0.001)                                  


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
for ii=0L, n_elements(field)-1 do begin
   
   current_run = field[ii]
   w     = where(data.run eq  current_run, N_total)
   x     = where(data.run eq  current_run and data.PQSOLOWZ gt 0.400, N_lowz)
   y     = where(data.run eq  current_run and data.PQSOMIDZ gt 0.424, N_midz)
   z     = where(data.run eq  current_run and data.PQSOHIZ  gt 0.400,  N_hiz)

   fields_total[ii] = N_total/area[ii]
   fields_lowz[ii]  = N_lowz/area[ii]
   fields_midz[ii]  = N_midz/area[ii]
   fields_hiz[ii]  = N_hiz/area[ii]
   
   ;; NOTE MIGHT - NOT SURE - BE BETTER IN JUST ONE 
   ;; FOR LOOP, WITH MANY WHERE's....

   ;; Total Number Counts in 17.75 􏰑 i < 22.45
   ;; in 0.25 mag bins...
;   for jj=0L, 19 do begin
 ;     bright_mag = 17.50+(jj*0.25)
 ;     faint_mag  = 17.50+((jj+1)*0.25)
   
;   w = where(data.PQSOLOWZ ge 0.2 and PSFMAG_g ge 20.00 and PSFMAG_g lt 21.00, N)            

    ;  endfor
  
   print, ii,  field[ii],  area[ii], N_total/area[ii], N_lowz/area[ii], N_midz/area[ii], N_hiz/area[ii]

endfor




charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.05
YTICKLEN  = 0.05

binsize=0.01
plothist, data.PQSO, $
          bin=binsize, $
          /ylog, $
          xrange=[0.0,1.0], $
          xstyle=1, ystyle=1, $
          xthick=xthick, ythick=ythick, $
          XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
          charsize=charsize, $
          charthick=charthick/1.8, $
          thick=thick,$ 
          xtitle='P(QSO,z)', $
          ytitle='log_10(No. of objects)'
plothist, data.PQSOLOWZ, thick= thick, bin=0.01, /ylog, /over, color=80
plothist, data.PQSOMIDZ, thick= thick, bin=0.01, /ylog, /over, color=160
plothist, data.PQSOHIZ, thick= thick, bin=0.01, /ylog, /over, color=240 

xyouts, 0.6, 300., 'PQSOLOWZ', charsize=charsize, charthick=charthick/1.8, color=80
xyouts, 0.6, 100., 'PQSOMIDZ', charsize=charsize, charthick=charthick/1.8, color=160
xyouts, 0.6,  33., 'PQSOHIZ', charsize=charsize, charthick=charthick/1.8, color=240


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
xmin = 17.70      ; i-band PSF mag
xmax = 22.50

;; y-ranges
ymin =  0.1       ;; To be close to 2SLAQ QSO, Croom09b
ymax = 50.0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
;set_plot, 'ps'
;device, filename='temp.eps', $
;        xsize=8.0, ysize=8.0, $
 ;;       xoffset=0.2, yoffset=0.2, $
   ;     /inches, /color, /encapsulated

plotsym, 0, 0.4, /fill
;plot, (u[w_stars] - g[w_stars]), $
;      (g[w_stars] - r[w_stars]), $
;      psym=8, $
 ;     position=[0.20, 0.20, 0.98, 0.98], $
  ;    xrange=[xmin, xmax], $
   ;   yrange=[ymin, ymax], $
    ;  xstyle=1, ystyle=1, $
     ; xthick=xthick, ythick=ythick, $
;      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
 ;     charsize=charsize, $
  ;    charthick=charthick, $
   ;   thick=thick,$ 
    ;  /nodata, $
     ; xtitle=' (u - g) ', $
      ;ytitle=' (g - r) ', $
      ;color=black

;xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

;device, /close
;set_plot, 'X'     

end
