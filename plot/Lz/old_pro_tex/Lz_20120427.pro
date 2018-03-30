;+
; NAME:
;   program_name
;
; PURPOSE:
;   My "ultimate", all-singing, all-dancing, L-z plot...
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

  print
  print
  print, '------------------------------------------------------'
  print, '------------------------------------------------------'
  print, '!!!!! You have to run the "red" routine before this !!'
  print
  print, ' IDL> red         '
  print, '   red, omega0 = 0.228, omegalambda = 0.772, h100=0.702    '
  print
  ;; (Komatsu et al. 2011) 
;  print, '  red, omega0=0.30, omegalambda=0.70, h100=0.70   '
  print, '------------------------------------------------------'
  print, '------------------------------------------------------'
  print


;;
;;
;;   C H O I C E S !!!!
;;
;;      D A T A
;;
choice_plot_boss = 'y'
;read, choice_plot_boss, PROMPT=' - Plot BOSS points?? y/n  '

choice_plot_boss21 = 'y'
;read, choice_plot_boss, PROMPT=' - Plot BOSS21 points?? y/n  '

choice_plot_sdss = 'y'
;read, choice_plot_sdss, PROMPT=' - Plot SDSS DR7Q points?? y/n  '

choice_plot_2slaq = 'n'
;read, choice_plot_2slaq, PROMPT=' - Plot 2SLAQ points?? y/n  '

choice_plot_ibandlimits = 'y'
;read, choice_plot_2slaq, PROMPT=' - Plot i-band limits points?? y/n  '

;;
;;      ``SIDE-ON'' HISTOGRAMS
;;
choice_plot_side_nofz = 'y'
;read, choice_plot_2slaq, PROMPT=' - Plot i-band limits points?? y/n  '

choice_plot_side_AbsI = 'y'
;read, choice_plot_2slaq, PROMPT=' - Plot i-band limits points?? y/n  '


;;
;;      G R I D 
;;
choice_plot_grid = 'n'
;read, choice_plot_2slaq, PROMPT=' - Plot i-band limits points?? y/n  '


print
print, 'choice_plot_boss:         ', choice_plot_boss 
print, 'choice_plot_boss21:       ', choice_plot_boss21

print, 'choice_plot_sdss:         ', choice_plot_sdss 
print, 'choice_plot_2slaq:        ', choice_plot_2slaq
print, 'choice_plot_ibandlimits;  ', choice_plot_ibandlimits 
print
print, 'choice_plot_side_nofz     ', choice_plot_side_nofz
print, 'choice_plot_side_AbsI     ', choice_plot_side_AbsI
print
print, 'choice_plot_grid          ',  choice_plot_ibandlimits 
print



;;
;;
;; D R 7 Q
;;
dr7q_full = mrdfits('../../data/dr7qso.fits', 1)
print
print, 'SDSS DR7Q read-in  ', n_elements(dr7q_full)

w_dr7q_uniform = where(dr7q_full.USELFLAG eq 1 and dr7q_full.imag gt 0.0, N_dr7q_uniform)  
dr7q_uni = dr7q_full[w_dr7q_uniform]
print, 'SDSS DR7Q UNIFORMs         ', n_elements(dr7q_uni);, N_dr7q_uniform
;;
;; 10,000 random DR7Q quasars...
;;
result = RANDOMU(100, 10000)
res    = long(result*float(n_elements(dr7q_uni)))
dr7q = dr7q_uni[[res]]
print, 'SDSS DR7Q UNIFORMs to plot ', n_elements(dr7q);, N_dr7q_uniform
print
print

;;
;; B O S S
;;
qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
xdspallfile = qsomaskdir+'/compfiles/xdspallmask_DR9.fits'
boss = mrdfits(xdspallfile,1)

;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41+2LL^42+2LL^43
;;         16492675464192
target_flag = 34084861508608LL
boss_qsos   = where( ((boss.boss_target1 and target_flag) ne 0) and $
                     boss.PSFMAG[3] gt 0. and $
                     boss.zWarning eq 0 and $
                     boss.z ge 2.20 and boss.z lt 5.50, N_boss_qsos)
boss_z    = boss[boss_qsos].z
boss_imag = boss[boss_qsos].PSFMAG[3]
print
print
print, 'N_BOSS_QSOs (0.02<z<5.50), zWarn=0, i-PSG mag >0.       ', N_boss_qsos
print, '(but only actually going to be using boss[0:10000]....) '

result = RANDOMU(100, 10000)
res    = long(result*float(N_boss_qsos))
boss_z    = boss_z[res]
boss_imag = boss_imag[res]
print
print

;; 
;;  B O S S  2 1
;;
data = mrdfits('/cos_pc19a_npr/BOSS/QLF/data/VAR_BOSS_boss21.fits',1)
print
print, ' (Full) boss21 READ-IN, ', n_elements(data)
print
print 
qsos_in = where((strtrim(data.TYPE) eq 'QSO'      or $
                strtrim(data.TYPE) eq 'QSO_BAL'   or $
                strtrim(data.TYPE) eq 'QSO_?'  ) and $
                data.zscan ge 1.00,   N_QSOs_in)                                         
print
print, 'No. of QSOs from BOSS 21  ', N_qsos_in, '  with minmax(zscan[qsos_in])  ', minmax(data[qsos_in].zscan)
print
boss21 = data[qsos_in]


;;
;; 2 S L A Q
;;
readcol, '../../data/2SLAQ_QSO_mini.cat', ra_2slaq, dec_2slaq, $
         u_2slaq,  g_2slaq, r_2slaq, i_2slaq, z_2slaq, $
         red_2slaq
w = where(red_2slaq ge 0.02, N)
i_2slaq   = i_2slaq[w]
red_2slaq = red_2slaq[w]
print
print, 'N_2SLAQ_QSOs ', N_elements(i_2slaq)
print
print


;;
;; Reading in the Richards06 k-correction 
;; Normalized at z=2. 
readcol, '../../data/Richards06_Table04.dat', kcor_redshift, kcor, /silent
print
print, 'Richards06_kcor.dat READ-IN', n_elements(kcor)
print



print, 'Doing DLUMS..... '
print
if choice_plot_boss eq 'y' then begin
   dlums_boss  = DLUMINOSITY(boss_z)  / 1e6
   boss_Abs_iM =  boss_imag  - (5 * alog10(dlums_boss))  - 25.00 - kcor(fix(boss_z/0.01))
endif
print, '   BOSS Absolute i-band Mags calculated '

if choice_plot_boss21 eq 'y' then begin
   dlums_boss21  = DLUMINOSITY(boss21.zscan) / 1e6
   boss21_Abs_iM =  boss21.psfmag[3]  - (5 * alog10(dlums_boss21))  - 25.00 - kcor(fix(boss21.zscan/0.01))
endif
print, '   BOSS21 Absolute i-band Mags calculated '

if choice_plot_sdss eq 'y' then begin
   dlums_dr7q  = DLUMINOSITY(dr7q.z)    / 1e6
   dr7q_Abs_iM =  dr7q.IMAG - (5 * alog10(dlums_dr7q))   - 25.00 - kcor(fix(dr7q.z/0.01))
endif
print, '   SDSS DR7Q Absolute i-band Mags calculated '

if choice_plot_2slaq eq 'y' then begin
   dlums_2slaq    = DLUMINOSITY(red_2slaq)    / 1e6
   twoslaq_Abs_iM =  i_2slaq - (5 * alog10(dlums_2slaq))   - 25.00 - kcor(fix(red_2slaq/0.01))
endif
print, '   2SLAQ QSOs Absolute i-band Mags calculated '


;; set an i-band limit as a guide...
if choice_plot_ibandlimits eq 'y' then begin
   i_band_mag = findgen(600)
   i_band_mag_bright = findgen(600)
   
   for ii=0ll, N_elements(i_band_mag)-1 do begin
      i_band_mag_bright[ii] = 18.00
      i_band_mag[ii]        = 22.00
   endfor
   i_band_mag_red = (findgen(600)/100)+0.001
   
   dlums_ibandmag     = DLUMINOSITY(i_band_mag_red)  / 1e6
   Abs_i_limit        = i_band_mag        - (5* alog10(dlums_ibandmag)) - 25.00 - kcor(fix(i_band_mag_red/0.01))
   Abs_i_limit_bright = i_band_mag_bright - (5* alog10(dlums_ibandmag)) - 25.00 - kcor(fix(i_band_mag_red/0.01))
endif
print
print, 'Absolute i-band Mag LIMITS calculated '
print
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

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03



;; x-ranges
xmin_iband = 15.1
xmax_iband = 23.0

;; y-ranges
ymin_iband = -0.01
ymax_iband =  5.50


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  z  vs.  i-mag 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_plot, 'ps'
device, filename='iband_redshift_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

!p.multi=0
plotsym, 0, 0.4, /fill
plot,  boss_imag,  boss_z, $
       psym=8, $
       position=[0.20,0.20,0.96,0.96], $
       xrange=[xmin_iband, xmax_iband], $
       yrange=[ymin_iband, ymax_iband], $
       xstyle=1, ystyle=1, $
       xthick=xthick, ythick=ythick, $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       charsize=charsize, charthick=charthick, thick=thick,$ 
       xtitle=' i-band mag ', ytitle=' redshift  ', $
       /nodata, $
       color=black

plotsym, 0, 0.4, /fill
if choice_plot_sdss eq 'y' then begin
   oplot, dr7q.IMAG, dr7q.z,  psym=8, color=black
   xyouts, 15.9, 5.0, 'SDSS DR7Q', charsize=charsize, charthick=charthick*1.4, color=black
endif

plotsym, 0, 0.33, /fill
if choice_plot_boss eq 'y' then begin
   oplot, boss_imag, boss_z,  psym=8, color=red
   xyouts, 15.9, 4.5, 'BOSS DR9', charsize=charsize, charthick=charthick*1.4, color=red
endif

plotsym, 8, 0.33
if choice_plot_boss21 eq 'y' then begin
   oplot, boss21.psfmag[3], boss21.zscan,  psym=8, color=red
   xyouts, 15.9, 4.0, 'BOSS 21', charsize=charsize, charthick=charthick*1.4, color=red
endif


device, /close
set_plot, 'X' 




;; x-ranges
xmin = -0.1
xmax = 5.5

;; y-ranges
ymin= -20.00

;ymax = -15.00
;ymax = min(dr7q_Abs_iM) -.125
ymax = -30.00
;ymax = max( boss_Abs_iM)
BL_corner_x = 0.22
BL_corner_y = 0.22
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  R E D S H I F T   vs.  A B S O L U T E   M A G 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_plot, 'ps'
device, filename='Lz_temp.eps', $
        xsize=7.0, ysize=7.0, $
;        xsize=14.0, ysize=14.0, $
;       xsize=21.0, ysize=13.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

leq_sign = String(108B)
!p.multi=[0,2,2]
!p.font=-1

talk=1.0 ;;1.6 good if .eps is 14x14...

plot,  boss_z,  boss_Abs_iM, $
       psym=8, $
       position=[BL_corner_x,BL_corner_y,0.98,0.98], $
       xrange=[xmin, xmax], yrange=[ymin, ymax], $
       xstyle=1, ystyle=1, $
       xthick=xthick, ythick=ythick, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       xtickformat="(A1)", ytickformat="(A1)", $
;       xtitle=' redshift ', $
;       ytitle='M!I i !N', $  
       charsize=charsize*talk, charthick=charthick*talk, thick=thick, $
       /nodata


plotsym, 0, 0.3, /fill
if choice_plot_boss eq 'y' then begin
   oplot, boss_z,  boss_Abs_iM, psym=8, color=red
;   oplot, boss_z,  boss_Abs_iM, psym=8, color=black
endif

plotsym, 8, 0.4
if choice_plot_boss eq 'y' then begin
   oplot, boss21.zscan,  boss21_Abs_iM, psym=8, color=blue
;   oplot, boss_z,  boss_Abs_iM, psym=8, color=black
endif

plotsym, 0, 0.3, /fill
if choice_plot_sdss eq 'y' then begin
   oplot, dr7q.z,  dr7q_Abs_iM,  psym=8, color=black
endif

plotsym, 0, 0.3, /fill
if choice_plot_2slaq eq 'y' then begin
   oplot, red_2slaq,  twoslaq_Abs_iM,  psym=8, color=blue
endif

choice_plot_ibandlimits = 'n'
if choice_plot_ibandlimits eq 'y' then begin
   oplot, i_band_mag_red, Abs_i_limit,        linestyle=0, thick=8.0, color=orange
   oplot, i_band_mag_red, Abs_i_limit_bright, linestyle=0, thick=8.0, color=orange
endif

;;
;; Labels
;;
;; 2 samples:
;xyouts, 2.7, -22.0, 'SDSS DR7Q',  charsize=charsize, charthick=charthick*1.8, color=black
;xyouts, 2.7, -21.0, 'BOSS DR9 ',  charsize=charsize, charthick=charthick*1.8, color=red

;; 3 samples:
xyouts, 2.6, -22.25, 'SDSS DR7Q',  charsize=charsize, charthick=charthick*1.8, color=black
xyouts, 2.6, -21.5, 'BOSS DR9 ',  charsize=charsize, charthick=charthick*1.8, color=red
xyouts, 2.6, -20.75, 'BOSS S82 var',  charsize=charsize, charthick=charthick*1.8, color=blue

;xyouts, 2.7, -21.0, '2SLAQ QSO',  charsize=charsize*2.0, charthick=charthick*2.8, color=blue

if choice_plot_ibandlimits eq 'y' then begin
   xyouts, 2.7, -22.5, 'i-band=17.0 ',  charsize=charsize, charthick=charthick*1.8, color=red
   xyouts, 2.7, -21.5, 'i-band=21.8 ',  charsize=charsize, charthick=charthick*1.8, color=black
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;; PUTTING OUT L-z GRID ON THE PLOTS....
;; 
;; Croom09b redshift bins...
redshift_grid=[0.11,0.40, 0.68, 1.06, 1.06, 1.44, 1.82, 2.20, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

mag_grid = 25

grid_x = fltarr(N_elements(redshift_grid),mag_grid)  
grid_y = fltarr(N_elements(redshift_grid),mag_grid)  

for ii=0L,N_elements(redshift_grid)-1  do begin
   grid_x[ii,*] = redshift_grid[ii]
   grid_y[ii,*] = grid_x[ii]
endfor
   
for jj=0L,mag_grid-1 do begin
   dummy_mag = -23.00-(0.30*jj)  
   grid_y[*,jj]= dummy_mag
endfor

if choice_plot_grid eq 'y' then begin
   for ii=0L,N_elements(redshift_grid)-1 do oplot, grid_x[ii,*], grid_y[ii,*], color=green, thick=4., linestyle=2
   for jj=0L,mag_grid-1                  do oplot, grid_x[*,jj], grid_y[*,jj], color=green, thick=4., linestyle=2
endif


charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.05
YTICKLEN  = 0.05
xcharsize = 2.8 
ycharsize = 1.4
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  REDSHIFT HISTOGRAM
;;
plothist, boss_z, bin=0.05,  $
          xhist, yhist, /noplot

talk=0.8
if choice_plot_side_nofz eq 'y' then begin
   plot,   yhist,  xhist,  $  ;; doesn't really matter if /nodata, but probably should be xhist, yhist...
           position=[BL_corner_x, 0.14, 0.98, BL_corner_y], $
           xrange=[xmin, xmax], yrange=[min(yhist),max(yhist)*1.05], $
           xstyle=1, ystyle=1, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
           xtitle=' redshift ', $
           xthick=4.2, ythick=4.2*talk, thick=4.2*talk, $
           yticks=2, $ ;ytickformat='(a1)', $
           /nodata, $
           xcharsize=xcharsize*talk, ycharsize=ycharsize*talk, charthick=6.2
   
;AXIS, YAXIS=1, YRANGE =[ymin, ymax], YSTYLE = 1

   if choice_plot_boss eq 'y' then begin
      plothist, boss_z,       bin=0.05, /over, thick =4.2*talk, color=red
;      plothist, boss_z,       bin=0.05, xhist, yhist, /noplot
      plothist, boss21.zscan, bin=0.05, /over, peak=max(yhist), thick =4.2*talk, color=blue, linesyle=2
   endif
   
   if choice_plot_sdss eq 'y' then begin
;      plothist,  dr7q_full.z, bin=0.05, /over, thick =4.2*talk, color=black, linestyle=1
      plothist,  dr7q.z, bin=0.05, /over, thick =4.2*talk, color=black
   endif
   
   if choice_plot_2slaq eq 'y' then begin
      plothist,  red_2slaq, xhist, yhist, bin=0.05, /over, thick =4.2*talk, color=blue
      oplot, xhist, yhist*5., thick=thick,  color=blue
   endif
endif



charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03
xcharsize = 2.4 
ycharsize = 1.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Absolute Magnitude histgram (``on the side'')
;;

plothist, boss_Abs_iM, $
          bin = 0.1,  $
          ;peak = 0.50, $  ;; If xrange is from 0.0 to 0.59, then need this in here...
          xhist, yhist, /noplot
talk=0.8
if choice_plot_side_AbsI eq 'y' then begin
   plot, yhist,  xhist,  $
         position=[0.12,BL_corner_x,BL_corner_y,0.98], $
         xrange=[0.0,max(yhist)*1.05], yrange=[ymin, ymax], $
;         xrange=[0.0,0.59], yrange=[0.00, max(yhist)*1.05], $
         xstyle=1, ystyle=1, XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
         xthick=4.2, ythick=4.2*talk, thick=4.2*talk, $
         xcharsize=ycharsize*talk, ycharsize=xcharsize*talk, charthick=4.2, $
         xticks=2,  xtickformat='(a1)', $
         ytickformat='(a1)', $
         ytitle='M!I i !N', $  
         /nodata
;AXIS, YAXIS=1, YRANGE =[ymin, ymax], YSTYLE = 1

;;http://cow.physics.wisc.edu/~craigm/idl/archive/msg03916.html
numticks = 4
labels = [' ', '-22', '-24', '-26', '-28', '-30'] 

xpos = Replicate(!x.Window[0] -0.02, numticks+1)
;ypos = !y.Window[0] + (!y.Window[1] - !y.Window[0]) *
;                                      Findgen(numticks + 1) /
;                                      numticks
;; if we're going from -20.00 to -31.00, across ypos=.98 .26, thent
;; that's 13.09 per 2 mags...
;ypos = !y.Window[0] + (13.09) * Findgen(numticks + 1) / numticks
ypos = [0.26, 0.404, 0.548, 0.692, 0.836, 0.980]
;; Position the labels:
FOR j=0, numticks DO BEGIN
   XYOutS, xpos[j], ypos[j]-0.06, labels[j], $
           Alignment=0.0, Orientation=90, /Normal, $
           charsize=ycharsize*1.4, charthick=4.2
endfor

;;   
;; B O S S
;;
if choice_plot_boss eq 'y' then begin
   for i=0L, n_elements(xhist) - 1L do begin
      if i eq 0L or i eq n_elements(xhist) - 1L then begin 
         oplot, [0,yhist[i]],         [xhist[i]+0.05, xhist[i]+0.05], thick = 4.2*talk, color = red
         oplot, [yhist[i], yhist[i]], [xhist[i]-0.05, xhist[i]+0.05], thick =4.2*talk,  color = red
      endif else begin
         oplot, [yhist[i-1], yhist[i]], [xhist[i]-0.05, xhist[i]-0.05], thick =4.2*talk, color = red
         oplot, [yhist[i],yhist[i]],    [xhist[i]-0.05, xhist[i]+0.05], thick =4.2*talk, color = red
      endelse
   endfor
endif 

if CHOICE_PLOT_BOSS21 eq 'y' then begin
   plothist,  boss21_Abs_iM, bin = 0.1, xhist, yhist, peak=427,   /noplot 
   ;; peak =427 comes from yhist value for;; plothist, boss_Abs_iM, bin = 0.1
   for i=0L, n_elements(xhist) - 1L do begin
      if i eq 0L or i eq n_elements(xhist) - 1L then begin 
         oplot, [0,yhist[i]],         [xhist[i]+0.05, xhist[i]+0.05], thick=4.2*talk, color=blue
         oplot, [yhist[i], yhist[i]], [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color=blue
      endif else begin
         oplot, [yhist[i-1], yhist[i]], [xhist[i]-0.05, xhist[i]-0.05], thick=4.2*talk, color= blue
         oplot, [yhist[i],yhist[i]],    [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color= blue
      endelse
   endfor
endif


;;
;; S D S S   D R 7 Q
;;
   if choice_plot_sdss eq 'y' then begin
;      plothist, dr7q_Abs_iM, bin = 0.1, peak = 0.50, xhist, yhist,    /noplot 
      plothist, dr7q_Abs_iM, bin = 0.1, xhist, yhist,    /noplot 
      for i=0L, n_elements(xhist) - 1L do begin
         if i eq 0L or i eq n_elements(xhist) - 1L then begin 
            oplot, [0,yhist[i]],         [xhist[i]+0.05, xhist[i]+0.05], thick=4.2*talk, color=black
            oplot, [yhist[i], yhist[i]], [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color=black
         endif else begin
            oplot, [yhist[i-1], yhist[i]], [xhist[i]-0.05, xhist[i]-0.05], thick=4.2*talk, color= black
            oplot, [yhist[i],yhist[i]],    [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color= black
         endelse
      endfor
   endif

   if choice_plot_sdss eq 'y' then begin
;      plothist, dr7q_Abs_iM, bin = 0.1, peak = 0.50, xhist, yhist,    /noplot 
      plothist, dr7q_Abs_iM, bin = 0.1, xhist, yhist,    /noplot 
      for i=0L, n_elements(xhist) - 1L do begin
         if i eq 0L or i eq n_elements(xhist) - 1L then begin 
            oplot, [0,yhist[i]],         [xhist[i]+0.05, xhist[i]+0.05], thick=4.2*talk, color=black
            oplot, [yhist[i], yhist[i]], [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color=black
         endif else begin
            oplot, [yhist[i-1], yhist[i]], [xhist[i]-0.05, xhist[i]-0.05], thick=4.2*talk, color= black
            oplot, [yhist[i],yhist[i]],    [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color= black
         endelse
      endfor
   endif

   
;;
;; 2 S L A Q   Q S O s
;;
   if choice_plot_2slaq eq 'y' then begin
      plothist, twoslaq_Abs_iM, bin = 0.1, peak = 0.50, xhist, yhist,    /noplot 
      for i=0L, n_elements(xhist) - 1L do begin
         if i eq 0L or i eq n_elements(xhist) - 1L then begin 
            oplot, [0,yhist[i]],         [xhist[i]+0.05, xhist[i]+0.05], thick=4.2*talk, color=blue
            oplot, [yhist[i], yhist[i]], [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color=blue
         endif else begin
            oplot, [yhist[i-1], yhist[i]], [xhist[i]-0.05, xhist[i]-0.05], thick=4.2*talk, color= blue
            oplot, [yhist[i],yhist[i]],    [xhist[i]-0.05, xhist[i]+0.05], thick=4.2*talk, color= blue
         endelse
      endfor
   endif  ;  choice_plot_2slaq eq 'y' then 
   
endif     ;  choice_plot_side_AbsI eq 'y' then begin



!p.multi=0

device, /close
set_plot, 'X'     


end
