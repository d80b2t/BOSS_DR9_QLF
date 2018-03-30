;+
; NAME:
;       qlf
;
; PURPOSE:
;       To calculate Quasar Luminosity fucntions
;
; CALLING SEQUENCE:
;       >  red, omega0=0.27, omegalambda=0.73, h100=0.719
;       >  .run qlf
;
; INPUTS:
;       An optical or infrared quasar redshift catalog.
;
; OUTPUTS:
;       various
;
; PROCEDURES CALLED:
;
; COMMENTS:
;       /usr/common/rsi/lib/general/LibAstro/ 
;-
print
print

;;
;; Setting up a flat, H0 = 71.9 kms^-1 Mpc^-1 cosmology
;; DLUMINOSITY will not work unless this is set-up...
;;
; red, omega0=0.27, omegalambda=0.73, h100=0.719
;;
;; Cosmology from R06:
print
print
print, '    red, omega0=0.30, omegalambda=0.70, h100=0.70    '
print
print
print
red, omega0=0.30, omegalambda=0.70, h100=0.70

readin = 'y'
read, readin, prompt='Are the variables already in memory?? y/n '
if readin eq 'n' then begin
   
   ;;
   ;;  S D S S    D R 3     (Richards et al. 2006)
   ;; 
   readcol, '../../data/SDSS_QSO_DR3.dat', $
            name_dr3, z_dr3, imag_dr3, Mi_dr3, del_gi_dr3, corr_dr3, $
            format='(a,f,f,f,f,f)', /silent
   print
   print, '  SDSS DR3 (from Richards et al. 2006) read-in', n_elements(z_dr3)
   print
   print
   
   ;;
   ;; Reading in Richards et al. (2006)
   ;; Table 1, the completeness corrections...
   readcol, '../../data/Richards06_Table01.dat', imag_comp, z_comp, point, radio, ext, /silent
   
   ;;
   ;; Reading in the Richards06 k-correction
   ;; Normalized at z=2. 
   readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
   print, '  Richards06_kcor.dat READ-IN', n_elements(kcor)
   print
   print
   
   ;;
   ;;
   ;;  2 S L A Q   Q S Os
   ;;
   ;;
   readcol, '../../data/2SLAQ_QSO_mini.cat', $
            ra_2SLAQ_all, dec_2SLAQ_all, $
            umag_2SLAQ_all, gmag_2SLAQ_all, rmag_2SLAQ_all, imag_2SLAQ_all, zmag_2SLAQ_all, $
            z_2SLAQ_all 
   print
   print, '2SLAQ QSOs READ-IN...'
   print
   w= where(z_2SLAQ_all ge 0.02, N_2SLAQ_usable)
   imag_2SLAQ = imag_2SLAQ_all[w]
   z_2SLAQ    = z_2SLAQ_all[w]
   print, 'No. of 2SLAQ QSOs that are being used... ', N_2SLAQ_usable

   ;;
   ;;  B O S S   Q U A S A R  spectroscopic information
   ;;
   ;; We either want:
   ;;    all the BOSS Stripe 82 data, w/o any completeness corrections, 
   ;; or the "XDCORE-SPALL" file 
   ;; or the boss21 area file...

   boss21 = mrdfits('/cos_pc19a_npr/BOSS/Stripe82/boss21/spall-v5_5_0_boss21_qsos.fits',1)
   print
   print, ' boss21 READ-IN, ', n_elements(boss21)
   print
   print 
   ;; 
   ;; These datasets should be made/sculpted some place else
   ;; e.g. $BOSSQSOMASK_DIR/trunk/data/
   qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
   xdspallfile = qsomaskdir+'/compfiles/xdspall.fits'
   data = mrdfits(xdspallfile,1)
      
   ;completeness_level = 0.50
   completeness_level = 0.70
   ;completeness_level = 0.95
   ;completeness_level = 0.99
   
   w_boss = where(data.poly_weight ge completeness_level and $
                  data.zWarning eq 0 and data.z ge 2.00, N_boss_qsos)
   ;; Here if completeness_level =0.700 and z ge 2.00, N should = 31,822

   ;;  O L D !!  Trying to reproduce June results...
;   data = mrdfits('../../../Stripe82/spAll-v5_4_14_S82.fits', 1)
;   w_boss = where( data.zWarning eq 0 and $
 ;                  (strtrim(data.OBJTYPE) eq 'QSO' or $
  ;                  strtrim(data.class)   eq 'QSO'), N_boss_qsos_wDups)

   ;;
   ;; THIS LINE SETS THE WHOLE THING UP!!
   ;;
   ;; For the full XDQSO sample...
   boss = data[w_boss] 
   ;; For Stripe 82
;   boss = mrdfits('../../../bossqsomask/trunk/data/spall/mini-spAll-v5_4_45_Stripe82_z2QSOs.fits',1)
   qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
   stripe82file = qsomaskdir+'spall/mini-spAll-v5_4_45_Stripe82_z2QSOs.fits'
   boss = mrdfits(stripe82file, 1)
   N_boss_qsos = n_elements(boss)
   
   print
   print
   print, '   Completeness Level set at  ', completeness_level
   print
   print, '   size(boss)  ', size(boss)
   print, '   N_boss_qsos ', N_boss_qsos
   print
   print

   
endif

if readin eq 'y' then begin
   print
   print
   print, '   size(boss)  ', size(boss)
   print, '   N_boss_qsos ', N_boss_qsos
   print
   print
endif

;;
;; BOSS COMPLETENESS CORRECTIONS FROM McGREER ET AL. (2011...)
;;
;readcol, '../../completeness/grid_fid_iPQSOMIDZ_sel.txt', $
;readcol, '../../completeness/grid_newcont_iPQSOMIDZ_sel.txt', $
readcol, '../../completeness/grid_dusty_iPQSOMIDZ_sel.txt', $
         i_start, i_end, z_start, z_end, selfn
print
print, '   Selection function READ-IN ', n_elements(selfn)
print
print 
;; 
;; So, instead of taking 
red_selfn   = 0.5*(z_start+z_end)
iband_selfn = 0.5*( i_start+i_end)

;; Running the DLUMINOSITY command from the red.pro routine
;; Returns LUMINOSITY DISTANCES (using the above cosmology)
;; Divide by 1e6 to get into Mpc
if readin eq 'n' then begin
   
   print, 'Doing DLUMS....'
   dlums_dr3    = DLUMINOSITY(z_dr3)    / 1e6
   ;dlums_comp   = DLUMINOSITY(z_comp)  / 1e6 
   dlums_boss   = DLUMINOSITY(boss.z)   / 1e6
   dlums_2slaq  = DLUMINOSITY(z_2slaq)  / 1e6
   dlums_boss21 = DLUMINOSITY(boss21.z) / 1e6
   
;; Selec. fnc.
;;   dlums_selfn  = DLUMINOSITY(red_selfn) / 1e6
   dlums_z_start  = DLUMINOSITY(z_start) / 1e6
   dlums_z_end  = DLUMINOSITY(z_end) / 1e6

   print, 'DLUMS calculated '
   print
   print
endif

;; Working out the ABSOLUTE MAGNITUDE for each object
;; Recall 1st year (UG!) notes: 
;; DIST_MOD = 5 log (D_L /10pc) which means if D_L is in Mpc:

;; Note for the K-correction:
;; From from Richards et al. (2006)
;;   "The sign convention of the K-correction, K(z), 
;;   is defined by Oke & Sandage (1968) as 
;;     m_intrinsic = m_observed - K(z). 
;; Also, Hogg+2002 has:
;;     m_R = M_Q + DM + K_QR where
;; where, you *observe* through R, abut want to know for emitted frame Q. 
;; 

 
Abs_iM_dr3    = imag_dr3         - (5 * alog10(dlums_dr3))    - 25.00 - kcor(fix(z_dr3/0.01)) 
;Abs_iM_comp  = iflux            - (5 * alog10(dlums_zphot))  - 25.00 - kcor(fix(zphot/0.01)) 
Abs_iM_boss   = boss.PSFMAG[3]   - (5 * alog10(dlums_boss))   - 25.00 - kcor(fix(boss.z/0.01)) 
Abs_iM_2slaq  = imag_2slaq       - (5 * alog10(dlums_2slaq))  - 25.00 - kcor(fix(z_2slaq/0.01)) 
Abs_iM_boss21 = boss21.PSFMAG[3] - (5 * alog10(dlums_boss21)) - 25.00 - kcor(fix(boss21.z/0.01)) 


;; Step up four Abs_iMs based on the combo of i_start/end and z_start end
;Abs_iM_selfn = iband_selfn    - (5 * alog10(dlums_selfn))  - 25.00 - kcor(fix(boss.z/0.01)) 
Abs_iM_istart_zstart = i_start  - (5 * alog10(dlums_z_start))  - 25.00 - kcor(fix(z_start/0.01)) 
Abs_iM_istart_zend   = i_start  - (5 * alog10(dlums_z_end))    - 25.00 - kcor(fix(z_end/0.01)) 
Abs_iM_iend_zstart   = i_end    - (5 * alog10(dlums_z_start))  - 25.00 - kcor(fix(z_start/0.01)) 
Abs_iM_iend_zend     = i_end    - (5 * alog10(dlums_z_end))    - 25.00 - kcor(fix(z_end/0.01)) 




print, 'Absolute i-band Mags calculated '
print
print

plot, boss.z, Abs_iM_boss, $
      ps=3, xrange=[1.8, 3.8], xstyle=1, yrange=[-23.0, -29.5], ystyle=1
oplot, z_start, Abs_iM_istart_zstart,color=100, ps=3
oplot, z_end, Abs_iM_iend_zend,color=220, ps=3


xthick    = 4.0 
ythick    = 4.0 
charsize  = 3.2
charthick = 6.0

diff_in_Mi_dr3s = (Mi_dr3 - Abs_iM_dr3)

plot, z_dr3, diff_in_Mi_dr3s, $
      ps=3,$
      xrange=[-0.2, 6.0], $
      xstyle=1, $
      ystyle=1, $
      xthick=xthick, $
      ythick=ythick, $
      charsize=charsize/1.4, $
      charthick=charthick/2.4, $
      xtitle=' redshift', $
      ytitle=' M_i (DR3)  -  M_i (NPR) '
      


;; Plot Redshift vs. i-band Abs. Mag as a sanity check
;; cf. Rig 17 Richards et al. (2006, AJ, 131, 2766)
set_plot, 'ps'
!p.multi=0
;device, filename='Abs_iMag_vs_zphot_Stripe82_temp.ps', $
device, filename='Abs_iMag_vs_zphot_xdspall_temp.ps', $
        /inches, /color, $
        xoffset=0.2, yoffset=0.2, $
        xsize=12, ysize=12

;plot, zphot, Abs_iM_zphot, $
plot, z_dr3, Abs_iM_dr3, $
      psym=3, $
      /nodata, $
      position=[0.22, 0.22, 0.98, 0.98], $
      xstyle=1, $
      ystyle=1, $
      xrange=[-0.2, 6.0], $
      yrange=[-17.0, -33.0], $
      xthick=xthick, $
      ythick=ythick, $
      charsize=charsize, $
      charthick=charthick, $
      xtitle='Redshift, z', $
      ytitle='M!Ii!N[z=0] '

loadct, 6
plotsym, 0, 0.3, /fill

oplot,  boss.z,  Abs_iM_boss, psy=8, color= red
;oplot, zphot_limit, Abs_Km_limit,  linestyle=0, thick=4, color=40
;oplot, zphot_limit, Abs_Km_bright, linestyle=0, thick=4, color=40
oplot,  boss21.z,  Abs_iM_boss21, psy=2, color= blue

;xyouts, 2.7, -23.0, 'K-band and i<20  ',  color=100, charsize=2.6, charthick=6
;xyouts, 2.7, -22.0, 'Objects w/ k-corr  ',  color=140, charsize=2.6, charthick=6
;xyouts, 2.7, -21.0, 'K-band limit = 18.4',  color=200, charsize=2.6, charthick=6
;xyouts, 2.7, -20.0, 'i-band limit = 21.3' ,  charsize=2.6, charthick=6

loadct, 0
device, /close
set_plot, 'X'
close, /all
print, 'Abs_iMag_vs_zphot_xdspall_temp.ps  made... '
print
print


;; THINK ONE WANTS TO DO THE COMPLETENESS CORRECTION OUTSIDE
;; THE ``QLF LOOP'', and essentially on an object-by-object basis...
print, '  size(boss) ', size(boss)
print
print

delta_IM  = 0.0362500 ;; (Am assuming each bin is supposed to be 0.0725 delta_m wide..)
delta_red = 0.025
count=0
boss_wgt = fltarr( N_boss_qsos)
for ii=0ll, N_boss_qsos-1 do begin
   count=count+1
;   w_selfn = where((Abs_iM_selfn-delta_IM) lt ABS_IM_BOSS[ii] and $
 ;                  (Abs_iM_selfn+delta_IM) gt ABS_IM_BOSS[ii] and $
  ;                 (red_selfn+delta_red)   ge boss[ii].z      and $
   ;                (red_selfn-delta_red)   lt boss[ii].z, N_selfn)
   w_selfn = where(boss[ii].psfmag[3] ge i_start and $
                   boss[ii].psfmag[3] lt i_end and $
                   boss[ii].z         ge z_start and $
                   boss[ii].z         lt z_end, N_selfn)

if N_selfn ne 1 then print, ' *********** ', ii, boss[ii].psfmag[3], boss[ii].z
   
   if boss[ii].z ge 2.00 then begin
      count=count+1
;     if N_selfn eq 0 then print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn
      if N_selfn eq 1 then begin
        ; print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn, selfn[w_selfn]
         boss_wgt[ii] = 1./selfn[w_selfn]
;         boss_wgt[ii] = 1. ;; IF e.g. Stripe 82
      endif
;      if N_selfn gt 1 then print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn, selfn[ w_selfn]
   endif
endfor


;;  A R E A S 
;;
;; Full Sky                       = 41,253   deg^2 = 1
;; Total Area Stripe 82 for BOSS  =    220.0 deg^2 = 0.00533295
;; 
;; print, ((180./!dpi)^2.)*4.*!dpi 

;;
;; The given (sector) areas, for different completenesses...
;; see e.g. bossqsomask/trunk/data/compfiles/Completeness_areas_temp.dat  
;;
;; So, if using XDQSO (as generated by the $BOSSQSOMASK_DIR .pro
;; routines AND setting zWarning = 0, then you need the "Overall"
;; completeness areas...

if completeness_level eq 0.700 then Area_boss= 1849.04 /41253.
Area_boss75 = 2623.98 /41253.
Area_boss95 = 1437.25/41253.

Area_Stp82 =  (220.0/41253.)
Area_dr3   =  (1622./41253.)
;; From Richards06 (p.4) ``The uncertainty on this area is of order 10
;; deg^2''
area_2slaq  = (191.1/41253.)
area_boss21 = ( 15.0/41253.)

;; The number to be actually used!!
;area_boss = Area_boss95
area_boss = Area_Stp82


;; From Richards06, Section 6.1
mag_bins = 0.30
red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

openw, 11, 'My_QLF_iband_dr3_temp.dat'
openw, 12, 'My_QLF_iband_boss_temp.dat'
openw, 13, 'My_QLF_iband_boss_wgt_temp.dat'
openw, 16, 'My_QLF_iband_2slaq_temp.dat'
openw, 17, 'My_QLF_iband_boss21_temp.dat'

openw, 14, 'My_QLF_iband_boss_narrowZ_temp.dat'
openw, 15, 'My_QLF_iband_boss_narrowZ_wgt_temp.dat'

print, 'ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs, mean_z, V_comov, Num_den, log_Num_den'

bright_i    = 15.00
faint_i_loz = 19.10
faint_i_hiz = 20.20

N_QSOs_boss_wgt = 0.0

;;
;; WORKING OUT THE QLF
;;
print
print, ' WORKING OUT THE QLF....'
print
for ii = 0L, 10 do begin 
   ;; Eeeek, argh...
   ;; R06 z-bins not evenly spaced, especially at high-z...
   redmin = red_bins[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 37 do begin 
                                Abs_iM_min = ((jj+1)*0.30) + (-30.90)
                                Abs_iM_max = (jj*0.30)     + (-30.90)
      ;; To get the bins at 0.00, 0.25, 0.50, 0.75...
                                ;Abs_iM_min = ((jj+1)*0.25)  + (-30.125)
;      Abs_iM_max = (jj*0.25)      + (-30.125)
      Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
      
      w_qlfbin_dr3 = where(Abs_IM_dr3 le Abs_iM_min and $
                           Abs_IM_dr3 gt Abs_iM_max and $
                           z_dr3 ge redmin and z_dr3 lt redmax, N_QSOs_dr3)
      
      w_qlfbin_2slaq = where(Abs_iM_2slaq lt Abs_iM_min and $
                             Abs_iM_2slaq gt Abs_iM_max and $
                             z_2slaq ge redmin and z_2slaq lt redmax, N_QSOs_2slaq)
      
      w_qlfbin_boss21 = where(Abs_iM_boss21 lt Abs_iM_min and $
                              Abs_iM_boss21 gt Abs_iM_max and $
                              boss21.z ge redmin and boss21.z lt redmax, N_QSOs_boss21)
      
      w_qlfbin_boss = where(Abs_iM_boss lt Abs_iM_min and $
                            Abs_iM_boss gt Abs_iM_max and $
                            boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)
      
      if N_QSOs_boss ne 0 then N_QSOs_boss_wgt = total(boss_wgt[w_qlfbin_boss])
      
;      w_selfn = where(Abs_iM_selfn lt Abs_iM_min and $
 ;                     Abs_iM_selfn gt Abs_iM_max and $
  ;                    red_selfn ge redmin and red_selfn lt redmax, N_selfn)

      ;;
      ;;  F o r   S D S S   D R 3  
      ;;
      if N_QSOs_dr3 gt 0 then begin
         ;mean_z =  mean(z_dr3[w_qlfbin_dr3])
         
         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)
         
         Num_den_sum = float(N_QSOs_dr3) / ((V_comov_top-V_comov_bottom)*area_dr3)
;         Num_den_sum = (float(N_QSOs_dr3)* (1./mean(corr_dr3[w_qlfbin_dr3])))  / ((V_comov_top-V_comov_bottom)*area_dr3)

         ;; Since we have 0.3 mag bins but want to display in whole
         ;; (1.0) Mag bin
         Num_den_sum_per_mag = Num_den_sum * (1.0 /mag_bins)
         log_Num_den_sum_dr3 = alog10(Num_den_sum_per_mag)

         error_dr3 = (10^log_Num_den_sum_dr3) /sqrt(N_QSOs_dr3)

         Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
         z_bin = (redmin + redmax) /2.
         ;print, ii, jj,  z_bin, Abs_mag_bin, N_QSOs_dr3
         printf, 11, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_dr3, N_QSOs_dr3, error_dr3, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
      endif

      ;;
      ;;  F o r   2 S L A Q   Q S Os
      ;;
      if N_QSOs_2slaq le 0 then printf, 16, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_2slaq,  999.999
      if N_QSOs_2slaq gt 0 then begin
         
         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)
         
         Num_den_sum = float(N_QSOs_2slaq) / ((V_comov_top-V_comov_bottom)*area_2slaq)
;         Num_den_sum = (float(N_QSOs_2slaq)* (1./mean(corr_dr3[w_qlfbin_dr3])))  / ((V_comov_top-V_comov_bottom)*area_dr3)
         
         ;; Since we have 0.3 mag bins but want to display in whole
         ;; (1.0) Mag bin
         Num_den_sum_per_mag = Num_den_sum * (1.0 /mag_bins)
         log_Num_den_sum_2slaq = alog10(Num_den_sum_per_mag)
         
         error_2slaq = (10^log_Num_den_sum_2slaq) /sqrt(N_QSOs_2slaq)
         
         Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
         z_bin = (redmin + redmax) /2.
         printf, 16, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_2slaq, N_QSOs_2slaq, error_2slaq, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
      endif

      ;;
      ;;  F o r   B O S S   Q S Os
      ;;
;      if N_QSOs_boss le 0 then printf, 12, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss,     999.999
      if N_QSOs_boss gt 0 then begin
         mean_z  = mean(boss[w_qlfbin_boss].z)

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum         = float(N_QSOs_boss)     / ((V_comov_top-V_comov_bottom)*Area_boss)
         Num_den_sum_wgt     = float(N_QSOs_boss_wgt) / ((V_comov_top-V_comov_bottom)*Area_boss)

         Num_den_sum_per_mag     = Num_den_sum     * (1.0 /mag_bins)
         Num_den_sum_per_mag_Wgt = Num_den_sum_wgt * (1.0 /mag_bins)

         log_Num_den_sum_boss     = alog10(Num_den_sum_per_mag)
         log_Num_den_sum_boss_wgt = alog10(Num_den_sum_per_mag_wgt)

         error_boss = (10^log_Num_den_sum_boss) /sqrt(N_QSOs_boss)

        ; Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
        ; z_bin = (redmin + redmax) /2.

       ;  if N_selfn ne 0 then print, ii,jj, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss, N_selfn, selfn[w_selfn]

         printf, 12, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss,     N_QSOs_boss, error_boss, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
         printf, 13, z_bin, Abs_mag_bin, Num_den_sum_wgt, log_Num_den_sum_boss_wgt, N_QSOs_boss_wgt, error_boss, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
      endif


      ;;
      ;;  F o r   B O S S   2 1    Q S Os
      ;;
      ;if N_QSOs_boss21 le 0 then printf, 16, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss,     999.999
      if N_QSOs_boss21 gt 0 then begin
         mean_z  = mean(boss[w_qlfbin_boss21].z)

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum         = float(N_QSOs_boss21)     / ((V_comov_top-V_comov_bottom)*Area_boss21)
         Num_den_sum_per_mag     = Num_den_sum     * (1.0 /mag_bins)

         log_Num_den_sum_boss21     = alog10(Num_den_sum_per_mag)

         error_boss21 = (10^log_Num_den_sum_boss21) /sqrt(N_QSOs_boss21)

         printf, 17, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss21, N_QSOs_boss21, error_boss21, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
      endif

   endfor
endfor
close, 11
close, 12
close, 13
close, 16
close, 17




mag_bins = 0.30
;red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]
red_bins = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]

print
print, ' WORKING OUT THE QLF "NARROW Z"....'
print
for ii = 0L, 10 do begin 
   ;; Eeeek, argh...
   ;; R06 z-bins not evenly spaced, especially at high-z...
   redmin = red_bins[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 37 do begin 
      Abs_iM_min = ((jj+1)*0.30) + (-30.90)
      Abs_iM_max = (jj*0.30)     + (-30.90)
      Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
      
      w_qlfbin_boss = where(Abs_iM_boss lt Abs_iM_min and $
                            Abs_iM_boss gt Abs_iM_max and $
                            boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)

      if N_QSOs_boss ne 0 then N_QSOs_boss_wgt = total(boss_wgt[w_qlfbin_boss])
      
      ;;
      ;; FOR BOSS QSOs
      ;;
      if N_QSOs_boss le 0 then printf, 14, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss, 999.999
      if N_QSOs_boss gt 0 then begin
         mean_z  = mean(boss[w_qlfbin_boss].z)

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum         = float(N_QSOs_boss)     / ((V_comov_top-V_comov_bottom)*Area_boss)
         Num_den_sum_wgt     = float(N_QSOs_boss_wgt) / ((V_comov_top-V_comov_bottom)*Area_boss)

         Num_den_sum_per_mag     = Num_den_sum     * (1.0 /mag_bins)
         Num_den_sum_per_mag_Wgt = Num_den_sum_wgt * (1.0 /mag_bins)

         log_Num_den_sum_boss     = alog10(Num_den_sum_per_mag)
         log_Num_den_sum_boss_wgt = alog10(Num_den_sum_per_mag_wgt)

         error_boss = (10^log_Num_den_sum_boss) /sqrt(N_QSOs_boss)

        ; Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
        ; z_bin = (redmin + redmax) /2.

       ;  if N_selfn ne 0 then print, ii,jj, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss, N_selfn, selfn[w_selfn]

         printf, 14, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss,     N_QSOs_boss, error_boss, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
         printf, 15, z_bin, Abs_mag_bin, Num_den_sum_wgt, log_Num_den_sum_boss_wgt, N_QSOs_boss_wgt
      endif

   endfor
endfor



print, 'QLF calculated...'

close, 14
close, 15



print
print
print, '    Area that went into calculation.... ', area_boss * ((180./!dpi)^2.)*4.*!dpi
print
print 

set_plot, 'X'
close, /all

end
