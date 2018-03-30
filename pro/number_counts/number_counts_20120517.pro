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
;; I/P datafiles...
;;
qsomaskdir  =filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")

;; 
;;  F U L L   spAll    (actually needed for Stripe 82 calcs)
;;
datafile = qsomaskdir+'/spall/mini-spAll-v5_4_45.fits'
;data = mrdfits(datafile, 1)

;; 
;; For  X D Q S O - C O R E    D R 9
;;
xdspallfile = qsomaskdir+'/compfiles/xdspallmask.fits'
;xdspall     = mrdfits(xdspallfile,1)
boss        = xdspall

completeness_level = 0.85
print
print, 'completeness_level', completeness_level
print

;; 
;;  ``boss 21''
;;
;; Notes: Hmmmm....
;;    really have to check the similarities/differences etc. 
;;    between .z, .z_person, .zscan
;;
boss21_full = mrdfits('/cos_pc19a_npr/BOSS/QLF/data/VAR_BOSS_boss21.fits',1)
qsos_in = where((strtrim(boss21_full.TYPE) eq 'QSO'     or $
                strtrim(boss21_full.TYPE) eq 'QSO_BAL' or $
                strtrim(boss21_full.TYPE) eq 'QSO_?') and $
                boss21_full.plug_dec le 1.25,                N_QSOs_in)                                         
print
print, 'No. of QSOs included in the analysis, ', N_qsos_in, '  with minmax(zscan[qsos_in])  ', minmax(boss21[qsos_in].zscan)
print
print
boss21 = boss21_full[qsos_in]



;;
;; Stripe 82 cuts on the spAll file...
;;
;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19 +2LL^40+2LL^41 +2LL^42+2LL^43+2ll^44
target_flag = 34084861508608LL
;; Arguably the most inclusive, most sensible one...
target_flag_ancil = 2LL^0+2LL^1+2ll^2+2LL^3+2LL^4+2LL^5      +2LL^7+2LL^8+2LL^9   

w_boss_S82 = where(data.ra ge 300.0 or data.ra lt 60.0 and $ 
                   (data.dec le 1.25 and data.dec ge -1.25) $
                   and ( ((data.boss_target1 and target_flag) ne 0) $
                         or ((data.ancillary_target2 and target_flag_ancil) ne 0) ) $ 
                   and   data.specprimary eq 1 and data.z_conf_person ge 3 and $
                   data.z_person ge 2.20 and data.z_person le 3.50, N_boss_s82)

boss_S82 = data[w_boss_S82]

;;
;; Areas
;;
area_DR9 = 2236.    ;; for completeness_level = 0.85
area_S82 =  220.
;;
;; The area for boss21 is almost rectangulare and would be ~15.5
;; deg^2, if so. I'll take the CY number of 14.5 deg^2 to be
;; consistent with the QLF measurements as well. 
;;
area_boss21 = 14.5


print
print, ' **** APPLYING REDDENING EXTICTION CORRECTION!!! **** '
print
print

imag_boss     = boss.PSFMAG[3]     ;- boss.EXTINCTION[3]
imag_boss_S82 = boss_S82.PSFMAG[3] ;- boss_S82.EXTINCTION[3]
imag_boss21   = boss21.PSFMAG[3]   ;- boss21.EXTINCTION[3]



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
openw, 10, 'BOSS_DR9_XDQSO_iband_numcounts_temp.dat'
openw, 11, 'BOSS_DR9_Stripe82_iband_numcounts_temp.dat'
openw, 12, 'BOSS_DR9_boss21_iband_numcounts_temp.dat'

for jj=0L, no_bins-1 do begin
   ;; Total Number Counts in 17.75 􏰑 i < 22.45 ;; in 0.25 mag bins...
   bright_mag  = bright_mag_limit+(jj*0.25)
   faint_mag   = bright_mag_limit+((jj+1)*0.25)
   mid_mag[jj] = (bright_mag+faint_mag)/2.
   
   ;;
   ;; For  X D Q S O  -  D R 9  
   ;;
   dr9 = where(boss.poly_weight ge completeness_level and $
               boss.z_person ge 2.20 and boss.z_person lt 3.50 and boss.z_conf_person ge 3,  $
               N_midz)
   dr9 = where(boss.poly_weight ge completeness_level and $
               boss.z_person ge 2.20 and boss.z_person lt 3.50 and boss.z_conf_person ge 3 and $
               imag_boss ge bright_mag and imag_boss lt faint_mag, $
               N_counts)
   dr9 = where(boss.poly_weight ge completeness_level and $
               boss.z_person ge 2.20 and boss.z_person lt 3.50 and boss.z_conf_person ge 3 and $
               imag_boss lt faint_mag, $
               N_counts_sum)
   
   printf, 10, bright_mag, faint_mag, mid_mag[jj], N_counts, N_counts/area_DR9, N_counts_sum, N_counts_sum/area_DR9
   
   
   ;;
   ;; For the 2.2 < z < 3.5 on the Stripe... 
   ;;
   s82 = where(boss_s82.z_person ge 2.20 and boss_s82.z_person lt 3.50 $
               and boss_s82.z_conf_person ge 3 and $
               imag_boss_S82 ge bright_mag and imag_boss_S82 lt faint_mag, $
               N_counts_S82)
   
   s82 = where(boss_s82.z_person ge 2.20 and boss_s82.z_person lt 3.50 and $
               boss_s82.z_conf_person ge 3 and $
               boss_s82.PSFMAG[3] lt faint_mag, $ 
               N_counts_S82_sum)
   
   printf, 11, bright_mag, faint_mag, mid_mag[jj], N_counts_S82, N_counts_S82/area_S82, N_counts_S82_sum, N_counts_S82_sum/area_
   
   w_outside = where(boss_s82.z_person ge 2.20 and boss_s82.z_person lt 3.50 $
                     and boss_s82.z_conf_person ge 3 and $
                     imag_boss_S82 le bright_mag_limit or imag_boss_S82 ge faint_mag_limit, N_out)
   
   
   mid_mag_limit[jj] = faint_mag-0.25

   no_count_bin_S82[jj]     = N_counts_S82
   no_count_bin_S82_err[jj] = sqrt(N_counts_S82)

   no_count_bin_S82_limit[jj] = N_counts_S82_limit


   ;;
   ;; For the 1.0 < z < 2.2  objects in boss21
   ;;
   ;; !! Take care since z_person != zscan here... 
   ;; (Why, I'm not sure!) 
   w_boss21_magbin = where(boss21.zscan ge 1.00 and $
;                           boss21.z_conf_person ge 3 and $  ;; can take this or leave it... hmmmmm.... 
                           imag_boss21 ge bright_mag and imag_boss21 lt faint_mag, $
                           N_counts_boss21)
;   w_boss21_magbin = where(boss21.z_person ge 1.00 and $ ;;   boss21.z_person lt 2.20 and $
;                           boss21.z_conf_person ge 3 and $
;                           boss21.PSFMAG[3] ge bright_mag and boss21[qsos_in].PSFMAG[3] lt faint_mag, $
 ;                          N_counts_boss21)

   print, jj, bright_mag, faint_mag, mid_mag[jj], N_counts_boss21, N_counts_boss21/area_boss21
   printf, 12, bright_mag, faint_mag, mid_mag[jj], N_counts_boss21, N_counts_boss21/area_boss21


;   print, jj, bright_mag, faint_mag, mid_mag[jj], no_count_bin_S82[jj],  no_count_bin_S82_limit[jj]
endfor
print
print
printf, 11, '# ', N_out, ' objects in S82 too bright or too faint...' 

;; Errors
no_count_bin_S82_err_hi = no_count_bin_S82+no_count_bin_S82_err
no_count_bin_S82_err_lo = no_count_bin_S82-no_count_bin_S82_err

close, 10
close, 11
close, 12

close, /all

end
