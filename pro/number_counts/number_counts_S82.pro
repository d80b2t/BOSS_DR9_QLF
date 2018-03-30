;+
; NAME:
;   number_counts_S82.pro
;
; PURPOSE:
;   This is just a fairly simple bit of code, that 
;   works out some QSO number counts for the 10k or
;   so QSOs on the Stripe, and those associated with 
;   boss21 
;
; CALLING SEQUENCE:
;    .run number_counts_s82
;
; INPUTS:
;   spall-v5_5_0_boss21  
;       ;; the v5_5_0 data for the 7 plates that were taken as part as
;       ;; the boss21 chunk...
;
; COMMENTS:
;
; USEFUL URLs:
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



data        = mrdfits('Stripe82_tenk_specPrim.fits',1)                                                             
data_boss21 = mrdfits('../../../data/spAll/spall-v5_5_0_boss21.fits', 1)

;------------------------------
; Only the Variable targets
;
w = where(((data_boss21.ancillary_target2 and 256L) NE 0), N)
;w = where(data_boss21.plug_dec le 1.25, N) ;; N=1776
data_boss21 = data_boss21[w]


faint_mag_limit  = 24.50
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

;------------------------------------------
; Full BOSS Stripe 82
;
area_Stripe82 = 220.0
openw, 10, 'BOSS_Stripe82_iband_numcounts.dat' 

;------------------------------------------
; BOSS boss21 (4 of 7 plates...)
;

area_boss21 = 5.96 ; 6.60442

openw, 11, 'BOSS_boss21_gband_numcounts.dat' 
openw, 12, 'BOSS_boss21_iband_numcounts.dat' 


printf,10, '# bright_mag, faint_mag, mid_mag[jj], no_count_bin_S82[jj],  no_count_bin_S82_limit[jj]'
printf,11, '# bright_mag, faint_mag, mid_mag[jj], no_count_bin_boss21[jj],  no_count_bin_boss21_limit[jj]'
printf,12, '# bright_mag, faint_mag, mid_mag[jj], no_count_bin_boss21[jj],  no_count_bin_boss21_limit[jj]'

for jj=0L, no_bins-1 do begin
   ;; Total Number Counts in 17.75 􏰑 i < 22.45 ;; in 0.25 mag bins...
   bright_mag  = bright_mag_limit+(jj*0.25)
   faint_mag   = bright_mag_limit+((jj+1)*0.25)
   mid_mag[jj] = (bright_mag+faint_mag)/2.
   
   ;------------------------------------------
   ; Full BOSS Stripe 82
   ;
   w = where(data.PSFMAG[3] ge bright_mag and data.psfmag[3] lt faint_mag $
             and data.z ge 2.20 and data.z le 3.50 and $
             data.zwarning eq 0, N)
   ;print, bright_mag, faint_mag, mid_mag[jj], N, N/area_Stripe82
   printf, 10, bright_mag, faint_mag, mid_mag[jj], N, N/area_Stripe82
   

   ;------------------------------------------
   ; BOSS boss21 (4 of 7 plates...)
   w = where(data_boss21.PSFMAG[1] ge bright_mag and data_boss21.psfmag[1] lt faint_mag $
             and data_boss21.z ge 0.4 and data_boss21.z le 2.1 and $
             data_boss21.zwarning eq 0, N)
   print, bright_mag, faint_mag, mid_mag[jj], N, N/area_boss21
   printf, 11, bright_mag, faint_mag, mid_mag[jj], N, N/area_boss21
   
   w = where(data_boss21.PSFMAG[3] ge bright_mag and data_boss21.psfmag[3] lt faint_mag $
             and data_boss21.z ge 0.3 and data_boss21.z le 2.2 and $
             data_boss21.zwarning eq 0, N)
   ;print, bright_mag, faint_mag, mid_mag[jj], N, N/area_boss21
   printf, 12, bright_mag, faint_mag, mid_mag[jj], N, N/area_boss21
   
endfor
close, 10
close, 11
close, 12

close, /all
print
print


end 

