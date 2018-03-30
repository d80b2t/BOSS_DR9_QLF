;+
; NAME:
;   number_counts_S82.pro
;
; PURPOSE:
;   This is just a simple bit of code, that works out some QSO number
;   counts for the 10k or so QSOs on the Stripe, and those associated
;   with the 2SLAQ QSOs 
;
; CALLING SEQUENCE:
;    .run number_counts_2slaq.pro 
;
; INPUTS:
;   2SLAQ_QSO_mini.cat
;       ;; the 54k object strong 2SLAQ QSO (input) catalog
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

print
print
readcol, '../../data/2SLAQ_QSO_mini.cat', $
         ra_2SLAQ_all, dec_2SLAQ_all, $
         umag_2SLAQ_all, gmag_2SLAQ_all, rmag_2SLAQ_all, imag_2SLAQ_all, zmag_2SLAQ_all, $
         z_2SLAQ_all 
print
print, '2SLAQ QSOs READ-IN...'
print
; w= where(z_2SLAQ_all ge 0.02, N_2SLAQ_usable)
; imag_2SLAQ = imag_2SLAQ_all[w]
; z_2SLAQ    = z_2SLAQ_all[w]
; print, 'No. of 2SLAQ QSOs that are being used... ', N_2SLAQ_usable


faint_mag_limit  = 22.50
bright_mag_limit = 15.50
mag_binsize = 0.25

range   = (faint_mag_limit - bright_mag_limit)
no_bins = (faint_mag_limit - bright_mag_limit)/mag_binsize

mid_mag                 = fltarr(no_bins)
mid_mag_limit           = fltarr(no_bins)


;; To match Fig 12 Richards et al. (2006)
low_z  = 0.30
high_z = 2.20

;------------------------------------------
; Full 2SLAQ QSOs area 
;
area_2slaq = 191.1
openw, 10, '2SLAQ_QSOs_iband_numcounts.dat' 


printf,10, '#          bright_mag, faint_mag, imag_2SLAQ, no_count_bin_2SLAQ, n_2slaq_imag'
for jj=0L, no_bins-1 do begin
   ;; Total Number Counts in 17.75 􏰑 i < 22.45 ;; in 0.25 mag bins...
   bright_mag  = bright_mag_limit+(jj*0.25)
   faint_mag   = bright_mag_limit+((jj+1)*0.25)
   mid_mag[jj] = (bright_mag+faint_mag)/2.
   
   ;------------------------------------------
   ; 2SLAQ QSOs
   ;
   w = where(imag_2SLAQ_all ge bright_mag and imag_2SLAQ_all lt faint_mag $
             and z_2SLAQ_all ge low_z and z_2SLAQ_all le high_z, N)

   ;print, bright_mag, faint_mag, mid_mag[jj], N, N/area_Stripe82
   printf, 10, bright_mag, faint_mag, mid_mag[jj], N, N/area_2slaq
   
endfor
close, 10

close, /all
print
print


end 

