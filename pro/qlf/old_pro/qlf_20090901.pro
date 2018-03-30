;;+
;; NAME:
;;       qlf
;;
;; PURPOSE:
;;       To calculate Quasar Luminosity fucntions
;;
;; CALLING SEQUENCE:
;;       >  red, omega0=0.27, omegalambda=0.73, h100=0.719
;;       >  .run qlf
;;
;; INPUTS:
;;       An optical or infrared quasar redshift catalog.
;;
;; OUTPUTS:
;;       various
;;
;; PROCEDURES CALLED:
;;
;; COMMENTS:
;;       /usr/common/rsi/lib/general/LibAstro/ 
;;-
print
print

;; Setting up a flat, H0 = 71.9 kms^-1 Mpc^-1 cosmology
;; DLUMINOSITY will not work unless this is set-up...
red, omega0=0.27, omegalambda=0.73, h100=0.719
print
print

readin = 'y'
read, readin, prompt='Are the variables already in memory?? y/n '
if readin eq 'n' then begin
   readcol, 'R09_catalog_reduced_Stp82.dat', $
            ra, dec, $
            redphot, redlow, redup, redprob, $
            u, g, r, iflux, z, $
            uerr, gerr, rerr, ierr, zerr 
   print, 'R09_catalog_reduced_Stp82.dat READ-IN', n_elements(ra)
   print
   
   ;; Reading in the Richards06 k-correction
   ;; Normalized at z=2. 
   readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
   print, 'Richards06_kcor.dat READ-IN', n_elements(kcor)
   print
   
   ;; Reading in a file that has the K-band and quasar redshift information
   ;;
   ;; This version of the Stripe82_matched_temp.dat file still
   ;; has 242 objects with redshifts = -1.0000
   
   ;readcol, 'Stripe82_matched.dat', raK, decK, ras, decs, g2, iflux, Kmag, zphot
   ;print, 'Stripe82_matched.dat READ-IN', n_elements(raK)
   
   readcol, 'lassdssmatchcat_mini_Stp82.dat', raK, decK, imag, imag_err, Kmag, Kmag_err, $
            zphot, zphot_err, zspec
   print, 'lassdssmatchcat_mini_Stp82.dat READ-IN', n_elements(raK)
   print
   print
endif
print
print


;; Running the DLUMINOSITY command from the red.pro routine
;; Returns LUMINOSITY DISTANCES (using the above cosmology)
;; Divide by 1e6 to get into Mpc

print, 'Doing DLUMS....'
dlums = DLUMINOSITY(redphot) / 1e6
dlumsK = DLUMINOSITY(zphot) / 1e6
print, 'DLUMS calculated '
print

;; Implementing the K-band FAINT  mag limit of 18.2 (or 18.4).
;; Implementing the K-band BRIGHT mag limit of 16.0
Kmag_limit    = 18.4
Kmag_bright   = 16.0

;; Setting up a "redshift" array that will
;; be used to produce Luminosity Distance values
zphot_limit   = (findgen(61))/10.

;; These luminosity distance values are then used
;; with the bright and faint Kmag limits for the 
;; e.g. Abs Mag vs. redshift plots.
dlumsK_limit =  DLUMINOSITY(zphot_limit) /1e6

;; Working out the ABSOLUTE MAGNITUDE for each object
;; Note DIST_MOD = 5 log (D_L /10pc) which means if D_L is in Mpc:
Abs_iM = iflux - (5 * alog10(dlums))  - 25.00 + kcor(fix(redphot/0.01))
print, 'Absolute i-band Mags calculated '
print


;; Picking out which objects are "good" objects for
;; K-band studies. 
w_Kmag      = where(Kmag gt 0.000, N_Kmag) ;; 11499 objecrs
w_good_Kmag = where(Kmag gt 0.000 and Kmag lt Kmag_limit, N_good_Kmag) ;; 9943 for i=18.2, 10900 for 18.4
w_Kmag_i20  = where(Kmag gt 0.000 and imag lt 20.0, N_Kmag_i20) ;; 5586 objects

;; Setting up the power-law k-correction for the K-bands
;; e.g. Richards et al. (2006)
alpha_nu     = -0.5 
redshift_arr = (findgen(61)/10.)
k_corr       = (-2.5*(1. + alpha_nu)) * (alog10(1 + redshift_arr))

;; Okay, this is where is starts to get slightly tricky, in that
;; to be consistent with R06, we should have our *continuum* 
;; k-correction equal zero at redshift 2, e.g. k_cont(z=2) = 0.000 
;; Thus since, 
;; k_corr = (-2.5*(1. + -0.5)) * (alog10(1 + 2)) = -0.596402
;; we have:
k_corr = k_corr + 0.596402

Abs_KM_total = Kmag[w_Kmag]      - (5 * alog10(dlumsK[w_Kmag]))      - 25.00 ;+ kcor(fix(zphot/0.01))
Abs_KM       = Kmag[w_good_Kmag] - (5 * alog10(dlumsK[w_good_Kmag])) - 25.00 ;+ k_corr(fix(zphot[w_good_Kmag]*10.))
Abs_KM_corr  = Kmag[w_good_Kmag] - (5 * alog10(dlumsK[w_good_Kmag])) - 25.00 + k_corr(fix(zphot[w_good_Kmag]*10.))
Abs_KM_i20   = Kmag[w_Kmag_i20]  - (5 * alog10(dlumsK[w_Kmag_i20])) - 25.00 ;+ kcor(fix(zphot/0.01))

Abs_Km_limit = Kmag_limit - (5 * alog10(dlumsK_limit)) - 25.00 ;+ k_corr(fix(zphot_limit*10.))
Abs_Km_bright = Kmag_bright - (5 * alog10(dlumsK_limit)) - 25.00 ;+ k_corr(fix(zphot[w_good_Kmag]*10.))

print, 'Abs_KMags calculated '
print, ' with Kmag_limit = ', Kmag_limit, ' and K_bright ', Kmag_bright
print
print

print, 'Kmags read-in', size(Kmag)



;; Plot Redshift vs. i-band Abs. Mag as a sanity check
;; cf. Rig 17 Richards et al. (2006, AJ, 131, 2766)
set_plot, 'ps'
!p.multi=0
device, filename='Abs_iMag_Kmag_vs_redphot_Stripe82_temp.ps', $
        /inches, /color, $
        xsize=12, ysize=12

plot, redphot, Abs_iM, $
      psym=3, $
      /nodata, $
      position=[0.22, 0.22, 0.98, 0.98], $
      xstyle=1, $
      ystyle=1, $
      xrange=[-0.2, 6.0], $
      yrange=[-17.0, -33.0], $
      xthick=4, $
      ythick=4, $
      charsize=3.2, $
      charthick=6, $
      xtitle='Redshift, z', $
      ytitle='M!Ii!N[z=0] '

loadct, 6
plotsym, 0, 0.3, /fill

oplot, zphot[w_Kmag],      Abs_KM_total,  psym=3, color=200
oplot, zphot[w_good_Kmag], Abs_KM,        psym=8, color=200
;oplot, zphot[w_good_Kmag], Abs_KM_corr,  psym=8, color=140
;oplot, zphot[w_Kmag_i20],  Abs_KM_i20,    psym=8, color=100

oplot, zphot_limit, Abs_Km_limit,  linestyle=0, thick=4, color=40
oplot, zphot_limit, Abs_Km_bright, linestyle=0, thick=4, color=40

xyouts, 2.7, -23.0, 'K-band and i<20  ',  color=100, charsize=2.6, charthick=6
xyouts, 2.7, -22.0, 'Objects w/ k-corr  ',  color=140, charsize=2.6, charthick=6
xyouts, 2.7, -21.0, 'K-band limit = 18.4',  color=200, charsize=2.6, charthick=6
xyouts, 2.7, -20.0, 'i-band limit = 21.3' ,  charsize=2.6, charthick=6


loadct, 0
device, /close
set_plot, 'X'
close, /all
print, 'Abs_iMag_vs_redphot_Stripe82_temp.ps made '
print
print




;; Area of Stripe 82 (Full Sky = 1 = 41,253 deg^2)
;; Area Stripe 82                    = 276.2 deg^2
;; Area Stripe 82 w/ K-band coverage = 276.2 deg^2
Area_Stp82      = 6.6950e-3
Area_Stp82_K    = 5.1348e-3  ;; probably a little less than this again...
;; THIS 2nd NUMBER STILL NEEDS TO BE DOUBLE    CHECKED!!!

openw, 11, 'My_QLF_iband_temp.dat'
openw, 12, 'My_QLF_Kband_temp.dat'
openw, 13, 'My_QLF_Kband_inter_temp.dat'
for ii = 0L, 15 do begin ;; redshift bins, delta_z=0.38 to z_max=5.70
   redmin  = ((ii)   * 0.38) - 0.08
   redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.

   for jj = 0L, 37 do begin 
      Abs_iM_min = ((jj+1)*0.30) + (-30.90)
      Abs_iM_max = (jj*0.30)     + (-30.90)
      
      Abs_KM_min = ((jj+1)*0.30) + (-32.90)
      Abs_KM_max = (jj*0.30)     + (-32.90)
      Abs_Kmag_bin = (Abs_KM_min + Abs_KM_max) /2.

      www = where(Abs_iM lt Abs_iM_min and Abs_iM gt Abs_iM_max $
                  and redphot ge redmin and redphot lt redmax, N_QSOs)
      
      wwK = where(Abs_KM lt Abs_KM_min and Abs_KM gt Abs_KM_max $
                  and zphot[w_good_Kmag] ge redmin and zphot[w_good_Kmag] lt redmax, N_QSO_Ks)

;; Be very careful here. 
;; Not including objects at i>20 my seem like a good idea due to
;; completeness issues, but what you end up of course doing is not
;; including the generally low-luminosity K-band objects, near the
;; faint limit, leading to "early" incompleteness at in fainter bins
;; at lower redshifts....
;; 
;      wwK = where(Abs_KM_i20 lt Abs_KM_min and Abs_KM_i20 gt Abs_KM_max $
;                 and zphot[w_Kmag_i20] ge redmin and zphot[w_Kmag_i20] lt redmax, N_QSO_Ks)


      printf, 13, ii, jj, z_bin, Abs_Kmag_bin, Abs_Km_limit[fix(z_bin*10.)], N_QSO_Ks  
      
      if N_QSO_Ks gt 0 then begin
         mean_z_K =  mean(zphot[wwK]) 

         ;; Divide by (1e6^3) to get into Mpc 
         V_comov_K  =  VCOMOVING(mean_z_K) / (1e6^3)
         Num_den_K  = (float(N_QSO_Ks))   / (V_comov_K * Area_Stp82_K)
         
         Num_den_sum_K = 0.00000
         for ll = 0L, N_QSO_Ks-1 do begin
            z_Kth_QSO     = zphot[wwK[ll]]
            V_comov_max_K =  VCOMOVING(z_Kth_QSO) / (1e6^3)
            Num_den_sum_K = (1.00/ (V_comov_max_K * Area_Stp82_K)) + Num_den_sum_K
         endfor

         Num_den_per_mag_K     = Num_den_K     * (1.0 /0.300)
         Num_den_sum_per_mag_K = Num_den_sum_K * (1.0 /0.300)

         log_Num_den_K     = alog10(Num_den_per_mag_K)
         log_Num_den_sum_K = alog10(Num_den_sum_per_mag_K)
;         print, ii, jj, redmin, redmax, Abs_KM_min, Abs_KM_max, N_QSO_Ks  
;                mean_z_K, V_comov, Num_den, log_Num_den

         z_bin = (redmin + redmax) /2.
         
         ;; Let's put in the K-band magnitude limit cut, since
         ;; we know we're *not* complete below this
         ;; (i.e. K_limit=18.4)
         ;; Mind the Abs. Mag convention, have to be more luminous
         ;; than mag limit, hence the lt in the logic.
         if Abs_KM_min lt Abs_Km_limit[fix(z_bin*10.)] then begin
            
            printf, 12, z_bin, Abs_Kmag_bin, log_Num_den_K, log_Num_den_sum_K, N_QSO_Ks
         endif
      endif
         

      if N_QSOs gt 0 then begin
         mean_z =  mean(redphot[www])
         
         ;; Divide by (1e6^3) to get into Mpc 
         V_comov =  VCOMOVING(mean_z) / (1e6^3)
         Num_den = (float(N_QSOs)) / (V_comov * Area_Stp82)

         ;; Alternatively you could implement the 1/Vmax here:
         ;;
         Num_den_sum = 0.00000
         for kk = 0L, N_QSOs-1 do begin
            z_Nth_QSO   = redphot[www[kk]]
            V_comov_max =  VCOMOVING(z_Nth_QSO) / (1e6^3)
            Num_den_sum = 1.00 / (V_comov_max * Area_Stp82) + Num_den_sum
         endfor
         
         ;; Since we have 0.3 mag bins but want to display in whole
         ;; (1.0) Mag bin
         Num_den_per_mag     = Num_den     * (1.0 /0.300)
         Num_den_sum_per_mag = Num_den_sum * (1.0 /0.300)
         
         log_Num_den     = alog10(Num_den_per_mag)
         log_Num_den_sum = alog10(Num_den_sum_per_mag)
         
;         print, ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs, $
;                mean_z, V_comov, Num_den, log_Num_den
         
         Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
         z_bin = (redmin + redmax) /2.
         
         printf, 11, z_bin, Abs_mag_bin, log_Num_den, log_Num_den_sum, N_QSOs
      endif

;PSU_IDL> print, mean(redphot[w])
;      2.01601
;PSU_IDL> print, VCOMOVING(2.01601)
;   5.9002174e+29
;PSU_IDL> print, VCOMOVING(2.01601) / (1e6^3)
;   5.9002175e+11
;PSU_IDL> print, 179. / 5.9002175e11
;  3.03379e-10
;PSU_IDL> print, 179. / (5.9002175e11*6.695e-3)
;  4.53142e-08
;PSU_IDL> print, alog10(4.53142e-8)  
;     -7.34377
;
   endfor
endfor
print, 'QLF calculated...'
close, 11
close, 12
close, 13







;; Starting to set up the 1/Vmax calculation
;;
; imin       = 14.327000
; imin_stp82 = 14.605000 
;; though only a VERY few, 10s, of objects are <14.6
;; so call it 14.5 for sake of argument.
;;
;;
;; However, for reasons still generally unknown, this
;; 1/Vmax calculation is still giving number densities
;; that are ~10-1000 times too small. (i.e. the Vmaxs 
;; themselves are too large).
;; 
imin = 14.50
imax = 21.30

;; Area of Stripe 82 (Full Sky = 1)
Area_Stp82 = 6.6950e-3

Dlum_min = 10^( (imin - Abs_iM -  25.00) / 5.00) 
Dlum_max = 10^( (imax - Abs_iM - 25.00) / 5.00) 
print, 'Dlum_min, Dlum_max calculated'
print

;; D_comov also called D_M in Hogg99
D_comov_min = Dlum_min / (1.00 + redphot)
D_comov_max = Dlum_max / (1.00 + redphot)
print, 'D_comov_min, D_comov_max calculated'
print


V_comov_min = ( (4d0 * !dpi)/3d0 * (D_comov_min^3)) * Area_Stp82 
V_comov_max = ( (4d0 * !dpi)/3d0 * (D_comov_max^3)) * Area_Stp82 
print, 'V_comovs done'
print

Vmax = V_comov_max - V_comov_min


;for i=0L, n_elements(redphot)-1 do begin
;   ;; To stop 'crazy' objects with M_i=-32 or the 
;   ;; redphot = -1 objects
;   if (Abs_iM[i] lt -15.00 and Abs_iM[i] gt -32.00) then begin
;      ;; Figure out which Abs_mag bin the object is in
;      binn = fix((Abs_m +34.00) *3.333)
;      ;; Do the 1/Vmax summation...
;      sum_one_over_Vmax[bin] = sum_one_over_Vmax[bin]  + (1.00/Vmax[i])
;   endif
;endfor









set_plot, 'X'
close, /all

end
