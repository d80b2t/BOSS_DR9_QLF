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

;; Setting up a flat, H0 = 71.9 kms^-1 Mpc^-1 cosmology
;; DLUMINOSITY will not work unless this is set-up...
;red, omega0=0.27, omegalambda=0.73, h100=0.719
;;
;; 
;; Cosmology from R06:
red, omega0=0.30, omegalambda=0.70, h100=0.70

print
print

readin = 'y'
read, readin, prompt='Are the variables already in memory?? y/n '
if readin eq 'n' then begin
   
   ;;
   ;;  S D S S    D R 3     (Richards et al. 2006)
   ;; 
   readcol, '../data/SDSS_QSO_DR3.dat', name_dr3, z_dr3, imag_dr3, Mi_dr3, del_gi_dr3, corr_dr3, format='(a,f,f,f,f,f)', /silent
   print
   print, ' SDSS DR3 (from Richards et al. 2006) read-in'
   print
   print

   ;;
   ;;  S D S S    D R 6 p Q          (Richards et al. 2009)
   ;;
   readcol, '../data/R09_catalog_reduced_Stp82.dat', $
            ra, dec, $
            zphot, zphotlow, zphothi, zphotprob, $
            u, g, r, iflux, z, $
            uerr, gerr, rerr, ierr, zerr 
   print, 'R09_catalog_reduced_Stp82.dat READ-IN', n_elements(ra)
   print
   print
   
   ;;
   ;; Reading in the Richards06 k-correction
   ;; Normalized at z=2. 
   readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
   print, 'Richards06_kcor.dat READ-IN', n_elements(kcor)
   print
   print
   
   ;;
   ;;  B O S S   Q U A S A R  spectroscopic information
   ;;
   data = mrdfits('../../Stripe82/spAll-v5_4_14_S82.fits', 1)

   ;;print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41
   ;;       3298535930880
   ;; print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19+2LL^40+2LL^41+2LL^42+2LL^43
   ;;     16492675464192LL
   target_flag = 16492675464192LL ;; 

   boss_qsos     = where( ((data.boss_target1 and target_flag) ne 0) $
                          and data.zWarning eq 0 and $
                          data.z ge 0.02 and data.z lt 6.10, N_boss_qsos)
   
   boss_qsos_hiz = where( ((data.boss_target1 and target_flag) ne 0) $
                          and data.zWarning eq 0 and $
                          data.z ge 2.20 and data.z lt 5.10, N_boss_qsos_hiz)

   w_boss = where( data.zWarning eq 0 and $
                   (strtrim(data.OBJTYPE) eq 'QSO' or $
                    strtrim(data.class) eq 'QSO'), N_boss_qsos)
   boss = data[w_boss]

   print
   print
   print, '  N_boss_qsos ', N_boss_qsos
   print, '  *** REMOVING DUPS!!! *** '
   print
   boss = boss[rem_dup(boss.ra)]
   print, '  size(boss) ', size(boss)
   print
endif
print


;; Running the DLUMINOSITY command from the red.pro routine
;; Returns LUMINOSITY DISTANCES (using the above cosmology)
;; Divide by 1e6 to get into Mpc

print, 'Doing DLUMS....'
dlums_dr3   = DLUMINOSITY(z_dr3)  / 1e6
;dlums_zphot = DLUMINOSITY(zphot)  / 1e6
;dlums_zphot = DLUMINOSITY(zphot)  / 1e6
dlums_boss  = DLUMINOSITY(boss.z) / 1e6

print, 'DLUMS calculated '
print
print



;;
;; Simple k-correction just for the time being...
;;
k_corr = -0.5
;; Setting up the power-law k-correction for the K-bands
;; e.g. Richards et al. (2006)
alpha_nu     = -0.5 
redshift_arr = (findgen(61)/10.)
;k_corr       = (-2.5*(1. + alpha_nu)) * (alog10(1 + redshift_arr))

;; Okay, this is where is starts to get slightly tricky, in that
;; to be consistent with R06, we should have our *continuum* 
;; k-correction equal zero at redshift 2, e.g. k_cont(z=2) = 0.000 
;; Thus since, 
;; k_corr = (-2.5*(1. + -0.5)) * (alog10(1 + 2)) = -0.596402
;; we have:
;
;k_corr = k_corr + 0.596402


;; Working out the ABSOLUTE MAGNITUDE for each object
;; Note DIST_MOD = 5 log (D_L /10pc) which means if D_L is in Mpc:
Abs_iM_dr3    = imag_dr3       - (5 * alog10(dlums_dr3))   - 25.00 - kcor(fix(z_dr3/0.01)) 
;Abs_iM_zphot  = iflux          - (5 * alog10(dlums_zphot)) - 25.00 - kcor(fix(zphot/0.01)) 
Abs_iM_boss   = boss.PSFMAG[3] - (5 * alog10(dlums_boss))  - 25.00 - kcor(fix(boss.z/0.01)) 

print, 'Absolute i-band Mags calculated '
print
print


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
device, filename='Abs_iMag_vs_zphot_Stripe82_temp.ps', $
        /inches, /color, $
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

;xyouts, 2.7, -23.0, 'K-band and i<20  ',  color=100, charsize=2.6, charthick=6
;xyouts, 2.7, -22.0, 'Objects w/ k-corr  ',  color=140, charsize=2.6, charthick=6
;xyouts, 2.7, -21.0, 'K-band limit = 18.4',  color=200, charsize=2.6, charthick=6
;xyouts, 2.7, -20.0, 'i-band limit = 21.3' ,  charsize=2.6, charthick=6

loadct, 0
device, /close
set_plot, 'X'
close, /all
print, 'Abs_iMag_vs_zphot_Stripe82_temp.ps made '
print
print



;; Area of Stripe 82 
;; Full Sky                 = 41,253 deg^2 = 1
;; Area Stripe 82 for BOSS  = 220.0  deg^2 = 0.00533295
;; 

Area_Stp82 =  (220.0/41253.)
Area_dr3   =  (1622./41253.)
;; THIS 2nd NUMBER STILL NEEDS TO BE PRECISION CHECKED!!!

;; From Richards06, Section 6.1
mag_bins = 0.30
red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

openw, 11, 'My_QLF_iband_dr3_temp.dat'
openw, 12, 'My_QLF_iband_boss_temp.dat'
print, 'ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs, mean_z, V_comov, Num_den, log_Num_den'

bright_i    = 15.00
faint_i_loz = 19.10
faint_i_hiz = 20.20

for ii = 0L, 10 do begin 
   ;; Eeeek, argh...
   ;; R06 z-bins not evenly spaced, especially at high-z...
   redmin = red_bins[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 37 do begin 
      Abs_iM_min = ((jj+1)*0.30) + (-30.90)
      Abs_iM_max = (jj*0.30)     + (-30.90)
      
      w_qlfbin_dr3 = where(Abs_IM_dr3 lt Abs_iM_min and $
                           Abs_IM_dr3 gt Abs_iM_max $
                           and z_dr3 ge redmin and z_dr3 lt redmax, N_QSOs_dr3)
      
      w_qlfbin_boss = where(Abs_iM_boss lt Abs_iM_min and $
                            Abs_iM_boss gt Abs_iM_max $
                            and boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)
      
;      w_qlfbin_zphot = where(Abs_iM_zphot lt Abs_iM_min and $
;                             Abs_iM_zphot gt Abs_iM_max $
;                             and zphot ge redmin and zphot lt redmax, N_QSOs_zphot)
      
;      print, ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs_dr3;, mean_z, V_comov, Num_den, log_Num_den
;      print, ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs_dr3;, mean_z, V_comov, Num_den, log_Num_den
      
      ;; FOR SDSS DR3
      if N_QSOs_dr3 gt 0 then begin
         ;mean_z =  mean(z_dr3[w_qlfbin_dr3])

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum = float(N_QSOs_dr3) / ((V_comov_top-V_comov_bottom)*area_dr3)

         ;; Alternatively you could implement the 1/Vmax here:
;         Num_den_sum = 0.00000
;         for kk = 0L, N_QSOs_dr3-1 do begin
;            z_Nth_QSO   = z_dr3[w_qlfbin_dr3[kk]]
            ;; works *reasonably* well...
;            V_comov_max =  VCOMOVING(z_Nth_QSO) / (1e6^3)   
            
;            if (z_dr3[w_qlfbin_dr3[kk]] le 2.90 and imag_dr3[w_qlfbin_dr3[kk]] le 19.10) then begin
 ;              DL_diffs = (10^((faint_i_loz-(Abs_iM_dr3[w_qlfbin_dr3[kk]])-25.0-kcor(fix(z_dr3[w_qlfbin_dr3[kk]]/0.01)))/5.0)) ; - $
  ;                        (10^((bright_i - (Abs_iM_dr3[w_qlfbin_dr3[kk]]) -25.0-kcor(fix(z_dr3[w_qlfbin_dr3[kk]]/0.01))) /5.0))
   ;            V_comov_max = (( DL_diffs / (1+  z_Nth_QSO))^3.)* ((4*!dpi)/3.)
    ;                            ;              print, z_dr3[w_qlfbin_dr3[kk]],  V_comov_max
     ;       endif else begin
      ;         DL_diffs = (10^((faint_i_hiz-(Abs_iM_dr3[w_qlfbin_dr3[kk]])-25.0-kcor(fix(z_dr3[w_qlfbin_dr3[kk]]/0.01)))/5.0)) ;- $
       ;                   (10^((bright_i - (Abs_iM_dr3[w_qlfbin_dr3[kk]]) -25.0-kcor(fix(z_dr3[w_qlfbin_dr3[kk]]/0.01))) /5.0))
        ;       V_comov_max = (( DL_diffs / (1+  z_Nth_QSO))^3.)* ((4*!dpi)/3.)
         ;   endelse
          ;
  
;            Num_den_sum = (1.00 / (V_comov_max * Area_dr3)) + Num_den_sum
 ;           if (ii eq 4 and jj eq 11) then begin
  ;             
   ;            print, kk, z_nth_QSO,  V_comov_max, (1.00 / (V_comov_max * Area_dr3)), Num_den_sum, $
    ;                  (Num_den_sum)*(1.0/mag_bins), (((Num_den_sum)*(1.0/mag_bins))/(10^(-6.87)))
     ;       endif
;         endfor
         
         ;; Since we have 0.3 mag bins but want to display in whole
         ;; (1.0) Mag bin
         Num_den_sum_per_mag = Num_den_sum * (1.0 /mag_bins)
         log_Num_den_sum_dr3 = alog10(Num_den_sum_per_mag)
         
         Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
         z_bin = (redmin + redmax) /2.
         printf, 11, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_dr3, N_QSOs_dr3
      endif
      
      ;;
      ;; FOR BOSS QSOs
      ;;
      if N_QSOs_boss gt 0 then begin
         mean_z  = mean(boss[w_qlfbin_boss].z)

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum = float(N_QSOs_boss) / ((V_comov_top-V_comov_bottom)*Area_Stp82)
         Num_den_sum_per_mag = Num_den_sum * (1.0 /mag_bins)
         
         log_Num_den_sum_boss = alog10(Num_den_sum_per_mag)
         
         Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
         z_bin = (redmin + redmax) /2.
         printf, 12, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss
      endif

   endfor
endfor
print, 'QLF calculated...'
close, 11
close, 12
close, 13



set_plot, 'X'
close, /all

end
