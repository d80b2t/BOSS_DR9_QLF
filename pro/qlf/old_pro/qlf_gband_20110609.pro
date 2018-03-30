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
   ;; Reading in the Richards06 k-correction
   ;; Normalized at z=2. 
   readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
   print, 'Richards06_kcor.dat READ-IN', n_elements(kcor)
   print
   print
   
   ;; **Relatively** good first guess...
   ;; Since I'm claiming (g-i)~0.2 for 2.2<z<3.0 
   ;; e.g. Fig. 11 Croom et al. (2009a)
   kcor = kcor + 0.2
   
   ;;
   ;;  B O S S   Q U A S A R  spectroscopic information
   ;;
   data = mrdfits('../../../Stripe82/spAll-v5_4_14_S82.fits', 1)
   ;data = mrdfits('../../../Stripe82/spAll-v5_4_35_anti-S82.fits', 1)
   ;data = mrdfits('../../../Stripe82/spAll-v5_4_35_boss12_CORE.fits',1)
   
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

;;
;; BOSS COMPLETENESS CORRECTIONS FROM McGREER ET AL. (2011...)
;;
readcol, '../../completeness/sims_200QSO_5000LOS_modelA/modelA_200QSO_5000LOS_gPQSO_sel.txt', $
         g_start, g_end, z_start, z_end, selfn
print
print, ' Selection function READ-IN ', n_elements(selfn)
print
print 
red_selfn   = 0.5*(z_start+z_end)
gband_selfn = 0.5*(g_start+g_end)


;; Running the DLUMINOSITY command from the red.pro routine
;; Returns LUMINOSITY DISTANCES (using the above cosmology)
;; Divide by 1e6 to get into Mpc
if readin eq 'n' then begin
   print, 'Doing DLUMS....'
   dlums_boss  = DLUMINOSITY(boss.z)    / 1e6
   dlums_selfn = DLUMINOSITY(red_selfn) / 1e6

   print, 'DLUMS calculated '
   print
   print
endif

;; Simple k-correction just for the time being...
;;
; k_corr = -0.5  ;; generally not used...
;; Setting up the power-law k-correction for the K-bands
;; e.g. Richards et al. (2006)
;alpha_nu     = -0.5 
;redshift_arr = (findgen(61)/10.)
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
;Abs_iM_dr3    = imag_dr3       - (5 * alog10(dlums_dr3))   - 25.00 - kcor(fix(z_dr3/0.01)) 
;Abs_iM_zphot  = iflux          - (5 * alog10(dlums_zphot)) - 25.00 - kcor(fix(zphot/0.01)) 
Abs_Mg_boss  = boss.PSFMAG[1] - (5 * alog10(dlums_boss))  - 25.00 - kcor(fix(boss.z/0.01)) 
Abs_Mg_selfn = gband_selfn    - (5 * alog10(dlums_selfn)) - 25.00 ;;- kcor(fix(boss.z/0.01)) 

print, 'Absolute g-band Mags calculated '
print
print


delta_Mg  = 0.05
delta_red = 0.025
count=0
boss_wgt = fltarr(N_boss_qsos)
for ii=0l, N_boss_qsos-1 do begin
   count=count+1
   w_selfn = where((Abs_Mg_selfn-delta_IM) lt ABS_Mg_BOSS[ii] and $
                   (Abs_Mg_selfn+delta_IM) gt ABS_Mg_BOSS[ii] and $
                   (red_selfn+delta_red)   ge boss[ii].z      and $
                   (red_selfn-delta_red)   lt boss[ii].z, N_selfn)
   
   if boss[ii].z ge 2.00 then begin
      count=count+1
;      if N_selfn eq 0 then print, ii, ABS_Mg_BOSS[ii], boss[ii].z, N_selfn
      if N_selfn eq 1 then begin
        ; print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn, selfn[w_selfn]
         boss_wgt[ii] = 1./selfn[w_selfn]
      endif
;      if N_selfn gt 1 then print, ii, ABS_Mg_BOSS[ii], boss[ii].z, N_selfn, selfn[w_selfn]
   endif
endfor





;; Area of Stripe 82 
;; Full Sky                 = 41,253 deg^2 = 1
;; Area Stripe 82 for BOSS  = 220.0  deg^2 = 0.00533295
;; 

Area_Stp82 =  (220.0/41253.)
Area_dr3   =  (1622./41253.)
;; THIS 2nd NUMBER STILL NEEDS TO BE PRECISION CHECKED!!!

;; Since we split the data up into 0.30mag
;; bins for counting...
mag_bins = 0.30

;; From Croom09 and Richards09 Section 6.1
red_bins = [0.40, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

openw, 11, 'My_QLF_gband_dr3_temp.dat'
openw, 12, 'My_QLF_gband_boss_temp.dat'
openw, 13, 'My_QLF_gband_boss_wgt_temp.dat'
;print, 'ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs, mean_z, V_comov, Num_den, log_Num_den'

bright_i    = 15.00
faint_i_loz = 19.10
faint_i_hiz = 20.20

for ii = 0L, 10 do begin 
   ;; Eeeek, argh...
   ;; R06 z-bins not evenly spaced, especially at high-z...
   redmin = red_bins[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 40 do begin 
      Abs_Mg_min = ((jj+1)*0.30) + (-30.90)
      Abs_Mg_max = (jj*0.30)     + (-30.90)
      
      w_qlfbin_boss = where(Abs_Mg_boss lt Abs_Mg_min and $
                            Abs_Mg_boss gt Abs_Mg_max $
                            and boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)
      
      if N_QSOs_boss ne 0 then N_QSOs_boss_wgt = total(boss_wgt[w_qlfbin_boss])

;      print, ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs_dr3;, mean_z, V_comov, Num_den, log_Num_den
;      print, ii, jj, redmin, redmax, Abs_iM_min, Abs_iM_max, N_QSOs_dr3;, mean_z, V_comov, Num_den, log_Num_den
      
      ;;
      ;; FOR BOSS QSOs
      ;;
      if N_QSOs_boss gt 0 then begin
         mean_z  = mean(boss[w_qlfbin_boss].z)

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum         = float(N_QSOs_boss) / ((V_comov_top-V_comov_bottom)*Area_Stp82)
         Num_den_sum_wgt     = float(N_QSOs_boss_wgt) / ((V_comov_top-V_comov_bottom)*Area_Stp82)

         Num_den_sum_per_mag = Num_den_sum * (1.0 /mag_bins)
         Num_den_sum_per_mag_Wgt = Num_den_sum_wgt * (1.0 /mag_bins)
         
         log_Num_den_sum_boss = alog10(Num_den_sum_per_mag)
         log_Num_den_sum_boss_wgt = alog10(Num_den_sum_per_mag_wgt)

         Abs_mag_bin = (Abs_Mg_min + Abs_Mg_max) /2.
         z_bin = mean_z 
         printf, 12, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss
         printf, 13, z_bin, Abs_mag_bin, Num_den_sum_wgt, log_Num_den_sum_boss_wgt, N_QSOs_boss_wgt
      endif

   endfor
endfor
close, 11
close, 12
close, 13

print, 'QLF calculated...'
print
print



set_plot, 'X'
close, /all

end
