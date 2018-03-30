;+
;
;    Just really scared to break the regular qlf_iband code!!
;
;- 

;;
;; Setting up a flat, H0 = 71.9 kms^-1 Mpc^-1 cosmology
;; DLUMINOSITY will not work unless this is set-up...
;;
;; red, omega0=0.27, omegalambda=0.73, h100=0.719
;;
;; Cosmology from R06:
print
print
print, '    red, omega0=0.30, omegalambda=0.70, h100=0.70    '
print
print
print
red, omega0=0.30, omegalambda=0.70, h100=0.70    


;;
;; Reading in the Richards06 k-correction
;; Normalized at z=2. 
;readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
readcol, 'kcorr_Miz2_mcgreer.dat', kcor_redshift, kcor, /silent
;print, '  Richards06_kcor.dat READ-IN', n_elements(kcor)
print, '  kcorr_Miz2_mcgreer.dat READ-IN', n_elements(kcor)
print
print

;;
;;  D A T A 
;;
data = mrdfits('/cos_pc19a_npr/BOSS/QLF/data/VAR_BOSS_boss21.fits',1)
print
print, ' (Full) boss21 READ-IN, ', n_elements(data)
print
print 

;; There are a whole bunch of TYPES from the Visual inspections, (see
;; VAR_BOSS.pro in /cos_pc19a_npr/BOSS/QLF/data). Here I'm going to be
;; liberal and include the QSO_?'s but will also check to see
;; what the effect is of not including them. 
;;
;; 
;N_QSO             1259 minmax(zscan[qso])            0.28100000       3.9700000
;N_QSO_BAL           60 minmax(zscan[qso_bal])         1.2830000       3.2640000;
;
;N_QSO_ques           106 minmax(zscan[qso_ques])        0.46000000       4.0290000
;
;N_star_ques          279 minmax(zscan[star_ques])       -1.0000000       7.0030000
;N_ques               102 minmax(zscan[ques])            -1.0000000       6.2800000
;N_Bad                315 minmax(zscan[bad])             -1.0000000       0.0000000
;N_star              1252 minmax(zscan[star])             0.0000000       0.0000000
;N_gal                 60 minmax(zscan[gal])            0.060000000      0.99300000

qsos_in = where(strtrim(data.TYPE) eq 'QSO'     or $
                strtrim(data.TYPE) eq 'QSO_BAL' or $
                strtrim(data.TYPE) eq 'QSO_?',   N_QSOs_in)                                         
boss21 = data[qsos_in]

qsos_in_ge1 = where(strtrim(data.TYPE) eq 'QSO'     or $
                    strtrim(data.TYPE) eq 'QSO_BAL' or $
                    strtrim(data.TYPE) eq 'QSO_?' and $
                    data.zscan ge 1.0, N_QSOs_in_ge1) 

qsos_in_ge1le2 = where(strtrim(data.TYPE) eq 'QSO'     or $
                       strtrim(data.TYPE) eq 'QSO_BAL' or $
                       strtrim(data.TYPE) eq 'QSO_?' and $
                       data.zscan ge 1.0 and data.zscan le 2.2, $
                       N_QSOs_in_ge1le2) 
print
print
print, 'No. of QSOs with TYPE eq QSO, QSO_BAL or QSO_?                   ', N_qsos_in,     ' with minmax(zscan[qsos_in]) ', minmax(data[qsos_in].zscan)
print, 'No. of QSOs with TYPE eq QSO, QSO_BAL or QSO_? and z>1           ', N_qsos_in_ge1, ' with minmax(zscan[qsos_in]) ', minmax(data[qsos_in_ge1].zscan)
print, 'No. of QSOs with TYPE eq QSO, QSO_BAL or QSO_? and z>1 and z<2.2 ', N_qsos_in_ge1le2, ' with minmax(zscan[qsos_in]) ', minmax(data[qsos_in_ge1le2].zscan)

print
print


;;
;; S E L E C T I O N    F U N C T I O N
;;
readcol, '/cos_pc19a_npr/BOSS/QLF/completeness/SelectionFunctionCompletness.txt', $
         gmag_selfn, selfn


;; CALCULATION THE D_LUMs...
dlums_boss21 = DLUMINOSITY(boss21.zscan) / 1e6

;; CALCULATION THE ABS MAGs...
Abs_iM_boss21 = boss21.PSFMAG[3] - (5 * alog10(dlums_boss21)) - 25.00 - kcor(fix(boss21.zscan/0.01)) 

print
print, 'Absolute i-band Mags calculated '
print
print, 'minmax(Abs_iM_boss21)  ', minmax(Abs_iM_boss21)
print 
print
print

w_large_AbsMag = where(Abs_iM_boss21 lt -30., N_large_AbsMag)

if N_large_AbsMag gt 0 then begin
   print
   print, 'minmax(Abs_iM_boss21[w_large_AbsMag]) ', minmax(Abs_iM_boss21[w_large_AbsMag])
   print 
   
   openw, 10, 'spAll_Quasars_noPSFMAGs_temp.dat'
   printf, 10, '# PLATE MJD FIBERID, RA, DEC, PLUG_RA, PLUG_DEC '
   for ii=0ll, N_large_AbsMag-1 do begin
      printf, 10, boss21[w_large_AbsMag[ii]].plate,  $
              boss21[w_large_AbsMag[ii]].mjd,  $
              boss21[w_large_AbsMag[ii]].fiberid, $
              boss21[w_large_AbsMag[ii]].ra, $
              boss21[w_large_AbsMag[ii]].dec, $
              boss21[w_large_AbsMag[ii]].plug_ra, $
              boss21[w_large_AbsMag[ii]].plug_dec, $
              format='(i5,2x,i6,2x,i4, f12.5, f12.5, 1x, f13.7, f13.7)'
   endfor
   close, 10
endif

w_sensible_AbsMag = where(Abs_iM_boss21 ge -30., N_sensible_AbsMag)


x_min = 0.00   ;; or min(boss21.zscan)
x_max = max(boss21.zscan)+0.5
   
y_min=max(Abs_iM_boss21[w_sensible_AbsMag])
y_max=min(Abs_iM_boss21[w_sensible_AbsMag])

plot, boss21.zscan, Abs_iM_boss21, $
      ps=3, $
      xrange=[x_min, x_max], xstyle=1, $
      yrange=[y_min, y_max], ystyle=1, $
      charsize=2.0, charthick=2.0, $
      xtitle='z, redshift', $
      ytitle='Absolute Magnitude, i-band'


count   = 0
counter = 0
boss21_wgt = fltarr(N_QSOS_IN)
delta_mag = (GMAG_SELFN[1] -GMAG_SELFN[0])/2.
for ii=0ll, N_QSOS_IN-1 do begin
   
   w_selfn = where(boss21[ii].psfmag[1] ge GMAG_SELFN-delta_mag and $
                   boss21[ii].psfmag[1] lt GMAG_SELFN+delta_mag, N_selfn)
   
   if N_selfn ne 1 then begin
;      print, ' *********** ', ii, boss21[ii].psfmag[1], boss21[ii].z, boss21[ii].type, format='(a,i6,d,d,2x,a)'
      counter = counter+1
   endif
   
   if N_selfn eq 1 then begin
     ; print, ii, ABS_IM_BOSS21[ii], boss21[ii].z, N_selfn, selfn[w_selfn]
      boss21_wgt[ii] = 1./selfn[w_selfn]
      count=count+1
   endif
endfor
print
print, 'No of objects that were too bright (g<18.625) or too faint (g>23.3750) for GMAG_SELFN ', counter
print, 'Total No. of objects used in the QLF calcs...', count
print
print

Area_boss21 = 14.5    /41253.

;; From Richards06, Section 6.1
mag_bins = 0.30
red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]


openw, 17, 'My_QLF_iband_boss21_wSelfn_temp.dat'
openw, 77, 'My_QLF_iband_boss21_wSelfn_formatted4paper_temp.dat'

print
print, ' WORKING OUT THE QLF....'
print
total_N_QSOs_boss21 =0
for ii = 0L, 10 do begin 
   redmin = red_bins[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 37 do begin 
      Abs_iM_min = ((jj+1) * mag_bins) + (-30.90)
      Abs_iM_max = ( jj   *  mag_bins) + (-30.90)
      Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
      
      w_qlfbin_boss21 = where(Abs_iM_boss21 lt Abs_iM_min and $
                              Abs_iM_boss21 gt Abs_iM_max and $
;                              boss21.z ge redmin and boss21.z lt redmax, N_QSOs_boss21)
                              boss21.zscan ge redmin and boss21.zscan lt redmax, N_QSOs_boss21)

      ;; Note to NPR: Using boss21.z     gives rise to My_QLF_iband_boss21_wSelfn_20120305.dat
      ;; Note to NPR: Using boss21.zscan gives rise to My_QLF_iband_boss21_wSelfn_20120828.day


      if N_QSOs_boss21 ne 0 then N_QSOs_boss21_wgt = total(boss21_wgt[w_qlfbin_boss21])
      
      ;;
      ;;  F o r   B O S S   2 1    Q S Os
      ;;
      ;if N_QSOs_boss21 le 0 then printf, 17, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss,     999.999
      if N_QSOs_boss21 gt 0 then begin
         mean_z      = mean(boss21[w_qlfbin_boss21].z)
         mean_AbsMag = mean(Abs_iM_boss21[w_qlfbin_boss21])
         
         V_comov_top    = VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom = VCOMOVING(redmin) / (1e6^3)
         
         Num_den_sum         = float(N_QSOs_boss21)   / ((V_comov_top-V_comov_bottom)*Area_boss21)
         Num_den_sum_per_mag = Num_den_sum     * (1.0 /mag_bins)
         
         log_Num_den_sum_boss21 = alog10(Num_den_sum_per_mag)
         
         error_boss21 = (10^log_Num_den_sum_boss21) /sqrt(N_QSOs_boss21)
         
         printf, 17, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss21, N_QSOs_boss21, error_boss21, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'

         printf, 77, mean_z,  ' & ', mean_AbsMag, ' & ', Abs_mag_bin,  ' & ', $
                 N_QSOs_boss21,  ' & ', log_Num_den_sum_boss21, ' & ', error_boss21*1e9, ' \\',  $
                 format='(f7.3,a, f9.3,a, f9.3,a, i7,a, f8.3,a, f8.3,a)'

         total_N_QSOs_boss21 = N_QSOs_boss21 + total_N_QSOs_boss21

      endif
   endfor
endfor
close, 17
close, 77

print, 'QLF calculated...'

print
print
print, '    No. of Quasars that went into calculation.... ', n_elements(boss21)
print
print, '    Area that went into calculation....           ', area_boss21 * ((180./!dpi)^2.)*4.*!dpi
print
print 

set_plot, 'X'
close, /all


end

