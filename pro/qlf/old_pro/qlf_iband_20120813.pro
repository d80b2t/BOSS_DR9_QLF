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
print, '    red, omega0=0.30, omegalambda=0.70, h100=0.70    '
print
print
red, omega0=0.30, omegalambda=0.70, h100=0.70
print
print

readin = 'y'
read, readin, prompt='Are the variables already in memory?? y/n '
print
print

pipe_or_person = 1
;read, pipe_or_person, prompt=' (0) Pipeline   or  (1) Person redshifts??  ' 
print
print

;completeness_level = 0.85
completeness_level = 82

;read, completeness_level, PROMPT=' - Whats the completeness level, 0.85 (DR9)  82 (S82) ?? '


if readin eq 'n' then begin
   
   ;;
   ;;  S D S S    D R 3    D A T A    (Richards et al. 2006)
   ;; 
   readcol, '../../data/Richards06_Table05.dat', $
            name_dr3, z_dr3, imag_dr3, Mi_dr3, del_gi_dr3, corr_dr3, $
            format='(a,f,f,f,f,f)', /silent
   print
   print, '  SDSS DR3 (aka Table 5 from Richards et al. 2006) read-in', n_elements(z_dr3)
   print
   print
   name_dr3_RA_hr  = strmid(name_dr3, 5, 2)                                    
   name_dr3_RA_min = strmid(name_dr3, 7, 2)                                    
   name_dr3_RA_sec = strmid(name_dr3, 9, 5)                                    

   name_dr3_Dec_hr  = strmid(name_dr3, 14, 3)                                    
   name_dr3_Dec_min = strmid(name_dr3, 17, 2)                                    
   name_dr3_Dec_sec = strmid(name_dr3, 19, 5)                                    

   dr3_RA  = ((((name_dr3_RA_sec/60.)+name_dr3_RA_min)/60.) +name_dr3_RA_hr )*15.
   dr3_Dec = (((name_dr3_Dec_sec/60.)+name_dr3_Dec_min)/60.)+name_dr3_Dec_hr
   
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
   ;;  B O S S   Q U A S A R  spectroscopic information
   ;;

   ;; We either want:
   ;;   all the BOSS Stripe 82 data, w/o any completeness corrections, 
   ;; or 
   ;;   the "XDCORE-SPALL" file 
   ;;       
   ;; BOSS21 is now dealt with separately (in  qlf_iband_boss21.pro)
   ;;
   
   ;; 
   ;; These datasets should be made/sculpted some place else
   ;; e.g. $BOSSQSOMASK_DIR/data/
   qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
;   xdspallfile = qsomaskdir+'/compfiles/xdspallmask_DR9.fits'  ;; 
   xdspallfile = qsomaskdir+'/compfiles/xdspallmask.fits'
   data = mrdfits(xdspallfile,1)
 
   ;; The completeness_level is now set above...
   print
   print, 'completeness_level', completeness_level
   print
   
   if pipe_or_person ne 1 then begin
      w_boss = where(data.poly_weight ge completeness_level and $
                     data.zWarning eq 0 and data.z ge 2.00, N_boss_qsos)
   ENDIF ELSE BEGIN 
      w_boss = where(data.poly_weight ge completeness_level and $
                     data.z_conf_person ge 3 and data.z_person ge 2.20 and data.z_person le 3.50, N_boss_qsos)
   ENDELSE

   ;; slightly cheeky, but will save a whole load of faff down the line...
   if pipe_or_person eq 1 then  data.z =  data.z_person
   
   
   ;; Dec 2011 note: Here if completeness_level =0.700 and z ge 2.00, N should = 31,822
   ;; Feb 2012 note: Here if completeness_level =0.850 and z_pipe   ge 2.00, N should = 24,348
   ;; Feb 2012 note: Here if completeness_level =0.850 and z_person ge 2.00, N should = 24,348

   ;;  O L D !!  Trying to reproduce June results...
;   data = mrdfits('../../../Stripe82/spAll-v5_4_14_S82.fits', 1)
;   w_boss = where( data.zWarning eq 0 and $
 ;                  (strtrim(data.OBJTYPE) eq 'QSO' or $
  ;                  strtrim(data.class)   eq 'QSO'), N_boss_qsos_wDups)

   ;;
   ;; THIS LINE SETS THE WHOLE THING UP!!
   ;;
   ;; For the full XDQSO sample...
   if completeness_level le 1.00 then boss = data[w_boss] 

   ;; 
   ;; For BOSS STRIPE 82 
   ;; 
   ;; Note this has been updated from a dataset that had 
   ;; 6250 z_pipe>2.00 quasars to a file that had 
   ;; 6239 z_pipe>2.00quasars to the current incarnation of the dataset that has
   ;; 6380 z_person>2.00 quasars. 
   ;;
;   boss = mrdfits('../../../bossqsomask/trunk/data/spall/mini-spAll-v5_4_45_Stripe82_z2QSOs.fits',1)
;   qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
   stripe82file = qsomaskdir+'spall/mini-spAll-v5_4_45_Stripe82_z2QSOs.fits'
   if completeness_level eq 82 then begin
      data = mrdfits(stripe82file, 1)
      ;; print, 2LL^10+2LL^11+2LL^12+2LL^13+2LL^14+2LL^15+2LL^16+2LL^17+2LL^18+2LL^19 +2LL^40+2LL^41 +2LL^42+2LL^43+2ll^44
      target_flag = 34084861508608LL
      ;; Arguably the most inclusive, most sensible one...
      target_flag_ancil = 2LL^0+2LL^1+2ll^2+2LL^3+2LL^4+2LL^5      +2LL^7+2LL^8+2LL^9   

      ;; 3.50 gives you 5476 NQs...
      ;; 4.00 gives you that last z-bin on the 10-panel plot...
      top_zbin = 4.00
      
      s82_qso_zpconf_zW23 = where(data.ra ge 300.0 or data.ra lt 60.0 and $ 
                                    (data.dec le 1.25 and data.dec ge -1.25) $
                                    and ( ((data.boss_target1 and target_flag) ne 0) $
                                          or ((data.ancillary_target2 and target_flag_ancil) ne 0) ) $ 
                                    and   data.specprimary eq 1 and data.z_conf_person ge 3 and $
                                    data.z_person ge 2.20 and data.z_person le top_zbin, N_s82_qso_zpconf_zW23)
      
      if pipe_or_person eq 1 then  data.z =  data.z_person

      boss = data[s82_qso_zpconf_zW23]
   endif

   N_boss_qsos = n_elements(boss)

   print
   print
   print
   print, '  ***************************************************'
   print, '  ***                                             ***'
   print, '  ***  Completeness Level set at  ', completeness_level
   print, '  ***                                             ***'
   print, '  ***************************************************'
   print
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
;; BOSS DR9 COMPLETENESS CORRECTIONS FROM McGREER ET AL. (2012...)
;;
;readcol, '../../completeness/grid_fid_iPQSOMIDZ_sel.txt', $
;readcol, '../../completeness/grid_newcont_iPQSOMIDZ_sel.txt', $
;readcol, '../../completeness/grid_dusty_iPQSOMIDZ_sel.txt', $
;selec_func_file = 'oldfiducial_grid_newlinesvrev_i_sel.txt'
;selec_func_file = 'oldfiducial_grid_newlinesv5revx1d_i_sel.txt'

selec_func_file = 'fiducial_linetweak_grid_i_sel.txt'
selec_func_path = '../../completeness/'+selec_func_file

readcol, selec_func_path,          i_start, i_end, z_start, z_end, selfn

print
print, '   Selection function READ-IN ', n_elements(selfn)
print
print 

;; Need to think about this really carefully!
;selfn_zero = where(selfn eq 0.00, N_selfn_zero)
selfn_zero = where(z_start ge 2.20 and z_end le 3.50 and selfn eq 0, N_selfn_zero)
selfn[selfn_zero] =1.000

plothist, selfn, $
          bin=0.02, /ylog, $
          charsize=2.0, charthick=2.2, thick=2.4, $
          xtitle='completeness', ytitle='no of bins', $
          title=selec_func_file

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
   dlums_boss   = DLUMINOSITY(boss.z)   / 1e6
   
;; Selec. fnc.
;;   dlums_selfn  = DLUMINOSITY(red_selfn) / 1e6
   dlums_z_start  = DLUMINOSITY(z_start) / 1e6
   dlums_z_end    = DLUMINOSITY(z_end) / 1e6

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

;; 
print
print, ' **** APPLYING REDDENING EXTICTION CORRECTION!!! **** '
print
print

imag_boss = boss.PSFMAG[3] - boss.EXTINCTION[3]

Abs_iM_dr3    = imag_dr3         - (5 * alog10(dlums_dr3))    - 25.00 - kcor(fix(z_dr3/0.01)) 
;Abs_iM_comp  = iflux            - (5 * alog10(dlums_zphot))  - 25.00 - kcor(fix(zphot/0.01)) 
;Abs_iM_boss   = boss.PSFMAG[3]   - (5 * alog10(dlums_boss))   - 25.00 - kcor(fix(boss.z/0.01)) 
Abs_iM_boss   = imag_boss        - (5 * alog10(dlums_boss))   - 25.00 - kcor(fix(boss.z/0.01)) 


;; Step up four Abs_iMs based on the combo of i_start/end and z_start end
;Abs_iM_selfn = iband_selfn    - (5 * alog10(dlums_selfn))  - 25.00 - kcor(fix(boss.z/0.01)) 
Abs_iM_istart_zstart = i_start  - (5 * alog10(dlums_z_start))  - 25.00 - kcor(fix(z_start/0.01)) 
Abs_iM_istart_zend   = i_start  - (5 * alog10(dlums_z_end))    - 25.00 - kcor(fix(z_end/0.01)) 
Abs_iM_iend_zstart   = i_end    - (5 * alog10(dlums_z_start))  - 25.00 - kcor(fix(z_start/0.01)) 
Abs_iM_iend_zend     = i_end    - (5 * alog10(dlums_z_end))    - 25.00 - kcor(fix(z_end/0.01)) 




print, 'Absolute i-band Mags calculated '
print
print

;plot, boss.z, Abs_iM_boss, $
 ;     ps=3, xrange=[1.8, 3.8], xstyle=1, yrange=[-23.0, -29.5], ystyle=1
;oplot, z_start, Abs_iM_istart_zstart,color=100, ps=3
;oplot, z_end, Abs_iM_iend_zend,color=220, ps=3


xthick    = 4.0 
ythick    = 4.0 
charsize  = 3.2
charthick = 6.0

;diff_in_Mi_dr3s = (Mi_dr3 - Abs_iM_dr3)
;plot, z_dr3, diff_in_Mi_dr3s, ps=3, xrange=[-0.2, 6.0], $
;      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
 ;     charsize=charsize/1.4, charthick=charthick/2.4, $
  ;    xtitle=' redshift', ytitle=' M_i (DR3)  -  M_i (NPR) '
   ;   

;; THINK ONE WANTS TO DO THE COMPLETENESS CORRECTION OUTSIDE
;; THE ``QLF LOOP'', and essentially on an object-by-object basis...
print, '  size(boss) ', size(boss)
print
print

;; "old values"... 
;delta_IM  = 0.0362500 ;; (Am assuming each bin is supposed to be 0.0725 delta_m wide..)
;delta_red = 0.025 

;; Newer values...
;delta_IM  = 0.0500 ;; (Am assuming each bin is supposed to be 0.100 delta_m wide..)
;delta_red = 0.025  ;; z_bins are/remain delta_z=0.05

count=0
counter=0
boss_wgt = fltarr( N_boss_qsos)
openw, 20, 'weightening_loop_temp.dat'
for ii=0ll, N_boss_qsos-1 do begin
   count=count+1
;   w_selfn = where((Abs_iM_selfn-delta_IM) lt ABS_IM_BOSS[ii] and $
 ;                  (Abs_iM_selfn+delta_IM) gt ABS_IM_BOSS[ii] and $
  ;                 (red_selfn+delta_red)   ge boss[ii].z      and $
   ;                (red_selfn-delta_red)   lt boss[ii].z, N_selfn)
   
   w_selfn = where(imag_boss[ii] ge i_start and $
                   imag_boss[ii] lt i_end and $
                   boss[ii].z         ge z_start and $
                   boss[ii].z         lt z_end, N_selfn)
   
;if N_selfn ne 1 then print, ' *********** ', ii, boss[ii].psfmag[3], boss[ii].z
   if N_selfn ne 1 then print, ' *********** ', ii, imag_boss[ii], boss[ii].z
   
   if boss[ii].z ge 2.20 then begin
      counter=counter+1
;     if N_selfn eq 0 then print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn
      if N_selfn eq 1 then begin
                                ; print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn, selfn[w_selfn]
         boss_wgt[ii] = 1./selfn[w_selfn]
;         boss_wgt[ii] = 1. ;; IF e.g. Stripe 82
      endif
;      if N_selfn gt 1 then print, ii, ABS_IM_BOSS[ii], boss[ii].z, N_selfn, selfn[ w_selfn]
   endif

if N_selfn ne 1 then printf, 20, ii, imag_boss[ii], boss[ii].z, ABS_IM_BOSS[ii], N_selfn, format='(i,d,d,d,i)'
if N_selfn eq 1 then printf, 20, ii, imag_boss[ii], boss[ii].z, ABS_IM_BOSS[ii], N_selfn, selfn[w_selfn], format='(i,d,d,d,i,d)'
endfor
close, 20

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

if completeness_level eq 0.000 then Area_boss = 3671.   /41253.
if completeness_level eq 0.750 then Area_boss = 2803.12 /41253.
if completeness_level eq 0.850 then Area_boss = 2236.38 /41253.
if completeness_level eq 0.950 then Area_boss = 1517.88 /41253.
if completeness_level eq 82    then Area_boss = 220.    /41253.

Area_Stp82 =  (220.0/41253.)
Area_dr3   =  (1622./41253.)
;; From Richards06 (p.4) ``The uncertainty on this area is of order 10
;; deg^2''
area_2slaq  = (191.1/41253.)  ;; now not used...

;; The number to be actually used!!
;area_boss = Area_boss95
;area_boss = Area_Stp82


;; From Richards06, Section 6.1
mag_bins = 0.30
red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

openw, 11, 'My_QLF_iband_dr3_temp.dat'
openw, 12, 'My_QLF_iband_boss_temp.dat'
openw, 13, 'My_QLF_iband_boss_wgt_temp.dat'
openw, 33, 'My_QLF_iband_boss_wgt_formatted4paper_temp.dat'

openw, 17, 'My_QLF_iband_boss_temp_and_wgt_temp.dat'

openw, 14, 'My_QLF_iband_boss_narrowZ_temp.dat'
openw, 15, 'My_QLF_iband_boss_narrowZ_wgt_temp.dat'
openw, 34, 'My_QLF_iband_boss_narrowZ_formatted4paper_temp.dat'
openw, 35, 'My_QLF_iband_boss_narrowZ_wgt_formatted4paper_temp.dat'

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
      Abs_iM_min = ((jj+1) * mag_bins) + (-30.90)
      Abs_iM_max = ( jj   *  mag_bins) + (-30.90)
      ;; To get the bins at 0.00, 0.25, 0.50, 0.75...
                                ;Abs_iM_min = ((jj+1)*0.25)  + (-30.125)
;      Abs_iM_max = (jj*0.25)      + (-30.125)
      Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
      
      w_qlfbin_dr3 = where(Abs_IM_dr3 le Abs_iM_min and $
                           Abs_IM_dr3 gt Abs_iM_max and $
                           z_dr3 ge redmin and z_dr3 lt redmax, N_QSOs_dr3)
      
      
      w_qlfbin_boss = where(Abs_iM_boss lt Abs_iM_min and $
                            Abs_iM_boss gt Abs_iM_max and $
                            boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)
      
      if N_QSOs_boss ne 0 then N_QSOs_boss_wgt = total(boss_wgt[w_qlfbin_boss])
      if N_QSOs_boss eq 0 then N_QSOs_boss_wgt = N_QSOs_boss
      
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
      ;;  F o r   B O S S   Q S Os
      ;;
;      if N_QSOs_boss le 0 then printf, 12, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss,     999.999
      if N_QSOs_boss gt 0 then begin
         mean_z      = mean(boss[w_qlfbin_boss].z)
         mean_Absmag = mean(Abs_iM_boss[w_qlfbin_boss])

         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)

         Num_den_sum         = float(N_QSOs_boss)     / ((V_comov_top-V_comov_bottom)*Area_boss)
         Num_den_sum_wgt     = float(N_QSOs_boss_wgt) / ((V_comov_top-V_comov_bottom)*Area_boss)

         Num_den_sum_per_mag     = Num_den_sum     * (1.0 /mag_bins)
         Num_den_sum_per_mag_Wgt = Num_den_sum_wgt * (1.0 /mag_bins)

         log_Num_den_sum_boss     = alog10(Num_den_sum_per_mag)
         log_Num_den_sum_boss_wgt = alog10(Num_den_sum_per_mag_wgt)

         error_boss     = (10^log_Num_den_sum_boss)     /sqrt(N_QSOs_boss)
         error_boss_wgt = (10^log_Num_den_sum_boss_wgt) /sqrt(N_QSOs_boss_wgt)

        ; Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
        ; z_bin = (redmin + redmax) /2.

       ;  if N_selfn ne 0 then print, ii,jj, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss, N_selfn, selfn[w_selfn]

         printf, 12, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss,     N_QSOs_boss, error_boss, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
         printf, 13, z_bin, Abs_mag_bin, Num_den_sum_wgt, log_Num_den_sum_boss_wgt, N_QSOs_boss_wgt, error_boss_wgt, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'

         ;printf, 33, mean_z, ' & ', Abs_mag_bin, ' & ', Num_den_sum_wgt, ' & ', $
         ;        log_Num_den_sum_boss_wgt, ' & ', N_QSOs_boss, ' & ', N_QSOs_boss_wgt, ' & ', error_boss_wgt, ' \\',$
         ;        format='(f9.5,a, f12.5,a, e16.6,a, f16.8,a, i8,a, i8,a, e16.6,a)'

         printf, 33, mean_z, ' & ', mean_AbsMag, ' & ', Abs_mag_bin, ' & ', $
                  N_QSOs_boss, ' & ', log_Num_den_sum_boss_wgt, ' & ', error_boss_wgt*1e9, ' \\',$
;                 format='(f9.5,a, f12.5,a, e16.6,a, f16.8,a, i8,a, e16.6,a)'
                 format='(f7.3,a, f9.3,a, f9.3,a, i7,a, f8.3,a, f8.3,a)'

         printf, 17, ii, jj, z_bin, Abs_mag_bin, N_QSOs_boss, N_QSOs_boss_wgt, Num_den_sum, Num_den_sum_wgt,  $
                 log_Num_den_sum_boss, log_Num_den_sum_boss_wgt, error_boss, error_boss_wgt, $
                 format='(i4,i4, f7.3, f10.3,  i8,i8,  e14.4,e14.4, f14.4,f14.4, e14.4,e14.4)'
      endif
   endfor
endfor
close, 11
close, 12
close, 13
close, 33
close, 17



mag_bins = 0.30
;red_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]
red_bins_narrow = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]

print
print, ' WORKING OUT THE QLF "NARROW Z"....'
print
for ii = 0L, 10 do begin 
   ;; Eeeek, argh...
   ;; R06 z-bins not evenly spaced, especially at high-z...
   redmin = red_bins_narrow[ii]        ;redmin  = ((ii)   * 0.38) - 0.08
   redmax = red_bins_narrow[ii+1]      ;redmax  = ((ii+1) * 0.38) - 0.08
   z_bin = (redmin + redmax) /2.
   
   for jj = 0L, 37 do begin 
      Abs_iM_min = ((jj+1)*mag_bins) + (-30.90)
      Abs_iM_max = (jj    *mag_bins) + (-30.90)
      Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
      
      w_qlfbin_boss = where(Abs_iM_boss lt Abs_iM_min and $
                            Abs_iM_boss gt Abs_iM_max and $
                            boss.z ge redmin and boss.z lt redmax, N_QSOs_boss)
      
      if N_QSOs_boss ne 0 then N_QSOs_boss_wgt = total(boss_wgt[w_qlfbin_boss])
      
      ;;
      ;; FOR BOSS QSOs
      ;;
      if N_QSOs_boss     le 0 then printf, 14, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss,     999.999
      if N_QSOs_boss_wgt le 0 then printf, 15, z_bin, Abs_mag_bin, -9999.999, -999.,     N_QSOs_boss_wgt, 999.999
      
      if N_QSOs_boss gt 0 then begin
         mean_z      = mean(boss[w_qlfbin_boss].z)
         mean_Absmag = mean(Abs_iM_boss[w_qlfbin_boss])
                  
         V_comov_top    =  VCOMOVING(redmax) / (1e6^3)
         V_comov_bottom =  VCOMOVING(redmin) / (1e6^3)
         
         Num_den_sum         = float(N_QSOs_boss)     / ((V_comov_top-V_comov_bottom)*Area_boss)
         Num_den_sum_wgt     = float(N_QSOs_boss_wgt) / ((V_comov_top-V_comov_bottom)*Area_boss)
         
         Num_den_sum_per_mag     = Num_den_sum     * (1.0 /mag_bins)
         Num_den_sum_per_mag_Wgt = Num_den_sum_wgt * (1.0 /mag_bins)
         
         log_Num_den_sum_boss     = alog10(Num_den_sum_per_mag)
         log_Num_den_sum_boss_wgt = alog10(Num_den_sum_per_mag_wgt)
         
         error_boss     = (10^log_Num_den_sum_boss)     / sqrt(N_QSOs_boss)
         error_boss_wgt = (10^log_Num_den_sum_boss_wgt) / sqrt(N_QSOs_boss_wgt)
         
        ; Abs_mag_bin = (Abs_iM_min + Abs_iM_max) /2.
        ; z_bin = (redmin + redmax) /2.

       ;  if N_selfn ne 0 then print, ii,jj, z_bin, Abs_mag_bin, Num_den_sum, log_Num_den_sum_boss, N_QSOs_boss, N_selfn, selfn[w_selfn]

         if ii eq 1 and N_QSOs_boss eq 1 then begin
            print,  w_qlfbin_boss
         endif
         
         ;; openw, 14, 'My_QLF_iband_boss_narrowZ_temp.dat'
         ;; openw, 15, 'My_QLF_iband_boss_narrowZ_wgt_temp.dat'
         ;; openw, 34, 'My_QLF_iband_boss_narrowZ_formatted4paper_temp.dat'
         ;; openw, 35, 'My_QLF_iband_boss_narrowZ_wgt_formatted4paper_temp.dat'
         
         printf, 14, z_bin, Abs_mag_bin, Num_den_sum,     log_Num_den_sum_boss,     N_QSOs_boss, error_boss, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'
         printf, 15, z_bin, Abs_mag_bin, Num_den_sum_wgt, log_Num_den_sum_boss_wgt, N_QSOs_boss_wgt, error_boss_wgt, $
                 format='(f9.5, f12.5, e16.6, f16.8, i8, e16.6)'

         printf, 34, mean_z, ' & ', mean_AbsMag, ' & ', Abs_mag_bin, ' & ', $
                  N_QSOs_boss, ' & ', log_Num_den_sum_boss, ' & ', error_boss*1e9, ' \\',$
;                 format='(f9.5,a, f12.5,a, e16.6,a, f16.8,a, i8,a, e16.6,a)'
                 format='(f7.3,a, f9.3,a, f9.3,a, i7,a, f8.3,a, f8.3,a)'

         printf, 35, mean_z, ' & ', mean_AbsMag, ' & ', Abs_mag_bin, ' & ', $
                  N_QSOs_boss, ' & ', log_Num_den_sum_boss_wgt, ' & ', error_boss_wgt*1e9, ' \\',$
;                 format='(f9.5,a, f12.5,a, e16.6,a, f16.8,a, i8,a, e16.6,a)'
                 format='(f7.3,a, f9.3,a, f9.3,a, i7,a, f8.3,a, f8.3,a)'


      endif

   endfor
endfor
close, 14
close, 15
close, 35

ww = [613, 4115, 678, 6175, 1110]


print, 'QLF calculated...'



print
print
print, '    No. of Quasars that went into calculation.... ', n_elements(boss)
print
print, '    Area that went into calculation....           ', area_boss * ((180./!dpi)^2.)*4.*!dpi
print
print 

set_plot, 'X'
close, /all

end
