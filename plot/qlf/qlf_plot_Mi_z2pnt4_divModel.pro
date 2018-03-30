;+
; NAME:
;    qlf_plot_Mi_z2pnt4
; 
; PURPOSE:
;    To plot Quasar Luminosity Functions
;    in M_i(z=~2.4)
;
; CALLING SEQUENCE:
;       .run qlf_plot_Mi_z2pnt4
;
; INPUTS:
;        Richards06_Table06.dat   - QLF from Richards06, z=2.01 bin only
;
; OUTPUTS:
;       .ps file
;
; COMMENTS:
;       /usr/common/rsi/lib/general/LibAstro/ 
;         
;-
print
print

red, omega0=0.30, omegalambda=0.70, h100=0.700

print
print
Mi2_offset = (-0.25)   ;; Conroy&White12 Eq. 6, by default
print, ' M_i(z=2)        = M_g(z=2) -  0.25  set as defaut....'
print, ' M_i(z=2) + 0.25 = M_g(z=2)          set as defaut....'
;read,  Mi2_offset, PROMPT=' - Value of  M_i(z=2) to M_g(z=2)  offset??  '
print


;; 
;;
;;  C H O I C E S ... 
;;
;;
print
print

choice_mydr3_points = 'n'
;read, choice_mydr3_points, PROMPT=' - Plot My DR3 points?? y/n  '

choice_BOSS_points   = 'y'
choice_S82_points    = 'n' 
choice_S82_errorbars = 'n'
choice_DR9_points    = 'n'
choice_DR9_errorbars = 'n'
choice_DR9_nowgt     = 'n'

;read, choice_BOSS_points, PROMPT=' - Plot BOSS points?? y/n  '
if choice_BOSS_points eq 'y' then begin
   choice_BOSS_line = 'y'
   ;read, choice_BOSS_line, PROMPT=' - Plot BOSS (data) line ?? y/n  '
   ;read, choice_S82_points, PROMPT=' - Plot S82 points?? y/n  '
   ;if choice_S82_points eq 'y' then read, choice_S82_errorbars, PROMPT=' - Plot S82 error bars?? y/n  '
   read, choice_DR9_points, PROMPT=' - Plot DR9 points?? y/n  '
   if choice_DR9_points eq 'y' then begin
     ;read, choice_DR9_nowgt,     PROMPT=' - Plot DR9 points with NO WGT?  y/n  '
      read, choice_DR9_errorbars, PROMPT=' - Plot DR9 error bars?? y/n  '
   endif
endif   

choice_how2plot = '2'
;read, choice_how2plot, PROMPT=' - Plot  (1) (log Phi_obs)/(log_Phi_model)  or  (2) log10(Phi/Phi)  ?  '


print
print
choice_fits = 'n'
read, choice_fits, PROMPT=' - Plot ``fits'' (e.g. PLE, LDDE, Crot09, HRH07) ?? y/n  '

choice_models_HRH07  = 'n'
if choice_fits eq 'y' then choice_models_HRH07  = 'y'
;read, choice_models_HRH07, PROMPT=' - Plot HRH07 models  ?? y/n  '
print
print

choice_models = 'n'
read, choice_models, PROMPT=' - Plot models (e.g. CW12, Milleniuum, Shen09, DeGraf12) ?? y/n  '
print
print


;; 
;;         D  A  T  A
;;
;;
;;  S D S S   D R 3
;;
readcol, '../../data/Richards06_Table06.dat', $
         z_R06, Mi_SDSS, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor
print
print, '../../data/Richards06_Table06.dat READ-IN', N_elements(z_R06)
print

R06_delta_up    = alog10((10^log_PhiR06)+(sigma_Phi*1e-9))
R06_delta_up = (log_PhiR06 -R06_delta_up)
R06_delta_down = alog10((10^log_PhiR06)-(sigma_Phi*1e-9))
R06_delta_down = (log_PhiR06 - R06_delta_down)

R06_err = alog10(sigma_Phi*1e-9)

;;
;;   2 S L A Q    Q S O
;;
readcol, '../../data/Croom09b_2SLAQ_QLF.dat', $
         Mg_2SLAQ, z_2SLAQ, NQ, log_phi_2SLAQ, log_phi_2SLAQ_down, log_phi_2SLAQ_up
;; Mi2_offset = (-0.25)
Mi_2SLAQ = Mg_2SLAQ + Mi2_offset 

;;
;;
;;   B O S S   D R 9
;;
;;
;;  N. B. ::  ALL   WITH correction/weighting (wgt)
;;
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_McGkcor.dat', $ 
;readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat', $   ;; used in the DR9 QLF paper!!
readcol, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor_rv1.dat', $   ;; 
;readcol, '../../pro/qlf/My_QLF_iband_boss_narrowZ_S82_McGkcor_rv1.dat', $ 
         z_bin_boss_wgt, Abs_mag_bin_boss, blah_boss_wgt, log_Phi_BOSS_wgt, raw_N_QSOs_boss_wgt, sigma_Phi_BOSS_wgt
print
print, '../../pro/qlf/My_QLF_iband_boss_wgt_McGkcor.dat  READ-IN', N_elements(z_bin_boss_wgt)
print

Mi_BOSS_wgt = Abs_mag_bin_boss 

boss_delta_up_wgt   = alog10((10^log_Phi_BOSS_wgt+sigma_Phi_BOSS_wgt))
boss_delta_up_wgt   = boss_delta_up_wgt - log_Phi_BOSS_wgt
boss_delta_down_wgt = alog10((10^log_Phi_BOSS_wgt-sigma_Phi_BOSS_wgt))
boss_delta_down_wgt = abs(boss_delta_down_wgt - log_Phi_BOSS_wgt)


;;
;;  B O S S   2 1
;;
;readcol, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_20120828.dat', $  
readcol, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_McGkcor.dat', $
         z_bin_boss21, Abs_mag_bin_boss21, blah_boss21, log_Phi_BOSS21, raw_N_QSOs_boss21, sigma_Phi_BOSS21
print
print, '../../pro/qlf/My_QLF_iband_boss21_wSelfn_20120305.dat  READ-IN', N_elements(z_bin_boss21)
print 

Mi_boss21 = Abs_mag_bin_boss21 

boss21_delta_up   = alog10((10^log_Phi_BOSS21 + sigma_Phi_BOSS21))
boss21_delta_up   = boss21_delta_up - log_Phi_BOSS21
boss21_delta_down = alog10((10^log_Phi_BOSS21 - sigma_Phi_BOSS21))
boss21_delta_down = abs(boss21_delta_down - log_Phi_BOSS21)


;;
;;   B O S S  + M M T 
;;
readcol, 'boss21MMT_z2pnt4.dat', $
         Mg_MMT, Nq_MMT, log_phi_MMT, sigma_phi_MMT				
print
print, 'BOSS+MMT data  READ-IN', N_elements(z_bin_boss_S82)
print
;; Mi2_offset = (-0.25)
Mi_MMT = Mg_MMT + Mi2_offset  


;;
;;   C O M B O  -  1 7 
;;
readcol, 'COMBO17_z2pnt4.dat', $
         M145_C17, MBband_C17, NQ_C17, log_Phi_C17, log_Phi_C17_npr
print
print, 'COMBO-17 data  READ-IN', N_elements(log_Phi_C17)
print
;; Again from e.g. Eqn. (6) of Conroy&White12.
Mi_C17 = MBband_C17 - 0.71

error_C17 = 10^(log_Phi_C17) / sqrt( NQ_C17) 
sigma_C17 =  log_Phi_C17 - alog10(error_C17)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;  M o d e l   f i t s...
;;
;;

;;
;;  From   C R O O M   e t  a l .   (2 0 0 9 b)
;;
if (choice_fits eq 'y') then begin
   readcol, '../../pro/models/Croom09b_PLE_toz5.dat', $
            ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   
   readcol, '../../pro/models/Croom09b_mLDDE_toz5.dat', $
            ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE

   readcol, '../../pro/models/Croom09b_LEDE_toz5.dat', $
            ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
   
   mag_PLE   = mag_PLE   + Mi2_offset
   mag_mLDDE = mag_mLDDE + Mi2_offset
   mag_LEDE  = mag_LEDE  + Mi2_offset   

   readcol, '../../pro/models/chi_sq_PLE_model_a_20120810.dat', $
            z_ple_model_a, Mi2_ple_model_a, phi_ple_model_a
   
   readcol, '../../pro/models/Croton09_qlf_model.dat', $
            z_Croton09, MbJ_Croton09, phi_Croton09
   Mi_Croton09 = MbJ_Croton09 - 0.71 
endif

if (choice_models eq 'y') then begin
   readcol, '../../pro/models/CW12/CW12_qsolf_z2.40.dat', $
            CW12_Mi2, log_Phi_CW12
   
   readcol, '../../pro/models/Fanidakis12/bjlf_z.2.2_nf12.dat', $
            Fanidakis_MbJ_2,  log_Fanidakis_Phi_total, Fanidakis_Phi 
   Fanidakis_Mi2 = Fanidakis_MbJ_2 - 0.71

   readcol, '../../pro/models/Shen09/model_A_z2.40.dat', $
             Log_LBol_Shen09, Mi2_Shen09, log_Phi_L, log_Phi_Shen09
   ;; density in units of Mpc^-3
   Phi_Shen09 = alog10(10^(log_Phi_Shen09) * (0.7^3))
   
   readcol, '../../pro/models/Marulli08/agnLF_best_2.42.dat', $
            L_bol_Maru08, log_phi_Maru08, log_phi_Maru08_min, log_Phi_Maru08_max

   ;; The conversion of Croom++ and Shen++ should work just fine
   ;;                               Mi (z = 2) = 90 − 2.5 log(Lbol /erg s−1 ),
   ;; For L_Q in Watts you have ::  M_i(z=2) = 72.5 - 2.5 log10(L_Q) 
   Mi2_Maru08  = 72.5 - 2.5*( (alog10(3.9e26)+L_bol_Maru08))
   Mi2_Maru08  = 72.5 - (2.5* alog10((10^L_bol_Maru08)*3.9d26))

endif


;;
;;  From   H  R  H   e t  a l .   (2 0 0 7) 
;;
if (choice_models_HRH07 eq 'y') then begin
   ;; 1 : observed luminosity in the band, in (log_{10}(L [erg/s]))
   ;; 2 : corresponding absolute monochromatic AB magnitude at the given frequency nu. 
   ;; 3 : corresponding observed monochromatic flux in (log_{10}(S_nu [milliJanskys])) 
   ;; 4 : corresponding bolometric luminosity (given the median BC
   ;;     corrections as a func. of Luminosity) 
   ;; 5 : comoving number density per unit log_{10}(luminosity) : 
   ;;     (make sure to correct by the appropriate factor to convert to 
   ;;     e.g. the number density per unit magnitude) 

   readcol, '/cos_pc19a_npr/BOSS/QLF/pro/models/HRH07/HRH07_Bband_z2.40.dat', $
            band_lum_HRH07, ABMag_HRH07, Monoflux_HRH07, bol_lum_HRH07, phi_HRH07_z2pnt4

   ;; Converting to Mi(z=2), assuming M_B ~ M_bJ
   Mi_HRH07 = ABMag_HRH07 - 0.71
   band_lum_sol = (10d^band_lum_HRH07)/(3.9e33)
   Bmag_HRH07   =  4.82 + ((-2.5)* alog10(band_lum_sol))
endif





;; 
;; Setting a sanity-check completeness limit. 
;; 
;; Assume (g-r)=0.11 and (g-i)=0.21 for 2.2<z<2.4 QSOs.
;;
boss_glimit = 22.00
boss_rlimit = 21.85
boss_ilimit = 21.80

red_bins = [0.49, 0.87, 1.25, 1.63, 2.01, 2.40, 2.80, 3.25, 3.75, 4.25, 4.75]

;; For BOSS in M_1450 at z=1.25, 2.40, 3.24, 4.25]
M_limit=[-21.413717, -23.154132, -23.948839,-24.642529] 

Mg_limit = M_limit 


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

charsize  = 2.6
charthick = 6.8
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.04
YTICKLEN  = 0.06

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Actually plotting stuff...
;;
;;
;;  Q L F  i-band  
;;   Richards+06 as a template...
;;

;; x-ranges
x_min = -16.001
;x_min = -18.250
x_max =  -31.50    ; -30.50

;; y-ranges
y_min = -1.00  ; -2.25
y_max =  1.00  ; 1.20  ;; Used to be -5.00 for the z=2.6 plot so take care here!!

plot_sym_size_R06  = 1.8
plot_sym_size_BOSS = 2.2

set_plot, 'ps'
!p.multi=0
device, filename='QLF_Mi_z2pnt4_divmodel_temp.ps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

boss_color=60
;;S82_color=160  ;; set below...

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  z  ~  2.4
;;
plot,  Mi_SDSS, log_PhiR06, $
       position=[0.22, 0.20, 0.96, 0.96], $
       xstyle=1, $
       ystyle=1, $
       xrange=[x_min, x_max], $
       yrange=[y_min, y_max], $
       xthick=xthick, $
       ythick=ythick, $
       charsize=charsize, $
       charthick=charthick, $
;       ytickformat='(a1)', $
       XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
       /nodata, $
       xtitle='!3M!Ii!N(z=2)', $
;       ytitle='log!I10!N !7U!3!Iobs!N / log!I10!N !7U!3!Imodel!N ', $
       ytitle='log!I10!N (!7U!3!Iobs!N/!7U!3!Imodel!N) '


loadct, clr_table
;;
;;   S  D  S  S    R 0 6
;;
w  =  where(z_R06 gt 2.20 and z_R06 lt 2.60) 
plotsym, 8, plot_sym_size_R06/1.2, /fill
;oplot,  Mi_SDSS[w], log_PhiR06[w], psym=8
Mi_SDSS_err = Mi_SDSS[w] - Mi_SDSS[w]

;divmodel    = Hopkins_full_mag(2.40, Mi_SDSS[w])
;divmodel     = best_fit_PLE(2.40, Mi_SDSS[w])
divmodel     = best_fit_linear_Phi(2.40, Mi_SDSS[w])

log_divmodel = alog10(divmodel)

oplot,  Mi_SDSS[w],  alog10((10^log_PhiR06[w]) /divmodel), $
        psym=8., color=black

ratio     =  alog10((10^log_PhiR06[w]) /divmodel)
diff_up   =  ratio - alog10((10^(log_phir06[w] + R06_delta_down[w])) / divmodel)
diff_down =  ratio - alog10((10^(log_phir06[w] + R06_delta_up[w])) / divmodel)

oploterror, Mi_SDSS[w], ratio, Mi_SDSS_err, diff_up,   $
            /hibar, errthick=12, psym=8, color=black, ERRcolor=black
oploterror, Mi_SDSS[w], ratio, Mi_SDSS_err, diff_down, $
            /lobar, errthick=12, psym=8, color=black, ERRcolor=black



;;
;;   2 S L A Q    Q S Os
;;
;;        Mg_2SLAQ, z_2SLAQ, NQ, log_phi_2SLAQ, log_phi_2SLAQ_down, log_phi_2SLAQ_up
w = where(z_2SLAQ ge 2.2 and z_2SLAQ lt 2.6, N)
loadct, clr_table
if choice_models eq 'y' or choice_fits eq 'y' then loadct, 0
color_2slaq = light_blue+32
sym_2slaq = 4

plotsym, sym_2slaq, plot_sym_size_R06*1.35, /fill
oplot,       Mi_2SLAQ[w], log_phi_2SLAQ[w], psym=8, color=black

;divmodel    = Hopkins_full_mag(2.40, Mi_2SLAQ[w])
;divmodel    = best_fit_PLE(2.40, Mi_2SLAQ[w])
divmodel     = best_fit_linear_Phi(2.40, Mi_2SLAQ[w])

log_divmodel = alog10(divmodel)

;oplot, Mi_2SLAQ[w], (log_phi_2SLAQ[w]/log_divmodel), psym=8, color=color_2slaq
oplot, Mi_2SLAQ[w], alog10(10^(log_phi_2SLAQ[w])/divmodel), psym=8, color=color_2slaq

ratio     =  alog10((10^log_Phi_2SLAQ[w]) /divmodel)
diff_up   =  ratio - alog10((10^(log_Phi_2SLAQ[w] +  log_phi_2SLAQ_up[w])) / divmodel)
diff_down =  ratio - alog10((10^(log_Phi_2SLAQ[w] +  log_phi_2SLAQ_down[w])) / divmodel)

plotsym, sym_2slaq, plot_sym_size_R06, /fill
Mi_2SLAQ_err = Mi_2SLAQ[w] - Mi_2SLAQ[w]
oploterror, Mi_2SLAQ[w], ratio, Mi_2SLAQ_err,  diff_up, $
            /hibar, errthick=12, psym=8, color=color_2slaq, ERRcolor=color_2slaq
oploterror, Mi_2SLAQ[w], ratio, Mi_2SLAQ_err, diff_down, $
            /lobar, errthick=12, psym=8, color=color_2slaq, errcolor=color_2slaq



;;
;;  B O S S + M M T 
;; 
;;   
PLOTSYM, 0, 1.8, /fill
loadct, clr_table
if choice_models eq 'y' or choice_fits eq 'y' then loadct, 0
oplot, Mi_MMT, log_phi_MMT, psym=8,  color = blue

;divmodel    = Hopkins_full_mag(2.40, Mi_MMT)
;divmodel     = best_fit_PLE(2.40, Mi_MMT)
divmodel     = best_fit_linear_Phi(2.40, Mi_MMT)

log_divmodel = alog10(divmodel)

;oplot, Mi_boss21[w21], (log_Phi_MMT/log_divmodel), color=blue, ps=8
oplot, Mi_MMT, alog10((10^log_Phi_MMT)/divmodel), color=blue, ps=8

ratio   =  alog10((10^log_Phi_MMT) /divmodel)
diff    =  ratio - alog10((10^(log_Phi_MMT +  sigma_Phi_MMT)) / divmodel)


Mi_MMT_err = Mi_MMT - Mi_MMT
oploterror, Mi_MMT, ratio, Mi_MMT_err, diff,  $
            /hibar, errthick=12, psym=8, color=blue, ERRcolor=blue
oploterror, Mi_MMT, ratio, Mi_MMT_err, diff, $
            /lobar, errthick=12, psym=8, color=blue, errcolor=blue

;;
;;  C O M B O - 1 7
;; 
c17_sym = 5
loadct, clr_table
clr_c17 = orange
if choice_models eq 'y' or choice_fits eq 'y' then loadct, 0
;PLOTSYM, c17_sym, 3.4, thick=10
PLOTSYM, c17_sym, 3.4, /fill
;oplot, Mi_C17, log_phi_C17, psym=8, color = orange

;divmodel    = Hopkins_full_mag(2.40, Mi_C17)
;divmodel     = best_fit_PLE(2.40, Mi_C17)
divmodel     = best_fit_linear_Phi(2.40, Mi_C17)

log_divmodel = alog10(divmodel)

PLOTSYM, c17_sym, 3.4, /fill
if choice_how2plot eq '1' then oplot, Mi_C17, (log_phi_C17/log_divmodel), psym=8, color = clr_c17
;if choice_how2plot eq '2' then 
;oplot, Mi_C17,    alog10((10^log_Phi_C17) /divmodel), psym=8, color = clr_c17

ratio  =  alog10((10^log_Phi_C17) /divmodel)
diff   =  ratio - alog10((10^(log_Phi_C17[w] + sigma_C17)) / divmodel)

Mi_C17_err = Mi_C17 - Mi_C17
;oploterror, Mi_C17, ratio, Mi_C17_err, diff,  $
 ;           /hibar, errthick=8, psym=8, color=clr_c17, ERRcolor=clr_c17
;oploterror, Mi_C17, ratio, Mi_C17_err, diff, $
 ;           /lobar, errthick=8, psym=8, color=clr_c17, errcolor=clr_c17


;;
;;    B O S S   D R 9
;;
plotsym, 0, plot_sym_size_BOSS, /fill
boss_color = red
if choice_models eq 'y' or choice_fits eq 'y' then begin
   loadct, 0
   boss_color = 60
endif

plotsym, 0, plot_sym_size_BOSS, /fill
if (choice_DR9_points eq 'y') then begin
   www = where(z_bin_boss_wgt gt 2.20 and z_bin_boss_wgt lt 2.60 and $
;   www = where(z_bin_boss_wgt gt 2.30 and z_bin_boss_wgt lt 2.40 and $
               log_Phi_boss_wgt gt -50 and Mi_BOSS_wgt le -24.2, N) 
   
;   divmodel = Hopkins_full_mag(2.40, Mi_BOSS_wgt[www])
   ;divmodel = best_fit_PLE(2.40, Mi_BOSS_wgt[www])
   divmodel = best_fit_linear_Phi(2.40, Mi_BOSS_wgt[www])

   oplot,  Mi_BOSS_wgt[www],  alog10((10^log_Phi_boss_wgt[www]) /divmodel), $
           psym=8., color=black
   
   ratio     = alog10((10^log_Phi_boss_wgt[www]) /divmodel)
   diff_up   = ratio - alog10((10^(log_Phi_boss_wgt[www] + boss_delta_up_wgt[www])) / divmodel)
   diff_down = ratio - alog10((10^(log_Phi_boss_wgt[www] + boss_delta_down_wgt[www])) / divmodel)
   
   if (choice_DR9_errorbars eq 'y') then begin
      Mi_boss_err = Mi_BOSS_wgt[www] - Mi_BOSS_wgt[www]
      oploterror, Mi_BOSS_wgt[www], ratio, Mi_BOSS_err, diff_up,   $
                  /hibar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
      oploterror, Mi_BOSS_wgt[www], ratio, Mi_BOSS_err, diff_down, $
                  /lobar, errthick=12, psym=8, color=boss_color, ERRcolor=boss_color
   endif 
endif



;;
;;  L A B E L S 
;;
if (choice_models eq 'n' and choice_fits eq 'n') then begin

x_xyouts = -16.55
y_xyouts =  -0.10  ;1.20 ;-0.50  ;; 1.00  ;;; or something aroudn -0.50???
x_off    =   1.3
y_off    =   0.12
offset   =   0.14

;loadct, 6
plotsym, 0, plot_sym_size_BOSS/1.2, /fill
legend, '', position=[x_xyouts, y_xyouts], box=0, psym=8, color=red, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-x_off, y_xyouts-y_off, 'BOSS DR9', $
        charsize=charsize, charthick =charthick*1.2, color=red

plotsym, 8, plot_sym_size_R06/1.2, /fill
legend, '', position=[x_xyouts, y_xyouts-offset], box=0, psym=8, color=black, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-x_off, y_xyouts-y_off-offset, 'SDSS', $
        charsize=charsize, charthick =charthick*1.2, color=black

;loadct, clr_table
;color_2slaq = light_blue+32
plotsym, sym_2slaq, plot_sym_size_R06, /fill
legend, '', position=[x_xyouts, y_xyouts-(2*offset)], box=0, psym=8, color=color_2slaq, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-x_off, y_xyouts-y_off-(2*offset), '2SLAQ', $
        charsize=charsize, charthick =charthick*1.2, color=color_2slaq

loadct, clr_table
plotsym, 0, 1.6, /fill
legend, '', position=[x_xyouts, y_xyouts-(3*offset)], box=0, psym=8, color=blue, $
        charsize=charsize/1.4, charthick =charthick*1.2
xyouts, x_xyouts-x_off, y_xyouts-y_off-(3*offset), 'boss21+MMT', $
        charsize=charsize, charthick =charthick*1.2, color=blue

loadct, clr_table
plotsym, 5, 1.6, /fill
;legend, '', position=[x_xyouts, y_xyouts-(4*offset)], box=0, psym=8, color=orange, $
;        charsize=charsize/1.4, charthick =charthick*1.2
;xyouts, x_xyouts-1.2, y_xyouts-0.25-(4*offset), 'COMBO-17', $
;        charsize=charsize, charthick =charthick*1.2, color=orange

endif

if (choice_fits eq 'y') then begin
   x_xyouts = -19.30    ;; to match the z~2.4 plot...
   y_xyouts = -6.40     ;; -6.20 if plotting 6 models;  -6.40 otherwise...
   x_off    = 0.50
   y_off    = 0.35     ;; 0.25
   offset = 0.60

   loadct, clr_table
   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(1*offset)], box=0, line=0, $
           color=purple, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(1*offset), 'PLE', $
           charsize=charsize, charthick =charthick*1.2, color=purple
   
   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(2*offset)], box=0, line=1, $
           color=light_blue, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(2*offset), 'Croton09', $
           charsize=charsize, charthick =charthick*1.2, color=light_blue

   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(3*offset)], box=0, line=3, $
           color=red-8, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(3*offset), 'mLDDE', $
           charsize=charsize, charthick =charthick*1.2, color=red-8

   loadct, clr_table
   if (choice_models_HRH07 eq 'y') then begin
      legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(4*offset)], box=0, line=0, $
              color=TURQUIOSE, thick =thick*2.4, charsize=1.2
      xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(4*offset), 'HRH07', $
              charsize=charsize, charthick =charthick*1.2, color=TURQUIOSE
   endif
   ;; loadct, 8 and color=light_blue produces that nice dark green...
endif


if (choice_models eq 'y') then begin
   x_xyouts = -18.75   ;; to match the z~2.4 plot...
   y_xyouts = -0.35     ;; to match the z~2.4 plot...
   x_off    = 0.50
   y_off    = 0.18     ;; 0.25
   offset = 0.24       ;; 0.40

   loadct, clr_table
   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(1*offset)], box=0, line=2, $
           color=red-8, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(1*offset), 'CW12', $
           charsize=charsize, charthick =charthick*1.2, color=red-8

   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(2*offset)], box=0, line=1, $
           color=light_blue, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off-(2*offset), 'Shen09', $
           charsize=charsize, charthick =charthick*1.2, color=light_blue

   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(3*offset)], box=0, line=0, $
           color=TURQUIOSE, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off -(3*offset), 'Marulli08', $
           charsize=charsize, charthick =charthick*1.2, color=TURQUIOSE

   legend, '' , position=[x_xyouts+offset, y_xyouts-0.05-(4*offset)], box=0, line=3, $
           color=purple, thick =thick*2.4, charsize=1.2
   xyouts, x_xyouts-(x_off*5), y_xyouts-y_off -(4*offset), 'Fanidakis12', $
           charsize=charsize, charthick =charthick*1.2, color=purple
endif

;;
;;   z  ~  2 . 4 
;;
x_xyouts = -18.25   ;; -26.75
y_xyouts =   0.75   ;;  -5.50
offset   =   0.00
xyouts, x_xyouts, y_xyouts+(offset), '!8z~2.4!3', charsize=2.8, charthick=8, color=black


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   F I T S     AND     M   O   D   E   L  S
;;
;;
if choice_fits eq 'y' then begin
   loadct, clr_table
   
   w_PLE   = where(z_PLE_model_a   gt 2.38 and z_ple_model_a lt 2.42, N_PLE)
   oplot, MI2_PLE_MODEL_A[w_PLE], alog10(Phi_PLE_model_a[w_PLE]), thick=12, color=purple, linestyle = 0
   
   w_Croton09 = where(z_Croton09 gt 2.38 and z_Croton09 lt 2.42, N_Croton09)
   oplot, Mi_Croton09[w_Croton09], alog10(phi_Croton09[w_Croton09]), thick=12, color=light_blue, linestyle = 1
   
   per_mag = alog10(2.5)
   oplot, Mi_HRH07, alog10(phi_HRH07_z2pnt4)-per_mag, thick=8, linestyle=0,color=turquiose
   
   loadct, 8
   ;w_PLE   = where(z_PLE   gt 2.35 and z_ple lt   2.45, N_PLE)
   ;oplot, mag_PLE[w_PLE], alog10(Phi_PLE[w_PLE]), thick=12, color=purple, linestyle = 0
   
   w_mLDDE = where(z_mLDDE gt 2.35 and z_mLDDE lt 2.45, N_mLDDE)
   oplot, mag_mLDDE[w_mLDDE], alog10(Phi_mLDDE[w_mLDDE]), thick=12, color=light_blue, linestyle = 3
   
    ;w_LEDE = where(z_LEDE   gt 2.35 and z_LEDE  lt 2.45, N_LEDE)
    ;oplot, mag_LEDE[w_LEDE], alog10(Phi_LEDE[w_LEDE]), thick=12, color=160, linestyle = 5
endif

if choice_models eq 'y' then begin
   loadct, clr_table
   thick=12
   ;; Conray & White (2012)
   oplot, CW12_Mi2, log_Phi_CW12, thick=thick, color=red-8, linestyle = 2
   
   ;; Shen (2009)
   oplot,  Mi2_Shen09, log_Phi_Shen09, thick=thick, color=light_blue, linestyle=1
   ;oplot,  Mi2_Shen09,     Phi_Shen09, thick=thick, color=light_blue, linestyle=2

   oplot,  Mi2_Maru08, log_phi_Maru08-1.0 , thick=thick,color=turquiose, linestyle=0

   ;; Fanidakis et al. (2012)
   oplot, Fanidakis_Mi2, Fanidakis_Phi, thick=thick, color=purple, linestyle=3

endif




loadct, 0
!p.multi=0

device, /close
set_plot, 'X'

close, /all

end
