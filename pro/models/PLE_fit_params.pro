;+
;
;
;
;
;-


;;
;;
;;  F I T S    T O    T H E    D A T A 
;;
;;
readcol, 'PLE_fit_params_S82.dat', redshift, Phi_Star, alpha, beta, chi_sq_min,	No_of_bins

;; 
;; SFR DENSITY DATA  FROM Bouwens_2010_ApJ_709_L133_Fig5
;;
readcol, '../../plot/evolution/Bouwens_2010_ApJ_709_L133_Fig5.dat', $
         z_Bouwens, rho_SFR, rho_SFR_wCorr 
;age_SFR = getage(z_Bouwens)


;;
;;
;;  M o d e l   f i t s...
;;
;;
print
choice_models = 'y'
;read, choice_models, PROMPT=' - Plot model lines?? y/n  '

choice_PLE_points    = 'y'
choice_LDDE_points   = 'n'
choice_LEDE_points   = 'n'
if (choice_models eq 'y') then begin
;   read, choice_PLE_points,  PROMPT=' - Plot PLE  model??  y/n  '
 ;  read, choice_LDDE_points, PROMPT=' - Plot LDDE model??  y/n  '
  ; read, choice_LEDE_points, PROMPT=' - Plot LEDE model??  y/n  '
   
   if choice_PLE_points eq 'y' then begin
      readcol, 'Croom09b_PLE_temp.dat', $
               ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   endif
   if choice_LDDE_points eq 'y' then begin
      readcol, 'Croom09b_mLDDE_temp.dat', $
               ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE
   endif
   if choice_LEDE_points eq 'y' then begin
      readcol, 'Croom09b_LEDE_temp.dat', $
               ii_LEDE, z_LEDE, jj_LEDE, mag_LEDE, e_d_LEDE, zc_LEDE, alpha, Phi_LEDE
   endif
endif

print
print
print, 'choice_models      ', choice_models
print, 'choice_PLE_points  ', choice_PLE_points    
print, 'choice_LDDE_points ', choice_LDDE_points   
print, 'choice_LEDE_points ', choice_LEDE_points  
print
print

Mstar0_g = -22.17
k1       = 1.46
k2       = -0.328

N_PLE_mag_bins = max(jj_PLE)+1
N_PLE_z_bins   = max(ii_PLE)+1

max_Phi_PLE       = fltarr(N_PLE_z_bins)
z_bin_for_max_Phi = fltarr(N_PLE_z_bins)
Mstar             = fltarr(N_PLE_z_bins)

for ii=0LL, N_PLE_z_bins-1 do begin
   w = where(ii eq ii_PLE, N)
   max_Phi_PLE[ii]       = max(Phi_PLE[w])
   z_bin_for_max_Phi[ii] = z_PLE[ii*N_PLE_mag_bins]

   Mstar[ii]  = Mstar0_g  - 2.5*(k1*z_PLE[ii*N_PLE_mag_bins] + k2*z_PLE[ii*N_PLE_mag_bins]*z_PLE[ii*N_PLE_mag_bins])
   print, ii, N,  max_Phi_PLE[ii], z_bin_for_max_Phi[ii], Mstar[ii]

endfor




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

charsize=3.0
charthick=6.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; x-ranges
xmin = 2.0   ;0.5
xmax = 4.0   

;; y-ranges
ymin = -7.5   ;;-9.5
ymax =  6.0   ;; 8.4


set_plot, 'ps'
device, filename='PLE_fit_params_temp.eps', $
        xsize=8.0, ysize=14.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 2, /fill
plot, redshift, Phi_Star, $
      psym=8, $
      position=[0.20, 0.10, 0.94, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $
      xtitle=' redshift ', $
      ytitle=' value ', $
      color=black

;; Just overplotting the star-formation rate of THE UNIVERSE!!!
;oplot, z_Bouwens, rho_SFR,       thick=8,  color=light_blue
;oplot, z_Bouwens, rho_SFR_wCorr, thick=8, color=yellow


oplot, redshift, (chi_sq_min/(No_of_bins-3.)), ps=8, color=green
oplot, redshift, beta,     ps=8, color=blue
oplot, redshift, alpha   , ps=8, color=red
oplot, redshift, Phi_Star, ps=8, color=black

plotsym, 0, 1.2, /fill
;oplot, z_bin_for_max_Phi, alog10(max_Phi_PLE), ps=2, color=orange
;oplot, z_bin_for_max_Phi, alog10(max_Phi_PLE)-1.5, ps=2, color=orange

;oplot, z_bin_for_max_Phi, MSTAR+18,            ps=8, color=turquiose


charsize_label  = charsize*1.4
charthick_label = charthick*1.4

xyouts_x = 2.30                      ;; 2.30 if xmin=2.00 ; 0.70 if xmin=0.50
;xyouts_y = [7.40, 6.40, 5.40, 4.40] ;; if y_minmax = -9.5, 8.4
xyouts_y = [5.10, 4.30, 3.50, 2.70]
xyouts, xyouts_x, xyouts_y[0], '!4v!E2!N!3/!4m!3',      charsize=charsize_label, charthick=charthick_label, color=green
xyouts, xyouts_x, xyouts_y[1], 'faint-end slope',  charsize=charsize_label, charthick=charthick_label, color=blue
xyouts, xyouts_x, xyouts_y[2], 'brigh-end slope',  charsize=charsize_label, charthick=charthick_label, color=red
xyouts, xyouts_x, xyouts_y[3], '!4u!E*!3',         charsize=charsize_label, charthick=charthick_label, color=black

;xyouts, xyouts_x, 3.40, 'max(!4u!E*!N!3!IPLE!N)-1.5',charsize=charsize_label, charthick=charthick_label, color=orange
;xyouts, xyouts_x, 2.40, '!3M!E*!N!3+18',     charsize=charsize_label, charthick=charthick_label, color=turquiose



device, /close
set_plot, 'X' 


end    
