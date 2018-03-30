
;+
;
;
;
;
;-

print
print
print, '     red, omega0=0.30, omegalambda=0.70, h100=0.70 '
print
print

readcol, '../../pro/models/PLE_fit_params_S82.dat', $
         redshift, Phi_Star, alpha, beta, chi_sq_min,	No_of_bins

;red, omega0=0.237, omegalambda=0.763, h100=0.73
;zz=[0.00, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]       
;tt=getage(zz)                                  
ttt = [13.7771  ,    8.84072    ,  6.11412     , 3.46345  ,    2.27406  ,    1.24499  ,   0.502578]

;;red, omega0=0.30, omegalambda=0.70, h100=0.70
zz=[-0.01, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]    
tt=getage(zz)   
tt=[13.7507, 8.42653, 5.75177, 3.22670, 2.11257, 1.15478, 0.465897]


; GETAGE(z)              - look-back time
; GETREDSHIFT(age)       - redshift for a given age

age = GETAGE(redshift)  

qso_redshift = [0.526,0.804,1.026,1.225,1.413,1.579, 1.745, 1.921, 2.131,2.475]
qso_lookback_time = getage(qso_redshift)
r_nought = qso_redshift * 6.0


;; 
;; SFR DENSITY DATA  FROM Bouwens_2010_ApJ_709_L133_Fig5
;;
readcol, 'Bouwens_2010_ApJ_709_L133_Fig5.dat', z_Bouwens, rho_SFR, rho_SFR_wCorr 
age_SFR = getage(z_Bouwens)


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
      readcol, '../../pro/models/Croom09b_PLE_temp.dat', $
               ii_PLE, z_PLE, jj_PLE, mag_PLE, k1_PLE, k2_PLE, alpha_PLE, Phi_PLE
   endif
   
   if choice_LDDE_points eq 'y' then begin
      readcol, '../../pro/models/Croom09b_mLDDE_temp.dat', $
               ii_mLDDE, z_mLDDE, jj_mLDDE, mag_mLDDE, e_d_mLDDE, zc_mLDDE, Phi_mLDDE
   endif
   
   if choice_LEDE_points eq 'y' then begin
      readcol, '../../pro/models/Croom09b_LEDE_temp.dat', $
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
k1 = 1.46
k2 = -0.328

;w_PLE_mag_bins = where(ii_PLE eq 0, N_PLE_mag_bins)
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
age_PLE_model = getage(z_bin_for_max_Phi)


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

charsize=2.2
charthick=6.8
thick=4.8
xthick=4
ythick=4
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; x-ranges
xmin = 13.7
xmax = 0.0

;; y-ranges
ymin = -10.5
ymax =  3.0


set_plot, 'ps'
device, filename='evolution_PLE_params_temp.eps', $
        xsize=12.0, ysize=6.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 1.2, /fill
plot, age, Phi_Star, $
      psym=8, $
;      position=[0.20, 0.20, 0.94, 0.90], $
      position=[0.12, 0.20, 0.96, 0.88], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      xstyle=8+1, ystyle=3, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $
      xtitle='Gyr ', $
      ytitle=' value ', $
      color=black

;;
;;
;; Just overplotting the star-formation rate of THE UNIVERSE!!!
;;
;;
oplot, age_SFR, rho_SFR,      thick=8,  color=light_blue
oplot, age_SFR, rho_SFR_wCorr-3.5, thick=8, color=yellow


oplot, age, (chi_sq_min/(No_of_bins-3.)), ps=8, color=green
oplot, age, beta,     ps=8, color=blue
oplot, age, alpha   , ps=8, color=red
oplot, age, Phi_Star, ps=8, color=black

oplot, age_PLE_model, alog10(max_Phi_PLE), ps=2, color=orange
oplot, age_PLE_model, MSTAR+16,            ps=8, color=turquiose


charsize_label  = charsize
charthick_label = charthick

xyouts, 12.30,  2.00, '!4v!E2!N!3/!4m!3', charsize=charsize_label, charthick=charthick_label, color=green
xyouts, 12.30,  1.00, 'faint-end slope',  charsize=charsize_label, charthick=charthick_label, color=blue
xyouts, 12.30,  0.00, 'brigh-end slope',  charsize=charsize_label, charthick=charthick_label, color=red
xyouts, 12.30, -1.00, '!4u!E*!3',         charsize=charsize_label, charthick=charthick_label, color=black
xyouts, 12.30, -2.00, 'max(!4u!E*!N!3!IPLE!N)',charsize=charsize_label, charthick=charthick_label, color=orange
xyouts, 12.30, -3.00, '!3M!E*!N!3+16.',     charsize=charsize_label, charthick=charthick_label, color=turquiose

xyouts, 7.0,  0.25, '!N!7q!I!3SFR!N',     charsize=charsize_label*1.4, charthick=charthick_label/1.2, color=light_blue
xyouts, 7.0, -1.00, '!N!7q!I!3SFR!N',     charsize=charsize_label*1.4, charthick=charthick_label/1.2, color=yellow
xyouts, 6.0, -1.10, 'w/Corr -3.5',     charsize=charsize_label/1.2, charthick=charthick_label/1.2, color=yellow



;
;zz=[0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
;zlab=['0','','1','2','3','5','10']
;tlab=['13.8', ' ', '6.1', '3.5', '2.3', '1.3', '0.5']
;;ttlab=['0','0.5','1','2','3','5','10']

;print, getage(zz)
 ;     13.7507      8.42653      5.75177      3.22670      2.11257      1.15478     0.465897
;zlab = getage(zz)
;zlab=['0.0', '0.5', '1.0', '2.0', '3.0', '5.0', '10.0' ]
;; zz= [-0.0200000     0.500000      1.00000      2.00000      3.00000      5.00000      10.0000]
zlab=['10', '5', '3', '2', '1.0', '0.50', '0.00', '0.00' ]

print
print, 'zz', zz

nticks=n_elements(zz)
tt = getage(zz)
zz = ['0.0', '0.5', '1.0', '2', '3', '5', '10']

print
print
print, 'zz', zz
print, 'zlab', zlab
print, 'tt', tt
print, 'nticks', nticks
print
print

;; zz = -0.0100000     0.500000      1.00000      2.00000      3.00000      5.00000      10.0000
;; gives *where* on the top axis the labels should go. 
;; needs to be converted to age, and flipped around. 

zz = ['0.0', '0.5', '1.0', '2', '3', '5', '10']

axis, xaxis=1,       xtickv=tt, $
      xtickname=zz, $
      xticks=nticks-1, $
      xthick=4.2, charthick=charthick, charsize=charsize, $ 
;      XTICKFORMAT='(F5.2)',  $
      xtitle='!8z!3, redshift' 
;      xtitle='!8 Age of the Universe (Gyr)',charsize=2


device, /close
set_plot, 'X' 


end    
