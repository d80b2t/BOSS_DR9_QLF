;+
;
;   H o w   w r o n g   i s   H R H 0 7,   **really...***
;
;-

plot_qlf_points = 'n' 

print
print
read, PROMPT='  Plot the Ross13, McGreer13, Willott10 data??  y/n  ', plot_qlf_points
print
print, ' plot_qlf_points   ', plot_qlf_points
print

;; Actual data used in e.g. Fig. 15 of Ross et al. (2013)
readcol, 's82_zbin_fit.txt', z_points, $
         log_phi, log_Phi_err_down, log_Phi_err_up, $
         Mstar,   Mstar_err_down,   Mstar_err_up, $
         alpha,   alpha_err_down,   alpha_err_up, $
         beta,    beta_err_down,    beta_err_up


;;
;;   M O D E L      P L E       z = 0.68 - 2.60 
;;   Line 1, Table 7, Palanque-Delabrouille et al. (2013) 
;;
z_bin_max   = 2.20 
z_bin_min   = 0.00 ;; The fit is only down to z=0.30
z_bin_width = 0.01
n_zfull     = ((z_bin_max - z_bin_min ) / z_bin_width) +1.
z_full      = (findgen(n_zfull)*z_bin_width)+z_bin_min

M_min       = -20.0
M_max       = -30.0
M_bin_width =   0.02
n_Mbin     = ((M_min - M_max) / M_bin_width) ;+1
M_bin_full  = (findgen(n_Mbin)*M_bin_width)+M_max

;; Fixed params...
alpha_PLE_fixed =  -1.41
 beta_PLE_fixed =  -3.33
     Mstar0_PLE =  -22.425  ;; Mg(z=0) = Mi(z=0) + 0.255 and Mg∗(0) = −22.17 from Croom++09, and PD++13. 
         k1_PLE =   0.02
         k2_PLE =  -0.33
   Phi_star_PLE =  10^(-5.84)
        z_pivot = 2.20
Mstar_pivot_PLE = -26.485 ;; Mg*(z_pivot) = -26.23
;; Want to cycle over "full" redshift range, in this case, 
;; z = 0.30 - 2.20
     z_PLE_PD = z_full
;; and want to cycle over wide Mag-range, say -20 to -30
     M_PLE_PD = M_BIN_FULL

       alpha_PLE_PD = fltarr(n_zfull)
        beta_PLE_PD = fltarr(n_zfull)
log_phi_star_PLE_PD = fltarr(n_zfull)
       Mstar_PLE_PD = fltarr(n_zfull)
       phi_PLE_PD13 = fltarr(n_zfull, n_Mbin)

for ii=0LL, n_zfull-1 do begin
       alpha_PLE_PD[ii] = alpha_PLE_fixed
        beta_PLE_PD[ii] = beta_PLE_fixed
log_phi_star_PLE_PD[ii] = alog10(Phi_star_PLE)
;       Mstar_PLE_PD[ii] = Mstar0_PLE - 2.5*(k1_PLE*z_PLE_PD[ii] + k2_PLE*z_PLE_PD[ii]*z_PLE_PD[ii])
       Mstar_PLE_PD[ii] = Mstar_pivot_PLE - 2.5*(k1_PLE*(z_PLE_PD[ii]-z_pivot) + k2_PLE*((z_PLE_PD[ii]-z_pivot)^2))
   for jj=0ll, n_Mbin-1 do begin
      dmag_PLE          = M_PLE_PD[jj] - Mstar_PLE_PD[ii]
      ple_denom_PLE     = (10.0^(0.4*(alpha_PLE_fixed+1)*(dmag_PLE))) + (10.0^(0.4*(beta_PLE_fixed+1)*(dmag_PLE))) 
      phi_PLE_PD13[ii,jj] = Phi_star_PLE / ple_denom_PLE
;;      print, ii, jj, z_PLE[ii],  M_PLE[jj], ple_denom_PLE, phi_PLE_R13[ii,jj]  
   endfor
endfor
print
print, '  PLE model from z=0.68 - 2.60 done (Line 1, Table 7,  Palanque-Delabrouille++13).... COMPUTED '
print


;;
;;   M O D E L      P L E       z = 0.30 - 2.20 
;;   Line 1, Table 8, Ross et al. (2013) 
;;
z_bin_max   = 2.20 
z_bin_min   = 0.00 ;; The fit is only down to z=0.30
z_bin_width = 0.01
n_zfull     = ((z_bin_max - z_bin_min ) / z_bin_width) +1.
z_full      = (findgen(n_zfull)*z_bin_width)+z_bin_min

M_min       = -20.0
M_max       = -30.0
M_bin_width =   0.02
n_Mbin     = ((M_min - M_max) / M_bin_width) ;+1
M_bin_full  = (findgen(n_Mbin)*M_bin_width)+M_max

;; Fixed params...
alpha_PLE_fixed =  -1.16
 beta_PLE_fixed =  -3.37
     Mstar0_PLE = -22.85
         k1_PLE =   1.241
         k2_PLE =  -0.249
   Phi_star_PLE =  10^(-5.96)
;; Want to cycle over "full" redshift range, in this case, 
;; z = 0.30 - 2.20
       z_PLE = z_full
;; and want to cycle over wide Mag-range, say -20 to -30
       M_PLE = M_BIN_FULL

      alpha_PLE  = fltarr(n_zfull)
       beta_PLE  = fltarr(n_zfull)
log_phi_star_PLE = fltarr(n_zfull)
       Mstar_PLE = fltarr(n_zfull)
     phi_PLE_R13 = fltarr(n_zfull, n_Mbin)

for ii=0LL, n_zfull-1 do begin
       alpha_PLE[ii] = alpha_PLE_fixed
        beta_PLE[ii] = beta_PLE_fixed
log_phi_star_PLE[ii] = -5.96
       Mstar_PLE[ii] = Mstar0_PLE - 2.5*(k1_PLE*z_PLE[ii] + k2_PLE*z_PLE[ii]*z_PLE[ii])
   for jj=0ll, n_Mbin-1 do begin
      dmag_PLE          = M_PLE[jj] - Mstar_PLE[ii]
      ple_denom_PLE     = (10.0^(0.4*(alpha_PLE_fixed+1)*(dmag_PLE))) + (10.0^(0.4*(beta_PLE_fixed+1)*(dmag_PLE))) 
      phi_PLE_R13[ii,jj] = Phi_star_PLE / ple_denom_PLE
;;      print, ii, jj, z_PLE[ii],  M_PLE[jj], ple_denom_PLE, phi_PLE_R13[ii,jj]  
   endfor
endfor
print
print, '  PLE model from z=0.30 - 2.20 done (Line 1, Table 8, Ross++13).... COMPUTED '
print

;;
;;   M O D E L     L E D E    F R O M    R O S S 1 3      z = 2.2 - 3.5
;;
z_bin_max   = 4.00 ;; Extending to z=6 (7) due to e.g. Fig. 19 of McGreer ++13, (Venemans discussion)
z_bin_min   = 2.20 ;; Extending to z=6 due to e.g. Fig. 19 of McGreer ++13
z_bin_width = 0.01
n_zfull     = ((z_bin_max - z_bin_min ) / z_bin_width) +1.
z_full      = (findgen(n_zfull)*z_bin_width)+z_bin_min

M_min       = -20.0
M_max       = -30.0
M_bin_width =   0.02
 n_Mbin     = ((M_min - M_max) / M_bin_width) ;+1
M_bin_full  = (findgen(n_Mbin)*M_bin_width)+M_max
 
;; ;;   Line 7, Table 8, Ross et al. (2013), fixed params...
   alpha_LEDE_fixed =  -1.29
    beta_LEDE_fixed =  -3.51
    Mstar2pnt2_LEDE = -26.57
            c1_LEDE =  -0.689
            c2_LEDE =  -0.809
Phi_star_2pnt2_LEDE =  10^(-5.93)
;; Want to cycle over "full" redshift range, in this case, 
       z_LEDE = z_full
;; and want to cycle over wide Mag-range, say -20 to -30
       M_LEDE = M_BIN_FULL

    alpha_LEDE = fltarr(n_zfull)
    beta_LEDE = fltarr(n_zfull)
       Mstar_LEDE = fltarr(n_zfull)
log_phi_star_LEDE = fltarr(n_zfull)
     phi_LEDE_R13 = fltarr(n_zfull, n_Mbin)

for ii=0LL, n_zfull-1 do begin
   alpha_LEDE[ii] = alpha_LEDE_fixed
    beta_LEDE[ii] = beta_LEDE_fixed

   ;; Eqn (11) of Ross13 or Eqn 2 of McGreer++13
   log_phi_star_LEDE[ii] = alog10(Phi_star_2pnt2_LEDE) + (c1_LEDE*( z_LEDE[ii]-2.2))
   ;; Eqn (12) of Ross13 or Eqn 2 of McGreer++13
   Mstar_LEDE[ii] = Mstar2pnt2_LEDE + (c2_LEDE*(z_LEDE[ii]-2.2))

   for jj=0ll, n_Mbin-1 do begin
      dmag_LEDE           = M_LEDE[jj] - Mstar_LEDE[ii]
      ple_denom_LEDE      = (10.0^(0.4*(alpha_LEDE_fixed+1)*(dmag_LEDE))) + (10.0^(0.4*(beta_LEDE_fixed+1)*(dmag_LEDE))) 
      phi_LEDE_R13[ii,jj] = 10*(log_phi_star_LEDE[ii]) / ple_denom_LEDE
;;      print, ii, jj, z_LEDE[ii],  M_LEDE[jj], ple_denom_LEDE, phi_LEDE_R13[ii,jj]  
   endfor
endfor
print
print, ' LEDE model from z= 2.20 - 4.00 done (Line 1, Table 8, Ross++13).... COMPUTED '
print


;;
;;   M O D E L    mod L E D E    F R O M    M c G R E E R 1 3   z =  4.0 - 6.0 (7.0)
;;
z_bin_max   = 7.00 ;; Extending to z=6 (7) due to e.g. Fig. 19 of McGreer ++13, (Venemans discussion)
z_bin_min   = 4.00 ;; Extending to z=6 due to e.g. Fig. 19 of McGreer ++13
z_bin_width = 0.01
n_zfull     = ((z_bin_max - z_bin_min ) / z_bin_width) +1.
z_full      = (findgen(n_zfull)*z_bin_width)+z_bin_min

M_min       = -20.0
M_max       = -30.0
M_bin_width =   0.02
 n_Mbin     = ((M_min - M_max) / M_bin_width) ;+1
M_bin_full  = (findgen(n_Mbin)*M_bin_width)+M_max
 
;; ;;   Line 3, Table 5, McGreer et al. (2013), fixed params...
   alpha_modLEDE_fixed =  -1.80
    beta_modLEDE_fixed =  -3.26
    Mstar2pnt2_modLEDE = -26.39 - 1.486   ;; M1450 to Mi
            c1_modLEDE =  -0.60
            c2_modLEDE =  -0.68
Phi_star_2pnt2_modLEDE =  10^(-5.93)
;; Want to cycle over "full" redshift range, in this case, 
       z_modLEDE = z_full
;; and want to cycle over wide Mag-range, say -20 to -30
       M_modLEDE = M_BIN_FULL

       alpha_modLEDE = fltarr(n_zfull)
        beta_modLEDE = fltarr(n_zfull)
       Mstar_modLEDE = fltarr(n_zfull)
log_phi_star_modLEDE = fltarr(n_zfull)
    phi_modLEDE_Mc13 = fltarr(n_zfull, n_Mbin)
    
for ii=0LL, n_zfull-1 do begin
   alpha_modLEDE[ii] = alpha_modLEDE_fixed
    beta_modLEDE[ii] = beta_modLEDE_fixed
   
   ;; Eqn (11) of Ross13 or Eqn 2 of McGreer++13
   log_phi_star_modLEDE[ii] = alog10(Phi_star_2pnt2_modLEDE) + (c1_modLEDE*( z_modLEDE[ii]-2.2))
   ;; Eqn (12) of Ross13 or Eqn 2 of McGreer++13
   Mstar_modLEDE[ii] = Mstar2pnt2_modLEDE + (c2_modLEDE*(z_modLEDE[ii]-2.2))
   ;; condition to keep  M*1450 > -27.0  i.e. M*i > - 28.486
   if Mstar_modLEDE[ii] lt -28.486 then Mstar_modLEDE[ii] = -28.486
   
   for jj=0ll, n_Mbin-1 do begin
      dmag_modLEDE            = M_modLEDE[jj] - Mstar_modLEDE[ii]
      ple_denom_modLEDE       = (10.0^(0.4*(alpha_modLEDE_fixed+1)*(dmag_modLEDE))) + (10.0^(0.4*(beta_modLEDE_fixed+1)*(dmag_modLEDE))) 
      phi_modLEDE_Mc13[ii,jj] = 10*(log_phi_star_modLEDE[ii]) / ple_denom_modLEDE
;;      print, ii, jj, z_LEDE[ii],  M_LEDE[jj], ple_denom_LEDE, phi_LEDE_R13[ii,jj]  
   endfor
endfor
print
print, ' modLEDE2 model from z= 4.00 - 7.00 done (Line 3, Table 5, McGreer++13).... COMPUTED '
print


;; 
;;  H R H 0 7     " F U L L " 
;;
z_bin_max   = 7.00 
z_bin_min   = 0.00
z_bin_width = 0.01
n_zfull     = ((z_bin_max - z_bin_min ) / z_bin_width) +1.
z_full      = (findgen(n_zfull)*z_bin_width)+z_bin_min

M_min       = -20.0
M_max       = -30.0
M_bin_width =   0.02
n_Mbin     = ((M_min - M_max) / M_bin_width) ;+1
M_bin_full  = (findgen(n_Mbin)*M_bin_width)+M_max

;; 
;; Line 6 ("Full"), Table 3, HRH07
;;     Px values from qlf_calculator.c PFH code!!  
;;
log_phi_star =  -4.8250643  ;; P0
log_Lstar0   =  13.035753   ;; P1, L_sol = 3.9x10^33 ergs s^-1 
       k_L_1 =   0.63150872 ;; P2
       k_L_2 = -11.763560   ;; P3
       k_L_3 = -14.249833   ;; P4
 gamma1_zero =   0.41698725 ;; P5
    k_gamma1 =  -0.62298947 ;; P6 
 gamma2_zero =   2.1744386  ;; P7  
  k_gamma2_1 =   1.4599393  ;; P8
  k_gamma2_2 =  -0.79280099 ;; P9
;;	P10=0.;P11=0.;P12=0.;P13=0.;P14=0.;}

;; Want to cycle over "full" redshift range, in this case, 
;; z = 0.30 - 2.20
  z_HRH07 = z_full
;; and want to cycle over wide Mag-range, say -20 to -30
  M_HRH07 = M_BIN_FULL
  
z_ref = 2.00
beta_min = 1.3
c1=6.24800 
k1=-0.370587
c2=8.99833 
k2=-0.0115970

Lum_sol = 3.839e33  ;; erg/s

 log_L_star      = fltarr(n_zfull)
     lband       = fltarr(n_zfull)
lband_ratio      = fltarr(n_zfull)
     Mstar_HRH07 = fltarr(n_zfull)
         gamma_1 = fltarr(n_zfull)
         gamma_2 = fltarr(n_zfull)
     alpha_HRH07 = fltarr(n_zfull)
      beta_HRH07 = fltarr(n_zfull)
log_phi_norm_HRH = fltarr(n_zfull)

for ii=0LL, n_zfull-1 do begin
   log_phi_norm_HRH[ii] = log_phi_star
   ;; Equation 10 from HRH07
   xsi = alog10((1. + z_HRH07[ii])/(1. + z_ref))
   ;; Equation 9 from HRH07
   log_L_star[ii]   = log_Lstar0 +  k_L_1*xsi +  k_L_2*xsi*xsi +  k_L_3*xsi*xsi*xsi 
   ;; faint-end slope, Equation 17 from HRH07
   gamma_1[ii]     = gamma1_zero  * (10^(xsi*k_gamma1))
   alpha_HRH07[ii] = (gamma_1[ii] +1.) *(-1.)
   ;; bright-end slope, Equation 19 from HRH07
   gamma_2[ii]    = (2.0 * gamma2_zero) * (   (10^(xsi* k_gamma2_1) + 10^(xsi* k_gamma2_2))^(-1.))
   beta_HRH07[ii] = (gamma_2[ii] +1.) *(-1.)
   if gamma_2[ii] lt beta_min then print, ii, gamma_2, beta_min
   ;; Equation 2 of HRH07 with 
   ;; (c1, k1, c2, k2) =  (6.25, -0.37, 9.00, -0.012) for
   lband_ratio[ii]= (c1 * (( 10^log_L_star[ii] / 1e10)^(k1))) + (c2 * (( 10^log_L_star[ii] / 1e10)^(k2)))
   ;; Convery L_bol (in log_L_sol) to L_Bband (also in log_Lsol)
   Mstar_HRH07[ii] = (-2.5) * (alog10(10^(log_L_star[ii]) / lband_ratio[ii]))
;   print, ii,  z_HRH07[ii], log_L_star[ii], Mstar_HRH07[ii], gamma_2[ii]
endfor

;plot, alog10(z_test_HRH+1), log_L_star, $
;plot, alog10(z_test_HRH+1), log_L_star, $
;      yrange=[11.3, 13.3], ystyl=1, $
;      xthick=xthick,     ythick=ythick, $
;      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
;      charsize=charsize, charthick=charthick,  thick=thick, $
;      xtitle='!6alog10(1+z_test)', ytitle='Log (L^*)'

print
;print, '  ``FULL'' model from z=0.00 - 7.00 (Line 6, Table 3, HRH07)....  COMPUTED '
print



;; Colour Table
;; http://ham.space.umn.edu/johnd/ct/ct-names.html
clr_table = 11
;; 11 is BLUE-RED
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
charthick = 4.8
thick     = 3.8
xthick    = 4.0
ythick    = 4.0
XTICKLEN  = 0.03
YTICKLEN  = 0.03


;; positions...
xpos_min = 0.20
xpos_max = 0.98
ypos_min = 0.20
ypos_max = 0.98

;; x-ranges
xmin = 0.0
xmax = 7.0

;; y-ranges
ymin = -19.
ymax = -28.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
;set_plot, 'ps'
;device, filename='temp.eps', $
;        xsize=8.0, ysize=8.0, $
;        xoffset=0.2, yoffset=0.2, $
;        /inches, /color, /encapsulated
;window, 1
plotsym, 0, 0.4, /fill
plot, z_PLE, Mstar_PLE, $
      psym=8, $
      position=[xpos_min, ypos_min, xpos_max, ypos_max], $
      xrange=[xmin, xmax], yrange=[ymin+5, ymax-5], $
      xstyle=1, ystyle=1, $
      xthick=xthick,     ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick,$ 
;      /nodata, $
      xtitle='!6!8z,!6 redshift', ytitle='!6M_star'
oplot, z_HRH07, Mstar_HRH07, thick=thick, linestyle=1



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='Mstar_redshift_temp.eps', $
        xsize=10.0, ysize=10.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

plotsym, 0, 0.4, /fill
plot, z_PLE, Mstar_PLE, $
      psym=8, $
      position=[xpos_min, ypos_min, xpos_max, ypos_max], $
      xrange=[xmin, xmax], yrange=[ymin, ymax], $
      xstyle=1, ystyle=1, $
      xthick=xthick,     ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick,$ 
;      /nodata, $
      xtitle='!6!8z,!6 redshift', ytitle='!6M_star'

oplot, z_HRH07, Mstar_HRH07+5., thick=thick, color=red
;xyouts, xmax*0.80, ymax*0.80, 'words',  charsize=charsize, charthick=charthick, color=color
;xyouts, xmax*0.80, ymax*0.80, No_of_words,  charsize=charsize, charthick=charthick, color=color

;legend, pos=[-22.2, -7.5], ' ',   box=0, thick=14, linestyle = 0, charsize=1.2
;xyouts,      -25.0, -7.8,  'PLE', charsize=2.2, charthick=8.

device, /close
set_plot, 'X'     


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; FULL N(z)'s
;;
set_plot, 'ps'
device, filename='HRH07_Fig8_redux_temp.eps', $
        xsize=14.0, ysize=8.5, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated
!p.multi=[0,3,2]
 ;; (.98-.14)/3 = .28  ==> .14, .42, .70, .98

;;
;; MIND YOURSELF!!
;z_test_HRH = alog10( z_test_HRH +1.)
;xmax = alog10(1+6.)
xmax = 7.0
linethick =    thick * 2.4
charsize  = charsize * 1.16

xposoff= 0.02
;;
;;  A L P H A    f a i n t - e n d    s l o p e...
;;
plot, z_HRH07, gamma_1, $
      position=[.08+xposoff, .60, .32+xposoff, .98], $
;      xrange=[xmin, xmax], yrange=[-0.3,1.15], $ ;; for the gamma_1 def
      xrange=[xmin, xmax], yrange=[-2.2,-0.40], $
      /nodata, $
      xtitle='z', ytitle='!7a, !6faint-end slope', $
      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick

if plot_qlf_points eq 'y' then begin
;; Actual data..
   plotsym, 0, 1., /fill
   oploterror, z_points, alpha, (alpha_err_up   - alpha),  psym=8, /HIBAR, errthick=linethick, color=16,  errcolor=24 
   oploterror, z_points, alpha, (alpha - alpha_err_down),  psym=8, /lobar, errthick=linethick, color=16, errcolor=24
;oploterror, z_points[14:15], alpha[14:15], (alpha_err_up[14:15]   - alpha[14:15]),  psym=8, /HIBAR, errthick=linethick, color=200,  errcolor=200 
;oploterror, z_points[14:15], alpha[14:15], (alpha[14:15] - alpha_err_down[14:15]),  psym=8, /lobar, errthick=linethick, color=200, errcolor=200
   plotsym, 0, 2., /fill
   oplot, z_points[14:15], alpha[14:15], psym=8, color=200
   
;; The various "best-fit" models...
   oplot, z_PLE_PD,  alpha_PLE_PD,  thick=linethick, color=blue
;;   PLE from Ross13 z=0.3-2.2
;   oplot, z_PLE,     alpha_PLE,     thick=linethick, color=blue  ;;
   oplot, z_LEDE,    alpha_LEDE,    thick=linethick, color=turquiose
   oplot, z_modLEDE, alpha_modLEDE, thick=linethick, color=red
endif 
oplot, z_HRH07,   alpha_HRH07,   thick=linethick, color=black

;;
;;  B E T A    b r i g h t - e n d    s l o p e...
;;
plot, z_HRH07, gamma_2, $
      position=[.40+xposoff, .60, .40+.24+xposoff, .98], $
;      xrange=[xmin, xmax], yrange=[1.0,2.65], $
      xrange=[xmin, xmax], yrange=[-1.6,-4.30], $
      /nodata, $
      xtitle='z', ytitle='!7b, !6bright-end slope', $
      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick

if plot_qlf_points eq 'y' then begin
   plotsym, 0, 1., /fill
   oploterror, z_points, beta, (beta_err_up   - beta),  psym=8, /HIBAR, errthick=linethick, color=16,  errcolor=24 
   oploterror, z_points, beta, (beta - beta_err_down),  psym=8, /lobar, errthick=linethick, color=16, errcolor=24
;oploterror, z_points[14:15], beta[14:15], (beta_err_up[14:15]   - beta[14:15]),  psym=8, /HIBAR, errthick=linethick, color=200,  errcolor=200 
;oploterror, z_points[14:15], beta[14:15], (beta[14:15] - beta_err_down[14:15]),  psym=8, /lobar, errthick=linethick, color=200, errcolor=200
   plotsym, 0, 2., /fill
   oplot, z_points[14:15], beta[14:15], psym=8, color=200
   
   oplot, z_PLE_PD,  beta_PLE_PD,  thick=linethick, color=blue
;   oplot, z_PLE,     beta_PLE,     thick=linethick, color=blue
   oplot, z_LEDE,    beta_LEDE,    thick=linethick, color=turquiose
   oplot, z_modLEDE, beta_modLEDE, thick=linethick, color=red
endif
oplot, z_HRH07,   beta_HRH07,   thick=linethick, color=black



plot, z_HRH07, log_phi_norm_HRH, $
;      position=[.72+xposoff, .60, .72+.24+xposoff, .98], $
      position=[.715+xposoff, .60, .72+.245+xposoff, .98], $
;      xrange=[xmin, xmax], yrange=[-5.7,-3.8], $
      xrange=[xmin, xmax], yrange=[-9.9,-4.5], $
;      xtitle='z', ytitle='!6Norm. log(phi)', $
      xtitle='z', ytitle='!6 log(!7u!6!15*!6!N)', $
      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick
if plot_qlf_points eq 'y' then begin
;; http://idlastro.gsfc.nasa.gov/ftp/pro/plot/oploterror.pro
   plotsym, 0, 1., /fill
   oploterror, z_points, log_phi, (log_Phi_err_up - log_phi),  psym=8, /hibar, errthick=linethick, color=16, errcolor=24
   oploterror, z_points, log_phi, (log_phi - log_Phi_err_down), psym=8, /lobar, errthick=linethick, color=16, errcolor=24
   plotsym, 0, 1.5, /fill
   oploterror, z_points[14:15], log_phi[14:15], (log_Phi_err_up[14:15] - log_phi[14:15]  ), psym=8, /HIBAR, errthick=linethick, color=200, errcolor=200 
   oploterror, z_points[14:15], log_phi[14:15], (log_phi[14:15] - log_Phi_err_down[14:15]), psym=8, /lobar, errthick=linethick, color=200, errcolor=200
   
   oplot, z_PLE_PD,  log_phi_star_PLE_PD,  thick=linethick, color=blue
;   oplot, z_PLE,     log_phi_star_PLE,     thick=linethick, color=blue
   oplot, z_LEDE,    log_phi_star_LEDE,    thick=linethick, color=turquiose
   oplot, z_modLEDE, log_phi_star_modLEDE, thick=linethick, color=red
endif
oplot, z_HRH07,   log_phi_norm_HRH,     thick=linethick, color=black


if plot_qlf_points eq 'y' then begin
;; LABELS
   xpos =   0.0
   xoff =   2.0
   yoff =  -0.80
   charsize_lab  = charsize /2.0
   charthick_lab = charthick/1.2

 ;  ypos =  -8.7
;   legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, psym=8, charsize=1.2, color=24
;   xyouts,      xpos+(0.3*xoff),  ypos+(0.5*yoff), '!6BOSS+MMT (Palanque++13) ', charsize=charsize_lab/1.1, charthick=charthick_lab, color=24
   ypos =  -8.7
   legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, psym=8, charsize=1.2, color=24
   xyouts,      xpos+(0.3*xoff),  ypos+(0.5*yoff), '!6BOSS DR9 (Ross++13) ', charsize=charsize_lab/1.1, charthick=charthick_lab, color=24
   ypos =  -9.2
   legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, psym=8, charsize=1.2, color=200
;xyouts,      xpos+(0.3*xoff),  ypos+(0.5*yoff), '!7z~!65 (McG13); !7z~!66 (Willott++10)', charsize=charsize_lab/1.1, charthick=charthick_lab , color=200
   xyouts,      xpos+(0.3*xoff),  ypos+(0.5*yoff), '!8z~!65 (McG13); !8z~!66 (Willott10)', charsize=charsize_lab/1.1, charthick=charthick_lab , color=200
endif


;;
;;
;;
plot, z_HRH07, log_L_star, $
      position=[0.08+xposoff,  .12, .32+xposoff, .50], $
      xrange=[xmin, xmax], yrange=[11.25,13.25], $
      /nodata, $
      xtitle='z', ytitle='!6Break Luminosity, log(L!13*!6!N)', $
      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick
oplot, z_HRH07, log_L_star, thick=linethick, color=black


plot, z_HRH07,  Mstar_HRH07, $
      position=[.40+xposoff,  .12, .40+.24+xposoff, .50], $
;      xrange=[xmin, xmax], yrange=[-19.0,-28.0], $
      xrange=[xmin, xmax], yrange=[-20.0,-31.0], $
      /nodata, $
      xtitle='z', ytitle='!6M!13*!6, Absolute Magnitude', $
      xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, charthick=charthick,  thick=thick
if plot_qlf_points eq 'y' then begin
   plotsym, 0, 1., /fill
   oploterror, z_points, MSTAR, (MSTAR_err_up   - MSTAR),  psym=8, /HIBAR, errthick=linethick, color=16,  errcolor=24 
   oploterror, z_points, MSTAR, (MSTAR - MSTAR_err_down),  psym=8, /lobar, errthick=linethick, color=16, errcolor=24
   plotsym, 0, 1.5, /fill
   oploterror, z_points[14:15], MSTAR[14:15], (MSTAR_err_up[14:15]   - MSTAR[14:15]),  psym=8, /HIBAR, errthick=linethick, color=200,  errcolor=200 
   oploterror, z_points[14:15], MSTAR[14:15], (MSTAR[14:15] - MSTAR_err_down[14:15]),  psym=8, /lobar, errthick=linethick, color=200, errcolor=200

   oplot, z_PLE_PD,  Mstar_PLE_PD,    thick=linethick, color=blue
;   oplot, z_PLE,     Mstar_PLE,       thick=linethick, color=blue
   oplot, z_LEDE,    Mstar_LEDE,      thick=linethick, color=turquiose
   oplot, z_modLEDE, Mstar_modLEDE,   thick=linethick, color=red
endif
oplot, z_HRH07,   Mstar_HRH07+3.7, thick=linethick, color=black


;; LABELS
xpos =   8.3
ypos = -31.0
xoff =   2.0
yoff =   0.80
charsize_lab  = charsize /2.0
charthick_lab = charthick/1.2

legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, linestyle = 0, charsize=1.2, color=black
xyouts,      xpos+(1.6*xoff),  ypos+(0.6*yoff), '"Full"', charsize=charsize_lab, charthick=charthick_lab 
xyouts,      xpos+(0.25*xoff), ypos+(2.0*yoff), 'Hopkins, Richards & Hernquist', charsize=charsize_lab, charthick=charthick_lab
xyouts,      xpos+(0.25*xoff), ypos+(3.1*yoff), '(2007)', charsize=charsize_lab, charthick=charthick_lab


if plot_qlf_points eq 'y' then begin
ypos = -31.0
ypos = -28.0
legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, linestyle = 0, charsize=1.2, color=blue
xyouts,      xpos+(1.6*xoff),  ypos+(0.6*yoff), '"PLE"', charsize=charsize_lab, charthick=charthick_lab, color=blue
xyouts,      xpos+(0.25*xoff), ypos+(2.0*yoff), 'Palanque-Delabrouille et al. (2013)', charsize=charsize_lab/1.1, charthick=charthick_lab, color=blue

;ypos = -28.0
;legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, linestyle = 0, charsize=1.2, color=blue
;xyouts,      xpos+(1.6*xoff),  ypos+(0.6*yoff), '"PLE"', charsize=charsize_lab, charthick=charthick_lab, color=blue
;xyouts,      xpos+(0.25*xoff), ypos+(2.0*yoff), 'Ross et al. (2013)', charsize=charsize_lab, charthick=charthick_lab, color=blue

ypos = -25.0
legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, linestyle = 0, charsize=1.2, color=turquiose
xyouts,      xpos+(1.6*xoff),  ypos+(0.6*yoff), '"LEDE"', charsize=charsize_lab, charthick=charthick_lab, color=turquiose
xyouts,      xpos+(0.25*xoff), ypos+(2.0*yoff), 'Ross et al. (2013)', charsize=charsize_lab, charthick=charthick_lab, color=turquiose

ypos = -22.0
legend, pos=[xpos,             ypos], ' ',      box=0, thick=12, linestyle = 0, charsize=1.2, color=red
xyouts,      xpos+(1.6*xoff),  ypos+(0.6*yoff), '"modLEDE2"', charsize=charsize_lab, charthick=charthick_lab, color=red
xyouts,      xpos+(0.25*xoff), ypos+(2.0*yoff), 'McGreer et al. (2013)', charsize=charsize_lab, charthick=charthick_lab, color=red
endif


;; http://ham.space.umn.edu/johnd/ct/ct-names.html
;loadct,0
;plot, z_test_HRH, log_L_star, $
;      xrange=[xmin, xmax], $ ; yrange=[-0.2,1.2], $
;      xtitle='z', ytitle='!6Norm. log(phi)', $
;      xstyle=1, ystyle=1, $
;      xthick=xthick,     ythick=ythick, $
;      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
;      charsize=charsize, charthick=charthick,  thick=thick, $
;      color=255
;loadct, clr_table

!p.multi=0

device, /close
set_plot, 'X'     





end


