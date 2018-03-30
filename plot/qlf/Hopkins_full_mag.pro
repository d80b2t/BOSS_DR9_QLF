

function hopkins_full_mag, z, mag

;+
;
; From Table 3 of HRH07
;
;-

readcol, '/cos_pc19a_npr/BOSS/QLF/pro/models/HRH07/HRH07_Bband_z2.00.dat', Bbol, AbMag_z20, flux, bol, phi_z20
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/models/HRH07/HRH07_Bband_z2.40.dat', Bbol, AbMag_z24, flux, bol, phi_z24
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/models/HRH07/HRH07_Bband_z3.25.dat', Bbol, AbMag_z32, flux, bol, phi_z32

Mi_HRH07 = ABMag_z20 - 0.71

;; Want to convert Mi2 input into MB:
;; Assume:  M_B ~ M_bJ
;;   Mi_HRH07 = ABMag_HRH07 - 0.71
;Bmag = mag_in + 0.71

;; Then want to convert Bmag to an "Mbol"

;;
P0=8.99833  & P1=6.24800 &  P2=-0.370587 & P3=-0.0115970

;lband = P0*pow(10.,P3*x) + P1*pow(10.,P2*x) ;

;z     = 2.00

z_ref = 2.00

log_phi_star  =  -5.325 
log_L_star_0  =  13.036

k_L_1 =   0.632
k_L_2 = -11.76
k_L_3 = -14.25

gamma_one_0 = 0.417
gamma_two_0 = 2.174

k_gamma_one   = -0.623
k_gamma_two_1 =  1.460
k_gamma_two_2 = -0.793

;M_star = 

xi =  alog10((1. + z)/(1. + z_ref))
;print
print, 'xi ', xi

log_L_star = (log_L_star_0) + (k_L_1*xi) + (k_L_2*xi*xi) + (k_L_3*xi*xi*xi)
;print
;print, 'log (L_star/L_sol)', log_L_star

gamma_one = gamma_one_0 * (10^(k_gamma_one * xi))
gamma_two = (2.0 *  gamma_two_0) * ((10^(xi*k_gamma_two_1)) + (10^(xi*k_gamma_two_2) ))^(-1)

;print, 'gamma_one', gamma_one
;print, 'gamma_two', gamma_two
;print

;print, ' (gamma_one+1.)*(-1.)  ',  (gamma_one+1.)*(-1.) 
;print, ' (gamma_two+1.)*(-1.)  ',  (gamma_two+1.)*(-1.)
;print

;alpha =  (gamma_one+1.)*(-1.)
;beta  =  (gamma_two+1.)*(-1.)

alpha = -1.2
beta  = -2.7

;print, 'alpha = -(gamma_one +1)' , alpha
;print, 'beta = -(gamma_two +1)' , beta
;print

;print, -2.5*alog10( 4.2370624d+46/3.9d33) + 4.72
;      -27.870001

L_star = (10d^(13.036))*3.9d33
;print, 'L_star, log_L_star ', L_star, alog10(L_star )
;print

;; M_bol_star - M_Bol_sun = -2.5 * alog10 (L_star/L_star)

;Mstar = (-2.5*alog10( L_star/3.9d33)) + 4.72
Mstar = -25.6

;print, 'Mstar  ', Mstar
;print

dmag = mag -  Mstar 
phi_demon = ( (10.0^(0.4*(alpha+1.)*(dmag))) + (10.0^(0.4*(beta+1.)*(dmag))) )

phi_star_dble_prime = 0.4*(10d^log_phi_star ) 

phi_model = phi_star_dble_prime / phi_demon

;; take into account HRH07 numbers are in 
;;   dphi/dlog_{10}(L)  [ Mpc^{-3} log_{10}(L)^{-1} ]
;; 
per_mag = alog10(2.5)

;plot, mag,       (phi_model), /ylog, xrange=[-15, -35]
;plot,  mag, alog10(phi_model)-per_mag, $
;       xrange=[-18.25, -30.5], yrange=[-9.2,-4.7], ps=2, xstyle=1, ystyle=1


;oplot, Mi_HRH07, alog10(phi_z20)-per_mag, color= 80, thick=3.;, ps=2
;oplot, Mi_HRH07, alog10(phi_z24)-per_mag, color=160, thick=3., ps=2
;oplot, Mi_HRH07, alog10(phi_z32)-per_mag, color=240, thick=3., ps=2

return, phi_model

end




