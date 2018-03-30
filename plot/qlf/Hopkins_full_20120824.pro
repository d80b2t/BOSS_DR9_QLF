

function hopkins_full, z, log_l_bol
;function hopkins_full, z, mag

;+
;
; From Table 3 of HRH07
;
;-

;; Want to convert Mi2 input into MB:
;; Assume:  M_B ~ M_bJ
;;   Mi_HRH07 = ABMag_HRH07 - 0.71
;Bmag = mag_in + 0.71

;; Then want to convert Bmag to an "Mbol"

;;
P0=8.99833  & P1=6.24800 &  P2=-0.370587 & P3=-0.0115970

;lband = P0*pow(10.,P3*x) + P1*pow(10.,P2*x) ;






z_ref = 2.00

log_phi_star = -4.825 
log_L_star_0 = 13.036

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
print
print, 'xi ', xi

log_L_star = (log_L_star_0) + (k_L_1*xi) + (k_L_2*xi*xi) + (k_L_3*xi*xi*xi)
print
print, 'log_L_star', log_L_star

gamma_one = gamma_one_0 * (10^(k_gamma_one * xi))
gamma_two = (2.0 *  gamma_two_0) * ((10^(xi*k_gamma_two_1)) + (10^(xi*k_gamma_two_2) ))^(-1)

print, 'gamma_one', gamma_one
print, 'gamma_two', gamma_two

alpha = (gamma_one+1)*(-1.)
beta  = (gamma_two+1)*(-1.)
print, 'alpha = -(gamma_one +1)' , alpha
print, 'beta = -(gamma_two +1)' , beta

log_L_star = log_L_star + alog10(3.9) + 33.

print
print, 'log_L_star', log_L_star
print

;L_star = 10d^(log_L_star)
;L_star = (log_L_star)

print
print, 'log_L_star', log_L_star , '  log_L_bol', log_L_bol
print

print,' (10^log_phi_star)  ',  (10^log_phi_star) 

x = log_L_bol - log_L_star

;phi_model = (10^log_phi_star) / ( ((L_bol / L_star)^gamma_one) +  ((L_bol / L_star)^gamma_two))
phi_model = 10^(log_phi_star -  alog10(10^(x*gamma_one) +  10^(x*gamma_two)))

;; 
;;  Bolometric corrections 
;;  e.g. Eqn. (2) of HRH07
;;

;; log_l_bol = (findgen(100)/12.)+8.

print, 'log_l_bol - 10. ', log_l_bol - 10.
lband = 0.

;xx =  log_l_bol - 10.
xx = (10^(log_l_bol - 10.)) / (3.9d33)
print, ' xx  ', xx
print

;; for -1, the B-band...
P0=8.99833 & P1=6.24800 & P2=-0.370587 & P3=-0.0115970

print, 'P0, P1, P2, P3, ', P0, P1, P2, P3

lband = (P0 * (10^(P3*xx))) + (P1*(10^(P2*xx)) )
print, 'lband ', lband

c1 =  6.25
k1 = -0.37
c2 =  9.00
k2 = -0.012
L_boller = (findgen(120)/20.)+8.6
bol_cor = (c1*((L_boller)^k1)) + (c2*((L_boller)^k2))

;;plot, log_l_bol, lband, xrange=[8.5, 14], /ylog

b_band = (10d^log_l_bol)/lband

print, 'b_band ',  alog10(b_band)

return, phi_model

end




