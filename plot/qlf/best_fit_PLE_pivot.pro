

function best_fit_PLE, z, mag

;+
;
; From our best-fit PLE model...
;
; mag will be an Absolute Mag from 
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


ACTUALLY!!!!

N.B. This is currently just a stand-in bit of code...



;; Best fit PLE across 0.4-2.2...
log_phi_star =  -6.11
alpha        =  -3.50
beta         =  -1.50
M_star_0     = -23.25 
k1           =   1.16
k2           =  -0.228


;; PLE "evolution"
M_star = M_star_0 - 2.5*(k1*z + k2*z*z)
;print

dmag = mag -  M_star 

phi_demon = ( (10.0^(0.4*(alpha+1.)*(dmag))) + (10.0^(0.4*(beta+1.)*(dmag))) )

phi_star            = (10d^log_phi_star )
phi_star_dble_prime = 0.4*(10d^log_phi_star ) 

phi_model = phi_star / phi_demon



return, phi_model

end




