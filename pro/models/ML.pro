;+
; NAME:
;   program_name
;
; PURPOSE:
;   Purpose here. 
;
; CALLING SEQUENCE:
;    program_name, [ option= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine ...
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   readspec
;   splog
;
; NOTES:
;
; REVISION HISTORY:
;   11-Jan-2011  v0.0.1     NPR
;-


;;
readcol, '../../data/Richards06_Table05.dat', $
         name, z_R06, i_R06, M_i_R06, delta_gi_R06, cor, $
         format='(a, f,f,f,f,f)'

;; Richards table 6, the BINNED QLF
readcol, '../../data/Richards06_Table06.dat', $
         z_R06_bin, M_i_R06_bin, log_Phi_R06, sigma_Phi_R06, fil, z_bar_R06, N_Q_R06, N_Qcor_R06

z_bins = [0.30, 0.68, 1.06, 1.44, 1.82, 2.20, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00]

N_sum=0.00
openw, 10, 'Richards06_Table6_NPR_temp.dat'
for ii=0,10 do begin
   for jj=0L, 30 do begin

      ;; starting faint and getting brighter...
;      mag_bin_min_R06 = -22.5-(jj*0.30)
 ;     mag_bin_max_R06 = -22.5-((jj+1)*0.30)

      ;; Start bright, getting fainter (to match R06, table 6)
      mag_bin_min_R06 = -31.80+(jj*0.30)
      mag_bin_max_R06 = -31.80+((jj+1)*0.30)

      mag_bin_R06     = 0.5*( mag_bin_min_R06 +  mag_bin_max_R06)
      
      w = where(z_R06   ge z_bins[ii] and z_R06 lt z_bins[ii+1] and $
                M_i_R06 ge  mag_bin_min_R06 and M_i_R06 lt mag_bin_max_R06, N)

;      mean_z = mean(z_R06[w])
      mean_z = 0.5* (z_bins[ii] +  z_bins[ii+1] )

      print, ii, jj, mean_z, mag_bin_R06, N

      if N gt 0 then printf, 10, mean_z, mag_bin_R06, N
      
      if N ge 3 then N_sum = N + N_sum
      
   endfor
endfor
close, 10
close, /all


;; Richards 06 Table 7 values...
; FIXED PARAMETERS
z_ref =   2.45
M_star = -26.00

; FREE PARAMETERS
log_Phi_star = -5.75
A_one   = 0.78                  ; +/-0.01
B_one   = 0.10
B_two   = 27.35
B_three = 19.27

;; Equation 5 of Richards et al. (2006)
;Phi(M,z

Phi_expt     = findgen(N_elements(z_R06_bin))
PHI_R06      = findgen(N_elements(z_R06_bin))
log_Phi_expt = findgen(N_elements(z_R06_bin))
chi_squared = 0.000
;; Looping over the binned QLF, working out a chi-squared's
for kk=0LL,N_elements(z_R06_bin)-1 do begin
   
;; Equation 8 of Richards et al. (2006)
   z = z_R06_bin[kk]
   M = M_i_R06_bin[kk]
   Phi_R06[kk]  = 10^(log_Phi_R06[kk])
   
   xi  = alog10( (1+z)/(1+z_ref))
   
;; Equation 7 of Richards et al. (2006)
   mu = M - (M_star + (B_one*xi) + (B_two*xi*xi) + (B_three*xi*xi*xi)) 
   
   
;; Equation 6 of Richards et al. (2006)
;; FIXED POWER LAW
   Phi_Star = 10^(log_Phi_star) 
   
   Phi_expt[kk] = Phi_star * (10^(A_one * mu) )
   
   log_Phi_expt[kk] = alog10(Phi_expt[kk])
   
   ;; Can I work out a chi_squared??
   ;; R06 has: 
   ;; "However, the 􏰢2 of this ML fit (as compared to the
   ;; binned QLF) is 394 for 69 degrees of freedom"
   ;;

;   chi_squared =  chi_squared + ( (( Phi_R06[kk] - Phi_expt[kk] )^2)  / Phi_expt[kk]    )
   chi_squared =  chi_squared + ( (( Phi_R06[kk] - Phi_expt[kk] )   / (SIGMA_PHI_R06[kk]*1e-9  ))^2)

;   chi_squared =  chi_squared + ( (( log_Phi_R06[kk] - log_Phi_expt[kk] )^2)  / log_Phi_expt[kk]    )

   print, kk, 10^(log_Phi_R06[kk]), Phi_expt[kk], z, z_ref, xi, mu, chi_squared
   
endfor 




end
