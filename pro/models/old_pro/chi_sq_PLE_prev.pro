;+
;
; Kinda assumes you've already run chi_sq_PLE...
;
;-

mag_bins = 30
mag_PLE  = fltarr(mag_bins)

alpha    = -1.30
beta     = -2.90
Phi_star = 10^(-5.92)

Phi_PLE_prev   = fltarr(mag_bins)

for jj=0L, mag_bins-1 do begin 
   mag_PLE[jj]      = min(Abs_mag_bin_boss_s82[w]) + (jj*0.30)  ;;   and put them on the same "absolute scale" 
   dmag             = mag_PLE[jj] - Mstar_g 
   ple_denom        = ((10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))))
   Phi_PLE_prev[jj] = Phi_star / ((10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag)))) 
   ;print, jj, mag_PLE[jj], dmag, ple_denom, Phi_PLE_prev[jj]
   print, ii, jj, alpha, beta, alog10(phi_star), mag_PLE[jj], ple_denom, Phi_PLE_prev[jj]
endfor 

oplot, mag_PLE, alog10(Phi_PLE_prev), thick=4, color=200, linestyle=2


end
