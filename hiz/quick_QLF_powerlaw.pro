


 faint_end_slope = -1.80
bright_end_slope = -3.26
break_magnitude  = -26.39
   normalisation = 10d^(-8.40)

;read,  faint_end_slope, PROMPT=' - Faint  end slope value?? '
;read, bright_end_slope, PROMPT=' - Bright end slope value?? '
;read,  break_magnitude, PROMPT=' - Break Magnitude?? '

   alpha = bright_end_slope
    beta = faint_end_slope
  M_star = break_magnitude
Phi_star = normalisation

;; No. of mag bins
mag_bins = 60

;; Set up phi...
mag = fltarr(mag_bins)
Phi = fltarr(mag_bins)

openw, 10, 'QLF_from_temp.dat'
printf, 10, '##' 
printf, 10, '## alpha  beta  M_star  Phi_star'
printf, 10,    alpha, beta, M_star, Phi_star,  $
        format='(f8.3, f8.3, f8.3, f8.3)' 
printf, 10, '## '

for jj=0ll, mag_bins-1 do begin
   
   mag[jj] = -32.0 + (jj*0.25)
   ;; => if mag_bins = 60, -32.00 < M_g < -17.25
   dmag   = mag[jj] - M_star 
   
   ple_denom =  ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )
   Phi[jj] = Phi_star / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

   print,      jj, mag[jj], phi[jj], alog10(Phi[jj])
   printf, 10, jj, mag[jj], phi[jj], alog10(Phi[jj])
endfor


plot, mag, phi, $
      /ylog, $
      xrange=[-22, -30], xstyle=1, $
      xtitle='!6magnitude', $
      ytitle='!6log!I10!N !7U!3(M!Ii!N(z=2)) [Mpc!E-3!N mag!E-1!N]'


end
