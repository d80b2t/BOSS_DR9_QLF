

readcol, 'My_QLF_iband_boss_wgt_formatted4paper_temp.dat', $
         redshift,a, mag,b, mmag,c, n, d,phi, $
         format='(d,a,d,a,d,a,d,a,d)'

;;Mi(z=2)= 72.5−2.5logLQ	(3)
Lbol = (mag-72.5)/(-2.5)

w = where(redshift ge 2.2 and redshift le 2.6, N)

;; trapezium rule...
trap = 0
for ii=0LL, N-2 do begin

   print, ii, trap, (10^Lbol[ii]-10^Lbol[ii+1])*(10^(phi[ii]) + 10^(phi[ii+1]))/2
   trap = trap + (10^Lbol[ii]-10^Lbol[ii+1])*(10^(phi[ii]) + 10^(phi[ii+1])/2.)

endfor

end



