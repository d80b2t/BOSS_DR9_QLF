


openw, 10, 'CW12_qsolf_allz_temp.dat'

readcol, 'CW12_qsolf_z0.50.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  0.50', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z1.00.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  1.00', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z1.50.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  1.50', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z2.00.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  2.00', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z2.40.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  2.40', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z2.80.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  2.80', mag[ii], log_phi[ii] 
endfor


readcol, 'CW12_qsolf_z3.25.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  3.25', mag[ii], log_phi[ii] 
endfor


readcol, 'CW12_qsolf_z3.75.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  3.75', mag[ii], log_phi[ii] 
endfor

readcol, 'CW12_qsolf_z4.75.dat', mag, log_phi
for ii=0ll, N_elements(mag)-1 do begin
   printf, 10, '  4.75', mag[ii], log_phi[ii] 
endfor


close, 10
close, /all
end
