


readcol, 'Richards06_Table06.dat', $
         z_R06, Mi_z2, log_PhiR06, sigma_Phi, fill, z_bar, N_Q, N_Qcor

openw, 10, 'Richards06_Table06_reformatted.dat'
for ii=0ll, n_elements( z_R06) -1 do begin
   
   printf, 10, z_R06[ii], ' & ', Mi_z2[ii], ' & ', Mi_z2[ii], ' & ', $
           N_Q[ii], ' & ', log_PhiR06[ii], ' & ', sigma_Phi[ii], ' \\',$
;                 format='(f9.5,a, f12.5,a, e16.6,a, f16.8,a, i8,a, e16.6,a)'
                 format='(f7.3,a, f9.3,a, f9.3,a, i7,a, f8.3,a, f8.3,a)'
endfor
close, 10

end
