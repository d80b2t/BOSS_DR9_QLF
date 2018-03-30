;+
;
; The extension of the Croom04 2QZ QLF model from 
; Croton, 2009, MNRAS, 394, 1109   (their Table 1). 
;
;-

fltarr_bins  = 150.
fltarr_bins2  = 60.

z_out         = (findgen(fltarr_bins))/25.
M_range_out   = (findgen(fltarr_bins2)*0.30)-31.05
phi_fit_a_out =  fltarr(fltarr_bins,fltarr_bins2)
Mstar_bJ_out  =  fltarr(fltarr_bins)

openw, 10, 'Croton09_model_temp.dat' 

for ii=0LL, n_elements(z_out)-1 do begin
   if z_out[ii] lt 3.0 then begin
      alpha_test     =  -3.31
      beta_test      =  -1.09
      Mstar0_bJ_test = -21.61
      k1_test        =   1.39   
      k2_test        =  -0.29  
      Phi_star_test  =   1.67e-6
   endif
   if z_out[ii] ge 3.0 then begin
      alpha_test     =  -3.31+(0.5*(z_out[ii] -3.))
      beta_test      =  -1.09
      Mstar0_bJ_test = -21.61
      k1_test        =   1.22   
      k2_test        =  -0.23  
      Phi_star_test  =   1.67e-6
   endif

      Mstar_bJ_out[ii]   = Mstar0_bJ_test  - 2.5*(k1_test*z_out[ii]  + k2_test*z_out[ii]*z_out[ii])

   for jj=0LL, n_elements(M_range_out)-1 do begin
      dmag_out      = M_range_out[jj] - Mstar_bJ_out[ii]
      ple_denom_test = (10.0^(0.4*(alpha_test+1)*(dmag_out))) + (10.0^(0.4*(beta_test+1)*(dmag_out))) 
      phi_fit_a_out[ii,jj] = Phi_star_test /ple_denom_test
      print, z_out[ii], M_range_out[jj], phi_fit_a_out[ii,jj]
      printf, 10, z_out[ii], M_range_out[jj], phi_fit_a_out[ii,jj]
   endfor
endfor
close, 10


end
