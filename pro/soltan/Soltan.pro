;;  Below I list a piece of Python code which gives the Croton09 fit to the QLF as a function of bJ magnitude and redshift.  This is a modified version of the Croom++04 fit.  I turned the densities into number per comoving (Mpc/h)^3 per magnitude from their per Mpc^3 per magnitude by dividing by hub**3 in the last line ... this of course assumes you've defined "hub" to be something useful like 0.72!  I was wondering whether we could overplot that line on your figure?  My by-eye comparisons show that Darren's modification does a pretty good job!

;;   lumfn_C09(mag, zz):
;;   The QSO luminosity function at redshift zz: dn/dM with magnitudes
;;   in the bJ system.
;;   Parameters for the LF fit from Croom++04, as modified by Croton09.
;;   This agrees with the C04 for z<3.
;;   """
;; 
;;

;; 
;; Copying straight from Darren's Code...
;;

;;  Redshift bins..
numz = 100
zz = findgen(numz)/10   

;; Virial Mass, 1e10 -> 1e15, in "alog10" steps ;-) 
mvir = 10.0^((findgen(60)/10)+10.0) 
numm = n_elements(mvir)

;Abs_Mag_model = fltarr(numm)
Abs_Mag_model = (findgen(numm)/4.)-32.

Mstar_g  = dblarr(numz)
Lstar_g  = dblarr(numz)
integral = dblarr(numz)

mass  = dblarr(numz, numm)
;mag   = dblarr(numz, numm)
Phi   = dblarr(numz, numm)
Phi_h = dblarr(numz, numm)

px = dblarr(numz)
py = dblarr(numz)

;; big G in units of (km s-1)^2 Mpc Msol^-1
G   =  4.3e-9   
;; H-nought
H0  = 100.0 
;; Something useful!  
hub = 0.70 


;print, alog10(1.67e-6)
;     -5.77728
;; MW's modified numbers
PhiStar   =   1.67e-6  ;;
Mstar0_bJ = -21.61d     ;; b_J magnitude value
Mstar0_g  = -22.24d     ;; converted to g-band (Croom09b

;;Convert from (log) bolometric luminosities in Watts to bJ magnitudes.
;;Uses Shen++09 bolometric calibration and Richards++06 conversion.
  ;  """
  ;  lgL  = lgLbol + 7.		# Convert from W to erg/s.
  ;  Miz2 = -2.5*lgL + 90.	# The i-band magnitude, k-corrected to z=2.
  ;  Mbj  = Miz2 + 0.71		# Finally convert to bJ.

Lstar0_g  = 10d^((Mstar0_g+(-90.))/(-2.5)-7.0)   ;; -7.0 puts it into Watts..


beta  =  -1.09
;; z< 3
;   k1,k2,a = 1.39,-0.29,-3.31
;;z>3
;   k1,k2,a = 1.22,-0.23,-3.31+0.5*(zz-3.0)

;; Croom et al. 2009 QLF
;; LEDE Model, Table 4
;PhiStar = 1.62e-6
;alpL = -3.48 & betL = -1.4
;Mstar0 = -22.24
;k1 = 1.23 & k2 = -0.206

for ii=0L, numz-1 do begin
   Hz = H0 * sqrt(0.3*(1.0+zz[ii])^3.0 + 0.7)
   ;Hz = H0 * sqrt(0.25*(1.0+zz[i])^3.0 + 0.75)
   
   if zz[ii] lt 3.0 then begin
      ;; Actual numbers from Croom09b,Table 4
      ;k1    =  1.23
      ;k2    = -0.206
      ;alpha = -3.48
      ;; Modified no.s from Croton/Martin
      k1    =  1.39d
      k2    = -0.29d
      alpha = -3.31d
   endif else begin
      ;; Croton modified values.... 
      ;; (no values for Croom09b above z=2.6...)
      k1    =  1.22d
      k2    = -0.23d
      alpha = -3.31d + (0.5d*(zz[ii]-3.0d))
   endelse
   
   ;; For both bJ and g-band...
   ;Mstar_bJ = Mstar0_bJ - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   Mstar_g[ii]  = Mstar0_g  - 2.5d*(k1*zz[ii] + k2*zz[ii]*zz[ii])
   Lstar_g[ii]  = 10d^((Mstar_g[ii]+(-90.d))/(-2.5d)-7.0d)
   
   for jj=0L, numm-1 do begin
      mag = -32.0 + (jj*0.25)
;      Abs_Mag_model[jj] = mag
      dmag   = mag - Mstar_g[ii] 

   ;   Lstar_g[ii]  = 10d^((Mstar_g[ii]+(-90.))/(-2.5)-7.0)
      Lum      = 10d^((mag+(-90.))/(-2.5)-7.0)
      dLum  = Lum - Lstar_g

;      Phi[ii, *] = PhiStar / ( (10.0^(0.4*(alpL+1)*(magQ-Mstar))) + (10.0^(0.4*(betL+1)*(magQ-Mstar))) )
      Phi[ii, jj] = PhiStar / ( (10.0^(0.4*(alpha+1)*(dmag))) + (10.0^(0.4*(beta+1)*(dmag))) )

;      print, ii, zz[ii], jj, mag, k1, k2, alpha, Phi[ii,jj]
     ; print, ii, zz[ii], jj, mag, Mstar_g, dmag,  Phi[ii,jj]  ;Lum, dLum

   endfor
   print, ii
   ;integral[ii] = Phi[ii,*]*6.8e14

  ; print, ii, zz[ii],  Mstar0_g, Mstar_g[ii], Lstar_g[ii], total(Phi[ii,*]), integral[ii]
endfor

Phi_sum = 0.0
redshift_bin =4
 lum_density = 0.0
for jj=0L, numm-1 do begin
      mag = -32.0 + (jj*0.25)
;;      Abs_Mag_model[jj] = mag
 ;     Lstar_g  = 10d^((Mstar0_g+(-90.))/(-2.5)-7.0)
;      Lum      = 10d^((mag+(-90.))/(-2.5)-7.0) ;; WATTS
      Lum      = 10d^((mag+(-90.))/(-2.5)) ;; ERG S^-1
 ;     dmag  = mag - Mstar_g 
  ;    dLum  = Lum - Lstar_g
;      print, mag, Mstar_g, dmag, Lum, dLum
      lum_density =   Lum*  phi[redshift_bin,jj] ;+ lum_density

   ;; TRAPEZIUM RULE!!
   ;;((phi[48,jj])-(phi[48,jj+1])/2.)*(0.250000)
   print, jj, -32.0 + (jj*0.25), Lum,  lum_density, phi[redshift_bin,jj], Phi_sum, phi[redshift_bin,jj]/ Phi_sum
   if jj lt numm-1 then  Phi_sum =  Phi_sum +  ((phi[redshift_bin,jj])-(phi[redshift_bin,jj+1])/2.)*(0.250000)
endfor



Phi_h = Phi / (hub^3)
print



end
