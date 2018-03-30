
;; From: 
;; http://cerebus.as.arizona.edu/~ioannis/research/red/
;; Set up cosmology...
;; (will need this for the Luminosity Distances...
;;
;; Not sure if this ever gets called properly, 
;; just type in at command-line...

 red, omega0=0.30, omegalambda=0.70, h100=0.70 


   readcol, '../../data/SDSS_QSO_DR3.dat', name_dr3, z_dr3, imag_dr3, Mi_dr3, del_gi_dr3, corr_dr3, format='(a,f,f,f,f,f)', /silent
   print
   print, ' SDSS DR3 (from Richards et al. 2006) read-in'
   print
   print
   
   ;;
   ;; Reading in Richards et al. (2006)
   ;; Table 1, the completeness corrections...
   readcol, '../../data/Richards06_Table01.dat', imag_comp, z_comp, point, radio, ext, /silent
   
   ;;
   ;; Reading in the Richards06 k-correction
   ;; Normalized at z=2. 
   readcol, 'Richards_2006_Table4.dat', kcor_redshift, kcor, /silent
   print, 'Richards06_kcor.dat READ-IN', n_elements(kcor)
   print
   print



;; Running the DLUMINOSITY command from the red.pro routine
;; Returns LUMINOSITY DISTANCES (using the above cosmology)
;; Divide by 1e6 to get into Mpc
   
   print, 'Doing DLUMS....'
   dlums_dr3    = DLUMINOSITY(z_dr3)  / 1e6

;; Working out the ABSOLUTE MAGNITUDE for each object
;; Recall 1st year (UG!) notes: 
;; DIST_MOD = 5 log (D_L /10pc) which means if D_L is in Mpc:

Abs_iM_dr3   = imag_dr3       - (5 * alog10(dlums_dr3))   - 25.00 - kcor(fix(z_dr3/0.01)) 


plot, z_dr3, Abs_iM_dr3, $
      psym=3, $
      /nodata, $
      position=[0.22, 0.22, 0.98, 0.98], $
      xstyle=1, $
      ystyle=1, $
      xrange=[-0.2, 6.0], $
      yrange=[-17.0, -33.0], $
      xthick=xthick, $
      ythick=ythick, $
      charsize=charsize, $
      charthick=charthick, $
      xtitle='Redshift, z', $
      ytitle='M!Ii!N[z=0] '


end
