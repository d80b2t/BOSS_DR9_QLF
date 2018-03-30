

function best_fit_linear_Phi, z, mag

;+
;
; From our best-fit PLE model...
;
; mag will be an Absolute Mag, in particular, 
; a M_i(z=2) measurement.  
;
;-

;;
;; Using Best fit IDMcG "linear Phi" model fit across 2.2-3.5
;; e.g. email from:
;;	From: 	Ian McGreer <imcgreer@as.arizona.edu>
;;	Subject: 	Re: super-fast update
;;	Date: 	October 15, 2012 9:17:52 AM PDT
;;	To: 	Nic Ross <npross@lbl.gov> 

;;log_phi_star =  -6.05
;Phi_star_0   =  -6.08   
;alpha        =  -1.42
;beta         =  -3.53
;M_star_0     = -26.70
;c1           =  -0.604
;c2           =  -0.678

;;
;;  From Table 8 from Ross et al. "rv1"   (aka the Version T-W-O... ;-)
;; 
     alpha = -1.30
      beta = -3.44
  M_star_0 = -26.51
        c1 = -0.684
        c2 = -0.843
Phi_star_0 = -5.92
  

;;
;; "pivot redshift''
;;
z_piv = 2.2

;;
;; Linear Phi "evolution"
;;
log_Phi_star = Phi_star_0 + ( c1 * (z - z_piv))
M_star       = M_star_0   + ( c2 * (z - z_piv))


dmag = mag -  M_star 

phi_demon = ( (10.0^(0.4*(alpha+1.)*(dmag))) + (10.0^(0.4*(beta+1.)*(dmag))) )

phi_star            = (10d^log_phi_star )
phi_star_dble_prime = 0.4*(10d^log_phi_star ) 

phi_model = phi_star / phi_demon



return, phi_model

end




