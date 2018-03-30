

readcol, '/cos_pc19a_npr/BOSS/QLF/plot/number_counts/nm_cuml_BOSS_R06_imag_20120518.dat', $
         BOSS_imag,  BOSS_N_DR9, BOSS_DR9_nm_sum, BOSS_DR9_nm_sum_corr,  BOSS_DR9_nm_sum_err_corr

;; can make very dataset specific...
mag_in = BOSS_imag
sigma  = BOSS_DR9_nm_sum_err_corr
counts = BOSS_DR9_nm_sum_corr

readcol, '/cos_pc19a_npr/BOSS/QLF/plot/number_counts/nm_cuml_boss21_imag_20120518.dat', $
         boss21_imag, boss21_N, BOSS_boss21_nm_sum, BOSS_boss21_nm_sum_err
w = where(boss21_N ne 0, N)

;; can make very dataset specific...
mag_in = boss21_imag[w]
sigma  = BOSS_boss21_nm_sum_err[w] 
counts = BOSS_boss21_nm_sum[w]



;; From Myers03...
;N_zero   = 10.0
;alpha_d  = 0.98
;beta_d   = 0.15
;m_nought = 19.10
;m_running = (findgen(100)/20.)+16.5 



;; alpha, the (bright-end in C09b) slope
alpha_bins = 250   ;; for realz...
alpha_bins = 100   ;; for realz...
;alpha_bins = 40     ;; for testing..
;alpha_bins = 10     ;; for testing..
alpha = fltarr(alpha_bins)

;; beta, the (faint-end slope in C09b)
beta_bins = 250    ;; for realz...
beta_bins = 100    ;; for realz...
;beta_bins = 40     ;; for testing..
;beta_bins = 10      ;; for testing...
beta = fltarr(beta_bins)

denom= fltarr(alpha_bins,beta_bins)

;; 
;; break magnitude...
;;
mag_bins = 40.
break_mag = fltarr(mag_bins)

;;
;; N_star, the normalization
;;
N_star_bins = 250
N_star_bins = 100
N_star = fltarr(N_star_bins)


alpha_norm = -2.00 
alpha_step = 0.05

beta_norm = -1.6 
beta_step = 0.1

mag_norm = 14.0
mag_step = 0.4

N_zero    = 0.2
N_step    = 0.02

;; some defaults:
; alpha_i = 5 
; beta_i  = 4
; mag_i = 15
; N_i = 50

chi_sq_min = 100000.
xxx_0min   = 100000.
for alpha_i = 0LL, alpha_bins-1 do begin
   alpha[alpha_i] = alpha_norm + (alpha_i*alpha_step)
   print, 'alpha[alpha_i] ', alpha[alpha_i]
  
   for beta_i = 0LL, beta_bins-1 do begin
      beta[beta_i] = beta_norm + (beta_i*beta_step)
      
      for mag_i = 0ll, mag_bins-1 do begin
         break_mag[mag_i] = mag_norm + (mag_i*mag_step)

;         print, alpha[alpha_i], beta[beta_i], break_mag[mag_i]
;         dmag    = break_mag[mag_i] - mag_in 
         dmag    = mag_in - break_mag[mag_i]
         denom = (10.0^((-1.)*alpha[alpha_i]*dmag)) + (10.0^((-1.)*beta[beta_i]*dmag))
         
         for N_i = 0L, N_star_bins-1 do begin
            N_star[N_i] = 10^(N_zero+(N_step*N_i))
            xxx = N_star[N_i]  / denom
            
;            if xxx[0] lt xxx_0min then  xxx_0min = xxx[0] 
            
            chi_sq_bin = ( (counts -  xxx)^2) $
                         /  (sigma^2)
            
            chi_sq = total(chi_sq_bin)
            


            if chi_sq lt chi_sq_min then begin
               print, N_i, alpha_i, beta_i, mag_i, chi_sq_min, chi_sq, $
                      alpha_norm + (alpha_i*alpha_step), beta_norm + (beta_i*beta_step), $
                       mag_norm + (mag_i*mag_step), N_star[N_i]
                     chi_sq_min  = chi_sq 
                     alpha_i_min = alpha_i
                     beta_i_min  = beta_i
                     N_i_min   = N_i
                     mag_i_min = mag_i

            endif
         endfor
      endfor 
   endfor
endfor

print
print
print, ' N_zero   = ', N_star[N_i_min]
print, ' alpha_d  = ', alpha_norm + (alpha_i_min*alpha_step)
print, ' beta_d   = ', beta_norm + (beta_i_min*beta_step)
print, ' m_nought = ', mag_norm + (mag_i_min*mag_step)
print
print

end
