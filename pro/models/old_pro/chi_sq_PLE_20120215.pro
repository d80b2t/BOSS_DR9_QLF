;+
;
; Okay, so this is just a fairly wee, fairly simple script that
; takes the "modified" LDDE model from Croom et al. (2009b), 
; and computes some chi_squared values for the BOSS DR9 dataset
; 
; I'm initially testing on 
;
;
;
;-


;;  THE   S T R I P E   8 2   DATA  - 6,250 Quasars strong
;readcol, '../qlf/My_QLF_iband_boss_narrowZ_S82.dat', $  ;; USE THIS FOR THE S82 plots!!!
readcol, '../qlf/My_QLF_iband_boss_narrowZ_S82_wExtinc_20120622.dat', $
         z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82

log_Phi_BOSS_S82 =  log_Num_den_boss_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)


;boss_ilimit = 21.80
;;red_bins = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]
;red_bins =     [2.40,       2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.90, 3.13, 3.35, 3.75, 4.00]
;dlums = DLUMINOSITY(red_bins) / 1e6
;Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ;+ kcor(fix(red_bins/0.01))
;red_bins =     [   2.40,   2.25,    2.35,    2.45,     2.55,   2.65,     2.75,    2.90,    3.13,   3.35,     3.75,    4.00]
Abs_iMag_limit =[-24.640, -24.470, -24.585, -24.695, -24.800, -24.901, -24.998, -25.137, -25.337, -25.514, -25.806, -25.973]
;; Note, Abs_iMag_limit[0] never used...

Abs_iMag_limit_bin = Abs_iMag_limit[10]
z_bin_range_min = 3.50
z_bin_range_max = 4.00

;ii = 45  ;; picking out zz[ii]=2.25 
;ii = 47  ;; picking out zz[ii]=2.35 
;ii = 49  ;; picking out zz[ii]=2.45
;ii = 51  ;; picking out zz[ii]=2.55
;ii = 53  ;; picking out zz[ii]=2.65
;ii = 55  ;; picking out zz[ii]=2.75
;ii = 58  ;; picking out zz[ii]=2.90
;ii = 63  ;; picking out zz[ii]=3.15
;ii = 67  ;; picking out zz[ii]=3.35
ii = 75  ;; picking out zz[ii]=3.75

;w = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max and raw_N_boss_QSOs_s82 ge 1, N)

w = where(z_bin_boss_s82 ge z_bin_range_min and z_bin_boss_s82 le z_bin_range_max $
          and raw_N_boss_QSOs_s82 gt 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit_bin , N_qlf_bins)

plotsym, 0, 1.6, /fill
plot, Abs_mag_bin_boss_s82[w] , log_Phi_BOSS_S82[w], psym=8 
Mi_boss_err_S82 = log_Num_den_boss_S82[w] - log_Num_den_boss_S82[w]
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_up_S82[w],  $
            /hibar, errthick=3., psym=8 ;, color=160, ERRcolor=160
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_down_S82[w], $   
            /lobar, errthick=3., psym=8 ;, color=160, errcolor=160



;; No. of redshift bins
zbins = 140. 
;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
zz = findgen(zbins)/20.    ;;  Redshift bins..

;; No. of mag bins
mag_bins = 30
mag_PLE = fltarr(mag_bins)

;; A, the normalization
phi_star_bins = 50
phi_star = fltarr(phi_star_bins)

;; alpha, the (bright-end in C09b) slope
alpha_bins = 100
alpha = fltarr(alpha_bins)

;; A, the normalization
beta_bins = 100
beta = fltarr(beta_bins)

;; Set up phi...
;Phi_LDDE = fltarr(zbins, mag_bins)
Phi_PLE = fltarr(mag_bins, phi_star_bins, alpha_bins, beta_bins)

print
print
print, zz[ii], z_bin_range_min , z_bin_range_max, Abs_iMag_limit_bin
print
print

;---------------------------------------------------------------
; 
;  P  L  E 
;
;   Good ol' PLE!! 
;   Section 6.1 from Croom et al. (2009b)
;   e.g. Table 2, Row 4, 

;    alpha = -3.33
;    beta = -1.41
Mstar0_g = -22.17
k1 = 1.46
k2 = -0.328
; Phi_star = 10^(-5.84)

Phi_star_norm = -7.27
alpha_norm = -5.50
beta_norm  = -3.00

alpha_step = 0.05
beta_step  = 0.05
Ai_step    = 0.05

openw, 9, 'Stripe82_PLE_temp.dat'


;for ii=0L, zbins-1 do begin
;; Equation (11) of Croom et al. (2009b)
Mstar_g  = Mstar0_g  - 2.5*(k1*zz[ii] + k2*zz[ii]*zz[ii])
    ;; => Mstar_g(z=2.25) = -26.2313
    ;; => Mstar_g(z=2.45) = -26.1905

for alpha_i = 0LL, alpha_bins-1 do begin
   ;alpha_i = 0
   ;alpha[alpha_i] = -3.33
   alpha[alpha_i] = alpha_norm + (alpha_i*alpha_step)
  ; alpha[alpha_i] = alpha_norm + (alpha_i*0.05)
   print, alpha_i

   for beta_i = 0LL, beta_bins-1 do begin
      ;beta[beta_i] = -1.41
      beta[beta_i] = beta_norm + (beta_i*beta_step)
      ;beta[beta_i] = beta_norm + (beta_i*0.05)


      for Ai = 0L, phi_star_bins-1 do begin
         Phi_star[Ai] = 10^(Phi_star_norm+(Ai*Ai_step))
         ;Phi_star[Ai] = 10^(Phi_star_norm+(Ai*0.05))
         
         for jj=0L, mag_bins-1 do begin
;      mag_PLE[jj] = -32.0 + (jj*0.25)      
;   mag_PLE[jj] = -32.0 + (jj*0.30)                       ;; to match the data's delta_mags
            mag_PLE[jj] = min(Abs_mag_bin_boss_s82[w]) + (jj*0.30) ;;   and put them on the same "absolute scale" 
                                                                                
            dmag    = mag_PLE[jj] - Mstar_g 
            
            ple_denom =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
            Phi_PLE[jj, Ai, alpha_i, beta_i] = Phi_star[Ai]  / ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
            
;            print, ii, jj, alpha_i, beta_i, Ai, zz[ii], alog10(phi_star[Ai]), mag_PLE[jj], ple_denom, Phi_PLE[jj, Ai, alpha_i, beta_i] 
;      printf, 9, ii, zz[ii], jj, mag_LDDE[jj], e_d[ii,jj], zc, Phi_LDDE[ii,jj], $
                                ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
         endfor
;         if (Ai mod 2) eq 1 then oplot, mag_PLE, alog10(Phi_PLE[*, Ai, alpha_i, beta_i]), color=100
      endfor
   endfor ;; beta
endfor;; alpha
close, 9

oplot, mag_PLE, alog10(Phi_PLE[*, 0,  0,  0]), color=40, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 10, 50, 0]), color=60, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 10, 99, 0]), color=80, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 35,  0, 49]), color=120, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 35, 50, 49]), color=160, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 35, 99, 49]), color=200, thick=2
oplot, mag_PLE, alog10(Phi_PLE[*, 24, 49, 49]), color=240, thick=2

openw, 11, 'chi_sq_values_PLE_inside_mag_loop_temp.dat'
printf, 11, '# Normalization, A  Abs i-mag     Phi_BOSS_S82, sigma_Phi_BOSS_S82,    Phi_LDDE_model,  chi_sq_bin, chi_sq '
printf, 11, '# (-9.57+(Ai*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]), sigma_Phi_BOSS_S82[w[kk]], Phi_LDDE[ii,kk, Ai], chi_sq_bin, chi_sq '

openw, 10, 'chi_sq_values_PLE_temp.dat'
printf, 10, '# Normalization A,     chi_sq_bin, chi_sq '
printf, 10, '# (-9.57+(Ai*0.01)),   chi_sq_bin, chi_sq '


;; http://en.wikipedia.org/wiki/Goodness_of_fit

chi_sq_min = 1000.
for alpha_i = 0LL, alpha_bins-1 do begin
   for beta_i = 0LL, beta_bins-1 do begin
      for Ai = 0L, phi_star_bins-1 do begin
         
         chi_sq=0         
         for kk=0L, N_qlf_bins-1 do begin
                                ;chi_sq_bin = ( log_Phi_BOSS_S82[w[kk]] -  alog10(Phi_LDDE[ii,kk]))^2 /  (alog10(sigma_Phi_BOSS_S82[w[kk]])^2)
            chi_sq_bin = ( 10^(log_Phi_BOSS_S82[w[kk]]) -  Phi_PLE[kk, Ai, alpha_i, beta_i])^2 $
                         /  (sigma_Phi_BOSS_S82[w[kk]]^2)
            chi_sq = chi_sq + chi_sq_bin
            
                                ;print, Ai, kk, Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]])*1e7, sigma_Phi_BOSS_S82[w[kk]]*1e7, $
                                ;       Phi_LDDE[ii,kk,Ai]*1e7, chi_sq_bin, chi_sq
                                ;print, Abs_mag_bin_boss_s82[w[kk]], log_Phi_BOSS_S82[w[kk]], alog10(sigma_Phi_BOSS_S82[w[kk]]), Phi_LDDE[ii,kk], $
                                ;      chi_sq_bin, chi_sq[kk]
            
;            printf, 11, (Phi_star_norm+(Ai*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]),  sigma_Phi_BOSS_S82[w[kk]], $
;;                    Phi_LDDE = fltarr(mag_bins, A_norm_bins, 1, beta_bins)
 ;                   Phi_PLE[kk, Ai, alpha_i, beta_i], $
  ;                  chi_sq_bin, chi_sq, $
                    format='(f9.4,2x, f9.4,2x, e,e,e, f, f)'
         endfor
         
         if chi_sq lt chi_sq_min then begin
            print,  Ai, alpha_i, beta_i, $
                    chi_sq_min, chi_sq,  $
                    (Phi_star_norm+(Ai*Ai_step)), alpha_norm +(alpha_i*alpha_step), beta_norm +(beta_i*beta_step) 

            chi_sq_min = chi_sq 
            alpha_i_min = alpha_i
            beta_i_min  = beta_i
            Ai_min      = Ai
         endif
         
      endfor
   endfor
endfor
print
print
;print,  'alpha_i, beta_i, Ai, chi_sq_min, chi_sq,  (-9.57+(Ai*0.01)),  -3.00 +(beta_i*0.02), -4.00 +(alpha_i*0.02)'
print, 'Ai_min, alpha_i_min, beta_i_min, chi_sq_min, chi_sq, (Phi_star_norm+(Ai*Ai_step)), (alpha_norm+(alpha_i_min*alpha_step)), (beta_norm+(beta_i_min*beta_step)) '
print,  Ai_min, alpha_i_min, beta_i_min, chi_sq_min, chi_sq, (Phi_star_norm+(Ai*Ai_step)), (alpha_norm+(alpha_i_min*alpha_step)), (beta_norm+(beta_i_min*beta_step)) 
print
print

close, 10
close, 11
close, /all

;; Current best fit....
oplot, mag_PLE, alog10(Phi_PLE[*,Ai_min,alpha_i_min,beta_i_min]), thick=4





end
