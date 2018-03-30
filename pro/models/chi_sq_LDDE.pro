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



red, omega0=0.30, omegalambda=0.70, h100=0.700


;;  THE   S T R I P E   8 2   DATA  - 6,250 Quasars strong
readcol, '../qlf/My_QLF_iband_boss_narrowZ_S82.dat', $  ;; USE THIS FOR THE S82 plots!!!
         z_bin_boss_s82, Abs_mag_bin_boss_s82, blah_boss_s82, log_Num_den_boss_s82, raw_N_boss_QSOs_s82, sigma_Phi_BOSS_S82

log_Phi_BOSS_S82 =  log_Num_den_boss_S82

boss_delta_up_S82   = alog10 ((10^log_Phi_BOSS_S82 + sigma_Phi_BOSS_S82))
boss_delta_up_S82   = boss_delta_up_S82       - log_Phi_BOSS_S82
boss_delta_down_S82 = alog10 ((10^log_Phi_BOSS_S82 - sigma_Phi_BOSS_S82))
boss_delta_down_S82 = abs(boss_delta_down_S82 - log_Phi_BOSS_S82)


;boss_ilimit = 21.80
;;red_bins = [2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 3.00, 3.25, 3.50, 4.00]
;red_bins = [2.40, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.90, 3.13, 3.35, 3.75, 4.00]
;dlums = DLUMINOSITY(red_bins) / 1e6
;Abs_iMag_limit = boss_ilimit - (5 * alog10(dlums))  - 25.00 ;+ kcor(fix(red_bins/0.01))
Abs_iMag_limit =[-24.640, -24.470, -24.585, -24.695, -24.800, -24.901,-24.998, -25.137, -25.337, -25.514, -25.806, -25.973]


w = where(z_bin_boss_s82 ge 2.40 and z_bin_boss_s82 le 2.50 $
          and raw_N_boss_QSOs_s82 ge 1 and  $
          Abs_mag_bin_boss_s82 le Abs_iMag_limit[3] , N_qlf_bins)


;; No. of redshift bins
zbins = 140. 
;; with 100 bins and / 10 with, 100 bins => 0.00<z<9.90 in  delta_z=0.10
;; with 100 bins and / 20 with, 100 bins => 0.00<z<4.95 in  delta_z=0.05
zz = findgen(zbins)/20.    ;;  Redshift bins..

;; No. of mag bins
mag_bins = 30
mag_LDDE = fltarr(mag_bins)

;; A, the normalization
A_norm_bins = 50
A = fltarr(A_norm_bins)

;; alpha, the (bright-end in C09b) slope
alpha_bins = 100
alpha = fltarr(alpha_bins)

;; A, the normalization
beta_bins = 100
beta = fltarr(beta_bins)

;; Set up phi...
;Phi_LDDE = fltarr(zbins, mag_bins)
Phi_LDDE = fltarr(mag_bins, A_norm_bins, alpha_bins, beta_bins)


;---------------------------------------------------------------
; 
;  L  D  D  E
;   Luminosity-dependent density evolution 
;   Section 6.2 from Croom et al. (2009b)

;   alpha =  -3.70
;    beta =  -2.34
M_star_g = -26.69 
  M_g_c  = -23.90
  gamma  =   0.68 
zc_zero  =   2.47
      p1 =   6.28
      p2 =  -2.85
;       A = 10^(-9.27) ;;        A = 10^(-9.27)  is our default... 

e_d = fltarr(zbins, mag_bins)

openw, 9, 'Stripe82_LDDE_temp.dat'

;
; "Modified Luminosity-dependent density evolution" (mLDDE)
;  where the modification is to the behaviour of the 
;  "e_d" parameter... 
;

ii =49  ;; picking out zz[ii]=2.45
;for ii=0L, zbins-1 do begin

for alpha_i = 0LL, alpha_bins-1 do begin
   ;alpha_i = 0
   ;alpha[alpha_i] = -3.70
   alpha[alpha_i] = -4.00 +(alpha_i*0.02)

   
   for beta_i = 0LL, beta_bins-1 do begin
      beta[beta_i] = -3.00 +(beta_i*0.02)
      
      for Ai = 0L, A_norm_bins-1 do begin
         A[Ai] = 10^(-9.57+(Ai*0.01))
         
         for jj=0L, mag_bins-1 do begin
;      mag_LDDE[jj] = -32.0 + (jj*0.25)      
;   mag_LDDE[jj] = -32.0 + (jj*0.30)                       ;; to match the data's delta_mags
            mag_LDDE[jj] = min(Abs_mag_bin_boss_s82[w]) + (jj*0.30) ;;   and put them on the same "absolute scale" 
                                                                                
            dmag    = mag_LDDE[jj] - M_star_g 
            ;; N.B. This is a different M_star_g than the Mstar_g from
            ;; above!!! Grrr.....
            dmag_c  = mag_LDDE[jj] - M_g_c 
            
            ;; Equnation (18) of Croom et al. (2009b)
            zc = zc_zero / ( 1. + (10^(0.4*gamma*dmag_c)))
            
            ;; Equnation (17) of Croom et al. (2009b)
            top_ed    = 2.*(1+zc)^(p1)
            bottom_ed = (((1+zz[ii])/(1+zc))^(-1.*p1)) + (((1+zz[ii])/(1+zc))^(-1.*p2))
            e_d(ii,jj) = top_ed / bottom_ed 
            
;      Abs_Mag_model[jj] = mag
            
            ple_denom =  ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
            Phi_LDDE[jj, Ai, alpha_i, beta_i] = (A[Ai] * e_d[ii,jj])  / ( (10.0^(0.4*(alpha[alpha_i]+1)*(dmag))) + (10.0^(0.4*(beta[beta_i]+1)*(dmag))) )
            
            print, ii, jj, alpha_i, beta_i, Ai, zz[ii], alog10(A[Ai]), mag_LDDE[jj], top_ed, bottom_ed, e_d[ii, jj],  (A[Ai] * e_d[ii,jj]), ple_denom
;      printf, 9, ii, zz[ii], jj, mag_LDDE[jj], e_d[ii,jj], zc, Phi_LDDE[ii,jj], $
                                ;             format='(i5, f9.5,   i5, f16.7, f16.7,      f16.5, e20.8)'
         endfor
      endfor
   endfor ;; beta
endfor;; alpha
close, 9

plotsym, 0, 1.6, /fill
plot,  Abs_mag_bin_boss_s82[w] , log_Phi_BOSS_S82[w], psym=8 
Mi_boss_err_S82 = log_Num_den_boss_S82[w] - log_Num_den_boss_S82[w]
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_up_S82[w],  $
            /hibar, errthick=3., psym=8;, color=160, ERRcolor=160
oploterror, Abs_Mag_bin_boss_S82[w], log_Phi_boss_S82[w], Mi_boss_err_S82, boss_delta_down_S82[w], $   
            /lobar, errthick=3., psym=8;, color=160, errcolor=160

oplot, mag_LDDE, alog10(Phi_LDDE[*,40,0, 0]), color=40, thick=2
oplot, mag_LDDE, alog10(Phi_LDDE[*,40,0, 40]), color=60, thick=2
oplot, mag_LDDE, alog10(Phi_LDDE[*,40,0,99]), color=80, thick=2
oplot, mag_LDDE, alog10(Phi_LDDE[*,40,99,0]), color=120, thick=2
oplot, mag_LDDE, alog10(Phi_LDDE[*,40,99,99]), color=160, thick=2
;oplot, mag_LDDE, alog10(Phi_LDDE[*,40,0,80]), color=200, thick=2
;oplot, mag_LDDE, alog10(Phi_LDDE[*,40,0,99]), color=240, thick=2

openw, 11, 'chi_sq_values_LDDE_inside_mag_loop_temp.dat'
printf, 11, '# Normalization, A  Abs i-mag     Phi_BOSS_S82, sigma_Phi_BOSS_S82,    Phi_LDDE_model,  chi_sq_bin, chi_sq '
printf, 11, '# (-9.57+(Ai*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]), sigma_Phi_BOSS_S82[w[kk]], Phi_LDDE[ii,kk, Ai], chi_sq_bin, chi_sq '

openw, 10, 'chi_sq_values_LDDE_temp.dat'
printf, 10, '# Normalization A,     chi_sq_bin, chi_sq '
printf, 10, '# (-9.57+(Ai*0.01)),   chi_sq_bin, chi_sq '


;; http://en.wikipedia.org/wiki/Goodness_of_fit
;chi_sq = fltarr(N_qlf_bins, A_norm_bins)
chi_sq_min = 1000.
for alpha_i = 0LL, alpha_bins-1 do begin
   for beta_i = 0LL, beta_bins-1 do begin
      for Ai = 0L, A_norm_bins-1 do begin
         
         
         chi_sq=0         
         for kk=0L, N_qlf_bins-1 do begin
                                ;chi_sq_bin = ( log_Phi_BOSS_S82[w[kk]] -  alog10(Phi_LDDE[ii,kk]))^2 /  (alog10(sigma_Phi_BOSS_S82[w[kk]])^2)
            chi_sq_bin = ( 10^(log_Phi_BOSS_S82[w[kk]]) -  Phi_LDDE[kk, Ai, alpha_i, beta_i])^2 $
                         /  (sigma_Phi_BOSS_S82[w[kk]]^2)
            chi_sq = chi_sq + chi_sq_bin
            
                                ;print, Ai, kk, Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]])*1e7, sigma_Phi_BOSS_S82[w[kk]]*1e7, $
                                ;       Phi_LDDE[ii,kk,Ai]*1e7, chi_sq_bin, chi_sq
                                ;print, Abs_mag_bin_boss_s82[w[kk]], log_Phi_BOSS_S82[w[kk]], alog10(sigma_Phi_BOSS_S82[w[kk]]), Phi_LDDE[ii,kk], $
                                ;      chi_sq_bin, chi_sq[kk]
            
            printf, 11, (-9.57+(Ai*0.01)), Abs_mag_bin_boss_s82[w[kk]], 10^(log_Phi_BOSS_S82[w[kk]]),  sigma_Phi_BOSS_S82[w[kk]], $
;                    Phi_LDDE = fltarr(mag_bins, A_norm_bins, 1, beta_bins)
                    Phi_LDDE[kk, Ai, alpha_i, beta_i], $
                    chi_sq_bin, chi_sq, $
                    format='(f9.4,2x, f9.4,2x, e,e,e, f, f)'
         endfor
         
         if chi_sq lt chi_sq_min then begin
            print,  alpha_i, beta_i, Ai, chi_sq_min, chi_sq,  (-9.57+(Ai*0.01)),  -3.00 +(beta_i*0.02), -4.00 +(alpha_i*0.02)
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
print,  'alpha_i_min, beta_i_min, Ai_min, chi_sq_min, chi_sq,  (-9.57+(Ai*0.01)),  -3.00 +(beta_i_min*0.02), -4.00 +(alpha_i_min*0.02)'
print,  alpha_i_min, beta_i_min, Ai_min, chi_sq_min, chi_sq,  (-9.57+(Ai_min*0.01)),  -3.00 +(beta_i_min*0.02), -4.00 +(alpha_i_min*0.02)
print
print

close, 10
close, 11
close, /all


oplot, mag_LDDE, alog10(Phi_LDDE[*,Ai_min,alpha_i_min,beta_i_min]), thick=4


end
