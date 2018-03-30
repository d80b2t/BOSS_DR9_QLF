

;;
;; 2SLAQ data from Scott
;;
;; Croom et al. (2009b), 399, 1755, Fig. 3. 
;;
;; 0.4 < z < 2.1 and Mg < −22.5 (or MbJ	< −22.5). 
;; The 2SLAQ points (filled circles) have been calculated 
;; using a z = 0 power-law K-correction with αν = −0.5 to match R05. 
;;
readcol, '../../pro/number_counts/nm_all_mg-22.5alphakcorr.out', $
         gmag_2SLAQ, n_2SLAQ_gmag,   dn_2SLAQ_gmag,  n_q_2SLAQ_gmag

;;
;; My 2SLAQ QSO (i-band) calculations...
;;
readcol, '../../pro/number_counts/2SLAQ_QSOs_iband_numcounts.dat', $
         bright_mag, faint_mag, imag_2SLAQ, no_count_bin_2SLAQ, n_2slaq_imag



;;
;; Richards et al. (2006)
;;
;; Quasar Number Counts (0.3 < z < 2.2)
;;
readcol, '../../pro/number_counts/Richards06_Table2_0pnt3z2pnt2.dat', $
         SDSS_mag, $
         SDSS_N_g, SDSS_N_gerr, SDSS_N_lt_g, SDSS_N_lt_gerr, N_Q_gmag, $
         SDSS_N_i, SDSS_N_ierr, SDSS_N_lt_i, SDSS_N_lt_ierr, N_Q_imag 

;; Richards et al. (2006)
;; Quasar Number Counts (3 < z < 5)
;;
readcol, '../../pro/number_counts/Richards06_Table3_3z5.dat', $
         SDSS_mag_hiz, $
         SDSS_N_i_hiz, SDSS_N_ierr_hiz, SDSS_N_lt_i_hiz, SDSS_N_lt_ierr_hiz, $
         N_Q_imag_hiz, N_Q_corr_imag_hiz

;;
;;  D A T A    F O R   B O S S   D R 9  
;;
;readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_Stripe82_numcounts.dat', $
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_Stripe82_iband_numcounts.dat', $
         bright_mag, faint_mag, BOSS_imag, BOSS_N_Q_S82, BOSS_num_count_S82

;;
;;  D A T A    F O R    S T R I P E  8 2
;;
;readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_Stripe82_numcounts.dat', $
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_DR9_Stripe82_iband_numcounts_temp.dat', $
         bright_mag, faint_mag, BOSS21_imag, BOSS21_N_Q_S82, BOSS21_num_count_iband

;;
;;  D A T A    F O R    ``b o s s  2 1 ''
;;
;readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_Stripe82_numcounts.dat', $
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_boss21_gband_numcounts.dat', $
         bright_mag, faint_mag, BOSS21_gmag, BOSS21_N_Q_S82, BOSS21_num_count_gband
readcol, '/cos_pc19a_npr/BOSS/QLF/pro/number_counts/BOSS_boss21_iband_numcounts.dat', $
         bright_mag, faint_mag, BOSS21_imag, BOSS21_N_Q_S82, BOSS21_num_count_iband




;; Colour Table
clr_table =13
loadct, clr_table

;; Colours for clr_table =13
black      =   0
purple     =  32
deep_blue  =  48
blue       =  64
light_blue =  80
turquiose  = 128
green      = 150
yellow     = 210
orange     = 232
red        = 254

charsize=2.6
charthick=4.8
thick=3.8
xthick=4
ythick=4
XTICKLEN  = 0.05
YTICKLEN  = 0.05



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  g-band N(M)'s
;;

set_plot, 'ps'
device, filename='Croom09_2SLAQ_Fig3_gmag_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

;; x-ranges
xmin = 18.0  ;; 18.0 in Croom09b
xmax = 24.0  ;; 22.0  in Croom09b

;; y-ranges
ymin = 0.30001 ;; 0.30001  in Croom09b
ymax = 19.9    ;; 19.9  in Croom09b

plotsym, 0,1.2, /fill
plot, gmag_2SLAQ, n_2SLAQ_gmag, $
      /ylog, $
      ps=2, $
      position=[0.24, 0.22, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      ystyle=1, xstyle=1, $
      xthick=xthick, ythick=ythick, $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
      /nodata, $
      xtitle=' g ', $
      ytitle='N(g) (deg!E-2!N [0.25 mag]!E-1!N)'

oplot, gmag_2SLAQ, n_2SLAQ_gmag, ps=8, color=turquiose+16.

oplot,  BOSS21_gmag, BOSS21_num_count_gband, $
       ps=8,$
       color=blue+8.

;   top = 10^(0.)
middle = 10^(0.0)
bottom = 10^(-0.2)

charthick=charthick*1.6
;xyouts, 20., top,    '2SLAQ 0.4<z<2.1', charsize=charsize, charthick=charthick, color=black
xyouts, 20., middle, '2SLAQ 0.4<z<2.1  ', charsize=charsize, charthick=charthick, color=turquiose+16.
xyouts, 20., bottom, 'boss21 0.4<z<2.1', charsize=charsize, charthick=charthick, color=blue+8.
charthick=charthick/1.6



device, /close
set_plot, 'X'     





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  i-band N(M)'s
;;
set_plot, 'ps'
device, filename='Richards_SDSS_Fig3_imag_temp.eps', $
        xsize=8.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated

;; x-ranges
xmin = 15.01   ;; 15.01 is nice...
xmax = 23.99   ;; 23.49 is okay

;; y-range
ymin = 0.00009   ;; 0.0002 is okay...
ymax = 40.0     ;; 10.0 is okay w/o boss21...

plotsym, 0,1, /fill
plot, SDSS_mag,  SDSS_N_i, $
      /ylog, $
      ps=2, $
      position=[0.26, 0.18, 0.96, 0.96], $
      xrange=[xmin, xmax], $
      yrange=[ymin, ymax], $
      ystyle=1, xstyle=1, $
      xthick=xthick, ythick=ythick, $
      YTICKNAME = ['10!E-4!N','0.001', '0.01', '0.1', '1', '10'], $
;      YTICKNAME = ['10', '1', '0.1', '0.01', '0.001'], $
;      YTICKS=5., $
      yminor=10., $
;      YTICKFORMAT='(i2, i1, f3.2, f4.3, f5.4)', $
      XTICKLEN=XTICKLEN, YTICKLEN=YTICKLEN, $
      charsize=charsize, $
      charthick=charthick, $
      thick=thick,$ 
;      /nodata, $
      xtitle=' i ', $
      ytitle='N(i) (deg!E-2!N [0.25 mag]!E-1!N)'

oplot, SDSS_mag_hiz,  SDSS_N_i_hiz, $
       ps=2,$
       color=red


;----------------------------------------------------
;  2 S L A Q   Q S Os 
;  Based on Croom09b, g-band and 0.3<z<2.1
;
;; Assuming (g-i) = 0.15ish...
;oplot, (gmag_2SLAQ-0.15), n_2SLAQ_gmag, ps=2,     color=green

;----------------------------------------------------
;  2 S L A Q   Q S Os 
;  Based on NPR, i-band, 0.3<z<2.20
;
oplot, imag_2SLAQ, n_2SLAQ_imag*4., ps=2,     color=green


oplot,  BOSS_imag, BOSS_num_count_S82, $
       ps=2,$
       color=deep_blue-8.

oplot,  BOSS21_imag, BOSS21_num_count_iband, $
       ps=8,$
       color=light_blue+8.

   top    = 10^(-2.5)
middle    = 10^(-2.9)
bottom    = 10^(-3.3)
below_btm = 10^(-3.7)

charthick=charthick*1.6
xyouts, 18., top,    'SDSS 0.3<z<2.2', charsize=charsize, charthick=charthick, color=black
xyouts, 18., middle, 'SDSS   3<z<5  ', charsize=charsize, charthick=charthick, color=red
xyouts, 18., bottom, 'BOSS 2.2<z<3.5', charsize=charsize, charthick=charthick, color=deep_blue-8.
xyouts, 18., below_btm, 'boss21 0.3<z<2.2', charsize=charsize, charthick=charthick, color=light_blue+8.

charthick=charthick/1.6

device, /close
set_plot, 'X'     

end
