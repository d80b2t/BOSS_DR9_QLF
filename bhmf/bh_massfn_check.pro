


;; 
;; Shankar_2009_ApJ_690_20 Table 3.
;;

readcol, 'Shankar_2009_ApJ_2009_690_20_Table3.dat', zShankar, logM_Shankar, $
         logP, logPnoCT, logPDCT, logPHRH, logPExt

;;
;; Output form MW's bh_massfn.py....
;;

readcol, 'bhmf_atzeq0.dat',   MBH0,          Phi0
readcol, 'bhmf_atzeq2.dat',   MBH2,          Phi2    

;; for z>3     
;;   k1,k2,a = 1.22,-0.23,-3.31+0.5*(zz-3.0)
readcol, 'bhmf_atzeq6.dat',       MBH6,          Phi6    
;;   k1,k2,a = 1.22,-0.23,-3.31+0.1*(zz-3.0)
readcol, 'bhmf_atzeq6_mod.dat',   MBH6_mod,          Phi6_mod    

!p.multi=[0,1,2]

thick     = 4.0
charsize  = 2.6
charthick = 4.0

;;
;;  z =  0.0  and  2.0
;;
plot, mbh2, phi2, $
      position=[0.12, 0.16, 0.52, 0.96], $
      thick=thick, $
      /xlog, /ylog, $
      xtitle='!3M!IBH!N', $
      ytitle='!3log !4u!3(M!IBH!N) [Mpc/h]!E3!N dex!E-1!N', $
      charsize=charsize, $
      xrange=[1e6, 1e10], $
      yrange=[1e-6, 1e-1], $
      xstyle=1, ystyle=1

;; redshift z=2.0 Shankar...
shankar2 = where(zshankar eq 2.000000, N)
oplot, 10^(LOGM_SHANKAR[shankar2]), 10^(LOGP[shankar2]), ps=2, color=200
oplot, 10^(LOGM_SHANKAR[shankar2]), 10^(LOGP[shankar2]), thick=thick, color=200


;; redshift z=0.0
oplot,  MBH0,          Phi0, thick=thick ;, color=100
;; redshift z=0.0 Shankar...
shankar0 = where(zshankar eq 0.0200000, N)
oplot, 10^(LOGM_SHANKAR[shankar0]), 10^(LOGP[shankar0]), ps=2, color=200
oplot, 10^(LOGM_SHANKAR[shankar0]), 10^(LOGP[shankar0]), thick=thick, color=200

xyouts, 2e7, 1e-2, 'Shankar09', charsize=charsize, charthick=charthick, color=200
xyouts, 4e6, 1e-2, 'z=0.0 ', charsize=charsize, charthick=charthick, color=200
;xyouts, 4e6, 3e-3, 'z=2.0', charsize=charsize, charthick=charthick, color=200

xyouts, 4e6, 3e-5, 'z=2.0', charsize=charsize, charthick=charthick

;;
;;  z =  6.0
;;
plot, mbh6, phi6, $
      position=[0.58, 0.16, 0.96, 0.96], $
      thick=4, $
      /xlog, /ylog, $
      charsize=charsize, $
      xrange=[1e7, 1e10], $
      yrange=[1e-10, 1e-4], $
      xtitle='!3M!IBH!N', $
;      ytitle='!3log !4u!3(M!IBH!N) [Mpc/h]!E3!N dex!E-1!N', $
      xstyle=1,       ystyle=1
oplot,  MBH6_mod,          Phi6_mod, thick=thick, color=100

shankar6 = where(zshankar eq 5.99000, N)
oplot, 10^(LOGM_SHANKAR[shankar6]), 10^(LOGP[shankar6]), ps=2, color=200
oplot, 10^(LOGM_SHANKAR[shankar6]), 10^(LOGP[shankar6]), thick=thick, color=200

xyouts, 2e8, 1e-5, 'z=6.0, Shankar09', charsize=charsize, charthick=charthick, color=200
xyouts, 2e7, 1e-7, 'z=6.0', charsize=charsize, charthick=charthick
;xyouts, 1e8, 1e-6, 'z=6.0', charsize=charsize, charthick=charthick

close, /all

!p.multi=0

end




