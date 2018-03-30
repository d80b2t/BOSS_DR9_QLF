
;; Fair bit of ``inheritance'' from: /cos_pc19a_npr/BOSS/Stripe82/BOSS_QSOs_ontheStripe.pro

;; The DR9 dataset...
data_full = mrdfits('/cos_pc19a_npr/BOSS/bossqsomask/trunk/data/spall/mini-spAll-v5_4_45.fits',1)
data = data_full


;; The XD-CORE targets from Stripe 82..
xdcore = mrdfits('/cos_pc19a_npr/BOSS/bossqsomask/trunk/data/xdcore/xdcore_Stripe82.fits',1)



target_flag = 34084861508608LL
high_z = 2.15

s82_qso_targs_specPrim = where(data.ra ge 300.0 or data.ra lt 60.0 and $ 
                               (data.dec le 1.25 and data.dec ge -1.25) $
                               and ((data.boss_target1 and target_flag) ne 0) and $
                               data.specprimary eq 1, N_s82_qso_targs_specPrim)

tenk = data[s82_qso_targs_specPrim]


spherematch, tenk.ra, tenk.dec, xdcore.ra, xdcore.dec, 2./3600., m10, m11, distxx 
;; REALLY CUTE WAY TO PICK OUT THOSE OBJECTS THAT *ARE* COINCIDENT
;; BETWEEN LISTS, AND THOSE THAT ARE *NOT* (!!!!) THANKS TO EHUFF!!

;; m10 => The 10k specPrim=1 QSO_targs on S82, crossed with the XDCORE objects...
res = histogram(m10)                                                                                    
;; XD and !XD... 
xd    = where(res eq 1 and tenk.zWarning eq 0 and tenk.z ge 2.20, N_xd)
notxd = where(res eq 0 and tenk.zWarning eq 0 and tenk.z ge 2.20, N_notxd)


end
