


   qsomaskdir=filepath(root=getenv("BOSSQSOMASK_DIR"), "data/")
   xdspallfile = qsomaskdir+'/compfiles/xdspallmask_DR9.fits'
   data = mrdfits(xdspallfile,1)

gcut       = where(data.PSFMAG[1] le 22.00, N_gcut)
rcut       = where(data.PSFMAG[2] le 21.85, N_rcut)
g_and_rcut = where(data.PSFMAG[1] le 22.0 OR data.PSFMAG[2] le 21.85,  N_g_and_rcut)
icut       = where(data.PSFMAG[3] le 21.80, N_icut)

print 
print, 'N_gcut      ', N_gcut
print, 'N_rcut      ', N_rcut
print, 'N_g_and_rcut', N_g_and_rcut
print, 'N_icut      ', N_icut
print
print

;data = data[icut]
;data = data


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


charsize  = 2.6
charthick = 8.8
thick     = 4.8
xthick    = 6.0
ythick    = 6.0
XTICKLEN  = 0.02
YTICKLEN  = 0.03


set_plot, 'ps'
device, filename='temp.ps', $
        xsize=10.0, ysize=8.0, $
        xoffset=0.2, yoffset=0.2, $
        /inches, /color, /encapsulated


mag_bin = 0.02

plothist, data.PSFMAG[1], bin=mag_bin, $
          position=[0.22,0.22,0.96,0.96], $
          color=green, $
          thick=thick,$
          xthick=xthick, ythick=ythick,XTICKLEN=XTICKLEN, YTICKLEN = YTICKLEN, $
          xrange=[17.5, 23.0], xstyle=1, $
;          yrange=[17.5, 24.0], ystyle=1, $
          charthick=charthick, charsize=charsize, $
          xtitle='PSF magnitude'

plothist, data.PSFMAG[2], bin=mag_bin, $
          color=red, thick=thick, /over
          
plothist, data.PSFMAG[3], bin=mag_bin, $
          color=purple, thick=thick, $
          /over


device, /close
set_plot, 'X'
close, /all


;; Plotting the relevant histograms...

w  = where(data.zWarning eq 0, N)             
ww = where(data.zWarning eq 0 and data.psfmag[3] le 21.80, NN)

plothist, data[w].z, bin=0.05, yrange=[0,4000], xrange=[0.5, 4.5]
plothist, data[ww].z, bin=0.05, color=200, /over  

end



