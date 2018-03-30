
little_h = 1.0
red, omega0=0.30, omegalambda=0.70, h100=little_h 

z_grid = INDGEN(2200)*0.0001
dl = DLUMINOSITY(z_grid)/1e6                   
M_lim =  17.77 - (5.*alog10(dl))-25.           

dl  = DLUMINOSITY(0.085)/1e6               
Mlim =  17.77 - (5.*alog10(dl))-25.        
print, Mlim                                

charsize  = 1.8
charthick = 1.6

window, 0, xs=650, ys=650

plot, z_grid, M_lim, yrange=[-17,-24], xrange=[0.00,0.22], xstyl=1,ystyle=1, $
      charsize=charsize, charthick=charthick

xyouts, 0.08, -19.0, '!6h=',            charsize=charsize, charthick=charthick
xyouts, 0.10, -19.0, little_h,          charsize=charsize, charthick=charthick
;xyouts, 0.12, -18.0, 'M_lim(z=0.085)=', charsize=charsize, charthick=charthick

little_h = 0.70
red, omega0=0.30, omegalambda=0.70, h100=.70                               
dl = DLUMINOSITY(z_grid)/1e6                                               
M_lim =  17.77 - (5.*alog10(dl))-25.                                       

dl  = DLUMINOSITY(0.085)/1e6               
Mlim =  17.77 - (5.*alog10(dl))-25.        
print, Mlim                                

oplot, z_grid, M_lim, color=222                                             

xyouts, 0.08, -18.0, '!6h=',  color=222,          charsize=charsize, charthick=charthick
xyouts, 0.10, -18.0, little_h,color=222,          charsize=charsize, charthick=charthick

print 
print 

GZ_data = mrdfits('/cos_pc19a_npr/programs/GZ/data/dr6_zoo_for_Nic.fits',1)

oplot, GZ_data.redshift, GZ_data.MR, ps=3, color=222
oplot, GZ_data.redshift, GZ_data.MR-(5*alog10(little_h)), ps=3 

print 
print, 'little h', little_h 
print 
print 

end
