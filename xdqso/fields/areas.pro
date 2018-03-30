;;
;
;
; Data created from:
;   
;SELECT fieldID, run,  rerun, camcol, field, nTotal, nObjects, nStars,
;quality, raMin, raMax, decMin, decMax, primaryArea 
;FROM field
;WHERE primaryArea > 0.00 and rerun = 301

data = mrdfits('fields.fits',1)

min_run = min(data.run)
max_run = max(data.run)

sum_area = 0.0
openw, 10, 'SDSS_DR8_field_areas.dat'
printf, 10, '# SDSS RUN,  RUN Area,  Cumulative Area  '
for ii=min_run,max_run do begin

   w = where(data.run eq ii, N)
   
   if N gt 0 then begin
      sum_area =  sum_area+        total(data[w].Primaryarea)
      print, ii,  total(data[w].Primaryarea), sum_area   
      printf, 10, ii,  total(data[w].Primaryarea), sum_area, $
              format='(i6, f16.8, f16.8)'
   endif

endfor
close, 10
close, /all

end
