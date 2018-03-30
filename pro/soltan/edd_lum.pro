

;; http://en.wikipedia.org/wiki/Eddington_luminosity

four_pi     = 4.*!dpi
G           =  6.67300e-11     ;;m^3 kg^-1 s^-2
mass_proton =  1.67262158e-27 ;;kilograms
c           = 3e8
sigma_T     = 6.6524586e-29  ;m2  6.6524586e-25 cm^2  =0.66524586\ldots~\textrm{barn} 
M_sol       = 2e30 ;m2  6.6524586e-25 cm^2  =0.66524586\ldots~\textrm{barn} 

;; Stefan–Boltzmann constant is given in SI by:[1]
;;
;;sigma = 5.670 400(40) \times 10^{-8}\ \textrm{W}\,\textrm{m}^{-2}\,\textrm{K}^{-4}. 
;;sigma = \frac{2\pi^5 k^4}{15c^2h^3}= 5.670 400 \times 10^{-8}\, \mathrm{J\, s^{-1}m^{-2}K^{-4}}, 


print, ' 4\pi   ', four_pi
print, '   G    ',  G
print, ' m_p    ', mass_proton
print, ' c      ', c
print, 'sigma_T ',6.6524586e-29

L_edd = (four_pi*G*mass_proton*c)/(sigma_T) 

print, ' L_edd =  ',  (four_pi*G*mass_proton*c)/(sigma_T) ,         '   M (in kg)'
print, ' L_edd =  ',  ((four_pi*G*mass_proton*c)/(sigma_T))*M_sol , '   M_sol'
print
print

end
