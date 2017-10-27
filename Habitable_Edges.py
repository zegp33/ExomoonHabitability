# Compute the habitable edge  (HE & Io) for an exomoon around a Binary star.
# by Zusi , zusi.gonzalez@udea.edu.co, Medellin, Colombia. Based in Heller&Zuluaga code.
# created 2017-09-19, last modification 2017-09-19
#-*- coding:utf-8 -*-
try:
 	from pylab import *
	pylabload = 'y'
except ImportError:
	print "pylab could not be imported."
	pylabload = 'n'

if pylabload == 'n':
	try:
		from numpy import *
	except ImportError:
		print "numpy could not be imported."

# 1) *** Constants
G        = 6.673 * 10**(-11.)                              # [m^3 / kg / s^2], Newton's Gravitational constant
h        = 6.626068 * 10.**(-34.)                          # [m^2 * kg / s], Planck constant
k_B      = 1.3806503 * 10.**(-23.)                         # [m^2 * kg / s^2 / K], Boltzmann constant
c        = 299792458.                                      # [m / s], speed of light
sigma_SB = 2. * pi**5. * k_B**4. / (15. * c**2. * h**3.)   # [kg / s^3 / K^4], Stefan-Boltzmann constant
AU       = 149.598 * 10**9.                                # [m], Astronomical Unit
R_sun    = 6.96 * 10**8.                                   # [m], solar radius
T_sun    = 5778.                                           # [K], solar surface or effective temperature
L_sun    = 4. * pi * R_sun**2. * sigma_SB * T_sun**4.      # [W], solar luminosity
R_jup    = (71492000. + 66854000.) / 2.                    # [m], mean of Jupiter's equatorial and polar radius
R_ear    = 6378. * 1000.                                   # [m], Earth radius
M_jup    = 1.8986 * 10**27.                                # [kg], Jupiter mass
M_ear    = 5.9736 * 10**24.                                # [kg], Earth mass
M_mar    = 6.4185 * 10.**23.                               # [kg], Mars' mass



# 2) *** Basic Parameters

#Star A  (4.4Gyr)
T_s=6269.50
R_s=1.682e+00 * R_sun
L_s=1.866e+00 * L_sun

#Star B (4.4 Gyr)
Ts2=5905.9
R_s2=1.012e+00 * R_sun
L_s2=1.085e+00 * L_sun

# Planet Kepler 1647b (Kostov et. al, 2016) & BHMcalc2 (Zuluaga,Mason, Cuartas)
alpha_p = 0.3
e_sp = 0.0581 
a_sp = 2.705 * AU   
Q=0.00744973882642177
Q_mean=6.54e19
Rho=1.00
Rho_mean=1.03e+03
M_p=0.983 *M_jup 
R_p=1.07994e+00 *R_jup

#BHZ Narrow (BHMcalc2 data)
inneredge = 2.071 * AU
outeredge = 3.610 * AU

# Hypotetical Exomoon
alpha_m_OPT = 0.4
alpha_m_IR = 0.05
tau_m = 638. #Tiempo de maxwell: parametro de marea, escala de tiempo de disipacion# ojo no es t maxwell es delta
k_2m = 0.3  #Love number according to a Ganymede with an ocean (20 km) "Tides in Astronomy & Astrophysics"
e_pm=0.01	#Excentricidad cte de la luna
A_pm = arange(0.9,30.001,0.05)    # continuous range to be tested for the Io and RG limits
M_m = 0.1 * M_ear                  # Mass moon
rmf = 0.42                       # with a rock-to-mass fraction of 68%
#imf = 0.42

#Radio de la luna (Fortney,2007) Eq (7) y (8)
try:
    R_m = ( (0.0912*imf + 0.1603)*(log10(M_m/M_ear))**2. + (0.3330*imf + 0.7387)*log10(M_m/M_ear)\
    + (0.4639*imf + 1.1193) ) * R_ear # [m], the moon's radius from Fortney et al. (2007), Eq. (7)
except NameError:
    R_m = ( (0.0592*rmf + 0.0975)*(log10(M_m/M_ear))**2. + (0.2337*rmf + 0.4938)*log10(M_m/M_ear)\
    + (0.3102*rmf + 0.7932) ) * R_ear # [m], the moon's radius from Fortney et al. (2007), Eq. (8)


# Runaway greenhouse 
a = 0.7344 
l = 2.425e6             # J/kg
r = 461.5               # J/kg K
p0 = 1e4                # Pa
tref = 273.13           # K
pref = 610.616          # Pa
k0 = 0.055
pstar = pref*exp(l/(r*tref))
kappa = 0.01            # m^2/kg
g_m = G * M_m/(R_m**2)
F_crit = a * sigma_SB * (  l / (r * log(pstar / sqrt(2.*p0*g_m/k0)  ) ) )**4.


# 3) *** Derived parameters


#Average flux of binary [W/m**2]
S_bin=(L_s+L_s2)/(2.*4.*pi*a_sp**2.)


#Equilibrium temperature of the planet [K]
Teq_p=(S_bin*(1-alpha_p)/sigma_SB)**(1./4)


# Internal temperature of the planet.

Tint_p=(((Q*Q_mean)/(4.*pi*sigma_SB))*((4.*pi*Rho*Rho_mean)/(3.*M_p))**(2./3.))**(1./4.)

#Temperatura del planeta
T_p = (Teq_p**4  + Tint_p**4)**(1./4) 

#Archivos de salida
file_out = open("Fglob_a_e_cte.dat", 'w')
file_out2 = open("h_m_a_e_cte.dat", 'w')
file_out3=open("fstar_Distance.dat",'w')
file_out4=open("frefle_planet_Distance.dat",'w')
file_out5=open("termal_energy_Distance.dat",'w')
file_out6=open("All_contibution_Fglob_e.dat",'w')


#***4)Calculos Habitable Edges
def HEexcentricidad():
    #1)Calculo de HE de acuerdo a la excentricidad
    file_out7 = open("HE_IO_e.dat", 'w')
    E_pm = arange(0.001, 0.1,0.001)	#Rango de Excentricidades
    #E_pm = array([0.001, 0.01,0.1]) #Para 3 excentrcidades dadas
    for e_pm in E_pm:
            #print e_pm
            beta = sqrt(1.-e_pm**2)
            f1 = 1. +(31./2.*e_pm**2.) + (255./8.*e_pm**4) + (185./16.*e_pm**6.) + (25./64.*e_pm**8.)
            f2 = 1. + (15./2.*e_pm**2.) +  (45./8.*e_pm**4.) +   (5./16.*e_pm**6.)
            f5 = 1. +  (3.*e_pm**2.)+  (3./8.*e_pm**4.)
            a_pm = A_pm * R_p
            Z_m = 3.*(G*G)*k_2m*(M_p*M_p)*(M_p+M_m)*((R_m**5.)/(a_pm**9.))*tau_m
            h_m = (Z_m/(beta**15.)*(f1-(f2*f2))/(f5))/(4.*pi*(R_m*R_m))
            F_glob = ((L_s+L_s2)/2.)*(1.-alpha_m_OPT)/(16.*pi*a_sp**2.*sqrt(1.-e_sp**2))*(1.+pi*R_p**2.*alpha_p/(2.*a_pm**2.))+R_p**2.*sigma_SB*T_p**4./a_pm**2.*(1.-alpha_m_IR)/4.+h_m

            # DA BRAIN a_HE = a_pm[where(F_glob < F_crit)[0][0]] a_Io
            = a_pm[where(h_m < 2)[0][0]] file_out7.write("%.5e %.5e
            %.3e \n" % (a_HE/R_p,a_Io/R_p,e_pm)) file_out7.write("\n")
            file_out7.close ()

"""
#2)HE Variando la masa del planeta y el radio.
MP=arange(200.0,1000.0,50.)*M_ear
file_out8 = open("HE_IO_Mp.dat", 'w')

for M_p in MP:
	#Radio del planeta
	R_p=((-8.935397162317*10**-7.*(M_p/M_ear)**2.)+(0.000732550054192*M_p/M_ear)+ 0.9195139359476)*R_jup #[m] Se obtiene de ajuste de datos de Tabla 4 de Fortney 2007. Para un plneta a 1 AU en 4.5 Gy.
	#R_p=1.07994e+00 *R_jup			#Radio K1647
	#print M_p/M_ear, R_p/R_jup
	# Internal temperature of the planet.
	Tint_p=(((Q*Q_mean)/(4.*pi*sigma_SB))*((4.*pi*Rho*Rho_mean)/(3.*M_p))**(2./3.))**(1./4.)
	#Temperatura del planeta
	T_p = (Teq_p**4  + Tint_p**4)**(1./4) 
	#HE & Io
	#Exomoon Habitability (Heller&Barnes,2013, pag 14)
	#print e_pm
	beta = sqrt(1.-e_pm**2)
	f1 = 1. +(31./2.*e_pm**2.) + (255./8.*e_pm**4) + (185./16.*e_pm**6.) + (25./64.*e_pm**8.)
	f2 = 1. + (15./2.*e_pm**2.) +  (45./8.*e_pm**4.) +   (5./16.*e_pm**6.)
	f5 = 1. +  (3.*e_pm**2.)+  (3./8.*e_pm**4.)
	a_pm = A_pm * R_p
	Z_m = 3.*(G*G)*k_2m*(M_p*M_p)*(M_p+M_m)*((R_m**5.)/(a_pm**9.))*tau_m
	h_m = (Z_m/(beta**15.)*(f1-(f2*f2))/(f5))/(4.*pi*(R_m*R_m))
	F_glob = ((L_s+L_s2)/2.)*(1.-alpha_m_OPT)/(16.*pi*a_sp**2.*sqrt(1.-e_sp**2))*(1.+pi*R_p**2.*alpha_p/(2.*a_pm**2.))+R_p**2.*sigma_SB*T_p**4./a_pm**2.*(1.-alpha_m_IR)/4.+h_m  

	# DA BRAIN
	a_HE = a_pm[where(F_glob < F_crit)[0][0]]
	a_Io = a_pm[where(h_m < 2)[0][0]]
	file_out8.write("%.5e	%.5e	%.3e \n" % (a_HE/R_jup,a_Io/R_jup,M_p/M_jup))
	file_out8.write("\n")
file_out8.close()

"""
#2)Contribuciones de Flujo Global
for a_pm in A_pm:
	beta = sqrt(1.-e_pm**2)
	f1 = 1. +(31./2.*e_pm**2.) + (255./8.*e_pm**4) + (185./16.*e_pm**6.) + (25./64.*e_pm**8.)
	f2 = 1. + (15./2.*e_pm**2.) +  (45./8.*e_pm**4.) +   (5./16.*e_pm**6.)
	f5 = 1. +  (3.*e_pm**2.)+  (3./8.*e_pm**4.)
	a_pm = A_pm * R_p
	Z_m = 3.*(G*G)*k_2m*(M_p*M_p)*(M_p+M_m)*((R_m**5.)/(a_pm**9.))*tau_m
	h_m = (Z_m/(beta**15.)*(f1-(f2*f2))/(f5))/(4.*pi*(R_m*R_m))
	F_glob = ((L_s+L_s2)/2.)*(1.-alpha_m_OPT)/(16.*pi*a_sp**2.*sqrt(1.-e_sp**2))*(1.+pi*R_p**2.*alpha_p/(2.*a_pm**2.))+R_p**2.*sigma_SB*T_p**4./a_pm**2.*(1.-alpha_m_IR)/4.+h_m 
	fstar=((L_s+L_s2)/2.)*(1.-alpha_m_OPT)/(16.*pi*a_sp**2.*sqrt(1.-e_sp**2))
	f_rp=((L_s+L_s2)/2.)*(1.-alpha_m_OPT)/(16.*pi*a_sp**2.*sqrt(1.-e_sp**2))*(pi*R_p**2.*alpha_p/(2.*a_pm**2.))
	f_t=R_p**2.*sigma_SB*T_p**4./a_pm**2.*(1.-alpha_m_IR)/4.
	# DA BRAIN
	a_HE = a_pm[where(F_glob < F_crit)[0][0]]
	a_Io = a_pm[where(h_m < 2)[0][0]]
	print a_HE/R_p,a_Io/R_p,F_crit
	
	file_out6.write("#F_glob		h_m[i]		 fstar		 f_rp		 f_t		a_pm/R_p\n")
	
	for i in range(500):
		file_out.write("%.5e	%.5e " % (F_glob[i],a_pm[i]/R_p))
		file_out.write("\n")
		file_out2.write("%.5e	%.5e " % (h_m[i],a_pm[i]/R_p))
		file_out2.write("\n")
		file_out3.write("%.5e	%.5e " % (fstar,a_pm[i]/R_p))
		file_out3.write("\n")
		file_out4.write("%.5e	%.5e " % (f_rp[i],a_pm[i]/R_p))
		file_out4.write("\n")
		file_out5.write("%.5e	%.5e " % (f_t[i],a_pm[i]/R_p))
		file_out5.write("\n")
		file_out6.write("%.5e	%.5e	%.5e	%.5e	%.5e	%.5e" % (F_glob[i], h_m[i], fstar, f_rp[i], f_t[i],a_pm[i]/R_p))
		file_out6.write("\n")
	exit(0)
	

file_out.close()
file_out2.close()
file_out3.close()
file_out4.close()
file_out5.close()
file_out6.close()

"""
print "\n\t\t*******Inicio\n\n\n\tPara una luna de una excentricidad de",e_pm,"y un radio de", R_m/R_ear,"R_ear sus Habitable Edges en radios planetarios son\n 1) HE", a_HE/R_p ,"\n 2) Io",a_Io/R_p

print "\t\t***Datos adicionales:\n"

print "\t*Temperatura interna del planeta [K]:\t",Tint_p
print "\n\t*Temperatura de equilibrio [K]:\t",Teq_p
print "\n\t*Albedo del planeta:\t", alpha_p
print "\n\t*Flujo promedio de la binaria (instantaneo) [W/m**2]:\t",S_bin
print "\n\t\t*******Fin\n\t\t*******El programa se ha ejecutado con exito.\n\n"
"""
print
