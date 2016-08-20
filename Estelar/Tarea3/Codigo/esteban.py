# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math as mt

#rho_c=3.17
# Constantes
G = 6.67*10**(-11) # Constante de gravitacion universal
c = 3*10**8 # Velocidad de la luz
m =  1.67*10**(-27) # Masa del neutron
M_sun = 1.989*10**(30) # Masa del Sol
rho_c = 2*10**17 # Densidad central de estrella de neutrones
h = 6.62*10**(-34) # Constante de Plank
n = float(rho_c/m) # Densidad de partículas
p_f = ((3*n/(8*np.pi))**(1.0/3))*h # Momentum de fermi
x = float(p_f/(m*c)) # Variable adimensional 
dr = 10.0 # 10 m de variacion
r = 1.0 # 1m radio inicial
P = ((np.pi*(m**4)*(c**5))/(3*h**3))*(np.sqrt((x**2)+1)*(2*(x**3)-3*x) + 3*mt.asinh(x))
e = ((np.pi*(m**4)*(c**5))/(h**3))*(np.sqrt((x**2)+1)*(2*(x**3)+x) - mt.asinh(x))
#e = 5.822*10**34
M = 0.0

mass = []
radius = []
pressure = []
edensity = []

while P > 0:
    
    mass.append(M)
    radius.append(r)
    pressure.append(P)
    edensity.append(e)
    
    dP_dr = -(G*((e+P)/(c**2))*(M + 4*mt.pi*(r**3)*P/(c**2)))/(r*(r - 2*G*M/(c**2)))
    dP_dx = (8*np.pi*(m**4)*(c**5)/(3*h**3))*((x**4)/(np.sqrt((x**2)+1)))
    dx = (dP_dr/dP_dx)*dr
    x += dx
    dm = ((4*np.pi*e*r**2)/(c**2))*dr
    M += dm
    P += (dP_dx)*dx
    de_dx = (8*np.pi*(m**4)*(c**5)/(h**3))*(((x**4)+(x**2))/(np.sqrt((x**2)+1)))
    e += (de_dx)*dx
    r += dr

Mass = np.array(mass)
Radius = np.array(radius)
Pressure = np.array(pressure)
e_density = np.array(edensity)

print Pressure

def politropa(dr):
    
    pmass = []
    ppressure = []
    pdensity = []
    pradius = []
    
    K = ((h**2)/(5*m))*(3.0/(8*np.pi))**(2.0/3)
    P = K*(rho_c/m)**(5.0/3)
    r = 1 # Radio inicial 
    p = rho_c
    M = 0
    
    while P > 0:
        
        dM = 4*np.pi*(r**2)*p*dr
        M += dM
        pmass.append(M)
        
        dP = -(G*M*p/(r**2))*dr
        P += dP
        
        if P < 0:
           ppressure.append(0)
           pdensity.append(0)
           pradius.append(r)
           
        else:
            
           ppressure.append(P) 
           p = ((P/K)*(m**(5.0/3)))**(3.0/5)
           pdensity.append(p)
           r += dr
           pradius.append(r)  
        
    return np.array(pmass)/M_sun, np.array(ppressure), np.array(pradius)/1000,np.array(pdensity)/rho_c    

pMass, pPressure, pRadius, pDensity = politropa(dr)

plt.plot(pRadius, pMass ,linestyle='-', color= 'blue', linewidth=1.2 )   

plt.plot(Radius/1000,Mass/M_sun,linestyle='-',color= 'red',linewidth = 1.2)    
plt.xlabel( r"$r (km)$",fontsize=14)
plt.ylabel(r"$\frac{M}{M_{\odot}}$",fontsize=14)
plt.show()





    
    
    
    
    
    
    
    
    
    





