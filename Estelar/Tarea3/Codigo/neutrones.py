import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import pylab as pl
import pylab
from pylab import *
from matplotlib import rc






def modelo(dr):
    pi=math.pi
    rho_c=5.9*10**17 #densidad central
    m_n= 1.68*10**(-27)  #masa neutron
    h=6.62*10**-34 #constante de planck
    c=3*10**8 #velocidad luz
    G=6.67*10**(-11) #constante universal
    p_f=((3/(8*pi))*(rho_c/m_n))**(1.0/3) * h #momento de fermi
    x=p_f/(m_n*c) #x
    M_sun = 1.989*10**(30)

    m=0.0 #masa inicial
    dm=0
    dP=0
    r=1 #radio inicial

#listas 

    lista_P=[]
    lista_E=[]
    lista_m=[]
    lista_r=[]    
     

    E=1.44*10**34 
    P=((pi*m_n**4*c**5)/(3*h**3))*((((x**2)+1)**0.5*(2*x**3-3*x)+ 3*math.asinh(x)))


    while P>0:
        
        print(dm)
        
        dP_dr=-((G*(E/c**2.0 + P/c**2.0))*(m + 4.0*pi*r**3*P/c**2.0))/(r*(r-2*G*m/c**2))

                
        dP_dx=(((pi*m_n**4*c**5)/(3*h**3))*((8*x**4)/(x**2+1)**0.5))
        
        dx=((dP_dr)/(dP_dx))*dr        
        
        dP=(((pi*m_n**4*c**5)/(3*h**3))*((8*x**4)/(x**2+1)**0.5))*dx
        dE_dx=(((pi*m_n**4*c**5)/h**3)*((8*(x**4+x**2))/(x**2+1)**0.5)) 
        
        
        
        dE=(((pi*m_n**4*c**5)/h**3)*((8*(x**4+x**2))/(x**2+1)**0.5))*dx
        
                
        dm=((4*pi*E*r**2)/(c**2))*dr
                
          
        x=x+dx
        m=m+dm
        P=P+dP
        r=r+dr
        E=E+dE
        lista_P.append(P)
        lista_m.append(m)
        lista_r.append(r)
        lista_E.append(E)
    
    
    
    return(np.array(lista_P),np.array(lista_m)/M_sun,np.array(lista_r),np.array(lista_E))
    
    
Presion=np.array(modelo(10)[0])
Masa=np.array(modelo(10)[1])
Radio=np.array(modelo(10)[2])
Energia=np.array(modelo(10)[3])
Masa_O=Masa/(2*10**30)


figure(0) #grafico M(ec)
plt.plot(Radio,Energia)
plt.ylabel(r"$E_{c}$")
plt.xlabel("R[m]")
plt.title("")
plt.legend()
plt.draw()

figure(1) #grafico R(ec)
plt.plot(Radio,Presion)
plt.ylabel("P[pascal]")
plt.xlabel("Radio[m]")
plt.title("")
plt.legend()
plt.draw()

figure(2) #grafico M(R)
plt.plot(Radio,Masa)
plt.ylabel(r"$\frac{M}{M_\odot}$")
plt.xlabel("Radio[m]")
plt.title("")
plt.legend()
plt.draw()


































