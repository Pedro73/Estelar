import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import pylab as pl
import pylab
from pylab import *
from matplotlib import rc






def modelo(dr,e): #la funcion pide dr y el e, que doy abajo al llamar la funcion
    pi=math.pi
    rho_c=5.9*10**17 #densidad central
    m_n= 1.67*10**(-27)  #masa neutron
    h=6.62*10**-34 #constante de planck
    c=3*10**8 #velocidad l/uz
    G=6.67*10**(-11) #constante universal
    M_sun = 1.989*10**(30)
    p_f=((3/(8*pi))*(rho_c/m_n))**(1.0/3) * h #momento de fermi
    x=float(((p_f/1.25)/(m_n*c))) #x

    


    m=0.0 #masa inicial
    dm=0
    dP=0
    r=1 #radio inicial

#listas 

    lista_P=[]
    lista_E=[]
    lista_m=[]
    lista_r=[]    
     

    E=e
    P=((pi*m_n**4*c**5)/(3*h**3))*((((x**2)+1)**0.5*(2*x**3-3*x)+ 3*math.asinh(x))) #p(x)
   # print(E)

    while P>0:
        print(P)

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
    
    return(max(lista_m)/M_sun, max(lista_r))##devuevlo la masa maxima normalizada y radio maximo 
ec=[] 
Masa=[]
Radio=[]
ei=0.1*10**34   
for x in range (0,999): #llenamos la lista ec con valores entre 0.01*10**34 a 1*10**35
    ec.append((x/100.0)*10**34)#o algo por ahi xd
    
energias=[]

for en in ec: #le sumamos las energias de ec a una energia base ei(ei no es el ec que use en la pregunta 2 ni nada) 
    energias.append(ei+en)# es solo para tener un valor inicial o algo asi, no creo que import realmente

for energia in energias: #rellenar las 2 listas masa y radio con valores de la funcion
    Masa.append(modelo(10,energia)[0])
    Radio.append(modelo(10,energia)[1])

print(max(Masa),"masa")
ayy=0
for dato in Masa:
    
    if dato==max(Masa):
        print(Radio[ayy])
    ayy=ayy+1





figure(0) #grafico M(ec)
plt.plot(energias,Masa)
plt.ylabel(r"$\frac{M}{M_\odot}$")
plt.xlabel(r"$E_{c}$")
plt.title("")
plt.legend()
plt.draw()

figure(1) #grafico R(ec)
plt.plot(energias,Radio)
plt.ylabel("Radio[m]")
plt.xlabel(r"$E_{c}$")
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







