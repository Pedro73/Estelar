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
    c=5.9*10**8 #velocidad luz
    G=6.67*10**(-11) #constante universal
    M_sun = 1.989*10**(30)
    p_f=((3/(8*pi))*(rho_c/m_n))**(1.0/3) * h #momento de fermi
    x=float(((p_f))/(m_n*c)) #x
    


    m=0.0 #masa inicial
    dm=0
    dP=0
    r=1 #radio inicial

#listas 

    lista_P=[]
    lista_E=[]
    lista_m=[]
    lista_r=[]    
     

    E=(((pi*m_n**4*c**5)/h**3)*((((x**2)+1)**0.5*(2*x**3+x)-math.asinh(x)))) 
    P=((pi*m_n**4*c**5)/(3*h**3))*((((x**2)+1)**0.5*(2*x**3-3*x)+ 3*math.asinh(x)))


    while P>0:

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
    
    
    return(np.array(lista_P),np.array(lista_m)/M_sun,np.array(lista_r),np.array(lista_E)/max(lista_E))
    
    
#Presion=np.array(modelo(10)[0])
Masa=np.array(modelo(10)[1])
Radio=np.array(modelo(10)[2])
Energia=np.array(modelo(10)[3])
#Masa_O=Masa/(2*10**30)



#plt.legend()
#plt.ylabel(r"$\frac{M}{M_\odot}$")

def politropa(dr):
    
    
    h = 6.62*10**-34  
    m =  1.67*10**(-27)    
    rho_c = 5.9*10**17    
    n = float(rho_c/m)
    p_f = ((3*n/(8*np.pi))**(1.0/3))*h
    pi=math.pi

    M_sun = 1.989*10**(30)

    G = 6.67*10**(-11)
    c = 3*10**8
    x=float(((p_f))/(m*c))
    e = (((pi*m**4*c**5)/h**3)*((((x**2)+1)**0.5*(2*x**3+x)-math.asinh(x))))
    pmass = []
    ppressure = []
    pdensity = []
    pradius = []
    
    K = ((8*np.pi*h**2)/(15*m))*((3./(m*8*np.pi*c**2))**(5./3))#((h**2)/(5*m))*(3.0/(8*np.pi))**(2.0/3)
    P = K*e**(5./3)
    r = 1 # Radio inicial 
    M = 0
    
    while P > 0:
        
        dM = (4*np.pi*(r**2)*e/(c**2))*dr
        M += dM
        pmass.append(M)
        
        dP = -(G*M*e/((r*c)**2))*dr
        P += dP
        
        if P < 0:
           ppressure.append(0)
           pdensity.append(0)
           pradius.append(r)
           
        else:
            
           ppressure.append(P) 
           e = (P/K)**(3./5)#((P/K)*(m**(5.0/3)))**(3.0/5)
           pdensity.append(e)
           r += dr
           pradius.append(r)  
        
    return np.array(pmass)/M_sun, np.array(ppressure), np.array(pradius),np.array(pdensity)/max(pdensity) 
pMass, pPressure, pRadius, pDensity = politropa(10)

#plt.plot(pRadius, pMass ,linestyle='-', color= 'blue', linewidth=1.2 )   
M_sun = 1.989*10**(30)
#plt.plot(Radio,Masa/M_sun,linestyle='-',color= 'red',linewidth = 1.2)    
#plt.xlabel( r"$r (m)$",fontsize=14)
#plt.ylabel(r"$\frac{M}{M_{\odot}}$",fontsize=14)
#plt.show()
figure(0)
plt.plot(pRadius,pDensity,color="red",label="Politropa")
plt.plot(Radio,Energia,label="Estrella de neutrones")
plt.ylabel("E/Ec")
plt.xlabel("Radio[m]")
plt.title("")
plt.legend()
plt.draw()

figure(1)
plt.plot(pMass,pDensity,color="red",label="Politropa")
plt.plot(Masa,Energia,label="Estrella de neutrones")
plt.ylabel("E/Ec")
plt.xlabel(r"$\frac{M}{M_\odot}$")
plt.title("")
plt.legend()
plt.draw()

figure(2)
plt.plot(pRadius,pMass,color="red",label="Politropa")
plt.plot(Radio,Masa,label="Estrella de neutrones")
plt.ylabel(r"$\frac{M}{M_\odot}$")
plt.xlabel("Radio [m]")
plt.title("")
plt.legend()
plt.draw()

















