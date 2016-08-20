import numpy as np
import matplotlib.pyplot as plt

def modelo(dr,n): 
    
    l_P = [] 
    l_rho = [] 
    l_m = [] 
    l_r = [] 
        
    rho = 1 
    
    nn=float(n)
    P = rho**((nn+1)/nn) 
    
    m=0 
    dm=0 
    dP=0 
    r=0 
        
    while P> 0:
    
        cons=3*rho*(r**2) 
        dm= cons*dr 
        m=m + dm
        l_m.append(m)
                
                
        if r>=0.001: 
            hidro = -(rho*m)/(r**2)        
    
        if r < 0.001:
            hidro = -r
            
            
        dP = dr*hidro
        P = P + dP 
        
        if P>=0:
            
            l_P.append(P)
            rho = (P)**(1/((nn+1)/nn)) 
            l_rho.append(rho) 
            r = r + dr 
            l_r.append(r)        
        
        
        if P < 0:
            
            l_P.append(0) 
            l_r.append(r)
            l_rho.append(0)
            
                        
        for i in range (0,len(l_r)): 
            a=l_r[i]
            
            if a>l_r[i-1]:
                maxr=l_r[i]
                
   
            
    return   np.array(l_r)/maxr , np.array(l_rho) 




r_32,rho_32 = modelo(0.001,1.5) 
#r_43,rho_43 = modelo(0.001,1.33333)
#r_53,rho_53 = modelo(0.001,1.66667)
#r_340,rho_340 = modelo(0.001,3.4)


r_s, rho_s = np.loadtxt("sol.dat",usecols=(0,3),unpack=True)


for i in range (0,len(rho_s)): 
    a=rho_s[i]            
    if a>rho_s[i-1]:
        maxrho=rho_s[i]
                
#plt.plot(r_s,rho_s/maxrho,"y",label="sol")
plt.plot(r_32,rho_32,"r",label="n=3/2")
#plt.plot(r_43,rho_43,"b",label="n=4/3")
#plt.plot(r_53,rho_53,"g",label="n=5/3")
#plt.plot(r_340,rho_340,"m",label="n=3.4")


plt.legend()
plt.ylabel(r"$\rho / \rho_c$")
plt.xlabel(r"$r / r_o$")
plt.title("$\Delta$=0.001")

plt.show()
