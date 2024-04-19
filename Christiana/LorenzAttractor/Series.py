# -*- coding: utf-8 -*-
"""
For visulising the derived series from Lincoln Nebraska etc
"""
import matplotlib.pyplot as plt
import numpy as np

def Constseries(n,dt):
    '''
    Series for constant dt case, n is the upper limit ie n=nt
    '''
    
    summ =[ (1+dt)**i for i in range(n)]
    summ = np.array(summ)*(dt**2)/2.
    series = np.cumsum(summ)
    return series
    
def case2Series(N,dt,phi):
    '''
    attempt at the case2 dt +gamma*n, not in use
    '''
       
    e_n=[]
    e_n2 =0
    
    for n in range(1,N+1):
        
        prod = 1
        e_n1 = ((dt+(n-1)*phi)**2)/2.
        for m in range(1,n):
            for i in range(1,m+1):
                prod*= 1+dt+(n-i)*phi
            e_n2 = prod*((dt+(n-m-1)*phi)**2)/2.
        if e_n1+e_n2 > 100:
            print n
            e_n.append(0)
        else:
            e_n.append(e_n1+e_n2)
        
    return e_n
            
    
    
    
    

    
    
def case1series(N,dt,a):
    '''
    attempt at the case 1
    '''
    
    e_n = []
    e_n2 =0
    
    for n in range(1,N+1):
        e_n1 = (a**(2*(n-1))*dt**2)/2.
        prod =1
        for m in range(1,n):
            for i in range(1,m+1):
                prod *=(1+a**(n-i)*dt)
            e_n2 = prod*(dt**2 * a**(n-m-1))/2.
        if e_n1+e_n2 > 100:
            print n
            e_n.append(0)
        else:
            e_n.append(e_n1+e_n2)
    return e_n
    
    
def variable(variable_dt):
    '''
    series from nebraska lincoln, takes a list of dts, nt long does the cumulative sum
    according to the series
    '''
    es =[]
    e_n0=0
    for dt_n in variable_dt:
        e_n1 = (1+dt_n)*e_n0 +dt_n**2/2.
        es.append(e_n1)
        e_n0 = e_n1
    return es
        
    
#check values match calculation  
print case2Series(3,1,2)   
    

    
    
#case 1
variable_dtEu_alpha = [(0.0008)*1.000647930908985**i for i in range(2500)]
variable_timeEu_alpha = np.cumsum(variable_dtEu_alpha)


const_timeEu = [0.002]*2500
#case2
variable_dtEu_phi = [(0.0010029)+i*(7.9799919968e-7) for i in range(2500)]
variable_timeEu_phi = np.cumsum(variable_dtEu_phi)

#y_const = Constseries(2500, 0.002)

#y_alpha = case1series(150,0.0008,1.000647930908985)
#y_phi = case2Series(150,0.0010029,7.9799919968e-7)

#get values
y_alpha = variable(variable_dtEu_alpha)
y_phi = variable(variable_dtEu_phi)
y_const = variable(const_timeEu)
n = range(2500)

plt.plot(n,y_const,'k')
plt.plot(n,y_alpha,'steelblue')
plt.plot(n,y_phi,'firebrick')
plt.xlabel('n')
plt.ylabel('Max magnitude of error at n, $|e^{n}|$')
plt.show()
''' 
yCase2 = case2Series(2500,0.0010029,7.9799919968e-7)
xCase2 = range(len(yCase2))

plt.plot(xCase2,yCase2)
'''







