# -*- coding: utf-8 -*-
"""
Order of accuracy
"""

import matplotlib.pyplot as plt
import numpy as np
import lorenzParams
import LorenzAttractor


def orderOfAccuracy():
    
   #Set up dictionary of input parameters from lorenzParams file
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 100.
   

   #use non-cumulative error measure

   dt_range = [0.001,0.002,0.004]
   #dt_range = [0.00025,0.0005,0.001]
   d_err_dts = []
   str_coords = [line.rstrip('\n') for line in open('IC_selected.txt')]
   coords = []

   for str_coord in str_coords:
       float_coords = [ float(x) for x in str_coord.split(',')]
       coords.append(float_coords)
   
   #loop over the number of initial conditions to use

       
   for i in range(0,1):
   
       #create random intial conditions between -15 and 15
       x_o = coords[i][0]
       y_o = coords[i][1]
       z_o = coords[i][2]
       
       for index,dt in enumerate(dt_range):
           #make sure they run for the same duration
           lorenzParam['nt_RK'] = int(round(0.12/dt))
           lorenzParam['nt_FB'] = int(round(0.12/dt))
           lorenzParam['nt_Eu'] = int(round(0.12/dt))
           
           print dt
           #Get values
           x_Eu,y_Eu,z_Eu,d_Eu= LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)
           x_r,y_r,z_r,d_r = LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
           x_FB,y_FB,z_FB,d_FB = LorenzAttractor.ForwardBack(x_o,y_o,z_o, lorenzParam,1.,dt)
           x_hr,y_hr,z_hr,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
           
           #Cut the higher res run so it is comparable to the lower
           RK_mark = int(RK_factor)
           d_RKA =  d_RKh[1::RK_mark]
           d_RKA = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKA))
      
           
           # calculate the error a
           d_diff_RK = abs(d_RKA - d_r)
           #d_diff_RK = np.cumsum(np.array(d_diff_RK)**2)
    
           
           #d_diff_RK = [np.sqrt(r/ind) for ind,r in enumerate(d_diff_RK)]
           #d_diff_RK[0]=0 #Ics are the same
           dError = d_diff_RK[-1]
           print dError
           d_RKA = np.array(d_RKA)

           d_err_dts.append(dError)
           '''
           #try different error measure
           error = abs(d_r[-1] - d_RKA[-1])
           d_err_dts[index].append(error)
           '''
   
   print d_err_dts
   #d_err_dts  = np.mean(np.array(d_err_dts),axis = 1)
   print len(d_err_dts)
   print d_err_dts
   num = d_err_dts[2]-d_err_dts[1]
   denom = d_err_dts[1] - d_err_dts[0]
   m = (1./np.log(2.0))*np.log(num/denom)
   print m, '******'
   #plt.loglog(dt_range,d_err_dts, linestyle = '--', marker ='x', color = 'firebrick')

   xlog = np.log(dt_range)
   y2log = np.log(d_err_dts)
   n2, intercept = np.polyfit(xlog,y2log,1)
   print n2
   
   
#orderOfAccuracy() 

######################For finding gamma and alpha##################
def findAlphaDt(nt, T):
    a = 1.0
    x=0.0
    dt = 0.00002
    while x<T:

        if x< (T -0.01):
            a+=0.00001
            
        elif x< (T - 0.0001):
            a+=0.000001
            
        elif x<(T - 0.0000001):
            a+=0.000000001
            
            
        elif x<(T - 0.000000001):
            a+=0.00000000001
            
        elif x< (T -0.00000000001):
            a+=0.0000000000001
            
        else:
            a+=0.000000000000001

                
            
            
        x = dt*((1-a**nt)/(1-a))
 
        #print x
        
        
        

    return dt, a,x

dtRK,aRK,x1 = findAlphaDt(4250.0,8.5)
dtEu,aEu,x2 = findAlphaDt(300.0,0.6)
dtFB,aFB,x3 = findAlphaDt(300.0,0.6)
'''

dtRK,aRK,x1 = findAlphaDt(5500.0,11.0)
dtEu,aEu,x2 = findAlphaDt(400.0,0.8)
dtFB,aFB,x3 = findAlphaDt(850.0,1.7)
'''
'''
dtRK,aRK,x1 = findAlphaDt(6500.0,13.0)
dtEu,aEu,x2 = findAlphaDt(500.0,1.0)
dtFB,aFB,x3 = findAlphaDt(1400.0,2.8)
'''

print dtRK, aRK

print dtEu , aEu

print dtFB, aFB

'''
def findGammaDt(nt,T):
    dt = 0.005
    gamma = 2*(T-dt)/(nt*(nt-1))
    
    return dt, gamma


'''
def findGammaDt(nt,T,dt):


    
    gamma = 2*(T-dt*nt)/(nt*(nt-1))
    dtMax = dt +(nt-1)*gamma

    return dt,gamma,dtMax

dt2RK,yRK,dttR = findGammaDt(4250.0,8.5,0.00002)
dt2Eu,yEu ,dttE= findGammaDt(300.0,0.6,0.00002)
dt2FB, yFB, dttF= findGammaDt(300.0,0.6,0.00002)

'''  

dt2RK,yRK,dttR = findGammaDt(5500.0,11.0,0.00002)
dt2Eu,yEu ,dttE= findGammaDt(400.0,0.8,0.00002)
dt2FB, yFB, dttF= findGammaDt(850.0,1.7,0.00002)
'''
'''
dt2RK,yRK,dttR = findGammaDt(6500.0,13.0,0.00002)
dt2Eu,yEu ,dttE= findGammaDt(500.0,1.0,0.00002)
dt2FB, yFB, dttF= findGammaDt(1400.0,2.8,0.00002)
'''

print dt2RK,yRK,dttR
print dt2Eu,yEu,dttE
print dt2FB, yFB,dttF













  