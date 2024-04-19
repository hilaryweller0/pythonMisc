# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:53:37 2016

@author: christiana
"""
import numpy as np

#                          Time-sorting functions
########################################################################

def find_nearest_index(vals, targets):
    """ Searches through vals for the target values, returning the index
    of the value in vals that is closest numerically to the target value.
    If more than one value in vals is equally close to the target, the index
    of the first value will be returned.
    
    In this context vals with be the high resolution time values and 
    the targets with be the variable reolution time values. The indicies
    returned will provide the incidicies of the d_RKh array that can be 
    comapred to the d array for the variable timestep
    """
    indicies =[]
    start = 0
    #Loop over all target values
    for j in range(len(targets)):
        #re-set indicators
        min_index = -1
        min_diff = None
        # Loop over all the values in the given list from the last index
        #can do this since we know these time arrays will be sequential
        for i in range(start,len(vals)):
            # Find the absolute difference between this value and the target
            diff = abs(vals[i] - targets[j])
            # If this is the first time, or if the difference is smaller than
            # the smallest difference found so far, remember the difference and
            # the index
            if min_diff is None or diff < min_diff:
                min_diff = diff
                min_index = i
        indicies.append(min_index)
        start = min_index
    return indicies 




def time_match(var_t, highres_t):
    var_t_indicies =[]
    highres_t_indicies =[]
    var_t = np.ndarray.tolist(var_t)
    #find t instances that occur in both
    matching_vals= set(var_t).intersection(highres_t)
    for t in matching_vals:
        var_t_indicies.append(var_t.index(t))
        highres_t_indicies.append(highres_t.index(t))
    return var_t_indicies,highres_t_indicies
            

#                      Functions for lorenz attractor time-stepping
########################################################################




def Euler(lorenzParam,factor,x_o,y_o,z_o,dt):
    #Access needed params
    nt = int(lorenzParam['nt_Eu']*factor)
    
    #set up dt list if input is constant
    if type(dt) == float:
        dt_value = dt/factor
        dt = [dt_value]*nt

    
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    #add in Ics
    x[0]=x_o
    y[0]=y_o
    z[0]=z_o
    #position
    r = np.zeros(nt+1)
    r[0] =np.sqrt(x_o**2 +y_o**2+z_o**2)

    for n in range(0,nt):
        dx,dy,dz = lorenz_derivatives(x[n],y[n],z[n],lorenzParam)
        x[n+1] = x[n] + dx*dt[n]
        y[n+1] = y[n] + dy*dt[n]
        z[n+1] = z[n] + dz*dt[n]
            
        #u = dx
        #v = dy
        #w = dz
        
        
        #V = np.sqrt(v**2 +u**2 +w**2)
        #r[n+1] = r[n] + V*dt[n]
        r[n+1] = np.sqrt(x[n+1]**2+y[n+1]**2+z[n+1]**2)
    return x,y,z,r


def lorenz_derivatives(x,y,z,lorenzParam):
    '''
    takes the x y and z values and params and solves each
    of the lorenz derivatives for use in the Runge Kutta update scheme
    '''
    sigma = lorenzParam['sigma']
    rho = lorenzParam['rho']
    beta = lorenzParam['beta']
    
    dx = -sigma*(x-y)
    dy = rho*x -y -x*z
    dz = x*y - beta*z
    
    
    return dx,dy,dz
    
    
    
    
def Euler_step(x,y,z,lorenzParam):
    '''
    not in use
    '''    
    
    dt = lorenzParam['dt']/4.
    dx,dy,dz = lorenz_derivatives(x,y,z,lorenzParam)
    x_new = x + dt*dx
    y_new = y +dt*dy
    z_new = z +dt*dz
    return x_new,y_new,z_new
    
def CN(lorenzParam):
    '''
    not in use
    '''
        
    #Access needed params
    nt = 10*lorenzParam['nt']
    dt = lorenzParam['dt']/10.
    x_o = lorenzParam['x_o']
    y_o = lorenzParam['y_o']
    z_o = lorenzParam['z_o']
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    x[0]=x_o
    y[0]=y_oC
    z[0]=z_o
    
    for n in range(0,nt):
        
        x_eul,y_eul,z_eul= Euler_step(x[n],y[n],z[n],lorenzParam)
        
        dx_n,dy_n,dz_n = lorenz_derivatives(x[n],y[n],z[n],lorenzParam)
        
        dx_n1,dy_n1,dz_n1 = lorenz_derivatives(x_eul,y_eul,z_eul,lorenzParam)
        
        x[n+1] = x[n] +(dx_n +dx_n1)*dt/2.
        y[n+1] = y[n] +(dy_n +dy_n1)*dt/2.
        z[n+1] = z[n] +(dz_n +dz_n1)*dt/2.
    
    return x,y,z
    
    

    
def RungeKutta(lorenzParam,factor,x_o,y_o,z_o,dt):                 
    '''
    Takes the input paramaters and preforms a 4 stage runge kutta method
    for the lorenz eqns returns the values of x,y and z at all timesteps
    including intiial conditions and position
    '''
    
    #Access needed params
    nt = int(lorenzParam['nt_RK']*factor)
    if type(dt) == float:
        dt_value = dt/factor
        dt = [dt_value]*nt

    
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    #Ics
    x[0]=x_o
    y[0]=y_o
    z[0]=z_o
    #position
    r = np.zeros(nt+1)
    #v = np.zeros(nt+1)
    #w = np.zeros(nt+1)
    r[0] =np.sqrt(x_o**2 +y_o**2+z_o**2)
    #v[0] =0

    
    for n in range(0,nt):
  

    
        #At each timestep we must calculate k(1-4) for x(k),y(l) and z(m) 
        k1,l1,m1 = lorenz_derivatives(x[n],y[n],z[n],lorenzParam)
        
        k2,l2,m2 = lorenz_derivatives(x[n]+(dt[n]/2.)*k1,y[n]+(dt[n]/2.)*l1,z[n]+(dt[n]/2.)*m1,\
                  lorenzParam)
                
        k3,l3,m3 = lorenz_derivatives(x[n]+(dt[n]/2.)*k2, y[n]+(dt[n]/2.)*l2,z[n]+(dt[n]/2.)*m2,\
                  lorenzParam)
                
        k4,l4,m4 = lorenz_derivatives(x[n]+dt[n]*k3,y[n]+dt[n]*l3,z[n]+dt[n]*m3,\
                 lorenzParam)
        
        #Update Eqns
        x[n+1] = x[n] + (dt[n]/6.)*(k1+2*k2+2*k3+k4)
        y[n+1] = y[n] + (dt[n]/6.)*(l1+2*l2+2*l3+l4)
        z[n+1] = z[n] + (dt[n]/6.)*(m1+2*m3+2*m3+m4)
        #Calculate the speed to find the distance 
        '''
        u= (k1+2*k2+2*k3+k4)/6.
        v= (l1+2*l2+2*l3+l4)/6.
        w= (m1+2*m3+2*m3+m4)/6.
 
        
        
        V = np.sqrt(v**2 +u**2 +w**2)
        '''
        #r[n+1] = r[n] + V*dt[n]
        r[n+1] = np.sqrt(x[n+1]**2+y[n+1]**2+z[n+1]**2)
        
    return x,y,z,r
    
    
def ForwardBack(x_o,y_o,z_o,lorenzParam,factor,dt):
    '''
    Backwards foreward Matssuno, takes ICs, dt as either a float or list and paramter dictionary
    factor can be used to alter the resolution
    
    '''
    nt = int(lorenzParam['nt_FB']*factor)
    if type(dt) == float:
        dt_value = dt/factor
        dt = [dt_value]*nt

    
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    #Ics
    x[0]=x_o
    y[0]=y_o
    z[0]=z_o
    r = np.zeros(nt+1)
    r[0] =np.sqrt(x_o**2 +y_o**2+z_o**2)
    #paramters
    b =lorenzParam['beta']
    rho = lorenzParam['rho']
    sigma = lorenzParam['sigma']


    for n in range(nt):
        dx=-sigma*(x[n]-y[n])

        x[n+1] = x[n] + dx*dt[n]         
        #u = dx     

        # alternate between calculating y or z first

        if nt%2 == 0:
            
            #Update y then z
            dy = rho*x[n+1] - y[n] - x[n+1]*z[n]
            y[n+1] = y[n] +dy*dt[n]
            #v = dy
            dz = x[n+1]*y[n+1] - b*z[n]
            z[n+1] = z[n] +dz*dt[n]
            

        else:
            #Update z then y
            dz = x[n+1]*y[n] - b*z[n]
            z[n+1] = z[n] +dz*dt[n]
            w = dz
            dy = rho*x[n+1] - y[n] - x[n+1]*z[n+1]
            y[n+1] = y[n] +dy*dt[n]
            v = dy
                
        #V = np.sqrt(v**2 +u**2 +w**2)
        #r[n+1] = r[n] + V*dt[n]
        r[n+1] = np.sqrt(x[n+1]**2+y[n+1]**2+z[n+1]**2)
                    
    
    
    return x,y,z,r
    