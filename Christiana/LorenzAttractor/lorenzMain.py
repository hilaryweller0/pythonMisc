# -*- coding: utf-8 -*-
from __future__ import division

"""
Main to exeute the solving of the lorenz equations 
-- need to implement the constant addition time step growth 

"""

import matplotlib.pyplot as plt
import LorenzAttractor
import LorenzPlot
import numpy as np
import csv


def main():
    
   #Set up dictionary of input parameters from lorenzParams file
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   #access constant time-step size, this is the same accross all schemes
   dt = float(lorenzParam['dt'])
   
   #access list of variable step sizes for each scheme
   dt_variableRK = lorenzParam['variable_dtRK']
   dt_variableEu = lorenzParam['variable_dtEu']
   dt_variableFB = lorenzParam['variable_dtFB']
   
   #time axes for variable timestep
   var_timeRK = lorenzParam['variable_timeRK']
   var_timeEu = lorenzParam['variable_timeEu']
   var_timeFB = lorenzParam['variable_timeFB']

   
   #Number of timesteps required for each scheme
   ntRK = lorenzParam['nt_RK']
   ntEu = lorenzParam['nt_Eu']
   ntFB = lorenzParam['nt_FB']
   
   #keys for plotting
   varTypeRK = lorenzParam['var_typeRK']
   varTypeEu = lorenzParam['var_typeEu']
   varTypeBF = lorenzParam['var_typeFB']
   

   


   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 100.
   
   #create the time-array for the high reslution of the same length as each schemes duration
   highres_dtRK = [dt/RK_factor]*int(ntRK*RK_factor)
   highres_dtEu = [dt/RK_factor]*int(ntEu*RK_factor)
   highres_dtFB = [dt/RK_factor]*int(ntFB*RK_factor) 
   
   highres_tRK = np.cumsum(highres_dtRK)
   highres_tEu = np.cumsum(highres_dtEu)
   highres_tFB = np.cumsum(highres_dtFB)

   

   
   #RK_dt= dt/RK_factor
   
   
   #Initialise list to store the cumulative rms pos errors (d_error)
   #for each set of initial conditions and one for the averages v = variable
   d_errorRK =[]
   d_errorRKv=[]
   
   av_d_errorsRK = [] 
   av_d_errorsRKv =[]
   
   d_errorEu =[]
   d_errorEuv =[]
   
   av_d_errorsEu =[]
   av_d_errorsEuv =[]
   
   d_errorBF=[]
   d_errorBFv =[]
   
   av_d_errorsBF =[]
   av_d_errorsBFv =[]

   #set up list of ICs from file for use
   str_coords = [line.rstrip('\n') for line in open('IC_selected.txt')]

   coords = []

   for str_coord in str_coords:
       float_coords = [ float(x) for x in str_coord.split(',')]
       coords.append(float_coords)
   
   

   
   
   
   #loop over the number of initial conditions to use
   for i in range(100):

       #get ICs       
       x_o = coords[i][0]
       y_o = coords[i][1]
       z_o = coords[i][2]

       '''
       x_o = -15 +30*np.random.random()
       y_o = -15 +30*np.random.random()
       z_o = -15 +30*np.random.random()
       '''
   
       #Get values
   
       #high res
       xRKh,yRKh,zRKh,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
       #RK then Rk variable
       xRK,yRK,zRK,d_RK= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
       xRKv,yRKv,zRKv,d_RKv= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt_variableRK)
       #Euler then euler variable
       x_eu,y_eu,z_eu,d_eu = LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)
       xeu,yeu,zeu,d_euv = LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt_variableEu)
       #FB then FB variable
       xFB,yFB,zFB,d_FB = LorenzAttractor.ForwardBack(x_o,y_o,z_o,lorenzParam,1.,dt)
       xFBv,yFBv,zFBv,d_FBv = LorenzAttractor.ForwardBack(x_o,y_o,z_o,lorenzParam,1.,dt_variableFB)


   
       #Cut the higher res run so it is comparable to the lower
       RK_mark = int(RK_factor)
       d_RKE = d_RKh[1:(ntEu)*RK_mark:RK_mark]
       d_RKFB = d_RKh[1:(ntFB)*RK_mark:RK_mark]
       d_RKA =  d_RKh[1::RK_mark]
       #attach the Ics back
       d_RKA = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKA))
       d_RKE = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKE))
       d_RKFB = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKFB))
   
       
      
       #find the nearest indicies to compare high res and variable time-stepping, only do once 
       if i ==0:
           print 'getting indicies'
           indiciesRK = LorenzAttractor.find_nearest_index(highres_tRK, var_timeRK)
           indiciesEu = LorenzAttractor.find_nearest_index(highres_tEu, var_timeEu)
           indiciesFB = LorenzAttractor.find_nearest_index(highres_tFB, var_timeFB)
           print 'got indicies'
       
       
       
       #Get the difference between the high res and the variable solutions

       d_diff_RKv = [abs(d_RKh[indiciesRK[g]] - d_RKv[g]) for g in range(len(indiciesRK))]
       d_diff_Euv = [abs(d_RKh[indiciesEu[h]] - d_euv[h]) for h in range(len(indiciesEu))]
       d_diff_FBv = [abs(d_RKh[indiciesFB[k]] - d_FBv[k]) for k in range(len(indiciesFB))]
       
       #caluclate the RMS in 2 stages
       d_diff_RKv = np.cumsum(np.array(d_diff_RKv)**2)
       d_diff_Euv = np.cumsum(np.array(d_diff_Euv)**2)
       d_diff_FBv = np.cumsum(np.array(d_diff_FBv)**2)
       
       d_diff_RKv = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_RKv)]
       d_diff_Euv = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_Euv)]
       d_diff_FBv = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_FBv)]
      
   

   
       #pos_diff_CN = np.sqrt((x_A-x)**2 +(y_A-y)**2 +(z_A-z)**2)
       #pos_diff_RK = np.sqrt((x_RKA-x)**2 +(y_RKA-y)**2 +(z_RKA-z)**2)
   
   
       #caluclate the RMS for the constant time-step runs
       d_diff_RK = abs(d_RKA - d_RK)
       
       d_diff_RK = np.cumsum(np.array(d_diff_RK)**2)
       d_diff_RK = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_RK)]
       
       d_diff_eu = abs(d_RKE - d_eu)
       d_diff_eu = np.cumsum(np.array(d_diff_eu)**2)
       d_diff_eu = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_eu)]
       
       d_diff_BF = abs(d_RKFB - d_FB)
       d_diff_BF = np.cumsum(np.array(d_diff_BF)**2)
       d_diff_BF = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_BF)]
       
       
       
       #append all 6 errors to their lists
       d_errorRK.append(d_diff_RK)
       d_errorEu.append(d_diff_eu)
       d_errorBF.append(d_diff_BF)
       d_errorRKv.append(d_diff_RKv)
       d_errorEuv.append(d_diff_Euv)
       d_errorBFv.append(d_diff_FBv)
       print i
   print 'doing calcs' 
   
   #taek the average of each error for corresponding times over all ICs
   av_d_errorsRK  = np.mean(np.array(d_errorRK),axis = 0)
   av_d_errorsEu  = np.mean(np.array(d_errorEu),axis = 0)
   av_d_errorsBF  = np.mean(np.array(d_errorBF),axis = 0)
   av_d_errorsRKv = np.mean(np.array(d_errorRKv),axis =0)
   av_d_errorsEuv = np.mean(np.array(d_errorEuv),axis =0)
   av_d_errorsBFv = np.mean(np.array(d_errorBFv),axis =0)
   
   print 'RK' ,varTypeRK, (av_d_errorsRK[-1] - av_d_errorsRKv[-1])
   print 'Eu' ,varTypeEu, (av_d_errorsEu[-1] - av_d_errorsEuv[-1])
   print 'FB' ,varTypeBF, (av_d_errorsBF[-1] - av_d_errorsBFv[-1])
   

   #error comparison plots
   LorenzPlot.errorCompar(av_d_errorsRK,av_d_errorsRKv,\
                           dt,ntRK,'RK4',dt_variableRK,varTypeRK)
                           
   LorenzPlot.errorCompar(av_d_errorsEu,av_d_errorsEuv,\
                           dt,ntEu,'EulerF',dt_variableEu,varTypeEu)

   LorenzPlot.errorCompar(av_d_errorsBF,av_d_errorsBFv,\
                           dt,ntFB,'Back_F',dt_variableFB,varTypeBF)
   #singular error plots                        
   #LorenzPlot.errorPlot(av_d_errorsRK, lorenzParam,'RK4','const',ntRK)
   #LorenzPlot.errorPlot(av_d_errorsEu, lorenzParam,'Eu','const',ntEu)
   #LorenzPlot.errorPlot(av_d_errorsBF, lorenzParam,'BF','const',ntFB)
   



 
   LorenzPlot.lorentzPlotting(xRKh[1:],yRKh[1:],zRKh[1:],lorenzParam,ntRK*RK_factor)



   
if __name__ == '__main__':
    main()

