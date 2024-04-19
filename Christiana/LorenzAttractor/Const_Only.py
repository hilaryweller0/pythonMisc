# -*- coding: utf-8 -*-
from __future__ import division
"""
Created on Tue Jun 28 17:37:27 2016

@author: christiana
"""

# -*- coding: utf-8 -*-
"""
Copy of main for constant time-step only: allows quicker run time for testing small chanegs
before implementation

"""

import matplotlib.pyplot as plt
import LorenzAttractor
import LorenzPlot
import numpy as np
import csv
from mpl_toolkits.mplot3d import Axes3D


def main():
   #Set up dictionary of input parameters from lorenzParams file
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   #access constant time-step size, this is the same accross all schemes
   dt = float(lorenzParam['dt'])

   
   RKtime = lorenzParam['const_timeRK']
   Eutime = lorenzParam['const_timeEu']
   FBtime = lorenzParam['const_timeFB']
   #Number of timesteps required for each scheme
   ntRK = lorenzParam['nt_RK']
   ntEu = lorenzParam['nt_Eu']
   ntFB = lorenzParam['nt_FB']

  

   


   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 100.
   highres_dtRK = [dt/RK_factor]*int(ntRK*RK_factor)
   highres_dtEu = [dt/RK_factor]*int(ntEu*RK_factor)
   highres_dtFB = [dt/RK_factor]*int(ntFB*RK_factor) 
   highres_tRK = np.cumsum(highres_dtRK)
   highres_tEu = np.cumsum(highres_dtEu)
   highres_tFB = np.cumsum(highres_dtFB)

   

   
   #RK_dt= dt/RK_factor
   
   
   #Initialise list to store the total-distance-travelled errors
   #for each set of initial conditions and one for the averages
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
   ICxs=[]
   ICys =[]
   ICzs =[]
   for str_coord in str_coords:
       float_coords = [ float(x) for x in str_coord.split(',')]
       coords.append(float_coords)
   print len(coords)
   

   oldstr_coords = [line.rstrip('\n') for line in open('ICs.txt')]
   old_coords = []
   for oldstr_coord in oldstr_coords:
       oldfloat_coords = [ float(x) for x in oldstr_coord.split(',')]
       old_coords.append(oldfloat_coords)
       

   Xs = [old_coords[h][0] for h in range(len(old_coords))]
   Ys = [old_coords[h][1] for h in range(len(old_coords))]
   Zs = [old_coords[h][2] for h in range(len(old_coords))]

   
   #loop over the number of initial conditions to use
   for i in range(100):
       
       #create random intial conditions from a pre-set solution
       #goes up to 135000 b/c I ignored first 5000 points of high res solution to help
       #ensure it was on the attractor
       '''
       IC_index = np.random.randint(0,135000)
       x_o = coords[IC_index][0]
       y_o = coords[IC_index][1]
       z_o = coords[IC_index][2]
       

      
       
       #ICs = str( x_o) + ',' +str(y_o)+ ',' + str(z_o) +'\n'
       '''

       
       x_o = coords[i][0]
       y_o = coords[i][1]
       z_o = coords[i][2]
       ICxs.append(x_o)
       ICys.append(y_o)
       ICzs.append(z_o)
       '''
       x_o = -15 +30*np.random.random()
       y_o = -15 +30*np.random.random()
       z_o = -15 +30*np.random.random()
       '''
   
       #Get values
       xRK,yRK,zRK,d_RK= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
       print 'fixed run'
       xRKh,yRKh,zRKh,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
       print 'high res run'
   
       print 'variable run'
       x_eu,y_eu,z_eu,d_eu = LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)
   
       xFB,yFB,zFB,d_FB = LorenzAttractor.ForwardBack(x_o,y_o,z_o,lorenzParam,1.,dt)
 


   
       #Cut the higher res run so it is comparable to the lower
       RK_mark = int(RK_factor)
       d_RKE = d_RKh[1:(ntEu)*RK_mark:RK_mark]
       d_RKFB = d_RKh[1:(ntFB)*RK_mark:RK_mark]
       d_RKA =  d_RKh[1::RK_mark]
       d_RKA = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKA))
       d_RKE = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKE))
       d_RKFB = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKFB))
   
       
      
         

   
       #pos_diff_CN = np.sqrt((x_A-x)**2 +(y_A-y)**2 +(z_A-z)**2)
       #pos_diff_RK = np.sqrt((x_RKA-x)**2 +(y_RKA-y)**2 +(z_RKA-z)**2)
       d_diff_RK = abs(d_RKA - np.array(d_RK))
       d_diff_RK = np.cumsum(np.array(d_diff_RK)**2)
       d_diff_RK = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_RK)]
       
       d_diff_eu = abs(d_RKE - np.array(d_eu))
       d_diff_eu = np.cumsum(np.array(d_diff_eu)**2)
       d_diff_eu = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_eu)]
       
       d_diff_BF = abs(d_RKFB - np.array(d_FB))
       d_diff_BF = np.cumsum(np.array(d_diff_BF)**2)
       d_diff_BF = [np.sqrt(r/(index+1)) for index,r in enumerate(d_diff_BF)]
       
       
       
       
       d_errorRK.append(d_diff_RK)
       d_errorEu.append(d_diff_eu)
       d_errorBF.append(d_diff_BF)
       

       print i
   print 'doing calcs' 
   
   
   av_d_errorsRK  = np.mean(np.array(d_errorRK),axis = 0)
   av_d_errorsEu  = np.mean(np.array(d_errorEu),axis = 0)
   av_d_errorsBF  = np.mean(np.array(d_errorBF),axis = 0)


   '''
   LorenzPlot.errorPlot(av_d_errorsRK, lorenzParam,'RK4','const',ntRK)
   LorenzPlot.errorPlot(av_d_errorsEu, lorenzParam,'Eu','const',ntEu)
   LorenzPlot.errorPlot(av_d_errorsBF, lorenzParam,'BF','const',ntFB)
   '''
   
   #LorenzPlot.ICplot(ICxs,ICys,ICzs,lorenzParam,Xs,Ys,Zs)
   
   
   for d in d_errorRK:
       plt.plot(RKtime, d[1:], c='lightsteelblue')
       
   plt.plot(RKtime, av_d_errorsRK[1:],c='k')
   plt.xlabel('duration, t')
   plt.ylabel('RMS error in position: RK4 wrt RK4')
   plt.axvline(x = 8, ls = '--',color = 'b')
   plt.axvline(x = 13, ls = '--', color = 'b')
   plt.savefig('Spaghetti_plot_longRunRK4.pdf')
   plt.show()
   
   
   for d in d_errorEu:
       plt.plot(Eutime, d[1:], c='lightsteelblue')
   plt.plot(Eutime, av_d_errorsEu[1:],c='k')
   plt.xlabel('duration, t')
   plt.ylabel('RMS error in position: Euler wrt RK4')
   plt.axvline(x = 0.5, ls = '--',color = 'b')
   plt.axvline(x = 0.9, ls = '--', color = 'b')
   plt.savefig('Spaghetti_plot_longRunEu.pdf')
   plt.show()
   
   
   for d in d_errorBF:
       plt.plot(FBtime, d[1:],c='lightsteelblue')
       
       
   plt.plot(FBtime, av_d_errorsBF[1:], c='k') 
   plt.xlabel('duration, t')
   plt.ylabel('RMS error in position: FB wrt RK4') 
   plt.axvline(x = 1.2, ls = '--',color = 'b')
   plt.axvline(x = 2.8, ls = '--', color = 'b')
   plt.savefig('Spaghetti_plot_longRunFB.pdf')
   plt.show()

   LorenzPlot.errorPlot(av_d_errorsRK, lorenzParam,'RK','Const',ntRK,8.5,13)
   LorenzPlot.errorPlot(av_d_errorsEu, lorenzParam,'Eu','Const',ntEu,0.6,1)
   LorenzPlot.errorPlot(av_d_errorsBF, lorenzParam,'BF','Const',ntFB,0.6,2.8)
   
 
   LorenzPlot.lorentzPlotting(xRKh[1:],yRKh[1:],zRKh[1:],lorenzParam,ntRK)


   
   
if __name__ == '__main__':
    main()

