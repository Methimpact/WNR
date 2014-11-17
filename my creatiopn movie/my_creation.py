





##########################################################################
## This program simulates the local action density in case of  SU(2)     #
##  in 2D with periodic boundary condition using "Heat bath algorithm".  #   
##                                                                       #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: 2014-Jun-25        Sub: Lattice gauge theory                    #
## By Dibakar Sigdel        Collection: Python-25-06 -03                 #
##########################################################################



import os
import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math
import pylab as pylab


# Delerations --------------------------------
alpha = 16.0
total_itr = 100
l1 = 50
l2 = 50
l1l2 = l1*l2
Dx = float(l1) + 0.5
Dy = float(l2) + 0.5
dt = 1.0
#-------------------------------------------

#########################
#Functions
#########################

def cold_start(r,L1,L2,N): 
   I = np.identity(N)
   CD = [[[I for x in range(L1)]for y in range(L2)]for z in range(r)]
   return CD
   


def periodic_boundary(r,U,L1,L2):
  for r in range(2):
    for s in range(L1):
         U[r][s].append(U[r][s][0])
  for r in range(2):
    U[r].append(U[r][0])
  for r in range(2):
    U[r].append(U[r][1])
  flat_U = U
  return flat_U
  
 
def horizontal_links(UU,L1,L2):
 UV = cold_start(2,L1+2,L2+2,2)
 for s  in range(1,L1+1):
   for t in range(0,L2):
       #LK = np.matrix(UU[0][s][t])
       LK1 = np.matrix(UU[1][s+1][t])
       LK2 = np.matrix(UU[0][s][t+1]).getH()
       LK3 = np.matrix(UU[1][s][t]).getH()
       LK4 = np.matrix(UU[1][s][t-1]).getH()
       LK5 = np.matrix(UU[0][s][t-1])
       LK6 = np.matrix(UU[1][s+1][t-1])
       TT1 = np.dot(LK1,np.dot(LK2,LK3))
       TT2 = np.dot(LK4,np.dot(LK5,LK6))
       CI = TT1 + TT2
       det = ((CI[0,0]*CI[1,1]) - (CI[0,1]*CI[1,0]))
       ssm = det.real
       Am = CI/(math.sqrt(ssm))
       AD = np.matrix(Am).getH()  
       count = 0
       while (count < 1):
           r1 = random.random()
           r2 = random.random()
           r3 = random.random()
           x1 = -(math.log(r1)/((ssm)*alpha))
           x2 = -(math.log(r2)/((ssm)*alpha))
           C = (math.cos(2*math.pi*r3))**2
           A = x1*C
           delta = x2+A
           r4 = random.random()
           if (r4**2) < (1-(0.5*delta)):
               a0 = (1- delta)
               count = count+1
           else: 
               count = 0
       aa = math.sqrt(1 - (a0**2))
       t1 = random.random()
       theta = math.acos(2*t1-1)
       t2 = random.random()
       xi = math.pi*(2*t2 -1)
       a1 = aa*math.sin(theta)*math.cos(xi)
       a2 = aa*math.sin(theta)*math.sin(xi)
       a3 = aa*math.cos(theta)
       XX  = np.matrix([[complex(a0,a3),complex(a2,a1)],[complex(-a2,a1),complex(a0,-a3)]])
       
       NU = np.dot(XX,AD) 
       UV[0][s][t] = NU 
 return UV 
  
def verticle_links(UU,L1,L2): 
 VV = cold_start(2,L1+2,L2+2,2)       
 for s  in range(1,L1+1):
   for t in range(0,L2):
       #LK = np.matrix(UU[1][s][t])
       LK1 = np.matrix(UU[0][s][t])
       LK2 = np.matrix(UU[1][s+1][t])
       LK3 = np.matrix(UU[0][s][t+1]).getH()
       LK4 = np.matrix(UU[0][s-1][t+1]).getH()
       LK5 = np.matrix(UU[1][s-1][t]).getH()
       LK6 = np.matrix(UU[0][s-1][t])
       TT1 = np.dot(LK1,np.dot(LK2,LK3))
       TT2 = np.dot(LK4,np.dot(LK5,LK6))
       CI = TT1 + TT2
       det = ((CI[0,0]*CI[1,1]) - (CI[0,1]*CI[1,0]))
       ssm = det.real
       Am = CI/(math.sqrt(ssm))
       AD = np.matrix(Am).getH()  
       count = 0
       while (count < 1):
           r1 = random.random()
           r2 = random.random()
           r3 = random.random()
           x1 = -(math.log(r1)/((ssm)*alpha))
           x2 = -(math.log(r2)/((ssm)*alpha))
           C = (math.cos(2*math.pi*r3))**2
           A = x1*C
           delta = x2+A
           r4 = random.random()
           if (r4**2) < (1-(0.5*delta)):
               a0 = (1- delta)
               count = count+1
           else: 
               count = 0
       aa = math.sqrt(1 - (a0**2))
       t1 = random.random()
       theta = math.acos(2*t1-1)
       t2 = random.random()
       xi = math.pi*(2*t2 -1)
       a1 = aa*math.sin(theta)*math.cos(xi)
       a2 = aa*math.sin(theta)*math.sin(xi)
       a3 = aa*math.cos(theta)
       XX  = np.matrix([[complex(a0,a3),complex(a2,a1)],[complex(-a2,a1),complex(a0,-a3)]])
       
       NU = np.dot(XX,AD) 
       VV[1][s][t] = NU  
       
 return VV   
 
def revers_boundary(UU,L1,L2):
     UU[0][0] = UU[0][L1]
     UU[1][0] = UU[1][L2]
     FU = UU
     return FU

def  ave_plaquette(UV,L1,L2):
 sum_trp = 0.0  
 for s  in range(L1):
   for t in range(L2):
       UP1 = np.matrix(UV[0][s][t])
       UP2 = np.matrix(UV[1][s+1][t])
       UP3 = np.matrix(UV[0][s][t+1]).getH()
       UP4 = np.matrix(UV[0][s][t]).getH()
       UP = np.dot(UP1,np.dot(UP2,np.dot(UP3,UP4)))
       
       trup = 0.5*np.trace(UP).real
       sum_trp = sum_trp + (trup/float(l1l2))
 return sum_trp     


  
def plaquetteCount(UV,L1,L2):
  plaquette_count = [0 for x in range(l1l2)]
  ct = 0
  for s  in range(L1):
     for t in range(L2):
          UP1 = np.matrix(UV[0][s][t])
          UP2 = np.matrix(UV[1][s+1][t])
          UP3 = np.matrix(UV[0][s][t+1]).getH()
          UP4 = np.matrix(UV[0][s][t]).getH()
          UP = np.dot(UP1,np.dot(UP2,np.dot(UP3,UP4)))
          trup = 0.5*np.trace(UP).real
          #print trup
          plaquette_count[ct] = trup
          ct = ct+1
  
  return  plaquette_count 


def plaquette_stor_room(UU,total_itr):
    ll =0
    PSR  = [[0.0 for p in range(l1l2)]for q in range(total_itr)] 
    while ll < total_itr:
             n_UU = periodic_boundary(2,UU,l1,l2)
             n_VV = horizontal_links(n_UU,l1,l2)
             n_NU = verticle_links(n_VV,l1,l2)
             n_VU = revers_boundary(n_NU,l1,l2)
             PSR[ll] = plaquetteCount(n_VU,l1,l2)
             #aaa = ave_plaquette(n_VU,l1,l2)
             #print aaa
             UU = n_VU[0:2][0:l1][0:l2]
             ll = ll +1
             
    return PSR         
    
def CreateMovie(PSR,plotter, numberOfFrames, fps=10):
	
	for i in range(numberOfFrames):
		plotter(PSR,i)
		fname = '_tmp%05d.png'%i
		pylab.savefig(fname)
		pylab.clf()
 
	os.system("rm movie.mp4")
	os.system("ffmpeg -r "+str(fps)+" -b 1800 -i _tmp%05d.png movie.mp4")
	os.system("rm _tmp*.png")


def plotfunction(PSR,frame):
        dt = 1.0
        values = PSR[frame]
        #print values
        pylab.axis()
        s_map = [(float(p),float(q)) for p in range(1,l1+1) for q in range(1,l2+1)]

        for site in range(l1l2):
             if values[site] < 1.0 :
                if values[site] > 0.95 :
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Indigo')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.95 :
                if values[site] > 0.9 :
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Purple')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.9 :
                if values[site] > 0.85 :
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='DarkViolet')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.85:
                if values > 0.8:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Magenta')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.8:
               if values > 0.75:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='MediumSlateBlue')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.75:
               if values > 0.70:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='RoyalBlue')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.7:
               if values > 0.65:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='DeepSkyBlue')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.65:
               if values > 0.60:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='DarkTurquoise')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.60:
               if values > 0.55:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Cyan')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.55:
               if values > 0.50:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Aquamarine')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.50:
               if values > 0.45:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='SpringGreen')
                  pylab.gca().add_patch(rtgl) 
             if values[site] < 0.45:
               if values > 0.40:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Lime')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.40:
               if values > 0.35:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='GreenYellow')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.35:
               if values > 0.30:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Yellow')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.30:
               if values > 0.25:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Gold')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.25:
               if values > 0.20:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Orange')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.20:
               if values > 0.15:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='OrangeRed')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.15:
               if values > 0.10:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Red')
                  pylab.gca().add_patch(rtgl)
             if values[site] < 0.10:
               if values > 0.05:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='Crimson')
                  pylab.gca().add_patch(rtgl) 
             if values[site] < 0.05:
               if values > 0.00:
                  rtgl = pylab.Rectangle(s_map[site], dt,dt, fc='DarkRed')
                  pylab.gca().add_patch(rtgl)              
                
            #pylab.show()
        #xframe = [0.5,Dx]
        #yframe = [0.5,Dy]
        #dt = [1.0 , 1.0]
        #xdraw = [1.5,1.5]
        #ydraw = [1.5,1.5]
        #ct = 0
        #while ct < l1:
        #    pylab.plot(xframe, xdraw, 'b')
        #    xdraw[0] = xdraw[0] + dt[0]
        #    xdraw[1] = xdraw[1] +dt[1]
        #    ct = ct+1
        #ct = 0
        #while ct < l2:
        #    pylab.plot(ydraw,yframe ,'b')
        #    ydraw[0] = ydraw[0] + dt[0]
        #    ydraw[1] = ydraw[1] +dt[1]
        #    ct = ct+1
        #
        pylab.axis('scaled')
        pylab.axis([0.5, Dx, 0.5, Dy])
        pylab.xticks([])
        pylab.yticks([])
        
        #pylab.close()
        



UU = cold_start(2,l1,l2,2)
PSR = plaquette_stor_room(UU,total_itr)
#print PSR
CreateMovie(PSR,plotfunction,total_itr)
