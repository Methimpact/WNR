
# -*- coding: utf-8 -*-

##########################################################################
## This MODULE generates the initial configuration of SU(2) matrix for   #
##  L by L lattice in two dimension.                                     #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: Mon Sep 29 15:07:50 2014        Sub: Lattice gauge theory       #
## By Dibakar Sigdel        Collection: Python-014-09-12                 #
##########################################################################





import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math
import cmath as cmath
import WNR






def fun(s):
              if s ==0:
                  fn = 1
              else:
                  fn = 0
              return fn 
      

def  CI(NB):
            IN = np.matrix([[complex(0,0) for k in range(NB)]for l in range(NB)])
            for k in range(NB):
                for l in range(NB):
                    if k == l:
                        IN[k,l] = complex(1,0)
            return IN  




class Start(object):
            def __init__(self,L,N):
                    self.L = L
                    self.N = N

            def cold_start(self): 
                I = np.matrix(np.identity(self.N))
                UU = [[[I for x in range(self.L)]for y in range(self.L)]for r in range(2)]
                return UU 



            def SU2(self):
                   ai = complex(0,1)
                   r = [random.random(),random.random(),\
                        random.random()] 
                   xi = (math.pi*(2*r[0]-1))
                   theta =0.5*(math.acos(2*r[1]-1))
                   phi = (math.pi*(2*r[2]-1))
                   a = [0.0 for l in range(2)]
                   a = [math.cos(theta)*(cmath.exp(ai*phi)),\
                   math.sin(theta)*(cmath.exp(ai*xi))]
                   su2 = []
                   su2 = np.matrix([[a[0],a[1]],\
                       [-a[1].conjugate(),a[0].conjugate()]])
                   return su2

       
            def su2tosun(self,s,t):
        
                    SUN = CI(self.N)
                    SU2 = self.SU2()
                    SUN[s,s] = SU2[0,0]
                    SUN[s,t] = SU2[0,1]
                    SUN[t,s] = SU2[1,0]
                    SUN[t,t] = SU2[1,1]
                    return SUN

                            
            def sun_gnr(self):
                    SUNM = CI(self.N)
                    s = 1
                    while s < self.N:
                        t = s+1
                        while t < self.N+1:
                            SUN = self.su2tosun(s-1,t-1)
                            SUNM = np.dot(SUNM,SUN)
                            t = t+1
                        s = s+1
                        ZSUN = SUNM
                        return ZSUN
                      
                    
            def hot_start(self):
                    I = np.matrix(np.identity(self.N))
                    UU = [[[I for x in range(self.L)]for y in range(self.L)]for z in range(2)]
         
                    for i in range (2):
                         for j in range(self.L):
                             for k in range(self.L):
                                 SUN = self.sun_gnr()     
                                 UU[i][j][k] = SUN
                    return UU
         



class Selector(object):
    
            def __init__(self,N):
                    self.N = N
        
            def count(self):
                    ct =  0
                    for k in range(self.N):
                        ct = ct + k
                    return ct


            def  select(self):
                    ct = 0
                    s = 0
                    P = [0 for k in range(self.count())]
                    Q = [0 for k in range(self.count())]
                    while s < self.N-1:
                        t = s+1
                        while t < self.N:
                            P[ct] = s
                            Q[ct] = t
                            ct = ct+1
                            t = t+1
                        s = s+1
                    return P,Q

  
            def extractw(self,W,ct):
                    WD =  CI(2)
                    P,Q = self.select()
                    s = P[ct]
                    t = Q[ct]
                    WD[0,0] = W[s,s]
                    WD[0,1] = W[s,t]
                    WD[1,0] = W[t,s]
                    WD[1,1] = W[t,t]
                    return WD
       

      
            def expandv(self,VD,ct):
                    V = CI(self.N)
                    P,Q = self.select()
                    s = P[ct]
                    t = Q[ct]
                    V[s,s] = VD[0,0]
                    V[s,t] = VD[0,1]
                    V[t,s] = VD[1,0]
                    V[t,t] = VD[1,1]
                    return V       


class Pendelton(object):
    
            def __init__(self,alpha,k):
                    self.alpha = alpha
                    self.k = k
        

            def su2(self,a0): 
                    aa = math.sqrt(1 - (a0**2))
    
                    t = [random.uniform(0.0,1.0),\
                        random.uniform(0.0,1.0)]
           
                    theta = math.acos(2.0*t[0]-1)
                    xi = 2.0*math.pi*(t[1])
    
                    a =  [ a0,\
                    aa*math.sin(theta)*math.cos(xi),\
                    aa*math.sin(theta)*math.sin(xi),\
                    aa*math.cos(theta)]

                    XX  = np.matrix([[complex(a[0],a[3]),complex(a[2],a[1])] \
                            ,[complex(-a[2],a[1]),complex(a[0],-a[3])]])
           
                    return XX 
         

               
            def pendlgnr(self):
                     count = 0
                     while (count < 1):
                         r = [ random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0)]
                 
                         x  = [-(math.log(r[0])/(self.k*self.alpha)),\
                               -(math.log(r[1])/(self.k*self.alpha))]
                 
                         C = (math.cos(2.0*math.pi*r[2]))**2
                         A = x[0]*C
                         delta = x[1]+A
           
                         if (r[3]**2) < (1-(0.5*delta)):
                             a0 = (1- delta)
                             count = count+1
                             XX = self.su2(a0)
                             return XX
                         else: 
                            count = 0
        
                  
         
                                   
class Update(object):
    
            def __init__(self,UU,L,N):
                     self.UU = UU
                     self.L = L
                     self.N = N
  
                                                  
            def staple(self,r,t,s):
                    if r ==0:
                        Q = [1,0,1,1,0,1]
                    elif r==1:
                        Q = [0,1,0,0,1,0]
                    LK = np.matrix(self.UU[r][t][s])
                    D = [ np.matrix(self.UU[Q[0]][(s+1)%self.L][(t-1) + (self.L*fun(t))]).getH(),\
                          np.matrix(self.UU[Q[1]][(t-1) + (self.L*fun(t))][s]).getH(),\
                          np.matrix(self.UU[Q[2]][s][(t-1) + (self.L*fun(t))]),\
                          np.matrix(self.UU[Q[3]][(s+1)%self.L][t]),\
                          np.matrix(self.UU[Q[4]][(t+1)%self.L][s]).getH(),\
                          np.matrix(self.UU[Q[5]][s][t]).getH()]
           
       
                    W = np.dot(D[0],np.dot(D[1],D[2])) \
                     + np.dot(D[3],np.dot(D[4],D[5]))
       
                    WW = np.dot(LK,W)
                    return WW
       

            def  findZk(self,W,ct):
                    Nn = Selector(self.N)
                    WD = Nn.extractw(W,ct)
                    X = WD[0,0] + (WD[1,1]).conjugate()
                    Y = (WD[0,1]).conjugate() - WD[1,0]
                    k = cmath.sqrt(abs(X)**2 + abs(Y)**2).real
                    x = X/k
                    y = Y/k
                    Z = np.matrix([[(x).conjugate(), - (y).conjugate()] ,[y,x]])
                    return k,Z
      
                                                       
            def link(self,r,t,s,alpha):
                    LK =  np.matrix(self.UU[r][t][s])
                    W =  self.staple(r,t,s)
                    Nn = Selector(self.N)
                   
                    V = [CI(self.N) for lt in range(Nn.count())]
                    ct = 0
                    while ct < Nn.count():
                        k,Z = self.findZk(W,ct)
                        XX = Pendelton(alpha,k).pendlgnr()
                        VD = np.dot(XX,Z)
                        V[ct] = Nn.expandv(VD,ct)
                        ct = ct+1
          
                    NU = CI(self.N)
                    for q in range(Nn.count()):   
                       NU = np.dot(NU,V[q])
                    NNU = np.dot(NU,LK)
       
                    self.UU[r][t][s] = NNU
       
                    return self.UU 
                    
                    
                    
                    
 
class Calculate(object):
            def __init__(self,UU,L,N):
                    self.UU = UU
                    self.L = L
                    self.N = N

            def plqt(self,s,t):
                    D  = [ np.matrix(self.UU[0][s][t]),\
                           np.matrix(self.UU[1][(t+1)%self.L][s]),\
                           np.matrix(self.UU[0][(s+1)%self.L][t]).getH(),\
                           np.matrix(self.UU[1][t][s]).getH()]
                    return D         


            def  avplqt(self):
                    sum_trp = 0.0  
                    for s  in range(self.L):
                        for t in range(self.L):
                            D = self.plqt(s,t)
                            UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                            trup = (1.0 - ((1.0/float(self.N))*np.trace(UP).real))
                            sum_trp = sum_trp + (trup/float(self.L*self.L))
                    return sum_trp  
 

            def wloop11(self,s,t):
                    D = self.plqt(s,t)       
                    UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                    wtr =  (1.0/float(self.N))*np.trace(UP).real
                    return wtr 
 

            def wilsonlp(self,K):
                    I = np.matrix(np.identity(self.N))
                    WKK = np.matrix(np.identity( self.N))
                    PD = [[I for k in range(K)]for p in range(4)]
                    DD = [I for k in range(4)]
                    for s  in range(K):
                        PD[0][s] = np.matrix(self.UU[0][0][s])
                    for s in range(K):
                        PD[1][s] = np.matrix(self.UU[1][K][s])
                    for s in range(K):
                        t = K-s-1
                        PD[2][s] = np.matrix(self.UU[0][K][t]).getH()
                    for s in range(K):
                        x = K-s-1
                        PD[3][s] = np.matrix(self.UU[1][0][x]).getH()
                    for r in range(4):
                        for k in range(K):
                            DD[r] = np.dot(DD[r],PD[r][k])
                    WKK = np.dot(DD[0],np.dot(DD[1],np.dot(DD[2],DD[3])))
                    wilp =  (1.0/float(self.N))*np.trace(WKK).real    
                    return wilp    
         
 
            def polyacove(self):
                    Tx = [np.matrix(np.identity(self.N))  for i in range(self.L)]
                    Ty = [np.matrix(np.identity(self.N)) for i in range(self.L)]
       
                    T1 = np.matrix(np.identity(self.N))
                    T2 = np.matrix(np.identity(self.N))
    
                    for i in range(self.L):
                        Tx[i] = np.matrix(self.UU[0][0][i])
                        Ty[i] = np.matrix(self.UU[1][0][i])
           
                    for k in range(self.L):
                        T1 = np.dot(T1,Tx[k])
                        T2 = np.dot(T2,Ty[k])
       
                    chit1 =  (1.0/float(self.N))*np.trace(T1).real
                    chit2 =  (1.0/float(self.N))*np.trace(T2).real
       
                    return chit1,chit2   
       
def Mean_Error(stor_w11):
        nt = len(stor_w11)
        ver = [0.0 for k in range(nt)] 
        mw11 = 0.0
         
        for k in range (nt):
            mw11 = mw11+stor_w11[k]/float(nt)
        for l in range (nt):
            ver[l] = (stor_w11[l]-mw11)**2
        s_error = math.sqrt(sum(ver)/nt**2)
        return  mw11, s_error
            
                    
                            
                                    
                                                    

###################################################################
#Declerations------------------------
# Input = [SUN: 10,Nmax: 16,l: 30, itr:100,sct:20,dlamda:0.25]
SUN = 10
Nmax = 16
l = 30
itr = 100
sct = 20
tct = itr-sct
#------------------------------------------
dlamda = 0.25

plot_dir = "/Users/dibakarsigdel/Dropbox/Plots/"  
data_dir = "/Users/dibakarsigdel/Dropbox/Data/"  
sun = 2
while sun < SUN+1:
    N = sun
    Nvalue = 1
    x = [0.0 for k in range(Nmax)]
    y = [0.0 for k in range(Nmax)]
    y_error = [0.0 for k in range(Nmax)]
    while Nvalue < Nmax+1:
            lamda = dlamda*Nvalue
            alpha = (N**2)/lamda
            xx =  (N*N)/alpha 
            stor_w11 = [0.0 for k in range(itr-sct)]
            
            ll = 1
            U = Start(l,N).cold_start()
            while (ll < itr+1): 
                    print sun,Nvalue,ll
                    for s in range(l):
                        for t in range(l):
                            for r in range(2):
                                U = Update(U,l,N).link(r,s,t,alpha)
                    w11 = Calculate(U,l,N).wloop11(15,15)
                    if ll > sct:
                        stor_w11[ll-sct-1] = w11
                    ll = ll+1
                    
            yy,yser = Mean_Error(stor_w11)                   
            
                   
            x[Nvalue-1] = xx
            y[Nvalue-1] = yy
            y_error[Nvalue-1] = yser
            
            plt.figure(1)
            plt.scatter(xx,yy)
            plt.xlabel('lambda')
            plt.ylabel('W11')
            plt.grid()
            plt.errorbar(xx,yy, yerr = yser,fmt='8')
            plt.savefig(plot_dir + 'plot01.png')
            plt.show()
            
            Nvalue = Nvalue+1
            
            
    st = str(N)
    WNR.wnr('su'+st+'.dat',[x,y,y_error]).writer()
    WNR.wnr(data_dir +'su'+st+'.dat',[x,y,y_error]).writer()
    
    sun = sun+1





















         
