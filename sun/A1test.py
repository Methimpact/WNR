


import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math
import cmath as cmath
#import grosswitt
#from scipy import linalg 



def fun(s):
              if s ==0:
                  fn = 1
              else:
                  fn = 0
              return fn 
      

def  CI(N):
            IN = np.matrix([[complex(0,0) for k in range(N)]for l in range(N)])
            for k in range(N):
                for l in range(N):
                    if k == l:
                        IN[k,l] = complex(1,0)
            return IN  






def cold_start(L,N): 
                I = np.matrix(np.identity(N))
                UU = [[[I for x in range(L)]for y in range(L)]for r in range(2)]
                return UU 



          



def su2(a0): 
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
         

               
def pendlgnr(alpha,det):
                     count = 0
                     failed = 0  
                     while (count < 1):
                         r = [ random.uniform(0.0,1.0),\
                               random.uniform(0.0,1.0),\
                               random.uniform(0.0,1.0),\
                               random.uniform(0.0,1.0)]
                 
                         x  = [-(math.log(r[0])/(det*alpha)),\
                               -(math.log(r[1])/(det*alpha))]
                 
                         C = (math.cos(2.0*math.pi*r[2]))**2
                         A = x[0]*C
                         delta = x[1]+A
           
                         if (r[3]**2) < (1-(0.5*delta)):
                             a0 = (1- delta)
                             count = count+1
                             XX = su2(a0)
                             return failed , XX
                         else: 
                            count = 0
                            failed = failed+1
        


    
       


def count(N):
                    ct =  0
                    for k in range(N):
                        ct = ct + k
                    return ct


def  select(N):
                    ct = 0
                    s = 0
                    P = [0 for k in range(count(N))]
                    Q = [0 for k in range(count(N))]
                    while s < N-1:
                        t = s+1
                        while t < N:
                            P[ct] = s
                            Q[ct] = t
                            ct = ct+1
                            t = t+1
                        s = s+1
                    return P,Q

  
def extractw(N,W,ct):
                    WD =  CI(2)
                    P,Q = select(N)
                    s = P[ct]
                    t = Q[ct]
                    WD[0,0] = W[s,s]
                    WD[0,1] = W[s,t]
                    WD[1,0] = W[t,s]
                    WD[1,1] = W[t,t]
                    return WD
       

      
def expandv(N,VD,ct):
                    V = CI(N)
                    P,Q = select(N)
                    s = P[ct]
                    t = Q[ct]
                    V[s,s] = VD[0,0]
                    V[s,t] = VD[0,1]
                    V[t,s] = VD[1,0]
                    V[t,t] = VD[1,1]
                    return V       


def Qq(r):
                    if r ==0:
                        QQ = [1,0,1,1,0,1]
                    elif r==1:
                        QQ = [0,1,0,0,1,0]
                    return QQ
         
                                   

                                     
def staple(U,L,r,t,s):
                    Q = Qq(r)
                  
                    #LK = np.matrix(self.UU[r][t][s])
                    D = [(U[Q[0]][(s+1)% L][(t-1) + (L*fun(t))]).getH(),\
                         (U[Q[1]][(t-1) + (L*fun(t))][s]).getH(),\
                         (U[Q[2]][s][(t-1) + (L*fun(t))]),\
                         (U[Q[3]][(s+1)%L][t]),\
                         (U[Q[4]][(t+1)%L][s]).getH(),\
                         (U[Q[5]][s][t]).getH()]
           
                    W = np.dot(D[0],np.dot(D[1],D[2])) \
                        + np.dot(D[3],np.dot(D[4],D[5]))
                        
                    LK = U[r][t][s]
                    WW = np.dot(LK,W)
                    return WW
       
def  findZk(WD,ct):
                   
                    X = (WD[0,0] + (WD[1,1]).conjugate())
                    Y = ((WD[0,1]).conjugate() - WD[1,0])
                    k = cmath.sqrt((X*X.conjugate()) + (Y*Y.conjugate())).real
                    x = X/(k)
                    y = Y/(k)
                    Z = np.matrix([[(x).conjugate(), - (y).conjugate()] ,[y,x]])
                    return k,Z
                    
                     
                                  
def link(U,L,N,r,t,s,alpha):
                    LK =  U[r][t][s]
                    W  =  staple(U,L,r,t,s)
                    
                   
                    V = [CI(N) for lt in range(count(N))]
                    ct = 0
                    while ct < count(N):
                        WD =  extractw(N,W,ct)
                        k,Z = findZk(WD,ct)
                        failed, XX = pendlgnr(alpha,k)
                        VD = np.dot(XX,Z)
                        V[ct] = expandv(N,VD,ct)
                        W = np.dot(V[ct],W)
                        ct = ct+1
                    
                    for q in range(count(N)):   
                       LK = np.dot(V[q],LK)
       
                    U[r][t][s] = LK
       
                    return failed, U
                    
                    
                    
                    
 


def plqt(U,L,s,t):
                    D  =  [(U[0][s][t]),\
                           (U[1][(t+1)%L][s]),\
                           (U[0][(s+1)%L][t]).getH(),\
                           (U[1][t][s]).getH()]
                    return D         


def  avplqt(U,L,N):
                    sum_trp = 0.0  
                    for s  in range(L):
                        for t in range(L):
                            D = plqt(U,L,s,t)
                            UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                            trup = (1.0 - ((1.0/float(N))*np.trace(UP).real))
                            sum_trp = sum_trp + (trup/float(L*L))
                    return sum_trp  
 

def wloop11(U,L,N,s,t):
                    D = plqt(U,L,s,t)       
                    UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                    wtr =  (1.0/float(N))*np.trace(UP).real
                    return wtr 



def thermalization(N,l,titr,alpha):
            NN = count(N)
            ll = 1
            U = cold_start(l,N)
            while (ll < titr+1): 
                        totfailed = 0
                       
                        for s in range(l):
                           for t in range(l):
                                for r in range(2):
                                        failed, U = link(U,l,N,r,s,t,alpha)
                                        totfailed = totfailed + failed
                                        
                        ttr = totfailed + (NN*l*l*2)
                        success = (NN*l*l*2/float(ttr))*100
                        avp = avplqt(U,l,N)
                        print avp, "success=", success, '%'
                        plt.figure(100)
                        plt.scatter(0.0,0.0)
                        plt.scatter(0.0,1.0)
                        plt.scatter(ll,avp)
                        plt.show()
                        ll = ll+1
             
            return  


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
            
            
            
                    
def exportstorw(N,l,titr,alpha):
            sitr = 50
            storw = [0.0 for k in range(titr-sitr)]
            ll = 1
            U = cold_start(l,N)
            while (ll < titr+1): 
                    
                    for s in range(l):
                        for t in range(l):
                            for r in range(2):
                              failed,  U = link(U,l,N,r,s,t,alpha)
                    w11 = wloop11(U,l,N,5,5)
                    
                    print alpha,w11
                    if ll > sitr:
                        storw[ll-sitr-1] = w11
                    ll = ll+1
             
            return  storw
class wnr(object):
    
            def __init__ (self,filename,ver):
                        self.filename = filename
                        self.ver = ver
        

            def writer(self):
                        DataOut = np.column_stack(self.ver)    
                        np.savetxt(self.filename,DataOut)
                        return


            def reader(self):
                        colno = len(self.ver)
                        for k in range (colno):
                            self.ver[k]  = np.loadtxt(self.filename, unpack = True, usecols = [k])
                        return self.ver
            

def erorbar(N,l,titr):
            plot_dir = "/Users/dibakarsigdel/Dropbox/Plots/"  
            data_dir = "/Users/dibakarsigdel/Dropbox/Data/" 
            Nmax = 14
            Nvalue = 14
            dlamda = 0.25
            x = [0.0 for k in range(Nmax)]
            y = [0.0 for k in range(Nmax)]
            y_error = [0.0 for k in range(Nmax)]
            while Nvalue < Nmax+1:
                    lamda =  dlamda*Nvalue
                    alpha = (2.0*N)/lamda
                    x[Nvalue-1] =  (2.0*N)/alpha
                    storw = exportstorw(N,l,titr,alpha) 
                    y[Nvalue-1],y_error[Nvalue-1] = Mean_Error(storw) 
                    #print x[Nvalue-1],y[Nvalue-1]
                    plt.figure(104)
                    plt.xlabel('lambda')
                    plt.ylabel('W11')
                    plt.grid()
                    plt.errorbar(x,y, yerr = y_error, fmt='8')
                    st = str(N)
                    plt.savefig(plot_dir + 'plotsu'+st+'.png')
                    plt.show()
                    wnr('su'+st+'.dat',[x,y,y_error]).writer()
                    wnr(data_dir +'su'+st+'.dat',[x,y,y_error]).writer()
                    Nvalue = Nvalue+1
            
            return  
                                        
                                     
      

def  w_checker(N,l,titr):
            lamda =  2.15
            alpha = (2.0*N)/lamda
            storw = exportstorw(N,l,titr,alpha) 
            y,y_error = Mean_Error(storw)
            plt.figure(104)
            plt.errorbar(lamda,y,yerr = y_error,fmt ='8')
            plt.show
            return









         
###################################################################           
#Declerations------------------------
N = 7
l = 20
titr = 300

#------------------------------------------
#thermalization(N,l,titr,10.0)
#w_checker(N,l,titr)
erorbar(N,l,titr)









         
