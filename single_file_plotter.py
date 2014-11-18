# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:27:13 2014

@author: dibakarsigdel
"""



import numpy as np
import matplotlib.pyplot as plt






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



def fx(x):
    if x > 2 or x == 2:
        fnx = 1/x
    elif x < 2:
        fnx = 1-x/4.0
    return fnx
    
def databox(N): 
        xx = []
        yy = []
        for k in range (N):
            xt = 0.025*k
            yt = fx(xt)
            xx.append(xt)
            yy.append(yt)
        return  xx,yy
    








#===========================================================================
x = [0.0 for k in range(30)]
y = [0.0 for k in range(30)]
er = [0.0 for k in range(30)]
ver = [x,y,er]
xx,yy = databox(300)
wnr('su2.dat',ver).reader()
plt.figure(102)
plt.grid()
plt.title('Theoretical versus Computational')
plt.xlabel('beta')
plt.ylabel('single plaquette wilson loop')
plt.errorbar(ver[0],ver[1],yerr = ver[2],fmt = '8',label = 'computational')
plt.plot(xx,yy,'-',label='theoritical')
plt.legend(loc='lower right')
plt.show()




