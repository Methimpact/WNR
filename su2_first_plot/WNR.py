

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











#===========================================================================
x = [0.0 for k in range(30)]
y = [0.0 for k in range(30)]
er = [0.0 for k in range(30)]
ver = [x,y,er]
vr = [x,y]
wnr('su2comp.dat',ver).reader()
wnr('su2wtp.dat',vr).reader()
plt.figure(102)
plt.grid()
plt.title('Theoretical vs Computational results')
plt.xlabel('beta')
plt.ylabel('single plaquette wilson loop')
plt.errorbar(ver[0],ver[1],yerr = ver[2],fmt = '8',label = 'computational')
plt.plot(vr[0],vr[1],'-',label='theoretical')
plt.legend(loc='lower right')
plt.show()




