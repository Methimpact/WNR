



import matplotlib.pyplot as plt
import math as math
from scipy.special import iv
import numpy as np




class wplot(object):
    
   def __init__(self,L,K):
        self.L = L
        self.K = K
        
    
   def cgcoff(self):

            deta = math.pi/float(100)
            II_sum = [] 
            II_int = []
            
            II_int = [[[(math.sin(deta*(dx+1))**2) \
                   *(math.sin((2)*deta*(dx+1))/math.sin(deta*(dx+1)))  \
                    *(math.sin((ri+1)*deta*(dx+1))/math.sin(deta*(dx+1)))  \
                       *(math.sin((rii+1)*deta*(dx+1))/math.sin(deta*(dx+1)))*deta \
                        for dx in range(100)]\
                         for rii in range(100)] \
                            for ri in range(100)] 
                        

            II_sum =  [[sum(II_int[ri][rii]) \
                  for rii in range(100)] \
                    for ri in range(100)] 
                
                
            for ri in range(100):
               for rii in range(100):
                   if II_sum[ri][rii] > 0.5:
                          II_sum[ri][rii] = 1.0
                   else:
                          II_sum[ri][rii] = 0.0              
            

            return II_sum 


   def wclc(self,DD,alpha):
       
                xx = ((self.L*self.L)-(self.K*self.K))
                dr = [r+1 for r in range(200)]
                beta_r = [0.0 for r in range(200)]
                beta_r = [((r+1)*iv(r+1,alpha))/(iv(1,alpha))
                         for r in range(200)]
             
                factor = [0.0 for rii in range(200)]
                factor = [(beta_r[rii]/float(dr[rii])) for rii in range(200)]


                multiplier = [[0.0 for rii in range(100)] \
                             for ri in range(100)] 
                    
                    
                multiplier = [[((DD[ri][rii])*dr[rii] \
                              *(1/float(2.0*dr[ri]))*(factor[rii])**(self.K*self.K))\
                               for rii in range(100)]  \
                                  for ri in range(100)] 
                      
    
              
                sum_multiplier = [0.0 for ri in range(100)] 
                sum_multiplier = [ sum(multiplier[ri]) \
                                 for ri in range(100)] 
                       
                         

                numerator = [0.0]
                numerator = [(((factor[ri])**xx) \
                              *sum_multiplier[ri]) \
                              for ri in range(100)] 
                   
                    
  
                sum_numerator = sum(numerator)
                denominator = []                 
                denominator = [ (factor[ri])**(self.L*self.L) for ri in range(100)]

                sum_denominator = sum(denominator)

                final_value = []
                final_value = (sum_numerator/sum_denominator) 
 
                return final_value

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





#==============================================

ll = 1
x = [0.0 for k in range(30)]
y = [0.0 for k in range(30)]
while ll < 30 + 1:
  dalpha = ll*1.0
  x[ll-1]= dalpha
  xyz = wplot(30,1)
  DD = xyz.cgcoff()
  y[ll-1] = xyz.wclc(DD,dalpha)
  print y[ll-1]
  ll = ll+1
wnr('su2wtp.dat',[x,y]).writer()
plt.figure(102)
plt.plot(x,y,'-')
plt.show()





