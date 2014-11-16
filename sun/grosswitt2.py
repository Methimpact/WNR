

import matplotlib.pyplot as plt
import CSV



N = 300
def fx(x):
    if x > 2 or x == 2:
        fnx = 1/x
    elif x < 2:
        fnx = 1-x/4.0
    return fnx
    
    
    
x = [0.0 for k in range(N)]
y = [0.0 for k in range(N)]

for k in range (N):
    xt = 0.025*k
    yt = fx(xt)
    x[k] = xt
    y[k] = yt
    
    
plt.figure(104)
plt.plot(x,y,'-')
plt.show()


U = [x,y]

CSV.CSV('wlamda.csv').csv_writer(U)
