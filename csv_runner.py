# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 17:12:56 2014

@author: nauraj
"""
import CSV_module
import matplotlib.pyplot as plt

nr =  100
nc =  21
nnr = 42
nnc = 21

Puredata = CSV_module.CSV('data.csv').csv_reader(nr,nc,nnr,nnc)
         
[SN ,Galaxy ,D_Mpc , U_B_0 ,B_v_0 ,M_V ,\
      log_LHalpha, RH_arcm, PA1, PA2, RH_kpc ,\
      b_by_a, Center, Center_RA, Dec_dms ,\
      log_SFRD,log_Mgas,RD_unc,RD_arcmin,RD_kpc, SFR] = Puredata         


nau_ver = [SN,Galaxy,D_Mpc]

CSV_module.CSV('nau.csv').csv_writer(nau_ver)

    
plt.figure(1)
plt.scatter(SN,D_Mpc)   
plt.show()
    
    
    
