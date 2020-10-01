
# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from operator import truediv
from scipy.optimize import fsolve
import math
from scipy.interpolate import interp1d

#%%

corr = np.array([-1., 0., 1.])

common_settings = {'output':'tCl, pCl, lCl, mPk',
                   'lensing':'yes',
                   'P_k_max_1/Mpc':100.0,
                   'n_s':0.9656
                   }

M = Class()

Pk1 = [] # P(k) in (Mpc/h)**3
Pk2 = [] # P(k) in (Mpc/h)**3
Pk3 = [] # P(k) in (Mpc/h)**3
Pk4 = []# P(k) in (Mpc/h)**3
kk = np.logspace(-2.5,2,1000) # k in h/Mpc


for i in range(4):
    M.set(common_settings)
    
    if i==0:
        M.set({'use_Beltran_cdi':'no','n_cdi':0.8963,'f_cdi':0.04911,'A_s':2.215e-9,
               'c_ad_cdi':corr[i],'n_ad_cdi':0.0})
    elif i==1:
        M.set({'use_Beltran_cdi':'no','n_cdi':2.0860097796,'f_cdi': 0.6406936,'A_s':2.215e-9,
               'c_ad_cdi':corr[i],'n_ad_cdi':0.0})
    elif i==2:
        M.set({'use_Beltran_cdi':'no','n_cdi':2.46,'f_cdi':0.1568,'A_s':2.215e-9,
               'c_ad_cdi':corr[i],'n_ad_cdi':0.0})
    elif i== 3:
        M.set({'A_s':2.215e-9})

#        M.set({'use_Beltran_cdi':'yes','n_cdi':3.0,'alpha_iso':0.1,'A_glob':3.14e-9,
#               'ellipse':corr[i],'n_ad_cdi':0.0})

    M.compute()
    derived = M.get_current_derived_parameters(['P_RR_1','P_RR_2', 'P_II_1','P_II_2'])
    
    a = 100.0*derived['P_II_1']/(derived['P_RR_1']+derived['P_II_1'])
    b = 100.0*derived['P_II_2']/(derived['P_RR_2']+derived['P_II_2'])
    

    # get P(k) at redhsift z=0
    h = M.h() # get reduced Hubble for conversions to 1/Mpc
    if i==0:
        for k in kk:
            Pk1.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print(r'For cos(Delta) = -1,  100 beta (k_1) = %f'%a)
        print(r'For cos(Delta) = -1, 100 beta (k_2) = %f'%b)
    elif i==1:
        for k in kk:
            Pk2.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print(r'For cos(Delta) = 0,  100 beta (k_1) = %f'%a)
        print(r'For cos(Delta) = 0,  100 beta (k_2) = %f'%b)
    elif i==2:
        for k in kk:
            Pk3.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        print(r'For cos(Delta) = +1,  100 beta (k_1) = %f'%a)
        print(r'For cos(Delta) = +1,  100 beta (k_2) =%f'%b)
    elif i==3:
        for k in kk:
            Pk4.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)
        
    M.struct_cleanup()
    M.empty()




#%%
# plot P(k)

plt.figure(1)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlim(kk[0],kk[-1])
plt.ylim(1e-2,1e5)
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$', fontsize=15)
plt.plot(kk,Pk1,'b--',label=r'$\mathrm{cos} \Delta= -1$')
plt.plot(kk,Pk2,'r--', label=r'$\mathrm{cos} \Delta=  0$')
plt.plot(kk,Pk3,'g--', label = r'$\mathrm{cos} \Delta= 1$')
plt.plot(kk,Pk4,'k-', label = r'$\Lambda \mathrm{CDM}$')
plt.title('Planck 2018 limits on AD+CDI models at 95 % C.L.', fontsize=13)
plt.legend(loc='best', fontsize=13)
plt.show()

#%%
plt.figure(2)
plt.xscale('log')
#plt.yscale('log')
plt.grid()
plt.xlim(kk[0],kk[-1])
#plt.ylim(1e-2,1e5)
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', fontsize=15)
plt.ylabel(r'$P_{(\mathrm{AD}+\mathrm{CDI})}/P_{(\Lambda\mathrm{CDM})}-1$', fontsize=15)
plt.plot(kk,map(truediv, list(np.array(Pk1) - np.array(Pk4)), Pk4),'b--',label=r'$\mathrm{cos} \Delta= -1$')
plt.plot(kk,map(truediv, list(np.array(Pk2) - np.array(Pk4)), Pk4),'r--', label=r'$\mathrm{cos} \Delta=  0$')
plt.plot(kk,map(truediv, list(np.array(Pk3) - np.array(Pk4)), Pk4),'g--', label = r'$\mathrm{cos} \Delta= 1$')
plt.title('Planck 2018 limits on AD+CDI models at 95 % C.L.', fontsize=13)
plt.legend(loc='best', fontsize=13)
plt.show()
#plt.savefig('pk.pdf')
