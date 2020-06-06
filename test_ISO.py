
# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from scipy.interpolate import interp1d
#%%

settings = {'output':'tCl, pCl, lCl, mPk',
                   'lensing':'yes',
                   'P_k_max_1/Mpc':3.0,
                   'n_s':0.9656,
                   'use_Beltran_cdi':'yes',
                   'n_cdi':3.0,
                   'alpha_iso':0.3,
                   'A_glob':3.14e-9,
                   'ellipse':0.2,
                   'n_ad_cdi':0.0,
                   }

M = Class()

M.set(settings)
M.compute()

derived = M.get_current_derived_parameters(['P_RR_1','P_RR_2', 'P_II_1','P_II_2','P_RI_1','P_RI_2'])


print('P_RR_1 = %e'%derived['P_RR_1'])
print('P_RR_2 = %e'%derived['P_RR_2'])
print('P_II_1 = %e'%derived['P_II_1'])
print('P_II_2 = %e'%derived['P_II_2'])
print('P_RI_1 = %e'%derived['P_RI_1'])
print('P_RI_2 = %e'%derived['P_RI_2'])
#%%
# get P(k) at redhsift z=0
kk = np.logspace(-4,np.log10(3),1000) # k in h/Mpc
Pk = [] # P(k) in (Mpc/h)**3
h = M.h() # get reduced Hubble for conversions to 1/Mpc
for k in kk:
    Pk.append(M.pk(k*h,0.)*h**3) # function .pk(k,z)

#%%
# plot P(k)

plt.figure(2)
plt.xscale('log');plt.yscale('log');plt.xlim(kk[0],kk[-1])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
plt.plot(kk,Pk,'b-')
plt.show()
#plt.savefig('pk.pdf')

M.struct_cleanup()
M.empty()
