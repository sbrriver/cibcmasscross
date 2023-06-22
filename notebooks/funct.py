# import dependencies 
import emcee

import matplotlib.pyplot as plt
import numpy as np
import os, sys, shutil
import pandas as pd
import sklearn.metrics as slm

import corner
pkdir = "/pscratch/sd/s/sbrisin/boss_cib_cmass/cibcmass/hmvec"
sys.path.insert(0,pkdir)

#import hmvec as hm
from hmvec import hmvec as hm 
from scipy.interpolate import interp1d

'''
easy use functions to be able to run chains without dying
'''

def easyinit(a,b):
    import sys
    import numpy as np
    pkdir = "/pscratch/sd/s/sbrisin/boss_cib_cmass/cibcmass/hmvec"
    sys.path.insert(0,pkdir)


    #import hmvec as hm
    from hmvec import hmvec as hm 
    zs = np.linspace(a,b,20)
    ms = np.geomspace(1e10,1e17,200)
    ks = np.geomspace(1e-4,100,1001)
    hcos = hm.HaloModel(zs,ks,ms=ms)
    hcos.set_cibParams('vierro')
    # a and b are the z range
    hcos.cib_params['alpha'] = 0.2
    hcos.cib_params['beta'] = 1.6
    hcos.cib_params['gamma'] = 1.7 # not in Viero, so using Planck13
    hcos.cib_params['delta'] = 2.4
    hcos.cib_params['Td_o'] = 20.7
    logmeff = hcos.cib_params['logM_eff'] 
    hcos.cib_params['var'] = 0.3
    l0 = hcos.cib_params['L_o'] 
    hcos.add_hod(name="CMASS",mthresh=10**12+zs*0.)
    return print(f'successfully initalized with z ={a} - {b}')

def cig_mod(a,b,x,y,c,d):
    def easyinit(a,b):
        pkdir = "/pscratch/sd/s/sbrisin/boss_cib_cmass/cibcmass/hmvec"
        sys.path.insert(0,pkdir)


        #import hmvec as hm
        from hmvec import hmvec as hm 
        zs = np.linspace(a,b,20)
        ms = np.geomspace(1e10,1e17,200)
        ks = np.geomspace(1e-4,100,1001)
        hcos = hm.HaloModel(zs,ks,ms=ms)
        hcos.set_cibParams('vierro')
        # a and b are the z range
        hcos.cib_params['alpha'] = 0.2
        hcos.cib_params['beta'] = 1.6
        hcos.cib_params['gamma'] = 1.7 # not in Viero, so using Planck13
        hcos.cib_params['delta'] = 2.4
        hcos.cib_params['Td_o'] = 20.7
        logmeff = hcos.cib_params['logM_eff'] 
        hcos.cib_params['var'] = 0.3
        l0 = hcos.cib_params['L_o'] 
        hcos.add_hod(name="CMASS",mthresh=10**12+zs*0.)
        return print(f'successfully initalized with z ={a} - {b}')
    easyinit(a,b)
    for i, logmeff in enumerate(np.linspace(x,y,20)):
        for j, l0 in enumerate(np.linspace(c,d,20)):
            hcos.cib_params['L_o'] = l0
            hcos.cib_params['logM_eff'] = logmeff
            # Power spectra for CIB x galaxies
            PgI_1h = hcos.get_power_1halo('CMASS', 'cib', nu_obs=np.array([cib_freq]))
            PgI_2h = hcos.get_power_2halo('CMASS', 'cib', nu_obs=np.array([cib_freq]))

            Cl_gI_1h = hcos.C_gI(ells, hcos.zs, hcos.ks, PgI_1h, gzs = zs, gdndz= np.ones_like(zs))
            Cl_gI_2h = hcos.C_gI(ells, hcos.zs, hcos.ks, PgI_2h, gzs = zs, gdndz= np.ones_like(zs))
            tot = Cl_gI_1h + Cl_gI_2h
    return tot