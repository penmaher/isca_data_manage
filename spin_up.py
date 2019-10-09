import xarray as xr
import os
import pdb
import numpy as np

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

def plot_equatorial_temp_profile(temp, rh, plev):

    #expet the model to take a year or so spinning up from initial conditions
    # and using gray radiation.

    t_zm = np.mean(temp[:,:,int(temp.shape[2]/2),:], axis=-1)

    fig = plt.figure()
    ax = plt.axes()

    ax.plot(t_zm[0,:], plev, label='First month')
    ax.plot(t_zm[12,:], plev,label='1 year')
    ax.plot(t_zm[24,:], plev, label='2 years')
    ax.plot(t_zm[36,:], plev, label='3 years')
    ax.plot(t_zm[48,:], plev, label='4 years')
    ax.plot(t_zm[60,:], plev, label='5 years')
    ax.plot(t_zm[-1,:], plev, label='last month of sim') 

    ax.set_ylim([1000,1])
    plt.legend()


    #get a feel for when temp structure has settles and then look at tropopause
    #lower strasphere water content.


    rh_zm = np.mean(rh[:,:,int(rh.shape[2]/2),:], axis=-1)

    fig = plt.figure()
    ax = plt.axes()

    ax.plot(rh_zm[0,:], plev, label='First month')
    ax.plot(rh_zm[12,:], plev,label='1 year')
    ax.plot(rh_zm[24,:], plev, label='2 years')
    ax.plot(rh_zm[36,:], plev, label='3 years')
    ax.plot(rh_zm[48,:], plev, label='4 years')
    ax.plot(rh_zm[60,:], plev, label='5 years')

    ax.plot(rh_zm[-1,:], plev, label='last month of sim') 

    plt.legend()
    ax.set_ylim([1000,1])

    plt.show()

def test_model_spin_up(temp, rh, plev):     

    plot_equatorial_temp_profile(temp, rh, plev)

    pdb.set_trace()




