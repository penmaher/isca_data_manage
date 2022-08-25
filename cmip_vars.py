import pdb
from matplotlib import pyplot as plt
import numpy as np

from eddies.eddy_terms import Eddy_Flux

latent_heat_cond  = 2.5e6 #J/kg
fill_value = 1e20 

# list of var name pdf
# https://pcmdi.llnl.gov/mips/cmip5/docs/standard_output.pdf?id=51

def area_weight_avg(data, lat, lat_axis):
    return np.average(data, weights=np.cos(np.radians(lat)), axis=lat_axis).mean()


def calc_eddy_terms(data):
    #For aquapplanets:
    # expect the zonal anomalies to be zero so no stationary terms


    var_info = {}
    var_info['ua'] = ['usvs','ucomp_vcomp',
                      'Transient eddy momenum flux','m^2/m^2']
    var_info['ta'] = ['vsts','vcomp_temp',
                      'Transient eddy heat flux','m.K/s']
    var_info['hus'] = ['vsqs','sphum_v',
                      'Transient eddy moisture flux','m.kg/s.kg']
  
    for var in ['ua','ta', 'hus']:
        EF = Eddy_Flux() 
        mmc, stat_eddy, trans_eddy = EF.get_flux_terms(data, 'va', var,
                                     data[var_info[var][1]])
        var_name = var_info[var][0]

        data = data.assign({var_name:trans_eddy})
        data[var_name].attrs['long_name'] = var_info[var][2]
        data[var_name].attrs['units']     = var_info[var][3]
        data[var_name].attrs['_FillValue']     = fill_value

    remove_list = ['ucomp_vcomp', 'sphum_v', 'vcomp_temp']
    data = data.drop(remove_list)

    return data

def get_zonal_climatology(data):

    for var in data.keys():
        #test if already zonal mean (ie eddy terms)

        if 'lon' in data[var].coords.keys():        
            data = data.assign({var:data[var].mean(dim=['time','lon'])})
        elif 'time' in data[var].coords.keys():        
            data = data.assign({var:data[var].mean(dim=['time'])})
        else:
            print(var, ' is already time mean and zonal mean')

    return data

class CMIP_compliant():

    def __init__(self):

        self.tmp = None

    def rename_variables_monthly(self, data):
        
        name_swap = {
                     #common fields
                     'temp':'ta','rh':'hur', 'sphum':'hus', 'height':'zg',
                     'precipitation':'pr', 't_surf':'ts', 
                     #winds
                     'ucomp':'ua','vcomp':'va','omega':'wap',
                     #wind stresses
                     'flux_u':'tauu','flux_v':'tauv',
                     #surface winds
                     'u_10m':'uas', 'v_10m':'vas',
                     #near surface fields
                     'temp_2m':'tas','sphum_2m':'huss',
                     #radiative properties
                     'flux_lhe':'hfls','flux_t':'hfss',
                    }

        data = data.rename(name_swap)
 
        data = self.prep_precip(data)        

        data = self.calc_lw_sw(data)

        data = self.remove_vars(data)

        return data

    def save_variables_daily(self, data, path, name):

        name_var = {'ta':'temp',  'zg':'height', 'hus':'sphum',
                    'ua':'ucomp', 'va':'vcomp',  'wap':'omega'}

        for var in ['ta', 'zg', 'hus', 'ua', 'va', 'wap']:
            print("Current variable is: ", var)
            #save 3d fields
            data_var = data[name_var[var]]
            data_var.name = var
            file_out = ('{0}{1}_{2}.nc').format(path, name, var)
            data_var.to_netcdf(file_out)
            data_var.close()
           

        #save 2d fields
        name_swap = {
                     #common fields
                     'precipitation':'pr',
                     #near surface fields
                     'temp_2m':'tas','sphum_2m':'huss',
                     #radiative properties
                     'flux_lhe':'hfls','flux_t':'hfss',
                    }
        
        data = data.rename(name_swap)
        data = self.prep_precip(data)   
        data = self.remove_vars_daily(data)
        file_out = ('{0}{1}.nc').format(path, name)
        data.to_netcdf(file_out)

    def prep_precip(self, data):

        #get the convective component
        data = data.assign({'prc':data['pr']-data['condensation_rain']})
        data['prc'].attrs['long_name'] = 'convective precipitation'
        data['prc'].attrs['units']     = data['pr'].units
        data['prc'].attrs['_FillValue']     = fill_value

        #add evaporation
        evap = data['hfls']/latent_heat_cond 

        #units [W /(m.m)]  / [J/kg] = W.kg /J.m.m = W.kg / W.s.m.m = kg/(m2.s)
        data = data.assign({'evspsbl':evap})
        data['evspsbl'].attrs['long_name'] = 'evaporation'
        data['evspsbl'].attrs['units']     = data['pr'].units
        data['evspsbl'].attrs['_FillValue']     = fill_value
        
        #now put in mm/day
        constant = 60.*60.*24.
        data['pr']      =  data['pr'] * constant
        data['prc']     =  data['prc'] * constant
        data['evspsbl'] =  data['evspsbl'] * constant

        return data

    def calc_lw_sw(self, data):

        #LW: net = up - down
        lwus = data['soc_surf_flux_lw_down']+data['soc_surf_flux_lw']      #rlds

        #SW: net = down - up
        swut = data['soc_toa_sw_down']-data['soc_toa_sw']                  #rsut

        #net = swds - swus = swds - swds*albedo 
        #sw_net = data['soc_surf_flux_sw_down']*(1-data['albedo'])            
        #swus = data['soc_surf_flux_sw_down']*data['albedo']                  #rsus
        swus = data['soc_surf_flux_sw_down'] - data['soc_surf_flux_sw']

        data = self.prep_rad_terms(data, lwus, swut, swus)

        test_rad = False
        if test_rad:
 
            rad_list_terms = {'lwut':data['rlut'], 'lwds':data['rlds'], 'lwus':lwus,
                          'swut':swut, 'swdt':data['rsdt'], 'swds':data['rsds'], 'swus':swus}

            self.test_rad_terms(rad_list_terms, data['lat'])

        return data

    def prep_rad_terms(self, data, lwus, swut, swus):

        rad_terms_rename  = {'soc_olr':'rlut', 'soc_surf_flux_lw_down':'rlds',
                             'soc_toa_sw_down':'rsdt','soc_surf_flux_sw_down':'rsds'}

        data = data.rename(rad_terms_rename)

        #remove net terms which I don't need
        rad_terms_remove = ['soc_surf_flux_lw','soc_toa_sw', 'soc_surf_flux_sw']
        data = data.drop(rad_terms_remove)

        #add rad terms to the array
        rad_list = ['rlus', 'rsut', 'rsus']
        data = data.assign({'rlus':lwus, 'rsut':swut, 'rsus':swus})

        rad_cmip_name = {'rlut': 'long wave up TOA (OLR) and long wave up surf', 
                         'rlds': 'long wave down surface',
                         'rlus': 'long wave up surface',
                         'rsut': 'short wave up TOA', 
                         'rsdt': 'short wave down TOA',
                         'rsds': 'short wave down surface',
                         'rsus': 'short wave up surface'}

        for i in rad_list:
            data[i].attrs['units']='watts/m2'  
            data[i].attrs['long_name'] = rad_cmip_name[i]
            data[i].attrs['_FillValue']     = fill_value

        return data

    def test_rad_terms(self, output, lat):

        lwut_zm = np.mean(output['lwut'][0,...], axis=-1)
        lwds_zm = np.mean(output['lwds'][0,...], axis=-1)
        lwus_zm = np.mean(output['lwus'][0,...], axis=-1)

        swut_zm = np.mean(output['swut'][0,...], axis=-1)
        swus_zm = np.mean(output['swus'][0,...], axis=-1)
        swdt_zm = np.mean(output['swdt'][0,...], axis=-1)
        swds_zm = np.mean(output['swds'][0,...], axis=-1)

        #test the area weighted average 
        lwut_avg = area_weight_avg(lwut_zm, lat, 0) #251 W/m2 earth is 240
        lwds_avg = area_weight_avg(lwds_zm, lat, 0) #293 W/m2 earth is 330
        lwus_avg = area_weight_avg(lwus_zm, lat, 0) #384 W/m2 earth is 400 
        swut_avg = area_weight_avg(swut_zm, lat, 0) #107 W/m2 earth is 100
        swus_avg = area_weight_avg(swus_zm, lat, 0) #99 W/m2 earth is 23       
        swdt_avg = area_weight_avg(swdt_zm, lat, 0) #341 W/m2 earth is 340
        swds_avg = area_weight_avg(swds_zm, lat, 0) #260 W/m2 earth is 160      

        #check conservation
        lwc = lwut_zm + lwds_zm - lwus_zm
        swa = swdt_zm - swut_zm - swds_zm - swus_zm  

        plt.plot(lat,lwut_zm,label='lwut')
        plt.plot(lat,lwds_zm,label='lwds')
        plt.plot(lat,lwus_zm,label='lwus')
        plt.plot(lat,swut_zm,label='swut')
        plt.plot(lat,swus_zm,label='swus')
        plt.plot(lat,swdt_zm,label='swdt')
        plt.plot(lat,swds_zm,label='swds')

        plt.legend()
        filename = 'rad_budget.eps'
        plt.savefig(filename)
        plt.show()
        plt.close()

        pdb.set_trace()

    def remove_vars(self, data):

        remove_list = [
                       #rad terms
                       'soc_tdt_lw', 'soc_tdt_sw','soc_tdt_rad', 'albedo',
                        #precip terms
                       'condensation_rain','zsurf',
                       #extras
                       'average_T1', 'average_T2', 'average_DT',
                       'time_bounds', 'nv', 
                      ]

        data = data.drop(remove_list)

        return data

    def remove_vars_daily(self, data):

        remove_list = [
                      'temp', 'height', 'sphum',
                      'ucomp', 'vcomp', 'omega',
                       #precip terms
                       'condensation_rain', 
                       'hfss', 'hfls',
                       #vertical interpolation terms
                       'zsurf','ps',
                       #extras
                       'average_T1', 'average_T2', 'average_DT',
                       'time_bounds', 'nv', 
                      ]
        data = data.drop(remove_list)
        data = data.drop(['lonb','latb','pfull'])


        return data

