import xarray as xr
import os
import pdb
import numpy as np
import datetime

from cmip_vars import CMIP_compliant, calc_eddy_terms, get_zonal_climatology 


#to run without warnings: python -W ignore data_prep.py

def print_month_update(month):

    if month % 20 == 0:
        print('    current month is ', month)

class CorrectSingleFiles():

    def __init__(self):

        #range of data
        self.run_start = 24  #start at 24  for a 2 year spin up
        self.run_end   = 384
        self.runfmt = 'run%04d' # ie 0001 is first run month

        self.filename_model_mm = 'atmos_monthly_interp_all.nc'
        self.filename_itcz_mm  = 'itcz-mip_control_monthly.nc'
        self.filename_itcz_mm_zm  = 'itcz-mip_zm_clim.nc'


        self.filename_model_dd = 'atmos_daily_interp_all.nc'
        self.filename_itcz_dd  = 'itcz-mip_control_daily'

    def set_data_path(self, run_flag):

        self.data_path = '{0}/{1}{2}/'.format(os.environ['GFDL_DATA'],  
                          'ITCZ-MIP_aqua_sst_soc_low_res_take3',run_flag)


    def open_file(self, file_in):

        return xr.open_dataset(file_in)

    def open_single_monthly_file(self, month, file_in):

        filename = self.data_path + self.runfmt % month + '/' + file_in
        data_in  =  xr.open_dataset(filename, decode_times=False)

        return data_in

    def open_multiple_files(self, file_in, run_flag, chunks):

        self.set_data_path(run_flag)

        file_list = [self.data_path + self.runfmt % month + '/' + file_in
                     for month in range(self.run_start, self.run_end)]

        #open open the files one layer at a time
        data_in = xr.open_mfdataset(file_list, chunks=chunks,decode_times=False)    


        return data_in
 
    def run_monthly(self, CC, month, run_flag):
        self.set_data_path(run_flag)

        data_mm = self.open_single_monthly_file(month, self.filename_model_mm)

        data_out = CC.rename_variables_monthly(data_mm)

        file_out = ('{0}{1}/{2}').format(
                    self.data_path,  self.runfmt % month, self.filename_itcz_mm)

        # output updated CMIP var names for each month in each run00XXX dir
        data_out.to_netcdf(file_out) 

        data_out.close()
        data_mm.close()

    def run_daily(self, CC, run_flag):

        data_dd = self.open_multiple_files(self.filename_model_dd,
                       run_flag, {'pfull': 1})
        CC.save_variables_daily(data_dd, self.data_path, self.filename_itcz_dd)
        data_dd.close()

    def monthly_zm_clima(self, run_flag):

        self.set_data_path(run_flag)

        filename = ('{0}{1}').format(self.data_path, self.filename_itcz_mm)
        data = self.open_file(filename)
        data_zm = get_zonal_climatology(data)

        file_out = ('{0}{1}').format(self.data_path, self.filename_itcz_mm_zm)
        data_zm.to_netcdf(file_out)

        data.close()
        data_zm.close()


def get_eddy_terms(CSF, run_flag):

    pdb.set_trace()
    data_mm = CSF.open_multiple_files(CSF.filename_itcz_mm,
                  run_flag, {'pfull': 1})

    data_out = calc_eddy_terms(data_mm)

    file_out = ('{0}{1}').format(CSF.data_path, CSF.filename_itcz_mm)

    data_out.to_netcdf(file_out)
    
    data_out.close()
    data_mm.close()

def main():     

    do_monthly       = False
    do_daily         = False

    prep_files       = False
    prep_eddy_fluxes = False
    make_zonal_clima = False

    CSF = CorrectSingleFiles() 
    CC = CMIP_compliant() 

    exp_list = [#'_m40', '_m20', '_zero', '_p20', '_p40', '',
                #'_m40_4CO2', '_m20_4CO2', '_zero_4CO2', '_p20_4CO2', 
                '_p40_4CO2']

    for run_flag in exp_list:
        print('Current exp is: ', run_flag)

        if do_monthly:
            if prep_files:
                for month in range(CSF.run_start,CSF.run_end):
                    #step 1. Using single files, make cmip compliant 
                    print_month_update(month)
                    CSF.run_monthly(CC, month, run_flag)

            elif prep_eddy_fluxes:
                #step 2. Need full time field for eddy flux terms
                get_eddy_terms(CSF, run_flag) #takes 10 minutes per exp
            elif make_zonal_clima:
                #step 3. Open .nc file and take the zonal mean.
                CSF.monthly_zm_clima(run_flag)
            else:
                print('Flag not known - what do you want to do?')
                pdb.set_trace()
        elif do_daily:
            CSF.run_daily(CC, run_flag) # this takes about 20 minutes to run
            print("Finished making daily data")
        else:
            print("No option set")
            pdb.set_trace()
    pdb.set_trace()

if __name__ == "__main__" : 

    main()
