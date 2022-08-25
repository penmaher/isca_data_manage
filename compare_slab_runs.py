import numpy as np
import pdb
from matplotlib import pyplot as plt
from data_prep import CorrectSingleFiles, print_month_update, get_eddy_terms
from cmip_vars import CMIP_compliant

class CompareRuns():
    
    def __init__(self):

        self.data_path = '/scratch/pm366/OutputIsca/'
        self.exp_list  = {'old':'ITCZ-MIP_aqua_sst_soc_low_res_take2/', 
                          'new':'ITCZ-MIP_aqua_sst_soc_low_res_take3/'} 

    def prep_new_data(self, CSF, CC):

        run_flag = '_take3'

        #step 1. Using single files, make cmip compliant 
        for month in range(CSF.run_start,CSF.run_end):
            print_month_update(month) #print every 20 months only
            CSF.run_monthly(CC, month, run_flag)

        #step 2. Need full time field for eddy flux terms
        get_eddy_terms(CSF, run_flag) #takes 10 minutes per exp

        #step 3. Open .nc file and take the zonal mean.
        CSF.monthly_zm_clima(run_flag)

    def testing(self, CSF):
        data = {}
        for opt in ['old','new']:
            filename = '{}{}itcz-mip_zm_clim.nc'.format(
                       self.data_path,self.exp_list[opt])

            data[opt] = CSF.open_file(filename)

        self.plot_surface_diff(data)
  
    def plot_surface_diff(self, data):

        lat = data['old'].lat
        for var in ['ts', 'prc']:
            diff_plot = data['new'][var] - data['old'][var]
            plt.plot(lat,data['new'][var] , label = var + ' new')
            plt.plot(lat,data['old'][var] , label = var + ' old')
            plt.plot(lat,diff_plot, label = var + ' new-old')
            plt.legend()
            plt.show()

        pdb.set_trace()

     
def main():

    CR  = CompareRuns()
    CSF = CorrectSingleFiles() 
    CC  = CMIP_compliant() 

    #for take 3 prep the files
    do_prep_new_data = False
    if do_prep_new_data:
        CR.prep_new_data(CSF, CC)
    else:
        CR.testing(CSF)
    pdb.set_trace()

if __name__ == "__main__" : 

    main()
