import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import cmocean
import geopandas as gpd
import matplotlib.transforms
import os
import urllib
import glob
from netCDF4 import Dataset
import datetime
import sys
# OceanICU framework functions
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU')
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU\Data_Loading')

import data_utils as du

import CCI_OC_SPATIAL_AV as oc
import cci_sstv2 as sst

def day_8_composite(outfile,res,start_yr,end_yr):
    log,lag = du.reg_grid(lat =res,lon=res)
    # oc.oc_cci_average_day('E:/Data/OC-CCI/v6.0/8day','E:/Data/OC-CCI/v6.0/8day/0.25',log = log,lag=lag,start_yr=start_yr,end_yr=end_yr)
    #sst.cci_sst_8day('E:/Data/SST-CCI/v301','E:/Data/SST-CCI/v301/8day',start_yr=1997,end_yr=2022)
    # sst.cci_sst_spatial_average('E:/Data/SST-CCI/v301/8day',out_loc='E:/Data/SST-CCI/v301/8day/0.25',start_yr=1997,end_yr=2022,log=log,lag=lag,monthly=False,v3=True)
    import construct_input_netcdf as cinp

    cinp.driver8day(outfile,[['OC-CCI','chlor_a','E:/Data/OC-CCI/v6.0/8day/0.25/%Y/ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0_0.25.nc',0],
        ['OC-CCI','chlor_a_log10_rmsd','E:/Data/OC-CCI/v6.0/8day/0.25/%Y/ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0_0.25.nc',0],
        ['OSISAF','sea_ice_fraction','E:/Data/SST-CCI/v301/8day/0.25/%Y/ESA_CCI_8DAY_SST_%Y%m%d_0.25_deg.nc',0]],lon=log,lat=lag,start_yr=start_yr,end_yr=end_yr)
