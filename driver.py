#!/usr/bin/env python3
import sys
import os
# OceanICU framework functions
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU')
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU\Data_Loading')

from argo_chla_extract import argo_extract
from argo_average import argo_average
from oc_cci_chl_argo import chl_argo_relationship#
from oc_cci_kringing import oc_cci_fill, output_single_chla
from timeseries_plotting import timeseries_plotting,climatology_plotting,plot_flag,plot_flag_l,plot_chla_scatter
import data_utils as du
out_loc = 'E:/SCOPE/Argo/v0-2pre1'
argo_file = 'F:/Data/ARGO/argo_bio-profile_index_07082025.txt' #File containing all the current bio-argo profiles (this needs to be downloaded periodically)
argo_loc = 'F:/Data/ARGO/data' # Local location to save files.
# These are the latitude band (North first) to extract the Argo Chl-a data for...
lat_bounds = [[-40,-90],[90,45]]
output_files = [os.path.join(out_loc,'csv','argo_chla_southernocean'),os.path.join(out_loc,'csv','argo_chla_arctic')]
res = 0.25
daystep = 1
start_yr = 1997
end_yr = 2025

# bathy_file_raw = f'F:/Data/Bathymetry/GEBCO_2023.nc'
# bathy_file = f'F:/Data/Bathymetry/{res}DEG_GEBCO_2023.nc'
bathy_file_raw = os.path.join(out_loc,'netcdf','bath.nc')
bathy_file = os.path.join(out_loc,'netcdf','bath.nc')
# ice_name = 'OSISAF_sea_ice_fraction'
ice_name = 'OSISAF_ice_conc'

du.makefolder(out_loc)
du.makefolder(os.path.join(out_loc,'csv'))
du.makefolder(os.path.join(out_loc,'plots'))
du.makefolder(os.path.join(out_loc,'netcdf'))
du.makefolder(os.path.join(out_loc,'relationships'))

occci_file = os.path.join(out_loc,'netcdf',f'oc-cci_chlor_a_'+str(res)+'deg.nc')
"""
Monthy runs
"""
# argo_extract(argo_file,argo_loc,lat_bounds,output_files,out_loc=out_loc)
# argo_average(res,start_yr,end_yr,output_files,out_loc=out_loc)
# argo_average(res,start_yr,end_yr,['C:/Users/df391/OneDrive - University of Exeter/Post_Doc_ESA_SCOPE/Python/SCOPE-ArgoChla/csv/insitudb_chla_V3'],skiprows=28,sep='\t',dateti = True,datecol = 'Date/Time',dateformat='%Y-%m-%dT%H:%M',
#     lat_col = 'Latitude',lon_col = 'Longitude',chla_col ='Chl a [mg/m**3] (High Performance Liquid Chrom...)',format='.tab',out_loc=out_loc)
# argo_average(res,start_yr,end_yr,['C:/Users/df391/OneDrive - University of Exeter/Post_Doc_ESA_SCOPE/Python/SCOPE-ArgoChla/csv/insitudb_chla_V3'],skiprows=28,sep='\t',dateti = True,datecol = 'Date/Time',dateformat='%Y-%m-%dT%H:%M',
#     lat_col = 'Latitude',lon_col = 'Longitude',chla_col ='Chl a [mg/m**3] (Chlorophyll a, fluorometric o...)',format='.tab',extra='_fluoro',out_loc=out_loc)
# chl_argo_relationship(res,start_yr,end_yr,output_files,occci_file,plot=True,area_wei=True,gebco_file=bathy_file_raw,gebco_out=bathy_file,land_mask=False,netcdf_loc=out_loc)
oc_cci_fill(occci_file,bathy_file,daystep,res,reset = True,ice_name=ice_name,out_loc=out_loc) # ,outloc ='E:/SCOPE/NN/Ford_et_al_SOM_chla/inputs/chla')
timeseries_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],output=os.path.join(out_loc,'plots',f'global_map_{str(res)}_deg.png'))
climatology_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],output=os.path.join(out_loc,'plots',f'global_climatology_{str(res)}_deg.png'))
plot_flag(occci_file,output=os.path.join(out_loc,'plots',f'flag_pixels_{str(res)}_deg.png'))
plot_flag_l(occci_file,output=os.path.join(out_loc,'plots',f'flagl_pixels_{str(res)}_deg.png'))
plot_chla_scatter(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_chla_V3_{str(res)}deg.nc'),output=os.path.join(out_loc,'plots',f'insitu_verification_{str(res)}_deg.png'))
plot_chla_scatter(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_chla_V3_{str(res)}deg_fluoro.nc'),output=os.path.join(out_loc,'plots',f'insitu_verification_{str(res)}_deg_fluoro.png'))

# output_single_chla(occci_file,outloc ='F:/OceanCarbon4Climate/Experimental_Dataset/inputs/chla')
# """
# 8 day runs
# """
# from oc_cci_8day_average import day_8_composite
# occci_file = 'netcdf/oc-cci_chlor_a_'+str(res)+'deg_8day.nc'
# ice_name = 'OSISAF_sea_ice_fraction'
# daystep = 8/30.5
# # day_8_composite(occci_file,res,start_yr,end_yr)
# oc_cci_fill(occci_file,bathy_file,daystep,res,reset = True,ice_name=ice_name)
