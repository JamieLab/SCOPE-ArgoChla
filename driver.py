#!/usr/bin/env python3
import sys
# OceanICU framework functions
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU')
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU\Data_Loading')

from argo_chla_extract import argo_extract
from argo_average import argo_average
from oc_cci_chl_argo import chl_argo_relationship#
from oc_cci_kringing import oc_cci_fill
from timeseries_plotting import timeseries_plotting,climatology_plotting,plot_flag,plot_flag_l,plot_chla_scatter

argo_file = 'F:/Data/ARGO/argo_bio-profile_index.txt' #File containing all the current bio-argo profiles (this needs to be downloaded periodically)
argo_loc = 'F:/Data/ARGO/data' # Local location to save files.
# These are the latitude band (North first) to extract the Argo Chl-a data for...
lat_bounds = [[-40,-90],[90,45]]
output_files = ['argo_chla_southernocean','argo_chla_arctic']
res = 0.25
daystep = 1
start_yr = 1997
end_yr = 2023
occci_file = 'netcdf/oc-cci_chlor_a_'+str(res)+'deg.nc'
bathy_file = f'F:/Data/Bathymetry/{res}DEG_GEBCO_2023.nc'
# ice_name = 'OSISAF_sea_ice_fraction'
ice_name = 'OSISAF_ice_conc'

"""
Monthy runs
"""
# daystep = 1
# argo_extract(argo_file,argo_loc,lat_bounds,output_files)
# argo_average(res,start_yr,end_yr,output_files)
# argo_average(res,start_yr,end_yr,['insitudb_chla_V3'],skiprows=28,sep='\t',dateti = True,datecol = 'Date/Time',dateformat='%Y-%m-%dT%H:%M',
#     lat_col = 'Latitude',lon_col = 'Longitude',chla_col ='Chl a [mg/m**3] (High Performance Liquid Chrom...)',format='.tab')
# argo_average(res,start_yr,end_yr,['insitudb_chla_V3'],skiprows=28,sep='\t',dateti = True,datecol = 'Date/Time',dateformat='%Y-%m-%dT%H:%M',
#     lat_col = 'Latitude',lon_col = 'Longitude',chla_col ='Chl a [mg/m**3] (Chlorophyll a, fluorometric o...)',format='.tab',extra='_fluoro')
# # chl_argo_relationship(res,start_yr,end_yr,output_files,occci_file,plot=True)
# oc_cci_fill(occci_file,bathy_file,daystep,res,reset = True,ice_name=ice_name)
# timeseries_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],output=f'plots/global_map_{str(res)}_deg.png')
# climatology_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],output=f'plots/global_climatology_{str(res)}_deg.png')
# plot_flag(occci_file,output=f'plots/flag_pixels_{str(res)}_deg.png')
# plot_flag_l(occci_file,output=f'plots/flagl_pixels_{str(res)}_deg.png')
plot_chla_scatter(occci_file,f'netcdf/insitudb_chla_V3_{str(res)}deg.nc',output=f'plots/insitu_verification_{str(res)}_deg.png')
plot_chla_scatter(occci_file,f'netcdf/insitudb_chla_V3_{str(res)}deg_fluoro.nc',output=f'plots/insitu_verification_{str(res)}_deg_fluoro.png')

# """
# 8 day runs
# """
# from oc_cci_8day_average import day_8_composite
# occci_file = 'netcdf/oc-cci_chlor_a_'+str(res)+'deg_8day.nc'
# ice_name = 'OSISAF_sea_ice_fraction'
# daystep = 8/30.5
# # day_8_composite(occci_file,res,start_yr,end_yr)
# oc_cci_fill(occci_file,bathy_file,daystep,res,reset = True,ice_name=ice_name)
