#!/usr/bin/env python3
import sys
import os
# OceanICU framework functions
sys.path.append('/home/RSims/UExP-FNN-U/code')
sys.path.append('/home/RSims/UExP-FNN-U/code/Data_Loading')

from argo_chla_extract import argo_extract
from argo_average import argo_average
from oc_cci_chl_argo import chl_argo_relationship, plot_spatial_relationship,download_occci_monthly,download_osisaf_daily
from oc_cci_kringing import oc_cci_fill, output_single_chla
from timeseries_plotting import timeseries_plotting,climatology_plotting,plot_flag,plot_flag_l,plot_chla_scatter
from argo_histogram import plot_histogram
from daily_argo_compare import daily_argo_match, plot_daily_argo,plot_res_comparision
from plot_manuscript_independent_validation import plotting_independent_wintertime, plotting_independent_cloud,plotting_independent_final_krig,plotting_independent_underice
import data_utils as du

version = 'v1-0-test'
loc = '/data/DanFord/SCOPE_Argo_chla'
out_loc = os.path.join(loc,version)
argo_file = os.path.join(loc,'inputs','argo','argo_bio-profile_index_15032026.txt') #File containing all the current bio-argo profiles (this needs to be downloaded periodically)
argo_loc = os.path.join(loc,'inputs','argo') # Local location to save files.

# These are the latitude band (North first) to extract the Argo Chl-a data for...
lat_bounds = [[-40,-90],[90,45]]
output_files = [os.path.join(out_loc,'csv','argo_chla_southernocean'),os.path.join(out_loc,'csv','argo_chla_arctic')]
res = 0.25
daystep = 1
start_yr = 1997
end_yr = 2025


bathy_file_raw = os.path.join(out_loc,'netcdf',f'bath_{res}.nc')
bathy_file = os.path.join(out_loc,'netcdf',f'bath_{res}.nc')
cci_land_file = os.path.join(loc,'inputs','cci_land','ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0')

ice_name = 'OSISAF_ice_conc'
insitu_chla = os.path.join(loc,'inputs','insitu_chla','insitudb_Chla_v2026')

du.makefolder(out_loc)
du.makefolder(os.path.join(out_loc,'csv'))
du.makefolder(os.path.join(out_loc,'plots'))
du.makefolder(os.path.join(out_loc,'netcdf'))
du.makefolder(os.path.join(out_loc,'relationships'))

occci_file = os.path.join(out_loc,'netcdf',f'oc-cci_chlor_a_'+str(res)+'deg.nc')
occci_daily_loc = 'E:/Data/OC-CCI/v6.0/daily'
"""
Monthy runs
"""
# argo_extract(argo_file,argo_loc,lat_bounds,output_files,out_loc=out_loc,start_yr=start_yr,end_yr=end_yr)
# plot_histogram(out_loc,output_files,start_yr,end_yr)

# plot_daily_argo(out_loc,output_files)
# daily_argo_match(out_loc,output_files,occci_daily_loc)

# argo_average(res,start_yr,end_yr,output_files,out_loc=out_loc)

# argo_average(res,start_yr,end_yr,[insitu_chla],skiprows=0,sep=',',dateti = True,datecol = 'time',dateformat='%Y-%m-%dT%H:%M:%SZ',
#     lat_col = 'lat',lon_col = 'lon',chla_col ='chla_hplc',format='.csv',out_loc=out_loc,quality_flag_col_include=False,quality_flag_col='QF Chl (1 if observation as unknown c...)',quality_flag_val=0)
# argo_average(res,start_yr,end_yr,[insitu_chla],skiprows=0,sep=',',dateti = True,datecol = 'time',dateformat='%Y-%m-%dT%H:%M:%SZ',
#     lat_col = 'lat',lon_col = 'lon',chla_col ='chla_fluor',format='.csv',extra='_fluoro',out_loc=out_loc,quality_flag_col_include=False,quality_flag_col='QF Chl (1 if observation as unknown c...)',quality_flag_val=0,
#     min_val_apply=True,min_val=0.01)

# download_occci_monthly(os.path.join(loc,'inputs','oc-cci_monthly'),start_yr,end_yr)
# download_osisaf_daily(os.path.join(loc,'inputs','osisaf_daily'),start_yr,end_yr)
# generate_occci(occci_file,res,start_yr,end_yr,area_wei=True,gebco_file=bathy_file_raw,gebco_out=bathy_file,land_mask=False,netcdf_loc=out_loc,chla_loc=os.path.join(loc,'inputs','oc-cci_monthly'),cci_land_file=os.path.join(loc,'inputs','cci_land','ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0'),osisaf_loc=os.path.join(loc,'inputs','osisaf_daily'),multiprocess=True)
# chl_argo_relationship(res,start_yr,end_yr,output_files,occci_file,plot=True,area_wei=True,gebco_file=bathy_file_raw,gebco_out=bathy_file,land_mask=False,netcdf_loc=out_loc,chla_loc=os.path.join(loc,'inputs','oc-cci_monthly'),cci_land_file=os.path.join(loc,'inputs','cci_land','ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0'),osisaf_loc=os.path.join(loc,'inputs','osisaf_daily'),multiprocess=True)
# plot_spatial_relationship(res,output_files,netcdf_loc=out_loc)
# plot_res_comparision(out_loc,output_files,res)
# oc_cci_fill(occci_file,bathy_file,daystep,res,reset =True,ice_name=ice_name,out_loc=out_loc,plot_debug=False,plot_debug_loc=os.path.join(out_loc,'plots','debug')) # ,outloc ='E:/SCOPE/NN/Ford_et_al_SOM_chla/inputs/chla')
# timeseries_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],start_yr,end_yr,output=os.path.join(out_loc,'plots',f'global_map_{str(res)}_deg.png'))
# climatology_plotting(occci_file,bathy_file,[[-60,0],[-55,0],[-50,0],[-45,0]],output=os.path.join(out_loc,'plots',f'global_climatology_{str(res)}_deg.png'))
# plot_flag(occci_file,output=os.path.join(out_loc,'plots',f'flag_pixels_{str(res)}_deg.png'),res=res)
# plot_flag_l(occci_file,output=os.path.join(out_loc,'plots',f'flagl_pixels_{str(res)}_deg.png'),res=res)
# plot_chla_scatter(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg.nc'),output=os.path.join(out_loc,'plots',f'insitu_verification_{str(res)}_deg.png'))
# plot_chla_scatter(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg_fluoro.nc'),output=os.path.join(out_loc,'plots',f'insitu_verification_{str(res)}_deg_fluoro.png'))
# #
# plotting_independent_wintertime(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg.nc'),os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg_fluoro.nc'),out_loc,res)
# plotting_independent_cloud(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg.nc'),os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg_fluoro.nc'),out_loc,res)
# plotting_independent_final_krig(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg.nc'),os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg_fluoro.nc'),out_loc,res)
# plotting_independent_underice(occci_file,os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg.nc'),os.path.join(out_loc,'netcdf',f'insitudb_Chla_v2026_{str(res)}deg_fluoro.nc'),out_loc,res)
#
# output_single_chla(occci_file,outloc ='F:/OceanCarbon4Climate/GCB2026/inputs/chla')


# """
# 8 day runs
# """#
version = 'v1-0-8daytest'
loc = '/data/DanFord/SCOPE_Argo_chla'
out_loc = os.path.join(loc,version)
argo_file = os.path.join(loc,'inputs','argo','argo_bio-profile_index_15032026.txt') #File containing all the current bio-argo profiles (this needs to be downloaded periodically)
argo_loc = os.path.join(loc,'inputs','argo') # Local location to save files.

# These are the latitude band (North first) to extract the Argo Chl-a data for...
lat_bounds = [[-40,-90],[90,45]]
output_files = [os.path.join(out_loc,'csv','argo_chla_southernocean'),os.path.join(out_loc,'csv','argo_chla_arctic')]
res = 0.25
start_yr = 1997
end_yr = 2025
days = 8


bathy_file_raw = os.path.join(out_loc,'netcdf',f'bath_{res}.nc')
bathy_file = os.path.join(out_loc,'netcdf',f'bath_{res}.nc')
cci_land_file = os.path.join(loc,'inputs','cci_land','ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0')

ice_name = 'OSISAF_ice_conc'
insitu_chla = os.path.join(loc,'inputs','insitu_chla','insitudb_Chla_v2026')

du.makefolder(out_loc)
du.makefolder(os.path.join(out_loc,'csv'))
du.makefolder(os.path.join(out_loc,'plots'))
du.makefolder(os.path.join(out_loc,'netcdf'))
du.makefolder(os.path.join(out_loc,'relationships'))

occci_file = os.path.join(out_loc,'netcdf',f'oc-cci_8day_chlor_a_'+str(res)+'deg.nc')
chla_8day_loc = os.path.join(loc,'inputs','oc-cci_8day')

from oc_cci_chl_argo import download_occci_8day
# download_occci_8day(chla_8day_loc,start_yr,end_yr)
from oc_cci_8day_average import generate_occci_8day
# generate_occci_8day(occci_file,res,start_yr,end_yr,area_wei=True,gebco_file=bathy_file_raw,gebco_out=bathy_file,land_mask=False,chla_loc=chla_8day_loc,
#     cci_land_file=cci_land_file,osisaf_loc=os.path.join(loc,'inputs','osisaf_daily'),file_struct='ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0.nc'
#     ,file_struct_out='ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0_<res>_deg.nc',process=8,days=days,osisaf_file=os.path.join('%m','*<hemi>*%Y%m%d*.nc'),
#     osisaf_file_out='OSISAF_<days>day_%Y%m%d_<hemi>.nc')

# argo_extract(argo_file,argo_loc,lat_bounds,output_files,out_loc=out_loc,start_yr=start_yr,end_yr=end_yr)
# argo_average(res,start_yr,end_yr,output_files,out_loc=out_loc,monthly=False,days=days)
chl_argo_relationship(res,start_yr,end_yr,output_files,occci_file,plot=True,netcdf_loc=out_loc,lags=16)
