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

import data_utils as du

import CCI_OC_SPATIAL_AV as oc
import OSISAF_download as OSI
import cci_sstv2 as sst

def occci_8day_ftp_download(download_loc):
    import ftplib
    mysite = "oceancolour.org"
    username = "oc-cci-data"
    password = "ELaiWai8ae"
    remote_dir = "/occci-v6.0/geographic/netcdf/8day/chlor_a"
    local_dir = ""
    pattern = ".*\.nc$"

    ftp = ftplib.FTP(mysite, username, password)
    download_ftp_tree(ftp, remote_dir, local_dir, pattern=pattern, overwrite=False, guess_by_extension=True)

def generate_occci_8day(occci_file,res,start_yr,end_yr,area_wei=False,gebco_file = False,gebco_out=False,land_mask=False,chla_loc='F:/Data/OC-CCI/monthly/chlor_a',cci_land_file='E:/Data/Land-CCI/ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0',
    osisaf_loc='F:/Data/OSISAF',file_struct='',file_struct_out='',process=1,days=8,osisaf_file='',osisaf_file_out=''):
    import CCI_OC_SPATIAL_AV as OC
    import Data_Loading.ESA_CCI_land as landcci
    lon,lat = du.reg_grid(lon=res,lat=res)
    res = np.round(np.abs(lon[0]-lon[1]),2)
    landcci.generate_land_cci(cci_land_file+'.nc',cci_land_file+'.tif',lon,lat,gebco_out)

    du.makefolder(os.path.join(chla_loc,str(res)+'DEG_weighted'))

    # OC.oc_cci_average_8day(chla_loc,os.path.join(chla_loc,str(res)+'DEG_weighted'),start_yr = start_yr,end_yr = end_yr,log=lon,lag=lat,
    #     area_wei=area_wei,gebco_file = gebco_file,gebco_out=gebco_out,land_mask=land_mask,file_struct=file_struct,file_struct_out=file_struct_out,process=process,days=days)

    # hemispheres=['nh','sh']
    # for hemi in hemispheres:
    #     #OSI.OSISAF_daily_av(osisaf_loc,os.path.join(osisaf_loc,str(days)+'_day'),start_yr=start_yr,end_yr=end_yr,file_struct=osisaf_file,file_struct_out=osisaf_file_out,days=days,process=process,hemi=hemi)
    #     OSI.OSISAF_spatial_average_daily(os.path.join(osisaf_loc,str(days)+'_day'),start_yr=start_yr,end_yr=end_yr,out_loc=os.path.join(osisaf_loc,str(days)+'_day',str(res)+'DEG'),log=lon,lag=lat,hemi=[hemi],processes=process,days=days,
    #         file_struct=osisaf_file_out)
    # OSI.OSISAF_merge_hemisphere_daily(data=os.path.join(osisaf_loc,str(days)+'_day',str(res)+'DEG'),start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,processes=process,days=days,file_struct=osisaf_file_out,bathy_file=gebco_out)

    import construct_input_netcdf as cinp
    cinp.driver8day(occci_file,[['OC-CCI','chlor_a',os.path.join(chla_loc,str(res)+'DEG_weighted','%Y',file_struct_out.replace('<res>',str(res))),0],
        ['OC-CCI','chlor_a_log10_rmsd',os.path.join(chla_loc,str(res)+'DEG_weighted','%Y',file_struct_out.replace('<res>',str(res))),0],
        ['OSISAF','ice_conc',os.path.join(osisaf_loc,str(days)+'_day',str(res)+'DEG','%Y',(osisaf_file_out.split('.')[0]+'_'+str(res)+'_deg.nc').replace('<hemi>','COM').replace('<days>',str(days))),0]]
            ,lon=lon,lat=lat,start_yr=start_yr,end_yr=end_yr,copts={'zlib':True})

# def day_8_composite(outfile,res,start_yr,end_yr,oc):
#     log,lag = du.reg_grid(lat =res,lon=res)
#     oc.oc_cci_average_day('E:/Data/OC-CCI/v6.0/8day','E:/Data/OC-CCI/v6.0/8day/0.25',log = log,lag=lag,start_yr=start_yr,end_yr=end_yr)
#     #sst.cci_sst_8day('E:/Data/SST-CCI/v301','E:/Data/SST-CCI/v301/8day',start_yr=1997,end_yr=2022)
#     # sst.cci_sst_spatial_average('E:/Data/SST-CCI/v301/8day',out_loc='E:/Data/SST-CCI/v301/8day/0.25',start_yr=1997,end_yr=2022,log=log,lag=lag,monthly=False,v3=True)
#     # import construct_input_netcdf as cinp
#     #
#     # cinp.driver8day(outfile,[['OC-CCI','chlor_a','E:/Data/OC-CCI/v6.0/8day/0.25/%Y/ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0_0.25.nc',0],
#     #     ['OC-CCI','chlor_a_log10_rmsd','E:/Data/OC-CCI/v6.0/8day/0.25/%Y/ESACCI-OC-L3S-CHLOR_A-MERGED-8D_DAILY_4km_GEO_PML_OCx-%Y%m%d-fv6.0_0.25.nc',0],
#     #     ['OSISAF','sea_ice_fraction','E:/Data/SST-CCI/v301/8day/0.25/%Y/ESA_CCI_8DAY_SST_%Y%m%d_0.25_deg.nc',0]],lon=log,lat=lag,start_yr=start_yr,end_yr=end_yr)
