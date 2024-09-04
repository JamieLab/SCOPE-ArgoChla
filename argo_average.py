#!/usr/bin/env python

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
res = 0.25
lon,lat = du.reg_grid(lat=res,lon=res)
res_h = res/2
start_yr = 1997
end_yr = 2022
t_len = (end_yr-start_yr+1)*12

files = ['argo_chla_southernocean','argo_chla_arctic']

for file in files:
    out_file = 'netcdf/'+file+'_'+str(res)+'deg.nc'
    data = np.genfromtxt('csv/'+file+'.csv', delimiter=',')

    chl = np.zeros((len(lon),len(lat),t_len))
    chl[:] = np.nan
    print(chl.shape)
    print(data)

    yr = start_yr
    mon=1
    i=0
    time = []
    while yr<=end_yr:
        print(str(yr) + ' - ' + str(mon))
        f = np.where((data[:,0] == yr) & (data[:,1] == mon))[0]
        if f.size!=0:
            print('Not Zeros')
            for j in range(0,len(lat)):
                for k in range(0,len(lon)):
                    g = np.where((data[f,3] > lat[j]-res_h) & (data[f,3] < lat[j]+res_h) & (data[f,4] > lon[k]-res_h) & (data[f,4] < lon[k]+res_h))
                    chl[k,j,i] = np.nanmean(np.log10(data[f[g],5]))

        time.append(datetime.datetime(yr,mon,15))
        i = i+1
        mon = mon+1
        if mon == 13:
            yr = yr+1
            mon=1
    print(i)
    print(len(time))
    print(time)
    ref = datetime.datetime(1950,1,1)
    time_o = np.zeros((t_len))
    print(time_o.shape)
    for i in range(0,t_len):
        #print(i)
        time_o[i] = (time[i] - ref).days

    outp = Dataset(out_file,'w',format='NETCDF4_CLASSIC')
    outp.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
    outp.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
    outp.createDimension('lon',lon.shape[0])
    outp.createDimension('lat',lat.shape[0])
    outp.createDimension('time',t_len)
    sst_o = outp.createVariable('chl','f4',('lon','lat','time'),zlib=False)
    sst_o[:] = chl
    sst_o.units = 'log10(mgm-3)'
    sst_o.standard_name = 'Chlorophyll-a concentration'

    lat_o = outp.createVariable('latitude','f4',('lat'))
    lat_o[:] = lat
    lat_o.units = 'Degrees'
    lat_o.standard_name = 'Latitude'
    lon_o = outp.createVariable('longitude','f4',('lon'))
    lon_o.units = 'Degrees'
    lon_o.standard_name = 'Longitude'
    lon_o[:] = lon

    lon_o = outp.createVariable('time','f4',('time'))
    lon_o.units = 'days since 1950-01-01'
    lon_o.standard_name = 'Time'
    lon_o[:] = time_o
    outp.close()
