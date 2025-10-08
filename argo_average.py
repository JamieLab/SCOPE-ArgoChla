#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
from netCDF4 import Dataset
import datetime
import data_utils as du
# res = 0.25


# start_yr = 1997
# end_yr = 2023


# files = ['argo_chla_southernocean','argo_chla_arctic']
def argo_average(res,start_yr,end_yr,files,year_col = '# Year',month_col = 'Month',lat_col = 'Latitude (deg N)',lon_col = 'Longitude (deg W)',chla_col = 'chlorphyll-a (mgm-3)'
    ,skiprows=0,sep = ',',format = '.csv',dateti = False,datecol = '',dateformat=False,extra='',out_loc = '',quality_flag_col_include=False,quality_flag_col='',quality_flag_val=0,min_val_apply=False,min_val = 0.01):
    lon,lat = du.reg_grid(lat=res,lon=res)
    res_h = res/2
    t_len = (end_yr-start_yr+1)*12
    for file in files:
        file2 = os.path.split(file)
        out_file = os.path.join(out_loc,'netcdf',file2[1]+'_'+str(res)+'deg'+extra+'.nc')
        # data = np.genfromtxt('csv/'+file+'.csv', delimiter=',')
        data = pd.read_table(file+format,sep = sep,skiprows=skiprows)
        print(data)
        if dateti:
            time = data[datecol]
            temp = np.zeros((np.array(time).size,3)); temp[:] = np.nan
            for i in range(time.size):
                tt = datetime.datetime.strptime(time[i],dateformat)
                temp[i,0] = tt.year
                temp[i,1] = tt.month
                temp[i,2] = tt.day
            data['# Year'] = temp[:,0]
            data['Month'] = temp[:,1]
            data['Day'] = temp[:,2]
        print(len(data))
        if quality_flag_col_include:
            data = np.transpose(np.vstack((np.array(data[year_col]),np.array(data[month_col]),np.array(data[lat_col]),np.array(data[lon_col]),np.array(data[chla_col]),np.array(data[quality_flag_col]) )))
            f = np.where(np.isnan(data[:,-1]) == 1)[0]
            data[f,-1] = 0

            f = np.where(data[:,-1] == quality_flag_val)[0]
            data = data[f,:]
        else:
            data = np.transpose(np.vstack((np.array(data[year_col]),np.array(data[month_col]),np.array(data[lat_col]),np.array(data[lon_col]),np.array(data[chla_col]))))
        if min_val_apply:
            f = np.where(data[:,4] > min_val)[0]
            data = data[f,:]
        print(len(data))

        chl = np.zeros((len(lon),len(lat),t_len))
        chl[:] = np.nan
        print(chl.shape)
        print(data)
        print(data.shape)

        yr = start_yr
        mon=1
        i=0
        time = []
        while yr<=end_yr:
            print(str(yr) + ' - ' + str(mon))
            f = np.where((data[:,0] == yr) & (data[:,1] == mon))[0]
            print(len(f))
            if f.size!=0:
                print('Not Zeros')
                for j in range(0,len(lat)):
                    for k in range(0,len(lon)):
                        g = np.where((data[f,2] > lat[j]-res_h) & (data[f,2] < lat[j]+res_h) & (data[f,3] > lon[k]-res_h) & (data[f,3] < lon[k]+res_h))
                        chl[k,j,i] = np.nanmean(np.log10(data[f[g],4]))

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
