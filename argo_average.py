#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
from netCDF4 import Dataset
import datetime
import data_utils as du



def argo_average(res,start_yr,end_yr,files,year_col = '# Year',month_col = 'Month',lat_col = 'Latitude (deg N)',lon_col = 'Longitude (deg W)',chla_col = 'chlorphyll-a (mgm-3)'
    ,skiprows=0,sep = ',',format = '.csv',dateti = False,datecol = '',dateformat=False,extra='',out_loc = '',quality_flag_col_include=False,quality_flag_col='',quality_flag_val=0,min_val_apply=False,min_val = 0.01,
    monthly=True, days=8,day_col='Day'):
    lon,lat = du.reg_grid(lat=res,lon=res)
    res_h = res/2

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
            data = np.transpose(np.vstack((np.array(data[year_col]),np.array(data[month_col]),np.array(data[day_col]),np.array(data[lat_col]),np.array(data[lon_col]),np.array(data[chla_col]),np.array(data[quality_flag_col]) )))
            f = np.where(np.isnan(data[:,-1]) == 1)[0]
            data[f,-1] = 0

            f = np.where(data[:,-1] == quality_flag_val)[0]
            data = data[f,:]
        else:
            data = np.transpose(np.vstack((np.array(data[year_col]),np.array(data[month_col]),np.array(data[day_col]),np.array(data[lat_col]),np.array(data[lon_col]),np.array(data[chla_col]))))
        if min_val_apply:
            f = np.where(data[:,5] > min_val)[0]
            data = data[f,:]
        print(len(data))

        if monthly:
            t_len = (end_yr-start_yr+1)*12
            chl = np.zeros((len(lon),len(lat),t_len))
            chl[:] = np.nan
            chl_std = np.copy(chl)
            chl_count= np.copy(chl)
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
                mean,count,std = aggregate_to_grid(data[f,3],data[f,4],np.log10(data[f,5]),lon,lat)
                # print(len(f))
                # if f.size!=0:
                #     print('Not Zeros')
                #     for j in range(0,len(lat)):
                #         for k in range(0,len(lon)):
                #             g = np.where((data[f,2] > lat[j]-res_h) & (data[f,2] < lat[j]+res_h) & (data[f,3] > lon[k]-res_h) & (data[f,3] < lon[k]+res_h))
                #             chl[k,j,i] = np.nanmean(np.log10(data[f[g],4]))
                chl[:,:,i] = np.transpose(mean)
                chl_count[:,:,i] = np.transpose(count)
                chl_std[:,:,i] = np.transpose(std)
                time.append(datetime.datetime(yr,mon,15))
                i = i+1
                mon = mon+1
                if mon == 13:
                    yr = yr+1
                    mon=1
        else:
            t_len = int((end_yr-start_yr+1)*np.ceil(365/days))
            print(t_len)
            chl = np.zeros((len(lon),len(lat),t_len))
            chl[:] = np.nan
            chl_std = np.copy(chl)
            chl_count= np.copy(chl)

            i = 0
            time = []
            for yr in range(start_yr,end_yr+1):
                print(yr)
                d = datetime.datetime(yr,1,1)

                while d.year == yr:
                    # Collect all the data that falls into the period of days.
                    t = 0
                    for j in range(0, days):
                        if (d+datetime.timedelta(days=j)).year == yr:
                            f = np.where((data[:,0] == (d+datetime.timedelta(days=j)).year) & (data[:,1] == (d+datetime.timedelta(days=j)).month) & (data[:,2] == (d+datetime.timedelta(days=j)).day))[0]
                            if t == 0:
                                f_temp = np.copy(f)
                                t = 1
                            else:
                                f_temp = np.concatenate((f_temp,f),axis=0)
                    mean,count,std = aggregate_to_grid(data[f,3],data[f,4],np.log10(data[f,5]),lon,lat)
                    chl[:,:,i] = np.transpose(mean)
                    chl_count[:,:,i] = np.transpose(count)
                    chl_std[:,:,i] = np.transpose(std)
                    time.append(d)
                    i = i+1
                    d = d + datetime.timedelta(days=days)

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
        sst_o = outp.createVariable('chl','f4',('lon','lat','time'),zlib=True)
        sst_o[:] = chl
        sst_o.units = 'log10(mgm-3)'
        sst_o.standard_name = 'Chlorophyll-a concentration'

        sst_o = outp.createVariable('chl_std','f4',('lon','lat','time'),zlib=True)
        sst_o[:] = chl_std
        sst_o.units = 'log10(mgm-3)'
        sst_o.standard_name = 'Standard Deviation Chlorophyll-a concentration'

        sst_o = outp.createVariable('chl_count','f4',('lon','lat','time'),zlib=True)
        sst_o[:] = chl_count
        sst_o.units = 'Count'
        sst_o.standard_name = 'Chlorophyll-a concentration number of observations'

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


def aggregate_to_grid(lat, lon, values,lon_g,lat_g):

    # Convert inputs to NumPy arrays
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    values = np.asarray(values)

    # Number of grid cells
    nlat = len(lat_g)
    nlon = len(lon_g)
    resolution = np.abs(lat_g[1]-lat_g[0])

    # Grid cell indices
    ilat = np.floor(
        (lat - lat_g[0]+resolution/2) / resolution
    ).astype(int)

    ilon = np.floor(
        (lon - lon_g[0]+resolution/2) / resolution
    ).astype(int)

    # Filter invalid and out-of-range points
    valid = (
        np.isfinite(lat) &
        np.isfinite(lon) &
        np.isfinite(values) &
        (ilat >= 0) & (ilat < nlat) &
        (ilon >= 0) & (ilon < nlon)
    )

    ilat = ilat[valid]
    ilon = ilon[valid]
    values = values[valid]

    # Flatten 2D indices into 1D indices
    flat_idx = ilat * nlon + ilon

    # Aggregate sums and counts
    sums = np.bincount(
        flat_idx, weights=values,
        minlength=nlat * nlon
    )

    counts = np.bincount(
        flat_idx,
        minlength=nlat * nlon
    )

    total_sq = np.bincount(
        flat_idx,
        weights=values**2,
        minlength=nlat * nlon
    )

    # Calculate mean
    mean = np.full(nlat * nlon, np.nan)
    np.divide(sums, counts, out=mean, where=counts > 0)

    variance = np.divide(
        total_sq, counts,
        where=counts > 0
    ) - mean**2

    std = np.sqrt(variance)
    # Reshape to 2D grid
    mean = mean.reshape(nlat, nlon)
    counts = counts.reshape(nlat, nlon)
    std = std.reshape(nlat,nlon)

    return mean, counts, std
