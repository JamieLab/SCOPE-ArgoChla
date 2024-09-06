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

def save_chla(file,chla,flag,flag_l,lat,lon,time,out_loc,time_r):
    outp = Dataset(file,'w',format='NETCDF4_CLASSIC')
    outp.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
    outp.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
    outp.createDimension('lon',lon.shape[0])
    outp.createDimension('lat',lat.shape[0])
    outp.createDimension('time',len(time))

    sst_o = outp.createVariable('chl','f4',('lon','lat','time'),zlib=False)
    sst_o[:] = chla
    sst_o.units = 'log10(mgm-3)'
    sst_o.standard_name = 'Chlorophyll-a concentration'

    sst_o = outp.createVariable('flag','f4',('lon','lat','time'),zlib=False)
    sst_o[:] = flag
    sst_o.standard_name = 'Chlorophyll-a concentration flag'

    sst_o = outp.createVariable('flag_l','f4',('lon','lat','time'),zlib=False)
    sst_o[:] = flag_l
    sst_o.standard_name = 'Chlorophyll-a concentration flag for number of months used in Argo approach'

    lat_o = outp.createVariable('latitude','f4',('lat'))
    lat_o[:] = lat
    lat_o.units = 'Degrees'
    lat_o.standard_name = 'Latitude'
    lon_o = outp.createVariable('longitude','f4',('lon'))
    lon_o.units = 'Degrees'
    lon_o.standard_name = 'Longitude'
    lon_o[:] = lon

    lon_o = outp.createVariable('time','f4',('time'))
    lon_o.units = 'days since 1970-01-15'
    lon_o.standard_name = 'Time'
    lon_o[:] = time
    outp.close()
    if outloc:
        for i in range(len(time)):
            d = datetime.datetime(int(time_r[i,0]),int(time_r[i,1]),15)
            du.makefolder(os.path.join(out_loc,d.strftime('%Y')))
            file_t = os.path.join(out_loc,d.strftime('%Y/%Y_%m_oc_cci_chla_filled.nc'))

            outp = Dataset(file_t,'w',format='NETCDF4_CLASSIC')
            outp.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
            outp.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
            outp.createDimension('lon',lon.shape[0])
            outp.createDimension('lat',lat.shape[0])

            sst_o = outp.createVariable('chl','f4',('lon','lat'),zlib=False)
            sst_o[:] = chla[:,:,i]
            sst_o.units = 'log10(mgm-3)'
            sst_o.standard_name = 'Chlorophyll-a concentration'

            sst_o = outp.createVariable('flag','f4',('lon','lat'),zlib=False)
            sst_o[:] = flag[:,:,i]
            sst_o.standard_name = 'Chlorophyll-a concentration flag'

            lat_o = outp.createVariable('latitude','f4',('lat'))
            lat_o[:] = lat
            lat_o.units = 'Degrees'
            lat_o.standard_name = 'Latitude'
            lon_o = outp.createVariable('longitude','f4',('lon'))
            lon_o.units = 'Degrees'
            lon_o.standard_name = 'Longitude'
            lon_o[:] = lon
            outp.close()
res = 1
files = ['argo_chla_southernocean','argo_chla_arctic']
argo_b = np.loadtxt('relationships/backwards_'+files[0]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
argo_f = np.loadtxt('relationships/forwards_'+files[0]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)

argo_b_arctic = np.loadtxt('relationships/backwards_'+files[1]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
argo_f_arctic = np.loadtxt('relationships/forwards_'+files[1]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)

occci_file = 'netcdf/oc-cci_chlor_a_'+str(res)+'deg.nc'
outfile = 'output/oc-cci_chla_corrected_'+str(res)+'deg.nc'
outloc = 'E:/SCOPE/NN/Ford_et_al_SOM_chla/inputs/chla'
#outloc = False
c = Dataset(occci_file,'r')
chla = np.array(c['OC-CCI_chlor_a'])
ice = np.array(c['OSISAF_ice_conc'])
lat = np.array(c['latitude'])
lon = np.array(c['longitude'])
time = np.array(c['time'])
c.close()

ref = datetime.datetime(1970,1,15)
time_r = np.zeros((len(time),2))
for i in range(0,len(time)):
    time_r[i,0] = (ref + datetime.timedelta(days=int(time[i]))).year
    time_r[i,1] = (ref + datetime.timedelta(days=int(time[i]))).month
#print(time_r)
f = np.where((time_r[:,0] == 1997) & (time_r[:,1] == 9))[0]
#print(f)
chla = chla[:,:,f[0]:]
ice = ice[:,:,f[0]:]
time_r = time_r[f[0]:,:]
time = time[f[0]:]


c = Dataset('D:/Data/Bathymetry/'+str(res)+'DEG_GEBCO_2023.nc','r')
ocean = np.array(c['ocean_proportion'])
c.close()

"""
Starting the flagging and filling of the data...
"""
# Set up the flag mask so that if a pixel is not dealt with it will have a flag of 0
# Aim is to have no zeros in the flag array :-)
flag = np.zeros((chla.shape))
flag_l = np.copy(flag)
# Set the actual OC-CCI observations to 2
flag[np.isnan(chla) == 0] = 2
# Set the land pixel flags to 1 based on the ocean proportion data generated at 1 deg from
# GEBCO data
flag[ocean==0] = 1
chla[ocean==0] = np.nan # Just setting the land pixels to nan in the chl-a array

#Now the under ice values - so if chl-a value is nan, and the ice coveraged from OSISAF is greater than 90%
# we fill with a set value. And set the flag to 3.
f = np.where((np.isnan(chla) == 1) & (ice >= 0.9))
chla[f] = np.log10(0.05)
flag[f] = 3

# Filling gaps in the subtropics due to cloud cover. This has a larger effect in the early years when SeaWiFS was
# the only satellite up. Mainly affects the Gulf of Guinea (equatorial atlantic), Arabian Sea and the Indonesian regions.
# We only apply these averaging apporaches between 40S and 40N as we have another approach to fill the wintertime chl-a values
# in the polar regions.
f = np.where((flag == 0) & ((lat[np.newaxis,:,np.newaxis] > -40) & (lat[np.newaxis,:,np.newaxis] < 40)))
print(len(f[0]))
for i in range(0,len(f[0])):
    if (f[0][i] == lon[0]) |  (f[0][i] == lon[-1]):
        print('Lon skip')
    else:
        ch = np.sum(flag[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,f[2][i]] == 2)
        print(ch)
        if ch > 6:
            chla[f[0][i],f[1][i],f[2][i]] = np.nanmean(chla[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,f[2][i]])
            flag[f[0][i],f[1][i],f[2][i]] = 4
        else:
            if ch > 4:
                chla[f[0][i],f[1][i],f[2][i]] = np.nanmean(chla[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,f[2][i]])
                flag[f[0][i],f[1][i],f[2][i]] = 5
            else:
                ch = np.sum(flag[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,f[2][i]-1:f[2][i]+2] == 2)
                print('Deep ch: ' + str(ch))
                if ch > 11:
                    chla[f[0][i],f[1][i],f[2][i]] = np.nanmean(chla[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,f[2][i]-1:f[2][i]+2])
                    flag[f[0][i],f[1][i],f[2][i]] = 6

f = np.where((flag[:,:,0] == 0) & ((lat[np.newaxis,:] > -40) & (lat[np.newaxis,:] < 40)))
for i in range(0,len(f[0])):
    if (f[0][i] == lon[0]) |  (f[0][i] == lon[-1]):
        print('Lon skip')
    else:
        chla[f[0][i],f[1][i],0] = np.nanmean(chla[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,0:2])
        flag[f[0][i],f[1][i],0] = 7

f = np.where((flag[:,:,-1] == 0) & ((lat[np.newaxis,:] > -40) & (lat[np.newaxis,:] < 40)))
for i in range(0,len(f[0])):
    if (f[0][i] == lon[0]) |  (f[0][i] == lon[-1]):
        print('Lon skip')
    else:
        chla[f[0][i],f[1][i],-1] = np.nanmean(chla[f[0][i]-1:f[0][i]+2,f[1][i]-1:f[1][i]+2,-2:])
        flag[f[0][i],f[1][i],-1] = 7

"""
Here we apply the wintertime chl-a relationships developed from the Bio-Argo profilers for the Arctic and Southern Ocean
seperately. This uses a relative relationship between the nearest in time OC-CCI observation (either backwards or forwards) and the Bio-Argo
profiler chl-a to construct the wintertime chl-a values.
"""
# Here we check for pixels that flag as 0 (not dealt with yet) and 4,5,6 (for some reason gone through previous averaging - should be fixed now...), and that aren't
# highly ice covered.
f = np.where(((flag == 0) | (flag == 4) | (flag == 5) | (flag == 6)) & ((lat[np.newaxis,:,np.newaxis] < -40)) & (ice <0.7))
print(f)

for i in range(len(f[0])):
    t = 0
    l = 1
    while t == 0:
        t2 = 0
        while t2 == 0:
            g = flag[f[0][i],f[1][i],f[2][i]-l]
            if g == 2:
                #print(l)
                if f[2][i]-l >0:
                    chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]-l]) * (1-np.abs(argo_b[l,1])))
                    flag[f[0][i],f[1][i],f[2][i]] = 8
                    flag_l[f[0][i],f[1][i],f[2][i]] = -l
                    l2 = l
                    t2 = 1

            l=l+1
            if l > 6:
                l2 = l
                t2 =1
        t2 = 0
        l = 1
        while t2 == 0:
            #l2 = l
            if f[2][i]+l <chla.shape[2]:
                g = flag[f[0][i],f[1][i],f[2][i]+l]
                if g == 2:
                    if l < l2:
                        chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]+l]) * (1-np.abs(argo_f[l,1])))
                        flag[f[0][i],f[1][i],f[2][i]] = 9
                        flag_l[f[0][i],f[1][i],f[2][i]] = l
                        t=1
                    t2 = 1
            l = l+1
            if l > 6:
                t2 = 1
        t = 1

"""
Artic version
"""
f = np.where(((flag == 0) | (flag == 4) | (flag == 5) | (flag == 6)) & ((lat[np.newaxis,:,np.newaxis] > 40)) & (ice <0.7))
print(f)

for i in range(len(f[0])):
    t = 0
    l = 1
    while t == 0:
        t2 = 0
        while t2 == 0:
            g = flag[f[0][i],f[1][i],f[2][i]-l]
            if g == 2:
                #print(l)
                if f[2][i]-l >0:
                    chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]-l]) * (1-np.abs(argo_b_arctic[l,1])))
                    flag[f[0][i],f[1][i],f[2][i]] = 8
                    flag_l[f[0][i],f[1][i],f[2][i]] = -l
                    l2 = l
                    t2 = 1
            l=l+1
            if l > 6:
                t2 =1
        t2 = 0
        l = 1
        while t2 == 0:
            #l2 = l
            if f[2][i]+l <chla.shape[2]:
                g = flag[f[0][i],f[1][i],f[2][i]+l]
                if g == 2:
                    if l < l2:
                        chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]+l]) * (1-np.abs(argo_f_arctic[l,1])))
                        flag[f[0][i],f[1][i],f[2][i]] = 9
                        flag_l[f[0][i],f[1][i],f[2][i]] = l
                        t=1
                    t2 = 1
            l = l+1
            if l > 6:
                t2 = 1
        t = 1

"""
Marginal Ice Zone - For this region we assume the chl-a will transition from the under ice value, to the "open" ocean values.
We use a spatial mean around the pixel, and if this still returns no value we set the value to the underice chl-a value. This
provides a smooth transition between underice and "open ocean" regions.
"""
chla2 = np.copy(chla)
f = np.where(((flag == 0)) & ((lat[np.newaxis,:,np.newaxis] < -40) | (lat[np.newaxis,:,np.newaxis] > 40)) & ((ice >=0.7) | (ice < 0.9)))
for i in range(len(f[0])):
    if f[0][i]<2:
        c = np.concatenate(((chla2[0 : f[0][i]+2 , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel(), (chla2[chla.shape[0]-f[0][i] : chla.shape[0] , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel()))
        c = np.nanmean(c)
    elif f[0][i] > chla.shape[0]-2:
        #print(2)
        #c = np.log10(3)
        # print(chla.shape[0]- f[0][i])
        # print((chla2[0 : chla.shape[0]- f[0][i] , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel().shape)
        # print((chla2[f[0][i]-2 : chla.shape[0] , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel().shape)
        c = np.concatenate(((chla2[0 : chla.shape[0]- f[0][i] , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel(), (chla2[f[0][i]-2 : chla.shape[0] , f[1][i]-2:f[1][i]+3 , f[2][i]]).ravel()))
        c = np.nanmean(c)
    elif f[1][i] <2:
        c = np.nanmean((chla2[f[0][i]-2:f[0][i]+3,0:f[1][i]+3,f[2][i]]).ravel())
    else:
        c = np.nanmean((chla2[f[0][i]-2:f[0][i]+3,f[1][i]-2:f[1][i]+3,f[2][i]]).ravel())
    #print(c)
    if np.isnan(c) == 1:
        chla[f[0][i],f[1][i],f[2][i]] = np.log10(0.05)
    else:
        chla[f[0][i],f[1][i],f[2][i]] = c
    flag[f[0][i],f[1][i],f[2][i]] = 10

"""
Final filling is a linear interpoaltion through time - We have been through many combinations of how to fill the data by this point.
So final fill for a pixel that hasn't been assigned a chl-a value is to look at the two nearest temporal points and linear interpolate
temporally between them. We do not care if they are an observation (flag==2) or a previous filled value by this stage...
"""
chla2 = np.copy(chla)
f = np.where((flag == 0))
for i in range(len(f[0])):
    t = 0
    b = 0
    while t == 0:
        b = b+1
        c = np.isnan(chla2[f[0][i],f[1][i],f[2][i]-b])

        if c == 0:
            t = 1
        if b > 5:
            t=1
            b=0
    t = 0
    a = 0
    while t == 0:
        a = a+1
        c = np.isnan(chla2[f[0][i],f[1][i],f[2][i]+a])
        if c == 0:
            t = 1
        if (a > 5) | (f[2][i]+a+1 == chla.shape[2]):
            t=1
            a=0
    if (a!=0) & (b!=0):
        ch = np.interp([0],[-b,a],[chla2[f[0][i],f[1][i],f[2][i]-b],chla2[f[0][i],f[1][i],f[2][i]+a]])
        print(ch)
        chla[f[0][i],f[1][i],f[2][i]] = ch
        flag[f[0][i],f[1][i],f[2][i]] = 11

save_chla(outfile,chla,flag,flag_l,lat,lon,time,outloc,time_r)
