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


def generate_occci(occci_file,res,lon,lat,start_yr,end_yr):
    import CCI_OC_SPATIAL_AV as OC
    OC.oc_cci_average('D:/Data/OC-CCI/monthly/chlor_a','D:/Data/OC-CCI/monthly/chlor_a/'+str(res)+'DEG',start_yr = start_yr,end_yr = end_yr,log=lon,lag=lat)
    import Data_Loading.gebco_resample as ge
    ge.gebco_resample('D:/Data/Bathymetry/GEBCO_2023.nc',lon,lat,save_loc = 'D:/Data/Bathymetry/'+str(res)+'DEG_GEBCO_2023.nc')

    from Data_Loading.OSISAF_download import OSISAF_spatial_average
    from Data_Loading.OSISAF_download import OSISAF_merge_hemisphere
    # OSISAF_spatial_average(data='D:/Data/OSISAF/monthly',out_loc='D:/Data/OSISAF/monthly/'+str(res)+'DEG',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'NH')
    # OSISAF_spatial_average(data='D:/Data/OSISAF/monthly',out_loc='D:/Data/OSISAF/monthly/'+str(res)+'DEG',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'SH')
    OSISAF_merge_hemisphere('D:/Data/OSISAF/monthly/'+str(res)+'DEG', 'D:/Data/Bathymetry/'+str(res)+'DEG_GEBCO_2023.nc',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat)

    import construct_input_netcdf as cinp
    vars = [['OC-CCI','chlor_a','D:/Data/OC-CCI/monthly/chlor_a/'+str(res)+'DEG/%Y/%Y_%m_*.nc',0],
            ['OSISAF','ice_conc','D:/Data/OSISAF/monthly/'+str(res)+'DEG/%Y/%Y%m*COM.nc',0]]

    cinp.driver(occci_file,vars,start_yr = start_yr,end_yr = end_yr,lon = lon,lat = lat,fill_clim=False)

res = 1
lon,lat = du.reg_grid(lat=res,lon=res)
start_yr = 1997
end_yr = 2022
files = ['argo_chla_southernocean','argo_chla_arctic']


occci_file = 'netcdf/oc-cci_chlor_a_'+str(res)+'deg.nc'
generate_occci(occci_file,res,lon,lat,start_yr,end_yr)
for file in files:
    app = file + '_' + str(res)

    argo_file = 'netcdf/'+file+'_'+str(res)+'deg.nc'


    c = Dataset(argo_file,'r')
    argo = 10**np.array(c['chl'])
    argo[argo<0.005] = np.nan
    c.close()

    c = Dataset(occci_file,'r')
    occci = 10**np.array(c['OC-CCI_chlor_a'])
    time = np.array(c['time'])
    c.close()

    # f = np.where(lat >= 40)[0]
    # argo = argo[:,f,:]
    # occci = occci[:,f,:]
    # lat = lat[f]
    print(argo.shape)
    print(occci.shape)

    f = np.where(np.isnan(argo) == 0)
    print(f)
    print(len(f[0]))

    arg = []
    oc = []
    dif = []
    for i in range(0,len(f[0])):
        l = 0
        t = 0
        while t == 0:
            oct = occci[f[0][i],f[1][i],f[2][i]-l]
            if np.isnan(oct) == 0:
                arg.append(argo[f[0][i],f[1][i],f[2][i]])
                oc.append(oct)
                dif.append(l)
                t = 1
            else:
                l = l+1
            if l == 12:
                t=1
    arg = np.array(arg); oc = np.array(oc); dif = np.array(dif)
    fig = plt.figure(figsize=(15,7))

    meds = np.zeros((10,3))
    for i in range(0,10):
        g = np.where((dif == i))
        vals = (np.log10(arg[g])-np.log10(oc[g]))
        if i == 0:
            med = np.nanmedian(vals)
        vals = (10**((np.log10(arg[g])-med)) -oc[g]) / oc[g]
        plt.boxplot(vals,positions=[i],labels = [f'{i} \n N = {len(g[0])} \n Med = {str(round(np.nanmedian(vals),2))} '])
        meds[i,0] = i
        meds[i,1] = np.nanmedian(vals)
        meds[i,2] = len(g[0])
    plt.ylim([-1,1])
    plt.ylabel('Percentage difference to OC-CCI chl-a')
    fig.savefig('plots/backwards_'+file+'_'+str(res)+'deg.png',dpi=300)
    np.savetxt('relationships/backwards_'+file+'_'+str(res)+'deg.csv',meds,delimiter=',',header='month_diff,% difference,number_samples')
    fig = plt.figure()
    c = [-3,1]
    g = np.where(dif == 0)
    print(med)
    plt.scatter(np.log10(arg[g]),np.log10(oc[g]),label='No Bias Correction')
    plt.scatter(np.log10(arg[g])-med,np.log10(oc[g]),label='Bias Corrected')
    plt.plot(c,c,'k-')
    plt.xlim(c)
    plt.ylim(c)
    plt.xlabel('Argo Chl-a (log$_{10}$(mgm$^{-3}$))')
    plt.ylabel('OCCCI Chl-a (log$_{10}$(mgm$^{-3}$))')
    plt.legend()
    fig.savefig('plots/compare_'+file+'_'+str(res)+'deg.png',dpi=300)

    arg = []
    oc = []
    dif = []
    for i in range(0,len(f[0])):
        l = 0
        t = 0
        while t == 0:
            if  f[2][i]+ l >= argo.shape[2]:
                t = 1
            else:
                oct = occci[f[0][i],f[1][i],f[2][i]+l]
                if np.isnan(oct) == 0:
                    arg.append(argo[f[0][i],f[1][i],f[2][i]])
                    oc.append(oct)
                    dif.append(l)
                    t = 1
                else:
                    l = l+1
                if l == 12:
                    t=1
    arg = np.array(arg); oc = np.array(oc); dif = np.array(dif)
    fig = plt.figure(figsize=(15,7))
    meds = np.zeros((10,3))
    for i in range(0,10):
        g = np.where((dif == i))
        # vals = (arg[g]-oc[g])/((oc[g]))
        # if i == 0:
        #     med = np.nanmedian(vals)
        vals = (10**((np.log10(arg[g])-med)) -oc[g]) / oc[g]
        plt.boxplot(vals,positions=[i],labels = [f'{i} \n N = {len(g[0])} \n Med = {str(round(np.nanmedian(vals),2))} '])
        meds[i,0] = i
        meds[i,1] = np.nanmedian(vals)
        meds[i,2] = len(g[0])
    plt.ylim([-1,1])
    plt.ylabel('Percentage difference to OC-CCI chl-a')
    np.savetxt('relationships/forwards_'+file+'_'+str(res)+'deg.csv',meds,delimiter=',',header='month_diff,% difference,number_samples')
    fig.savefig('plots/forwards_'+file+'_'+str(res)+'deg.png',dpi=300)
    #plt.show()
