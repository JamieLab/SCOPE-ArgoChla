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
import weight_stats as ws
import data_utils as du

let = ['a','b','c','d','e','f','g']

def generate_occci(occci_file,res,lon,lat,start_yr,end_yr):
    import CCI_OC_SPATIAL_AV as OC
    OC.oc_cci_average('F:/Data/OC-CCI/monthly/chlor_a','F:/Data/OC-CCI/monthly/chlor_a/'+str(res)+'DEG',start_yr = start_yr,end_yr = end_yr,log=lon,lag=lat)
    import Data_Loading.gebco_resample as ge
    ge.gebco_resample('F:/Data/Bathymetry/GEBCO_2023.nc',lon,lat,save_loc = 'F:/Data/Bathymetry/'+str(res)+'DEG_GEBCO_2023.nc')

    # from Data_Loading.OSISAF_download import OSISAF_spatial_average
    # from Data_Loading.OSISAF_download import OSISAF_merge_hemisphere
    # OSISAF_spatial_average(data='F:/Data/OSISAF/monthly',out_loc='F:/Data/OSISAF/monthly/'+str(res)+'DEG',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'NH')
    # OSISAF_spatial_average(data='F:/Data/OSISAF/monthly',out_loc='F:/Data/OSISAF/monthly/'+str(res)+'DEG',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'SH')
    # OSISAF_merge_hemisphere('F:/Data/OSISAF/monthly/'+str(res)+'DEG', 'F:/Data/Bathymetry/'+str(res)+'DEG_GEBCO_2023.nc',start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat)

    import construct_input_netcdf as cinp
    vars = [['OC-CCI','chlor_a','F:/Data/OC-CCI/monthly/chlor_a/'+str(res)+'DEG/%Y/%Y_%m_*.nc',0],
        ['OC-CCI','chlor_a_log10_rmsd','F:/Data/OC-CCI/monthly/chlor_a/'+str(res)+'DEG/%Y/%Y_%m_*.nc',0],
            ['OSISAF','ice_conc','F:/Data/OSISAF/monthly/'+str(res)+'DEG/%Y/%Y%m*COM.nc',0]]

    cinp.driver(occci_file,vars,start_yr = start_yr,end_yr = end_yr,lon = lon,lat = lat,fill_clim=False)

# res = 0.25

# start_yr = 1997
# end_yr = 2023
# files = ['argo_chla_southernocean','argo_chla_arctic']

def chl_argo_relationship(res,start_yr,end_yr,files,occci_file,plot=False,lags = 9):
    if plot:
        widths = 0.75
        font = {'weight' : 'normal',
                'size'   :14}
        matplotlib.rc('font', **font)
        fig = plt.figure(figsize=(28,7*len(files)))
        gs = GridSpec(len(files),2, figure=fig, wspace=0.2,hspace=0.2,bottom=0.1,top=0.95,left=0.05,right=0.95)
        axs = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(len(files))]
        flatList = [element for innerList in axs for element in innerList]
        axs = flatList
        p=0
    lon,lat = du.reg_grid(lat=res,lon=res)
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
                if l == lags+1:
                    t=1
        arg = np.array(arg); oc = np.array(oc); dif = np.array(dif)
        # fig = plt.figure(figsize=(18,7))

        meds = np.zeros((lags+1,4))
        for i in range(0,lags+1):
            g = np.where((dif == i))
            vals = (np.log10(arg[g])-np.log10(oc[g]))
            if i == 0:
                med = np.nanmedian(vals)

            vals = (10**((np.log10(arg[g])-med)) -oc[g]) / oc[g]

            if plot:
                axs[p].boxplot(vals*100,positions=[i],widths=widths,labels = [f'{i} \n N = {len(g[0])} \n Med = {str(int(round(np.nanmedian(vals*100),0)))} '])
            meds[i,0] = i
            meds[i,1] = np.nanmedian(vals)
            meds[i,2] = np.nanmedian(np.abs(vals - meds[i,1])) * 1.4826 # Converting median absolute deviation to a robust standard deviation equivalent
            meds[i,3] = len(g[0])
        if plot:
            axs[p].set_ylim([-100,100])
            axs[p].set_ylabel('Percentage difference to OC-CCI chl-a')
            axs[p].text(0.93,0.95,f'('+let[p]+')',transform=axs[p].transAxes,va='top',fontweight='bold',fontsize = 25)
            p=p+1
        # fig.savefig('plots/backwards_'+file+'_'+str(res)+'deg.png',dpi=300)
        np.savetxt('relationships/backwards_'+file+'_'+str(res)+'deg.csv',meds,delimiter=',',header='month_diff,% difference,% difference uncertainty,number_samples')

        # fig, ax = plt.subplots(1,1,figsize=(7,7))
        # unit = 'log$_{10}$(mgm$^{-3}$)'
        # c = np.array([-3,2])
        # g = np.where(dif == 0)
        # print(med)
        # ax.scatter(np.log10(arg[g]),np.log10(oc[g]),label='No Bias Correction',s=6,color='b',zorder=3)
        # stats=ws.unweighted_stats(np.log10(arg[g]),np.log10(oc[g]),file)
        # #plt.scatter(np.log10(arg[g])-med,np.log10(oc[g]),label='Bias Corrected',s=6,color='r')
        # ax.plot(c,c,'k-')
        # ax.set_xlim(c); ax.set_ylim(c)
        # ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
        # ax.set_xlabel('Argo Chl-a (log$_{10}$(mgm$^{-3}$))')
        # ax.set_ylabel('OCCCI Chl-a (log$_{10}$(mgm$^{-3}$))')
        # rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
        # ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
        # ax.text(0.5,0.3,f'Unweighted Stats\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top')
        # ax.grid()
        # ax.set_title(file)

        #plt.legend()
        # fig.savefig('plots/compare_'+file+'_'+str(res)+'deg.png',dpi=300)

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
                    if l == lags+1:
                        t=1
        arg = np.array(arg); oc = np.array(oc); dif = np.array(dif)
        # if plot:
        # fig = plt.figure(figsize=(18,7))
        meds = np.zeros((lags+1,4))
        for i in range(0,lags+1):
            g = np.where((dif == i))
            # vals = (arg[g]-oc[g])/((oc[g]))
            # if i == 0:
            #     med = np.nanmedian(vals)
            vals = (10**((np.log10(arg[g])-med)) -oc[g]) / oc[g]
            if plot:
                axs[p].boxplot(vals*100,positions=[i],widths=widths,labels = [f'{i} \n N = {len(g[0])} \n Med = {str(int(round(np.nanmedian(vals*100),0)))} '])
            meds[i,0] = i
            meds[i,1] = np.nanmedian(vals)
            meds[i,2] = np.nanmedian(np.abs(vals - meds[i,1])) * 1.4826 # Converting median absolute deviation to a robust standard deviation equivalent
            meds[i,3] = len(g[0])
            print(meds[i,2])
        if plot:
            axs[p].set_ylim([-100,100])
            axs[p].text(0.93,0.97,f'('+let[p]+')',transform=axs[p].transAxes,va='top',fontweight='bold',fontsize = 25)
            #axs[p].set_ylabel('Percentage difference to OC-CCI chl-a')
            p=p+1


        np.savetxt('relationships/forwards_'+file+'_'+str(res)+'deg.csv',meds,delimiter=',',header='month_diff,% difference,% difference uncertainty,number_samples')
        #plt.show()
    if plot:
        fig.savefig('plots/relationships_'+str(res)+'deg.png',dpi=300)
