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
import weight_stats as ws

occci_loc = 'E:/Data/OC-CCI/v6.0/daily'
files = ['argo_chla_southernocean','argo_chla_arctic']
generate = False
plot = True

if generate:
    t = 0
    for file in files:
        data = np.genfromtxt('csv/'+file+'.csv', delimiter=',')
        matched_chl = np.zeros((len(data),1))
        matched_chl[:] = np.nan
        for i in range(len(data)):
            yr = int(data[i,0])
            mon = int(data[i,1])
            day = int(data[i,2])
            dat = datetime.datetime(yr,mon,day)
            lat = data[i,3]
            lon = data[i,4]

            oc_file = os.path.join(occci_loc,str(yr),'ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-'+dat.strftime('%Y%m%d')+'-fv6.0.nc')
            if glob.glob(oc_file):
                if t == 0:
                    c=Dataset(oc_file,'r')
                    oc_lat = np.array(c['lat'])
                    oc_lon = np.array(c['lon'])
                    c.close()
                    t=1

                lat_dif = np.abs(lat - oc_lat)
                lon_dif = np.abs(lon - oc_lon)

                f = np.where(lat_dif == np.min(lat_dif))[0][0]
                g = np.where(lon_dif == np.min(lon_dif))[0][0]
                print(f)
                print(g)
                c = Dataset(oc_file,'r')
                oc_chl = np.squeeze(np.array(c.variables['chlor_a'][0,f-1:f+2,g-1:g+2]))
                oc_chl[oc_chl>10000] = np.nan
                print(oc_chl)
                oc_chl = 10**np.nanmean(np.log10(np.ravel(oc_chl)))
                print(oc_chl)

                c.close()
                matched_chl[i] = oc_chl
        data = np.concatenate((data,matched_chl),axis=1)
        np.savetxt('csv/'+file+'_daily.csv',data,delimiter=',')

if plot:
    font = {'weight' : 'normal',
            'size'   :16}
    matplotlib.rc('font', **font)
    worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))
    unit = 'log$_{10}$(mgm$^{-3}$)'
    col = 2
    row = int(len(files))
    fig = plt.figure(figsize=(21,7*row))
    gs = GridSpec(row,col, figure=fig, wspace=0.2,hspace=0.25,bottom=0.1,top=0.9,left=0.07,right=0.95)
    axs = [[fig.add_subplot(gs[i, j]) for j in range(col)] for i in range(row)]
    flatList = [element for innerList in axs for element in innerList]
    axs = flatList
    i=0
    for file in files:
        data = np.genfromtxt('csv/'+file+'_daily.csv', delimiter=',')
        data[np.log10(data[:,5])<-2.2,5] = np.nan
        c = np.array([-3,1])

        axs[i].scatter(np.log10(data[:,5]),np.log10(data[:,7]),s=6,zorder=4,color='b')
        axs[i].plot(c,c,'k-')
        axs[i].set_title(file)
        stats = ws.unweighted_stats(np.log10(data[:,5]),np.log10(data[:,7]),file)
        print(stats)
        axs[i].plot(c,c*stats['slope']+stats['intercept'],'b--')
        #axs[i].scatter(np.log10(data[:,5])+stats['med_rel_bias'],np.log10(data[:,7]),s=6,color='r',label = 'Bias corrected')
        axs[i].grid()
        rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
        ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
        axs[i].text(0.55,0.35,f'Unweighted Stats\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=axs[i].transAxes,va='top')
        axs[i].set_xlim(c); axs[i].set_ylim(c);
        axs[i].set_xlabel('Argo Chl-a (log$_{10}$(mgm$^{-3}$))')
        axs[i].set_ylabel('OC-CCI Chl-a (log$_{10}$(mgm$^{-3}$))')

        i= i+1
        worldmap.plot(color='lightgrey',ax=axs[i])
        vmin = -0.5; vmax = 0.5
        a = axs[i].scatter(data[:,4],data[:,3],c=np.log10(data[:,5])-np.log10(data[:,7]),vmin = vmin,vmax=vmax,cmap = cmocean.cm.balance)
        cbar = plt.colorbar(a)
        cbar.set_label('Argo - OC-CCI ('+unit+')')
        i=i+1
    plt.suptitle('Daily OC-CCI matchups')
    fig.savefig('plots/argo_daily_validation.png',dpi=300)
