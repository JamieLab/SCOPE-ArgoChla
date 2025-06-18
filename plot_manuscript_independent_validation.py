#!/usr/bin/env python3
"""
Created by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 10/2024
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import patches
import cmocean
import datetime
import geopandas as gpd
import weight_stats as ws
import matplotlib.transforms
font = {'weight' : 'normal',
        'size'   :14}
matplotlib.rc('font', **font)
sz = 8

unit = 'log$_{10}$(mg m$^{-3}$)'

res = 0.25
file_occci = 'E:/SCOPE/Argo/netcdf/oc-cci_chlor_a_'+str(res)+'deg.nc'
file_chl = f'E:/SCOPE/Argo/netcdf/insitudb_chla_V3_{str(res)}deg.nc'
output = f'plots/insitu_verification_{str(res)}_deg_manuscript.png'
c = Dataset(file_occci,'r')
lat = np.array(c['latitude'])
chla = np.array(c['chl_filled'])
flags = np.array(c['chl_flag'])
c.close()

c = Dataset(file_chl,'r')
chl_in = np.array(c['chl'])
c.close()
c = np.array([-3,1.5])
# f = np.where((lat < bounds[0]) & (lat > bounds[1]))[0]
# print(f)
fig = plt.figure(figsize=(14,14))
gs = GridSpec(2,2, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
ax = fig.add_subplot(gs[1,0])
# chla_t = chla[:,f,:]
# flags_t = flags[:,f,:]
# chl_in = chl_in[:,f,:]
f = np.where(((flags == 5) | (flags == 6)) & (np.isnan(chl_in) == 0))
ax.scatter(chl_in[f],chla[f],zorder=4,color='b',label='Southern Hemisphere',s=sz)
stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
ax.text(0.45,0.3,f'Southern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
ax.plot(c,c,'k-')
ax.set_xlim(c); ax.set_ylim(c);
# ax.legend(loc=3)
ax.set_xlabel('Wintertime in situ Chlorophyll-a (log$_{10}$(mg m$^{-3}$))')
ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (log$_{10}$(mg m$^{-3}$))')
ax.text(0.03,0.95,f'(c)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

ax = fig.add_subplot(gs[0,0])
f = np.where(((flags == 7) | (flags == 8)) & (np.isnan(chl_in) == 0))
ax.scatter(chl_in[f],chla[f],zorder=4,color='b',label='Northern Hemisphere',s=sz)
stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
ax.text(0.45,0.3,f'Northern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)

ax.plot(c,c,'k-')
ax.set_xlim(c); ax.set_ylim(c);
# ax.legend(loc=3)
ax.set_xlabel('Wintertime in situ Chlorophyll-a (log$_{10}$(mg m$^{-3}$))')
ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (log$_{10}$(mg m$^{-3}$))')
ax.text(0.03,0.95,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

file_chl = f'E:/SCOPE/Argo/netcdf/insitudb_chla_V3_{str(res)}deg_fluoro.nc'
c = Dataset(file_occci,'r')
lat = np.array(c['latitude'])
chla = np.array(c['chl_filled'])
flags = np.array(c['chl_flag'])
c.close()

c = Dataset(file_chl,'r')
chl_in = np.array(c['chl'])
c.close()
c = np.array([-3,1.5])
ax = fig.add_subplot(gs[1,1])
# chla_t = chla[:,f,:]
# flags_t = flags[:,f,:]
# chl_in = chl_in[:,f,:]
f = np.where(((flags == 5) | (flags == 6)) & (np.isnan(chl_in) == 0))
ax.scatter(chl_in[f],chla[f],zorder=4,color='b',label='Southern Hemisphere',s=sz)
stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
ax.text(0.45,0.3,f'Southern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
ax.plot(c,c,'k-')
ax.set_xlim(c); ax.set_ylim(c);
# ax.legend(loc=3)
ax.set_xlabel('Wintertime in situ Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
ax.text(0.03,0.95,f'(d)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
ax = fig.add_subplot(gs[0,1])
f = np.where(((flags == 7) | (flags == 8)) & (np.isnan(chl_in) == 0))
ax.scatter(chl_in[f],chla[f],zorder=4,color='b',label='Northern Hemisphere',s=sz)
stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
ax.text(0.45,0.3,f'Northern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)

ax.plot(c,c,'k-')
ax.set_xlim(c); ax.set_ylim(c);
# ax.legend(loc=3)
ax.set_xlabel('Wintertime in situ Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
ax.text(0.03,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

fig.savefig(output,dpi=300)
