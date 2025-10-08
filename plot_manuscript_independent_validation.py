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
import os
font = {'weight' : 'normal',
        'size'   :14}
matplotlib.rc('font', **font)
sz = 8



def plotting_independent_wintertime(file_occci,hplc_file,fluoro_file,output_loc,res):
    unit = 'log$_{10}$(mgm$^{-3}$)'
    c2 = Dataset(file_occci,'r')
    lat = np.array(c2['latitude'])
    chla = np.array(c2['chl_filled'])
    flags = np.array(c2['chl_flag'])
    c2.close()

    c2 = Dataset(hplc_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()
    c = np.array([-3,3.0])

    fig = plt.figure(figsize=(14,14))
    gs = GridSpec(2,2, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
    ax = fig.add_subplot(gs[1,0])

    f = np.where(((flags == 5) | (flags == 6)) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',label='Southern Hemisphere',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Southern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('Wintertime in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Wintertime gapfilled Chlorophyll-a ((mg m$^{-3}$)')
    ax.text(0.03,0.95,f'(c)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    ax = fig.add_subplot(gs[0,0])
    f = np.where(((flags == 7) | (flags == 8)) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',label='Northern Hemisphere',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Northern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('Wintertime in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (mg m$^{-3}$)')
    ax.text(0.03,0.95,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)


    c2 = Dataset(fluoro_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()

    ax = fig.add_subplot(gs[1,1])

    f = np.where(((flags == 5) | (flags == 6)) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',label='Southern Hemisphere',s=sz)

    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Southern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('Wintertime in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (mg m$^{-3}$)')

    ax.text(0.03,0.95,f'(d)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    ax = fig.add_subplot(gs[0,1])
    f = np.where(((flags == 7) | (flags == 8)) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',label='Northern Hemisphere',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Northern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.set_xlabel('Wintertime in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (mg m$^{-3}$)')
    ax.grid()
    ax.text(0.03,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    fig.savefig(os.path.join(output_loc,'plots',f'manuscript_wintertime_independent_{res}.png'),dpi=300)

def plotting_independent_cloud(file_occci,hplc_file,fluoro_file,output_loc,res):
    unit = 'log$_{10}$(mgm$^{-3}$)'
    c2 = Dataset(file_occci,'r')
    lat = np.array(c2['latitude'])
    chla = np.array(c2['chl_filled'])
    flags = np.array(c2['chl_flag'])
    c2.close()

    c2 = Dataset(hplc_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()
    c = np.array([-3,3.0])

    fig = plt.figure(figsize=(14,7))
    gs = GridSpec(1,2, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
    ax = fig.add_subplot(gs[0,0])

    f = np.where((flags == 3) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Cloud filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Cloud gapfilled Chlorophyll-a (mg m$^{-3}$)')
    ax.text(0.03,0.95,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)



    c2 = Dataset(fluoro_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()

    ax = fig.add_subplot(gs[0,1])

    f = np.where((flags == 3) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)

    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Cloud filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Cloud gapfilled Chlorophyll-a (mg m$^{-3}$)')

    ax.text(0.03,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    fig.savefig(os.path.join(output_loc,'plots',f'manuscript_cloud_independent_{res}.png'),dpi=300)

def plotting_independent_final_krig(file_occci,hplc_file,fluoro_file,output_loc,res):
    unit = 'log$_{10}$(mgm$^{-3}$)'
    c2 = Dataset(file_occci,'r')
    lat = np.array(c2['latitude'])
    chla = np.array(c2['chl_filled'])
    flags = np.array(c2['chl_flag'])
    c2.close()

    c2 = Dataset(hplc_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()
    c = np.array([-3,3.0])

    fig = plt.figure(figsize=(14,7))
    gs = GridSpec(1,2, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
    ax = fig.add_subplot(gs[0,0])

    f = np.where((flags == 9) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Final Kriging filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Final Kriging gapfilled Chlorophyll-a (mg m$^{-3}$)')
    ax.text(0.03,0.95,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)



    c2 = Dataset(fluoro_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()

    ax = fig.add_subplot(gs[0,1])

    f = np.where((flags == 9) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)

    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Final Kriging filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Final Kriging gapfilled Chlorophyll-a (mg m$^{-3}$)')

    ax.text(0.03,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    fig.savefig(os.path.join(output_loc,'plots',f'manuscript_final_kriging_independent_{res}.png'),dpi=300)

def plotting_independent_underice(file_occci,hplc_file,fluoro_file,output_loc,res):
    unit = 'log$_{10}$(mgm$^{-3}$)'
    c2 = Dataset(file_occci,'r')
    lat = np.array(c2['latitude'])
    chla = np.array(c2['chl_filled'])
    flags = np.array(c2['chl_flag'])
    c2.close()

    c2 = Dataset(hplc_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()
    c = np.array([-3,3.0])

    fig = plt.figure(figsize=(14,7))
    gs = GridSpec(1,2, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
    ax = fig.add_subplot(gs[0,0])

    f = np.where((flags == 4) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Under Ice filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Under Ice gapfilled Chlorophyll-a (mg m$^{-3}$)')
    ax.text(0.03,0.95,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)



    c2 = Dataset(fluoro_file,'r')
    chl_in = np.array(c2['chl'])
    c2.close()

    ax = fig.add_subplot(gs[0,1])

    f = np.where((flags == 4) & (np.isnan(chl_in) == 0))
    ax.scatter(10**chl_in[f],10**chla[f],zorder=4,color='b',s=sz)

    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(10**c,10**(c*stats['slope']+stats['intercept']),'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.45,0.3,f'Under Ice filled\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)
    ax.plot(10**c,10**c,'k-')
    ax.set_xlim(10**c); ax.set_ylim(10**c);
    ax.set_yscale('log'); ax.set_xscale('log')
    # ax.legend(loc=3)
    ax.grid()
    ax.set_xlabel('in situ Chlorophyll-a (mg m$^{-3}$)')
    ax.set_ylabel('Under Ice gapfilled Chlorophyll-a (mg m$^{-3}$)')

    ax.text(0.03,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)

    fig.savefig(os.path.join(output_loc,'plots',f'manuscript_under_ice_independent_{res}.png'),dpi=300)
