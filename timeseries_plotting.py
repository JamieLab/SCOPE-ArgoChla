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
import matplotlib.colors as colors
import data_utils as du


import matplotlib.transforms
font = {'weight' : 'normal',
        'size'   : 22}
matplotlib.rc('font', **font)

def extract_bar(lon,lat,data,latg,long,ax,maxis,time,fig,yticks=True,label = 'Year',unc=False,unc_data=False,start_yr=1997,end_yr=2024,bottom=False):
    maxis.scatter(lon,lat,color='r',s=16*7)
    f2 = np.abs(lat - latg)
    g2 = np.abs(lon - long)
    f = np.where(np.min(f2) == f2)[0]
    g = np.where(np.min(g2) == g2)[0]

    ax.plot(time,10**data[g[0],f[0],:],'k')
    if unc:
        ax.fill_between(time,10**(data[g[0],f[0],:] - unc_data[g[0],f[0],:]),10**(data[g[0],f[0],:] + unc_data[g[0],f[0],:]),alpha = 0.6,color='k',linewidth=0)
    ax.set_ylim(10**np.array([-2,0.5]))
    ax.set_xlabel(label)
    ax.set_xlim([start_yr,end_yr+1])
    if yticks:
        ax.set_ylabel('Chlorophyll-a (mg m$^{-3}$)')
    if not yticks:
        ax.tick_params(labelleft=False)
    ax.grid()
    ax.set_yscale('log')
    lim = ax.get_ylim()
    if bottom:
        lim = lim[0]
        pad = -3
    else:
        lim = lim[1]
        pad=0

    arrow = patches.ConnectionPatch(
    [time[int(data.shape[2]/2)]+pad,lim],
    [lon,lat],
    coordsA=ax.transData,
    coordsB=maxis.transData,
    # Default shrink parameter is 0 so can be omitted
    color="black",
    arrowstyle="-|>",  # "normal" arrow
    mutation_scale=30,  # controls arrow head size
    linewidth=3,
    )
    fig.patches.append(arrow)
    return f[0],g[0]

def timeseries_plotting(file,bathy_file,loc,start_yr,end_yr,ref_time = datetime.datetime(1970,1,15),output='plots/global_map.png'):
    font = {'weight' : 'normal',
            'size'   : 22}
    matplotlib.rc('font', **font)
    c = Dataset(file,'r')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    chla = np.array(c['chl_filled'])
    unc = np.array(c['chl_filled_unc'])
    chla2 =  np.array(c['OC-CCI_chlor_a'])
    time = np.array(c['time'])
    c.close()

    c = Dataset(bathy_file,'r')
    ocean = np.array(c['ocean_proportion'])
    c.close()

    ocean = np.repeat(ocean[:, :, np.newaxis], chla.shape[2], axis=2)
    chla[ocean == 0] = np.nan


    worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))
    time_a = []
    for i in range(len(time)):
        time_a.append(ref_time + datetime.timedelta(days=int(time[i])))
    for i in range(len(time)):
        time_a[i] = time_a[i].year + ((time_a[i].month-1)/12)

    fig = plt.figure(figsize=(25,25))
    gs = GridSpec(1,1, figure=fig, wspace=0.2,hspace=0.2,bottom=0.33,top=0.67,left=0.1,right=1)
    ax1 = fig.add_subplot(gs[0,0]);
    me = 10**np.transpose(np.nanmean(chla,axis=2))
    pc = ax1.pcolor(lon,lat,me,cmap = cmocean.cm.algae,norm=colors.LogNorm(vmin=10**-2, vmax=10**1))
    ax1.text(0.92,1.05,f'(c)',transform=ax1.transAxes,va='top',fontweight='bold',fontsize = 35)
    #
    # ax = fig.add_subplot([0.76,0.05,0.15,0.25])
    # f,g = extract_bar(100,-55,chla,lat,lon,ax,ax1,time_a,fig,yticks=False,unc='True',unc_data=unc,start_yr=start_yr,end_yr=end_yr)
    # ax.text(0.85,0.95,f'(f)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    # ax.plot(time_a,10**chla2[g,f,:],'b')

    ax = fig.add_subplot([0.53,0.73,0.44,0.25])
    f,g = extract_bar(-25,60,chla,lat,lon,ax,ax1,time_a,fig,yticks=False,unc='True',unc_data=unc,start_yr=start_yr,end_yr=end_yr,bottom=True)
    ax.text(0.92,1.05,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**chla2[g,f,:],'b')

    ax = fig.add_subplot([0.06,0.73,0.44,0.25])
    f,g = extract_bar(-150,55,chla,lat,lon,ax,ax1,time_a,fig,yticks=True,unc='True',unc_data=unc,start_yr=start_yr,end_yr=end_yr,bottom=True)
    ax.text(0.92,1.05,f'(a)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**chla2[g,f,:],'b')
    #
    ax = fig.add_subplot([0.53,0.05,0.44,0.25])
    f,g = extract_bar(45,-60,chla,lat,lon,ax,ax1,time_a,fig,yticks=False,unc='True',unc_data=unc,start_yr=start_yr,end_yr=end_yr)
    ax.text(0.92,1.05,f'(e)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**chla2[g,f,:],'b')
    #
    ax = fig.add_subplot([0.06,0.05,0.44,0.25])
    f,g =extract_bar(-50,-55,chla,lat,lon,ax,ax1,time_a,fig,unc='True',unc_data=unc,start_yr=start_yr,end_yr=end_yr)
    ax.text(0.92,1.05,f'(d)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**chla2[g,f,:],'b')

    cbar = plt.colorbar(pc,ax=ax1)
    cbar.set_label('Chlorophyll-a (mg m$^{-3}$)');
    worldmap.plot(color="lightgrey", ax=ax1)
    # pc.set_clim([0,0.2])
    fig.savefig(output,dpi=300)

def climatology_plotting(file,bathy_file,loc,ref_time = datetime.datetime(1970,1,15),output='plots/global_climatology_map.png'):

    c = Dataset(file,'r')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    chla = np.array(c['chl_filled'])
    chla2 = np.array(c['OC-CCI_chlor_a'])
    time = np.array(c['time'])
    c.close()

    c = Dataset(bathy_file,'r')
    ocean = np.array(c['ocean_proportion'])
    c.close()

    ocean = np.repeat(ocean[:, :, np.newaxis], chla.shape[2], axis=2)
    chla[ocean == 0] = np.nan


    worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))
    time_a = []
    for i in range(len(time)):
        time_a.append(ref_time + datetime.timedelta(days=int(time[i])))
    for i in range(len(time)):
        time_a[i] = time_a[i].month
    time_a = np.array(time_a)

    clim = np.zeros((chla.shape[0],chla.shape[1],12))
    clim_std = np.copy(clim)
    clim2 = np.copy(clim)
    clim2_std = np.copy(clim)
    for i in range(1,13):
        f = np.where(time_a == i)
        clim[:,:,i-1] = np.nanmean(chla[:,:,f[0]],axis=2)
        clim_std[:,:,i-1] = np.nanstd(chla[:,:,f[0]],axis=2)

        clim2[:,:,i-1] = np.nanmean(chla2[:,:,f[0]],axis=2)
        clim2_std[:,:,i-1] = np.nanstd(chla2[:,:,f[0]],axis=2)


    time_a = np.array(range(1,13))
    fig = plt.figure(figsize=(25,16))
    gs = GridSpec(1,1, figure=fig, wspace=0.2,hspace=0.2,bottom=0.32,top=0.95,left=0.1,right=1)
    ax1 = fig.add_subplot(gs[0,0]);
    me = 10**np.transpose(np.nanmean(chla,axis=2))
    pc = ax1.pcolor(lon,lat,me,cmap = cmocean.cm.algae,norm=colors.LogNorm(vmin=10**-2, vmax=10**1))
    ax1.text(0.03,1.05,f'(a)',transform=ax1.transAxes,va='top',fontweight='bold',fontsize = 35)

    ax = fig.add_subplot([0.76,0.05,0.15,0.25])
    f,g = extract_bar(100,-55,clim,lat,lon,ax,ax1,time_a,fig,yticks=False,label = 'Month',unc=True,unc_data = clim_std)
    ax.text(0.85,0.95,f'(f)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    set_month(ax)
    width = 3
    ax.plot(time_a,10**clim2[g,f,:],'b-',linewidth=width)
    ax.plot(time_a,10**(np.ones((len(time_a)))*np.log10(0.3)),'b--')

    ax = fig.add_subplot([0.59,0.05,0.15,0.25])
    f,g = extract_bar(45,-60,clim,lat,lon,ax,ax1,time_a,fig,yticks=False,label = 'Month',unc=True,unc_data = clim_std)
    ax.text(0.85,0.95,f'(e)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    set_month(ax)
    ax.plot(time_a,10**clim2[g,f,:],'b-',linewidth=width)
    ax.plot(time_a,10**(np.ones((len(time_a)))*np.log10(0.3)),'b--')

    ax = fig.add_subplot([0.42,0.05,0.15,0.25])
    f,g = extract_bar(-25,60,clim,lat,lon,ax,ax1,time_a,fig,yticks=False,label = 'Month',unc=True,unc_data = clim_std)
    ax.text(0.85,0.95,f'(d)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**clim2[g,f,:],'b-',linewidth=width)
    set_month(ax)
    ax.plot(time_a,10**(np.ones((len(time_a)))*np.log10(0.3)),'b--')

    ax = fig.add_subplot([0.25,0.05,0.15,0.25])
    f,g = extract_bar(-50,-55,clim,lat,lon,ax,ax1,time_a,fig,yticks=False,label = 'Month',unc=True,unc_data = clim_std)
    ax.text(0.85,0.95,f'(c)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**clim2[g,f,:],'b-',linewidth=width)
    set_month(ax)
    ax.plot(time_a,10**(np.ones((len(time_a)))*np.log10(0.3)),'b--')

    ax = fig.add_subplot([0.07,0.05,0.15,0.25])
    f,g = extract_bar(-150,55,clim,lat,lon,ax,ax1,time_a,fig,label = 'Month',unc=True,unc_data = clim_std)
    ax.text(0.85,0.95,f'(b)',transform=ax.transAxes,va='top',fontweight='bold',fontsize = 25)
    ax.plot(time_a,10**clim2[g,f,:],'b-',linewidth=width)
    set_month(ax)
    ax.plot(time_a,10**(np.ones((len(time_a)))*np.log10(0.3)),'b--')

    cbar = plt.colorbar(pc,ax=ax1)
    cbar.set_label('Chlorophyll-a (mg m$^{-3}$)');
    worldmap.plot(color="lightgrey", ax=ax1)
    # pc.set_clim([0,0.2])
    fig.savefig(output,dpi=300)

def set_month(ax):
    ax.set_xticks([1,7,12])
    ax.set_xticklabels(['Jan','Jul','Dec'])
    ax.set_xlim([1,12])

def plot_flag(file,ref_time = datetime.datetime(1970,1,15),output='plots/flag_pixels.png',res=0.25):

    font = {'weight' : 'normal',
            'size'   : 14}
    matplotlib.rc('font', **font)
    cols = ['#332288','#44AA99','#882255','#DDCC77', '#117733', '#88CCEE','#999933','#CC6677']
    lab = ['OC-CCI Observation','Cloud Kriging','Under Ice','Southern Ocean Backwards','Southern Ocean Forwards','Arctic Backwards','Arctic Forwards','Final Kriging']
    c = Dataset(file,'r')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    flag = np.array(c['chl_flag'])
    time = np.array(c['time'])
    c.close()
    area = np.transpose(du.area_grid(lon,lat,res))
    time_a = []
    for i in range(len(time)):
        time_a.append(ref_time + datetime.timedelta(days=int(time[i])))
    for i in range(len(time)):
        time_a[i] = time_a[i].year + ((time_a[i].month-1)/12)

    fig = plt.figure(figsize=(15,7))
    gs = GridSpec(1,1, figure=fig, wspace=0.25,hspace=0.15,bottom=0.1,top=0.95,left=0.1,right=0.95)
    ax = fig.add_subplot(gs[0,0])

    for i in range(len(time)):
        print(time_a)
        flag_t = flag[:,:,i]
        tot = np.sum(np.ravel(area[np.where(flag_t != 1)]))
        #print(tot)
        bottom = 0
        for j in range(2,10,1):
            f = np.sum(np.ravel(area[np.where(flag_t == j)]))
            if i == 1:
                ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[j-2],label = lab[j-2],width = 1/12)
            else:
                ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[j-2],width = 1/12)
            bottom = bottom + (f/tot)*100
    ax.legend()
    ax.set_xlim([np.min(time_a),np.max(time_a)])
    ax.set_ylim([0,100])
    ax.set_ylabel('Percentage area contribution to flagging (%)')
    ax.set_xlabel('Year')
    fig.savefig(output,dpi=300)

def plot_flag_l(file,ref_time = datetime.datetime(1970,1,15),output = 'plots/flagl_pixels.png',res=0.25):
    font = {'weight' : 'normal',
            'size'   : 14}
    matplotlib.rc('font', **font)
    cols = ['#332288','#44AA99','#882255','#DDCC77', '#117733', '#88CCEE','#999933','#CC6677','k']
    #lab = ['OC-CCI Observation','Cloud Kriging','Under Ice','Southern Ocean Backwards','Southern Ocean Forwards','Arctic Backwards','Arctic Forwards','Final Kriging']
    c = Dataset(file,'r')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    flag = np.array(c['flag_l'])
    time = np.array(c['time'])
    c.close()
    area = np.transpose(du.area_grid(lon,lat,res))
    time_a = []
    for i in range(len(time)):
        time_a.append(ref_time + datetime.timedelta(days=int(time[i])))
    for i in range(len(time)):
        time_a[i] = time_a[i].year + ((time_a[i].month-1)/12)

    fig = plt.figure(figsize=(15,7))
    gs = GridSpec(1,1, figure=fig, wspace=0.25,hspace=0.15,bottom=0.1,top=0.95,left=0.1,right=0.95)
    ax = fig.add_subplot(gs[0,0])

    for i in range(len(time)):
        #print(time_a)
        flag_t = flag[:,:,i]
        tot = np.sum(np.ravel(area[np.where(flag_t != 0)]))
        #print(tot)
        bottom = 0
        t = 0
        for j in [-8,-6,-4,-2, 1,3,5,7]:#range(-8,9,2):
            #print(t)
            #print(j)
            f = np.sum(np.ravel(area[np.where((flag_t == j) | (flag_t == j+1))]))
            if i == 0:
                ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[t],label = str(j) +' + ' +str(j+1),width = 1/12)
            else:
                ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[t],width = 1/12)
            bottom = bottom + (f/tot)*100
            t=t+1

        f = np.sum(np.ravel(area[np.where((flag_t == 9) | (flag_t == -9))]))
        if i == 0:
            ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[t],label = str(9) +' + ' +str(-9),width = 1/12)
        else:
            ax.bar(time_a[i],(f/tot)*100,bottom=bottom,color=cols[t],width = 1/12)
        bottom = bottom + (f/tot)*100
        t=t+1
    ax.legend()
    ax.set_xlim([np.min(time_a),np.max(time_a)])
    ax.set_ylim([0,100])
    ax.set_ylabel('Percentage area contribution to filling (%)')
    ax.set_xlabel('Year')
    fig.savefig(output,dpi=300)

def plot_chla_scatter(file_occci,file_chl,output = 'plots/insitu_verification.png'):
    font = {'weight' : 'normal',
            'size'   :12}
    matplotlib.rc('font', **font)
    sz = 8

    unit = 'log$_{10}$(mgm$^{-3}$)'
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
    fig = plt.figure(figsize=(7,7))
    gs = GridSpec(1,1, figure=fig, wspace=0.25,hspace=0.25,bottom=0.1,top=0.95,left=0.13,right=0.95)
    ax = fig.add_subplot(gs[0,0])
    # chla_t = chla[:,f,:]
    # flags_t = flags[:,f,:]
    # chl_in = chl_in[:,f,:]
    f = np.where(((flags == 5) | (flags == 6)) & (np.isnan(chl_in) == 0))
    ax.scatter(chl_in[f],chla[f],zorder=4,color='b',label='Southern Hemisphere',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(c,c*stats['slope']+stats['intercept'],'b--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.50,0.25,f'Southern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)


    f = np.where(((flags == 7) | (flags == 8)) & (np.isnan(chl_in) == 0))
    ax.scatter(chl_in[f],chla[f],zorder=4,color='r',label='Northern Hemisphere',s=sz)
    stats = ws.unweighted_stats(chl_in[f],chla[f],'a')
    ax.plot(c,c*stats['slope']+stats['intercept'],'r--')
    rmsd = '%.2f' %np.round(stats['rmsd'],2); bias = '%.2f' %np.round(stats['med_rel_bias'],2); sl = '%.2f' %np.round(stats['slope'],2);
    ip = '%.2f' %np.round(stats['intercept'],2); n = stats['n']
    ax.text(0.03,0.97,f'Northern Hemisphere\nRMSD = {rmsd} {unit}\nBias = {bias} {unit}\nSlope = {sl}\nIntercept = {ip}\nN = {n}',transform=ax.transAxes,va='top',fontsize=14)

    ax.plot(c,c,'k-')
    ax.set_xlim(c); ax.set_ylim(c);
    ax.legend(loc=3)
    ax.set_xlabel('Wintertime in situ Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
    ax.set_ylabel('Wintertime gapfilled Chlorophyll-a (log$_{10}$(mgm$^{-3}$))')
    fig.savefig(output,dpi=300)
