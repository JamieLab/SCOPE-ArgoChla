#!/usr/bin/env python3
"""
Created by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 08/2023
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter,FFMpegWriter
import os
import cmocean
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import geopandas as gpd
import matplotlib as mpl

mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Users\\df391\\OneDrive - University of Exeter\\Python\ffmpeg\\bin\\ffmpeg.exe'

worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))

def animate(i):
    fig.clear()
    gs = GridSpec(2,1, figure=fig, wspace=0.2,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.98)
    ax = fig.add_subplot(gs[0,0])
    worldmap.plot(color="lightgrey", ax=ax)
    c = ax.pcolor(lon,lat,np.transpose(chla[:,:,i]),cmap=cmocean.cm.algae,vmin=-2,vmax=1)
    c2 = fig.colorbar(c,ax=ax)
    c2.set_label('Chl a (log$_{10}$(mg m$^{-3}$))')

    ax = fig.add_subplot(gs[1,0])
    worldmap.plot(color="lightgrey", ax=ax)
    cmdict = cmocean.cm.thermal(9)
    c = ax.pcolor(lon,lat,np.transpose(flag[:,:,i]),cmap='tab10',vmin=2,vmax=11)
    c2 = fig.colorbar(c,ax=ax)
    c2.set_label('Chl a Flag ')

c = Dataset('oc-cci_chla_corrected.nc','r')
lon = np.array(c['longitude'])
lat = np.array(c['latitude'])
chla = np.array(c['chl'])
flag = np.array(c['flag'])
flag[flag == 1] = np.nan
c.close()
fig = plt.figure(figsize=(15,15))
ani = FuncAnimation(fig, animate, interval=40, blit=False, repeat=True,frames=chla.shape[2])
ani.save('animated.mp4', dpi=300, writer=FFMpegWriter(fps=3))
