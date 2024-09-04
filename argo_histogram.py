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
font = {'weight' : 'normal',
        'size'   :14}
matplotlib.rc('font', **font)
worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))
files = ['argo_chl.csv','argo_chl_arctic.csv']
t = 0
fig = plt.figure(figsize=(15,15))
gs = GridSpec(2,2, figure=fig, wspace=0.2,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.98)

ax2 = fig.add_subplot(gs[1,0])
ax4 = fig.add_subplot(gs[1,1])
ax = fig.add_subplot(gs[0,:])
l = []; l2 = []
for file in files:
    data = np.genfromtxt(file, delimiter=',')
    a = ax.scatter(data[:,4],data[:,3],c=np.log10(data[:,5]),vmin=-1,vmax=1)
    l.append(data[:,0]+0.5)
    l2.append(data[:,1])
cbar = plt.colorbar(a)
cbar.set_label('Chl a (log$_{10}$(mg m$^{-3}$))')
worldmap.plot(color="lightgrey", ax=ax)
ax2.hist(l,stacked=False,bins=list(range(2008,2023)))
ax2.set_xticks(np.array(range(2008,2023,2))+0.5,list(range(2008,2023,2)))
ax2.set_xlabel('Year')
ax2.set_ylabel('Profile Frequency')
ax2.legend(['Southern Ocean','Arctic'])
ax4.hist(l2,bins=list(range(1,14)))
ax4.set_xticks(np.array(range(1,13,1))+0.5,list(range(1,13,1)))
ax4.set_xlabel('Month')
ax4.set_ylabel('Profile Frequency')
fig.savefig('plots/argo_histrogram.png',dpi=300)
