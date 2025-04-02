#!/usr/bin/env python3

"""
Code for extracting delayed mode chl-a from Bio-Argo profilers in a set latitude region (or globally...)
This code requires the argo_bio-profile_index.txt (list of all avaiable profiles) off the Argo FTP server.
This file is then interrogated for all the profiles with a CHL-A as a parameter, and fits within the
latitude band definition. The script then runs through the profile list and downloads the file using the FTP
protocol, and extracts the chl-a information. It then saves this information as a text file (.csv), with year, month,
day, lat, lon, chl-a (in the top 20m). This csv file is then given to the argo_average.py script to average these data
onto a monthly grid (.nc file).
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import cmocean
import geopandas as gpd
import matplotlib.transforms
import os
import glob
from netCDF4 import Dataset
import datetime
from ftplib import FTP

def makefolder(fold):
    if not os.path.exists(fold):
        os.makedirs(fold)

def checkfileexist(file):
    #print(file)
    g = glob.glob(file)
    #print(g)
    if not g:
        return False
    else:
        return True

def argo_extract(argo_file,argo_loc, lat_bounds,output_files,plot = False):

    ref = datetime.datetime(1950,1,1,0,0,0) # This is the reference time within the Argo files
    ftp_loc = 'ftp.ifremer.fr' # FTP Location
    loc = '/ifremer/argo/dac' # Location within the FTP that the Argo files are
    #Make the local argo folder
    makefolder(argo_loc)
    # Load the argo profiles list
    for k in range(len(lat_bounds)):
        argo_list = pd.read_table(argo_file,delimiter=',',skiprows=8)

        print(len(argo_list))

        # Split out the data within our region.
        f = np.where((argo_list['latitude'] <= lat_bounds[k][0]) & (argo_list['latitude']>= lat_bounds[k][1]))[0]
        argo_list = argo_list[argo_list['latitude']<= lat_bounds[k][0]]
        argo_list = argo_list[argo_list['latitude']>= lat_bounds[k][1]]
        print(len(argo_list))

        # And find profiles that contain CHLA (as thats the parameter we want..)
        argo_list = argo_list[argo_list['parameters'].str.contains("CHLA")]
        print(len(argo_list))

        print(argo_list)

        chl = []
        for i in range(0,len(argo_list)):
            print(str(i) + ': ' +argo_list['file'].iloc[i])
            s = argo_list['file'].iloc[i].split('/')
            makefolder(os.path.join(argo_loc,s[0]))
            makefolder(os.path.join(argo_loc,s[0],s[1]))
            makefolder(os.path.join(argo_loc,s[0],s[1],s[2]))
            if (checkfileexist(os.path.join(argo_loc,argo_list['file'].iloc[i])) == False) & (argo_list['file'].iloc[i].split('/')[-1][1] != 'R'): # Second term to only get delayed (so corrected) files....
                print('Downloading: ' + argo_list['file'].iloc[i])
                ftp = FTP(ftp_loc)
                ftp.login()
                l = argo_list['file'].iloc[i].split('/')
                #print(loc+'/'+'/'.join(l))
                ftp.cwd(loc+'/'+'/'.join(l[0:-1]))
                #print(l[-1])
                ftp.retrbinary("RETR " + loc+'/'+'/'.join(l) ,open(os.path.join(argo_loc,argo_list['file'].iloc[i]), 'wb').write)
                ftp.close()

            if (checkfileexist(os.path.join(argo_loc,argo_list['file'].iloc[i])) == True) & (argo_list['file'].iloc[i].split('/')[-1][1] != 'R'):
                c = Dataset(os.path.join(argo_loc,argo_list['file'].iloc[i]))
                #print(c.variables.keys())
                if 'CHLA_ADJUSTED' in c.variables.keys():
                    print('CHLA_ADJUSTED is variable!')

                    depth = np.array(c['PRES'])
                    ch = np.array(c['CHLA_ADJUSTED']).astype(np.float64)
                    qf = np.array(c['CHLA_ADJUSTED_QC'])
                    qf2 = np.zeros(qf.shape)
                    #print(qf2.shape)
                    for i in range(qf.shape[0]):
                        for j in range(qf.shape[1]):
                            try:
                                qf2[i,j] = int(qf[i,j])
                            except:
                                qf2[i,j] = 9

                    ch[qf2>1] = np.nan

                    qf2[np.isnan(ch) == 1] = np.nan

                    cp = np.sum(np.isnan(ch)==False,axis=1)
                    print(cp)
                    f = np.where(cp != 0)[0]
                    print(f)
                    if len(f)>0:
                        print(f)
                        ch = ch[f,:]
                        depth = depth[f,:]
                        qf2 = qf2[f,:]
                        lat = np.array(c['LATITUDE'])[f]
                        lon = np.array(c['LONGITUDE'])[f]
                        date = np.array(c['JULD'])[f]
                        print(cp)
                        flo_date = ref + datetime.timedelta(days=int(date))
                        c.close()
                        #ch[ch>30] = np.nan

                        f = np.where(depth < 20)
                        ch = np.nanmean(ch[f])
                        qf2 = np.nanmean(qf2[f])
                        print([flo_date.year,flo_date.month,flo_date.day,lat[0],lon[0],ch,qf2])
                        if (flo_date.year<2024) & (lat <90) & (lon < 180) & (np.isnan(ch) == 0) & (ch>0.001) :
                            chl.append([flo_date.year,flo_date.month,flo_date.day,lat[0],lon[0],ch,qf2])
                    else:
                        c.close()
            else:
                print('Skipping...')

        #print(chl)
        chl = np.array(chl)
        np.savetxt('csv/'+output_files[k]+".csv", chl, delimiter=",",header='Year,Month,Day,Latitude (deg N),Longitude (deg W),chlorphyll-a (mgm-3),Chlorophyll-a quality flag')
        if plot:
            worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))

            fig = plt.figure(figsize=(15,15))
            gs = GridSpec(2,3, figure=fig, wspace=0.2,hspace=0.25,bottom=0.07,top=0.98,left=0.07,right=0.95)
            ax = fig.add_subplot(gs[0, :])
            worldmap.plot(color='lightgrey',ax=ax)
            a = ax.scatter(chl[:,4],chl[:,3],c = np.log10(chl[:,5]),vmin=-1.5,vmax=0)
            plt.colorbar(a)
            ax = fig.add_subplot(gs[1, 0])
            plt.hist(chl[:,0])
            ax = fig.add_subplot(gs[1, 1])
            plt.hist(chl[:,1])
            ax = fig.add_subplot(gs[1, 2])
            plt.hist(np.log10(chl[:,5]),np.arange(-3.5,1,0.05))
            fig.savefig('plots/'+output_files[k]+'.jpg',dpi=300)
