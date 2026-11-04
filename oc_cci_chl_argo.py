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

let = ['a','b','c','d','e','f','g','h','i','j']

def download_occci_monthly(loc,start_yr,end_yr):
    from ftplib import FTP

    ftp_loc = 'oceancolour.org'
    username = 'oc-cci-data'
    password = 'ELaiWai8ae'
    ftp_fold = '/occci-v6.0/geographic/netcdf/monthly/chlor_a'
    ftp = FTP(ftp_loc)
    ftp.login(username,password)
    du.makefolder(loc)
    for i in range(start_yr,end_yr+1):
        du.makefolder(os.path.join(loc,str(i)))
        ftp.cwd(ftp_fold+'/'+str(i))
        filenames = ftp.nlst()
        for filename in filenames:
            print(filename)
            local_filename = os.path.join(loc,str(i), filename)
            file = open(local_filename, 'wb')
            ftp.retrbinary('RETR '+ filename, file.write)
            file.close()

    ftp.quit()

def download_occci_8day(loc,start_yr,end_yr):
    from ftplib import FTP

    ftp_loc = 'oceancolour.org'
    username = 'oc-cci-data'
    password = 'ELaiWai8ae'
    ftp_fold = '/occci-v6.0/geographic/netcdf/8day/chlor_a'
    ftp = FTP(ftp_loc)
    ftp.login(username,password)
    du.makefolder(loc)
    for i in range(start_yr,end_yr+1):
        du.makefolder(os.path.join(loc,str(i)))
        ftp.cwd(ftp_fold+'/'+str(i))
        filenames = ftp.nlst()
        for filename in filenames:
            print(filename)
            local_filename = os.path.join(loc,str(i), filename)
            file = open(local_filename, 'wb')
            ftp.retrbinary('RETR '+ filename, file.write)
            file.close()

    ftp.quit()

def download_osisaf_daily(loc,start_yr,end_yr):
    from ftplib import FTP
    du.makefolder(loc)

    ftp_loc = 'osisaf.met.no'

    loc1 = '/reprocessed/ice/conc/v3p1' # This is the main climate data record (OSI-430) covering 1978 to 2020
    loc2 = '/reprocessed/ice/conc-cont-reproc/v3p0' # This is the interim climate data record (OSI-430-a1) covering 2021 to August 2025
    loc3 = '/reprocessed/ice/conc-cont-reproc-amsr/v3p0' # This is the recommended itnerim climate record extension as OSI-430-a1 is currently have date ingestion issues

    switch_date = [[1,2021],[9,2025]]

    ftp = FTP(ftp_loc)
    ftp.login()
    du.makefolder(loc)
    t = 0
    for i in range(start_yr,end_yr+1):
        for j in range(1,13):
            if (i==switch_date[0][1]) & (j == switch_date[0][0]):
                t=1
            elif (i==switch_date[1][1]) & (j == switch_date[1][0]):
                t=2

            if t == 0:
                ftp_fold = loc1
            elif t==1:
                ftp_fold=loc2
            else:
                ftp_fold=loc3
            du.makefolder(os.path.join(loc,str(i),du.numstr(j)))
            ftp.cwd(ftp_fold+'/'+str(i)+'/'+du.numstr(j))
            filenames = ftp.nlst()
            for filename in filenames:
                print(filename)
                local_filename = os.path.join(loc,str(i),du.numstr(j), filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+ filename, file.write)
                file.close()
    ftp.quit()

def generate_occci(occci_file,res,start_yr,end_yr,area_wei=False,gebco_file = False,gebco_out=False,land_mask=False,chla_loc='F:/Data/OC-CCI/monthly/chlor_a',cci_land_file='E:/Data/Land-CCI/ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0',osisaf_loc='F:/Data/OSISAF',multiprocess=False):
    lon,lat = du.reg_grid(lon=res,lat=res)
    import CCI_OC_SPATIAL_AV as OC
    import Data_Loading.ESA_CCI_land as landcci
    landcci.generate_land_cci(cci_land_file+'.nc',cci_land_file+'.tif',lon,lat,gebco_out)

    du.makefolder(os.path.join(chla_loc,str(res)+'DEG_weighted'))
    if multiprocess:
        OC.multiprocess_oc_cci_average(chla_loc,os.path.join(chla_loc,str(res)+'DEG_weighted'),start_yr = start_yr,end_yr = end_yr,log=lon,lag=lat,area_wei=area_wei,gebco_file = gebco_file,gebco_out=gebco_out,land_mask=land_mask)
    else:
        OC.oc_cci_average(chla_loc,os.path.join(chla_loc,str(res)+'DEG_weighted'),start_yr = start_yr,end_yr = end_yr,log=lon,lag=lat,area_wei=area_wei,gebco_file = gebco_file,gebco_out=gebco_out,land_mask=land_mask)

    from Data_Loading.OSISAF_download import OSISAF_spatial_average, OSISAF_merge_hemisphere, OSISAF_monthly_av, multiprocess_OSISAF_monthly_av, multiprocess_OSISAF_spatial_average, multiprocess_OSISAF_merge_hemisphere
    if multiprocess:
        multiprocess_OSISAF_monthly_av(osisaf_loc,start_yr=start_yr,end_yr=end_yr)
        multiprocess_OSISAF_spatial_average(data=os.path.join(osisaf_loc,'monthly'),out_loc=os.path.join(osisaf_loc,'monthly',str(res)+'DEG'),start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'nh')
        multiprocess_OSISAF_spatial_average(data=os.path.join(osisaf_loc,'monthly'),out_loc=os.path.join(osisaf_loc,'monthly',str(res)+'DEG'),start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'sh')
        multiprocess_OSISAF_merge_hemisphere(os.path.join(osisaf_loc,'monthly',str(res)+'DEG'), gebco_out,start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat)
    else:
        OSISAF_monthly_av(osisaf_loc,start_yr=start_yr,end_yr=end_yr)
        OSISAF_spatial_average(data=os.path.join(osisaf_loc,'monthly'),out_loc=os.path.join(osisaf_loc,'monthly',str(res)+'DEG'),start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'NH')
        OSISAF_spatial_average(data=os.path.join(osisaf_loc,'monthly'),out_loc=os.path.join(osisaf_loc,'monthly',str(res)+'DEG'),start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat,hemi = 'SH')
        OSISAF_merge_hemisphere(os.path.join(osisaf_loc,'monthly',str(res)+'DEG'), gebco_out,start_yr=start_yr,end_yr=end_yr,log=lon,lag=lat)

    import construct_input_netcdf as cinp
    vars = [['OC-CCI','chlor_a',os.path.join(chla_loc,str(res)+'DEG_weighted','%Y','%Y_%m_*.nc'),0],
        ['OC-CCI','chlor_a_log10_rmsd',os.path.join(chla_loc,str(res)+'DEG_weighted','%Y','%Y_%m_*.nc'),0],
        ['OSISAF','ice_conc',os.path.join(osisaf_loc,'monthly',str(res)+'DEG','%Y','%Y%m*COM.nc'),0]]

    cinp.driver(occci_file,vars,start_yr = start_yr,end_yr = end_yr,lon = lon,lat = lat,fill_clim=False)
    c = Dataset(occci_file,'a')
    c.variables['OC-CCI_chlor_a'].file_location = vars[0][2]
    c.variables['OC-CCI_chlor_a_log10_rmsd'].file_location = vars[1][2]
    c.variables['OSISAF_ice_conc'].file_location = vars[2][2]
    c.close()

def chl_argo_relationship(res,start_yr,end_yr,files,occci_file,plot=False,lags = 9,netcdf_loc = 'E:/SCOPE/Argo'):
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
    #generate_occci(occci_file,res,lon,lat,start_yr,end_yr,area_wei=area_wei,gebco_file = gebco_file,gebco_out=gebco_out,land_mask=land_mask,chla_loc=chla_loc,cci_land_file=cci_land_file,osisaf_loc=osisaf_loc,multiprocess=multiprocess)
    long,latg = np.meshgrid(lon,lat)
    long = np.transpose(long)
    latg = np.transpose(latg)
    for file in files:
        app = file + '_' + str(res)
        file2 = os.path.split(file)
        argo_file = os.path.join(netcdf_loc,'netcdf',file2[1]+'_'+str(res)+'deg.nc')


        c = Dataset(argo_file,'r')
        argo = 10**np.array(c['chl'])
        # argo[argo<0.005] = np.nan
        c.close()

        c = Dataset(occci_file,'r')
        occci = 10**np.array(c['OC-CCI_chlor_a'])
        time = np.array(c['time'])
        c.close()


        print(argo.shape)
        print(occci.shape)

        f = np.where(np.isnan(argo) == 0)
        print(f)
        print(len(f[0]))

        arg = []
        oc = []
        dif = []
        latgg = []
        longg = []
        for i in range(0,len(f[0])):
            l = 0
            t = 0
            while t == 0:
                oct = occci[f[0][i],f[1][i],f[2][i]-l]
                if np.isnan(oct) == 0:
                    arg.append(argo[f[0][i],f[1][i],f[2][i]])
                    oc.append(oct)
                    dif.append(l)
                    longg.append(long[f[0][i],f[1][i]])
                    latgg.append(latg[f[0][i],f[1][i]])
                    t = 1
                else:
                    l = l+1
                if l == lags+1:
                    t=1
        arg = np.array(arg); oc = np.array(oc); dif = np.array(dif); latgg = np.array(latgg); longg = np.array(longg)
        # fig = plt.figure(figsize=(18,7))

        meds = np.zeros((lags+1,4))
        for i in range(0,lags+1):
            g = np.where((dif == i))
            print(len(g[0]))

            if i == 0:
                vals = arg[g] / oc[g]
                med = np.nanmedian(vals)

            vals = ((arg[g]/med) -oc[g]) / oc[g]

            if plot:
                if len(g[0]) != 0:
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
        np.savetxt(os.path.join(netcdf_loc,'relationships',f'backwards_'+file2[1]+'_'+str(res)+'deg.csv'),meds,delimiter=',',header='month_diff,% difference,% difference uncertainty,number_samples')
        output_data = np.transpose(np.vstack((arg/med,oc,dif,latgg,longg)))
        np.savetxt(os.path.join(netcdf_loc,'relationships',f'backwards_'+file2[1]+'_'+str(res)+'deg_underlyingdata.csv'),output_data,delimiter=',',header='Argo_chl-a,OCCCI_chl-a,time_lag_difference,latitude,longitude')
        # np.savetxt('relationships/median_'+file+'_'+str(res)+'deg.csv',np.array(med),delimiter=',',header='median')
        print(file)
        print('Bias Correction: ' + str(med))

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
        latgg = []
        longg = []
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
                        longg.append(long[f[0][i],f[1][i]])
                        latgg.append(latg[f[0][i],f[1][i]])
                        t = 1
                    else:
                        l = l+1
                    if l == lags+1:
                        t=1
        arg = np.array(arg); oc = np.array(oc); dif = np.array(dif); latgg = np.array(latgg); longg = np.array(longg)
        # if plot:
        # fig = plt.figure(figsize=(18,7))
        meds = np.zeros((lags+1,4))
        for i in range(0,lags+1):
            g = np.where((dif == i))
            # vals = (arg[g]-oc[g])/((oc[g]))
            # if i == 0:
            #     med = np.nanmedian(vals)
            vals = ((arg[g]/med) -oc[g]) / oc[g]
            if plot:
                if len(g[0]) != 0:
                    axs[p].boxplot(vals*100,positions=[i],widths=widths,labels = [f'{i} \n N = {len(g[0])} \n Med = {str(int(round(np.nanmedian(vals*100),0)))} '])
            meds[i,0] = i
            meds[i,1] = np.nanmedian(vals)
            meds[i,2] = np.nanmedian(np.abs(vals - meds[i,1])) * 1.4826 # Converting median absolute deviation to a robust standard deviation equivalent
            meds[i,3] = len(g[0])
            # print(meds[i,2])
        if plot:
            axs[p].set_ylim([-100,100])
            axs[p].text(0.93,0.97,f'('+let[p]+')',transform=axs[p].transAxes,va='top',fontweight='bold',fontsize = 25)
            #axs[p].set_ylabel('Percentage difference to OC-CCI chl-a')
            p=p+1


        np.savetxt(os.path.join(netcdf_loc,'relationships',f'forwards_'+file2[1]+'_'+str(res)+'deg.csv'),meds,delimiter=',',header='month_diff,% difference,% difference uncertainty,number_samples')
        output_data = np.transpose(np.vstack((arg/med,oc,dif,latgg,longg)))
        np.savetxt(os.path.join(netcdf_loc,'relationships',f'forwards_'+file2[1]+'_'+str(res)+'deg_underlyingdata.csv'),output_data,delimiter=',',header='Argo_chl-a,OCCCI_chl-a,time_lag_difference,latitude,longitude')
        #plt.show()
    if plot:
        fig.savefig(os.path.join(netcdf_loc,'plots',f'relationships_'+str(res)+'deg.png'),dpi=300)

def plot_spatial_relationship(res,files,lags = 9,netcdf_loc = 'E:/SCOPE/Argo'):
    worldmap = gpd.read_file(gpd.datasets.get_path("ne_50m_land"))
    cmap = cmocean.cm.balance
    cmap = cmocean.tools.crop_by_percent(cmap, 10, which='both', N=None)
    ##Forwards - Spring
    font = {'weight' : 'normal',
            'size'   :22}
    matplotlib.rc('font', **font)
    fig = plt.figure(figsize=(28,5*7))
    gs = GridSpec(5,2, figure=fig, wspace=0.1,hspace=0.2,bottom=0.05,top=0.95,left=0.03,right=0.92)
    axs = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(5)]
    flatList = [element for innerList in axs for element in innerList]
    axs = flatList
    cbar = fig.add_axes([0.92,0.3,0.01,0.4])
    for file in files:
        file2 = os.path.split(file)
        data = np.loadtxt(os.path.join(netcdf_loc,'relationships',f'forwards_'+file2[1]+'_'+str(res)+'deg_underlyingdata.csv'),delimiter=',',skiprows=1)
        print(data.shape)
        for i in range(lags+1):
            print(i)
            f = np.where(data[:,2] == i)[0]
            a = axs[i].scatter(data[f,4],data[f,3],c=((data[f,0] - data[f,1])/data[f,1])*100,vmin=-100,vmax=100,cmap=cmap)
            worldmap.plot(ax=axs[i],color='lightgrey')
            axs[i].set_title('Month lag ' + str(i))
            axs[i].text(0.93,0.97,f'('+let[i]+')',transform=axs[i].transAxes,va='top',fontweight='bold',fontsize = 25)
    cba = fig.colorbar(a,cax=cbar)
    cba.set_label('Percentage difference to OC-CCI chl-a')
    fig.savefig(os.path.join(netcdf_loc,'plots',f'spatial_fowards_difference_'+str(res)+'deg.png'))
    plt.close(fig)

    ##Backwards - Autumn

    font = {'weight' : 'normal',
            'size'   :22}
    matplotlib.rc('font', **font)
    fig = plt.figure(figsize=(28,5*7))
    gs = GridSpec(5,2, figure=fig, wspace=0.1,hspace=0.2,bottom=0.05,top=0.95,left=0.03,right=0.92)
    axs = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(5)]
    flatList = [element for innerList in axs for element in innerList]
    axs = flatList
    cbar = fig.add_axes([0.92,0.3,0.01,0.4])
    for file in files:
        file2 = os.path.split(file)
        data = np.loadtxt(os.path.join(netcdf_loc,'relationships',f'backwards_'+file2[1]+'_'+str(res)+'deg_underlyingdata.csv'),delimiter=',',skiprows=1)
        print(data.shape)
        for i in range(lags+1):
            print(i)
            f = np.where(data[:,2] == i)[0]
            a = axs[i].scatter(data[f,4],data[f,3],c=((data[f,0] - data[f,1])/data[f,1])*100,vmin=-100,vmax=100,cmap=cmap)
            worldmap.plot(ax=axs[i],color='lightgrey')
            axs[i].set_title('Month lag ' + str(i))
            axs[i].text(0.93,0.97,f'('+let[i]+')',transform=axs[i].transAxes,va='top',fontweight='bold',fontsize = 25)
    cba = fig.colorbar(a,cax=cbar)
    cba.set_label('Percentage difference to OC-CCI chl-a')
    fig.savefig(os.path.join(netcdf_loc,'plots',f'spatial_backwards_difference_'+str(res)+'deg.png'))
    plt.close(fig)
