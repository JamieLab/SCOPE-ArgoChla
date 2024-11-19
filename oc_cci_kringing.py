import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import urllib
import glob
from netCDF4 import Dataset
import datetime
import sys
import skgstat as skg
import scipy.interpolate

sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU')
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU\Data_Loading')

import data_utils as du

def save_chla(chla,flag,flag_l,lat,lon,time,out_loc,time_r):
    # outp = Dataset(file,'w',format='NETCDF4_CLASSIC')
    # outp.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
    # outp.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
    # outp.createDimension('lon',lon.shape[0])
    # outp.createDimension('lat',lat.shape[0])
    # outp.createDimension('time',len(time))
    #
    # sst_o = outp.createVariable('chl','f4',('lon','lat','time'),zlib=False)
    # sst_o[:] = chla
    # sst_o.units = 'log10(mgm-3)'
    # sst_o.standard_name = 'Chlorophyll-a concentration'
    #
    # sst_o = outp.createVariable('flag','f4',('lon','lat','time'),zlib=False)
    # sst_o[:] = flag
    # sst_o.standard_name = 'Chlorophyll-a concentration flag'
    #
    # sst_o = outp.createVariable('flag_l','f4',('lon','lat','time'),zlib=False)
    # sst_o[:] = flag_l
    # sst_o.standard_name = 'Chlorophyll-a concentration flag for number of months used in Argo approach'
    #
    # lat_o = outp.createVariable('latitude','f4',('lat'))
    # lat_o[:] = lat
    # lat_o.units = 'Degrees'
    # lat_o.standard_name = 'Latitude'
    # lon_o = outp.createVariable('longitude','f4',('lon'))
    # lon_o.units = 'Degrees'
    # lon_o.standard_name = 'Longitude'
    # lon_o[:] = lon
    #
    # lon_o = outp.createVariable('time','f4',('time'))
    # lon_o.units = 'days since 1970-01-15'
    # lon_o.standard_name = 'Time'
    # lon_o[:] = time
    # outp.close()
    if outloc:
        for i in range(len(time)):
            d = datetime.datetime(int(time_r[i,0]),int(time_r[i,1]),15)
            du.makefolder(os.path.join(out_loc,d.strftime('%Y')))
            file_t = os.path.join(out_loc,d.strftime('%Y/%Y_%m_oc_cci_chla_filled.nc'))

            outp = Dataset(file_t,'w',format='NETCDF4_CLASSIC')
            outp.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
            outp.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
            outp.createDimension('lon',lon.shape[0])
            outp.createDimension('lat',lat.shape[0])

            sst_o = outp.createVariable('chl','f4',('lon','lat'),zlib=False)
            sst_o[:] = chla[:,:,i]
            sst_o.units = 'log10(mgm-3)'
            sst_o.standard_name = 'Chlorophyll-a concentration'

            sst_o = outp.createVariable('flag','f4',('lon','lat'),zlib=False)
            sst_o[:] = flag[:,:,i]
            sst_o.standard_name = 'Chlorophyll-a concentration flag'

            lat_o = outp.createVariable('latitude','f4',('lat'))
            lat_o[:] = lat
            lat_o.units = 'Degrees'
            lat_o.standard_name = 'Latitude'
            lon_o = outp.createVariable('longitude','f4',('lon'))
            lon_o.units = 'Degrees'
            lon_o.standard_name = 'Longitude'
            lon_o[:] = lon
            outp.close()

# plot = False
# reset = False
# daystep = 8 / 30.5
# ref_step = 30.5
# ice = True
# ice_name = 'OSISAF_sea_ice_fraction'#'OSISAF_ice_conc'
# res = 0.25
# outloc = False#'E:/SCOPE/NN/Ford_et_al_SOM_chla/inputs/chla'
# # file = f'netcdf/oc-cci_chlor_a_{res}deg.nc'
# file = f'E:/SCOPE/Argo/8day/out.nc'
# bathy_file = f'F:/Data/Bathymetry/{res}DEG_GEBCO_2023.nc'

def unc_montecarlo(chla,chla_unc,per,per_unc,ens = 1000):
    chla_r = chla + np.random.normal(0,chla_unc,ens)
    per_r = per + np.random.normal(0,per_unc,ens)
    fun =  np.log10((10**chla_r) * np.abs((1-np.abs(per_r))))
    out = np.std(fun)
    return out

def oc_cci_fill(file, bathy_file,daystep,res,reset = False,ice_name = 'OSISAF_sea_ice_fraction',plot=False,ref_step = 30.5,outloc = False,ice = True):
    ref = datetime.datetime(1970,1,15)

    c = Dataset(bathy_file,'r')
    ocean = np.transpose(np.array(c['ocean_proportion']))
    c.close()
    ocean_sum = np.zeros((ocean.shape[0]))
    for i in range(ocean.shape[0]):
        f = np.where(ocean[i,:]>0)
        if len(f[0])>0:
            ocean_sum[i] = len(f[0])
        else:
            ocean_sum[i] = np.nan
    ocean_sum[ocean_sum<1] = np.nan
    print(ocean_sum)
    #
    c = Dataset(file,'a')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    time = np.array(c['time'])

    time_r = np.zeros((len(time),2))
    for i in range(0,len(time)):
        time_r[i,0] = (ref + datetime.timedelta(days=int(time[i]))).year
        time_r[i,1] = (ref + datetime.timedelta(days=int(time[i]))).month
    #
    if reset:
        cci = c['OC-CCI_chlor_a'].shape
        print(cci)
        chla = np.zeros((cci))
        chla[:] = np.nan

        if 'chl_filled' in c.variables.keys():
            c['chl_filled'][:] = chla
        else:
            sst_o = c.createVariable('chl_filled','f4',('longitude','latitude','time'),zlib=False)
            sst_o[:] = chla
            sst_o.units = 'log10(mgm-3)'
            sst_o.standard_name = 'Gap filled Chlorophyll-a concentration'

        if 'chl_filled_unc' in c.variables.keys():
            c['chl_filled_unc'][:] = chla
        else:
            sst_o = c.createVariable('chl_filled_unc','f4',('longitude','latitude','time'),zlib=False)
            sst_o[:] = chla
            sst_o.units = 'log10(mgm-3)'
            sst_o.standard_name = 'Gapfilled Root-mean-square-difference of log10-transformed chlorophyll-a concentration in seawater.'

        if 'chl_flag' in c.variables.keys():
            c['chl_flag'][:] = chla
        else:
            sst_o = c.createVariable('chl_flag','f4',('longitude','latitude','time'),zlib=False)
            sst_o[:] = chla
            sst_o.standard_name = 'Gap filled Chlorophyll-a concentration flag'

        if 'chl_process' in c.variables.keys():
            c['chl_process'][:] = np.zeros((cci[2]))
        else:
            sst_o = c.createVariable('chl_process','f4',('time'),zlib=False)
            sst_o[:] = np.zeros((cci[2]))
            sst_o.standard_name = 'Process indicator'
    c.close()



    c = Dataset(file,'r')
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    cci = c['OC-CCI_chlor_a'].shape
    #cci_unc = c['OC-CCI_chlor_a_log10_rmsd']
    c.close()
    for i in range(0,cci[2]):
        print('Timestep: ' + str(i) + ' of ' + str(cci[2]))
        c = Dataset(file,'a')
        process_flag = c['chl_process'][i]
        if process_flag == 0:
            print('Processing flag = 0')
            cci_data = np.transpose(np.array(c['OC-CCI_chlor_a'][:,:,i]))
            cci_unc = np.transpose(np.array(c['OC-CCI_chlor_a_log10_rmsd'][:,:,i]))
            flag = np.transpose(np.array(c['chl_flag'])[:,:,i])
            if ice:
                ice_data = np.transpose(np.array(c[ice_name])[:,:,i])
            c.close()
            f = np.where(np.isnan(cci_data) == 0)
            flag[f] = 2
            f = np.where(ocean == 0)
            flag[f] = 1

            cci_sum = np.zeros((ocean.shape[0]))
            for j in range(ocean.shape[0]):
                f = np.where(np.isnan(cci_data[j,:])==0)
                #print(f)
                if len(f[0])>0:
                    cci_sum[j] = len(f[0])

            total = cci_sum/ocean_sum
            print(np.nanmean(total))
            if np.nanmean(total) > 0.3:
                print('Total')
                if plot:
                    plt.figure()
                    plt.plot(lat,(ocean_sum))
                    plt.plot(lat,cci_sum)
                    plt.figure()
                    plt.plot(lat,total)
                    plt.plot([-90,90],[0.1,0.1],'k--')
                    plt.show()
                f = np.where(total>0.1)[0]

                gp = [f[0],f[-1]]
                print(gp)
                g = np.squeeze(np.where((lat<=lat[gp[1]]) & (lat>=lat[gp[0]])))
                lat_t = lat[g]
                # lon = lon[g]
                cci_data_t = cci_data[g,:]
                cci_unc_t = cci_unc[g,:]
                # cci = cci[:,g]
                ocean_t = ocean[g,:]
                flag_t = flag[g,:]
                # ocean = ocean[:,g]
                # plt.figure()
                # plt.pcolor(lon,lat,cci)
                # plt.show()
                #
                lon_g,lat_g = np.meshgrid(lon,lat_t)
                print(lon_g.shape)
                print(lat_g.shape)
                f = np.where((np.isnan(cci_data_t)==0) & (ocean_t > 0))
                print(len(f[0]))
                l =int(cci[1]*0.06)
                vlon = np.ravel(lon_g[f])[0::l]
                vlat = np.ravel(lat_g[f])[0::l]
                #
                vchl = np.ravel(cci_data_t[f])
                vunc = np.ravel(cci_unc_t[f])
                print(vchl.size)
                vchl = vchl[0::l]
                vunc = vunc[0::l]
                print(vchl.size)

                V = skg.Variogram(np.column_stack((vlon,vlat)), vchl, n_lags=70,model='exponential', normalize=False,maxlag=100,fit_method='lm')
                if plot:
                    V.plot(show=True)
                    wait = input("Press Enter to continue.")
                ok = skg.OrdinaryKriging(V, min_points=3, max_points=6, mode='exact')
                f = np.where((np.isnan(cci_data_t)==1) & (ocean_t > 0))
                print(len(f[0]))
                field = ok.transform(lon_g[f].flatten(), lat_g[f].flatten())
                cci_data_t[f] = field
                flag_t[f] = 3

                ok = skg.OrdinaryKriging(V, min_points=3, max_points=6, mode='exact',values = vunc)
                f = np.where((np.isnan(cci_unc_t)==1) & (ocean_t > 0))
                field = ok.transform(lon_g[f].flatten(), lat_g[f].flatten())
                cci_unc_t[f] = field

                cci_data[g,:] = cci_data_t
                cci_unc[g,:] = cci_unc_t
                flag[g,:] = flag_t
                if plot:
                    fig, axes = plt.subplots(2, 1, figsize=(21,7))
                    axes[0].pcolor(lon,lat,cci_data)

                    axes[1].pcolor(lon,lat,flag)
                    plt.show()

                if ice:
                    cci_data[(ice_data >= 0.9) & (flag != 2)] = np.log10(0.1)
                    cci_unc[(ice_data >= 0.9) & (flag != 2)] = 0.4
                    flag[(ice_data >= 0.9) & (flag != 2)] = 4

                c = Dataset(file,'a')
                c.variables['chl_filled'][:,:,i] = np.transpose(cci_data)
                c.variables['chl_filled_unc'][:,:,i] = np.transpose(cci_unc)
                c.variables['chl_flag'][:,:,i] = np.transpose(flag)
                c.variables['chl_process'][i] = 1
                #c.sync()
                c.close()
        else:
            print('Processing flag == 1')
            c.close()


    c = Dataset(file,'r')
    flag = np.array(c['chl_flag'])
    chla = np.array(c['chl_filled'])
    unc = np.array(c['chl_filled_unc'])
    ice = np.array(c[ice_name])
    lat = np.array(c['latitude'])
    lon = np.array(c['longitude'])
    c.close()

    k = np.zeros((chla.shape[2]))
    for i in range(chla.shape[2]):
        a = np.where(np.isnan(chla[:,:,i])==0)
        print(a)
        if len(a[0]) > 0:
            k[i] = 1
    print(k)

    flag_l = np.zeros((chla.shape))

    files = ['argo_chla_southernocean','argo_chla_arctic']
    argo_b_data = np.loadtxt('relationships/backwards_'+files[0]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
    argo_b = scipy.interpolate.interp1d(argo_b_data[:,0],argo_b_data[:,1])
    argo_b_unc = scipy.interpolate.interp1d(argo_b_data[:,0],argo_b_data[:,2])

    argo_f_data = np.loadtxt('relationships/forwards_'+files[0]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
    argo_f = scipy.interpolate.interp1d(argo_f_data[:,0],argo_f_data[:,1])
    argo_f_unc = scipy.interpolate.interp1d(argo_f_data[:,0],argo_f_data[:,2])

    argo_b_arctic_data = np.loadtxt('relationships/backwards_'+files[1]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
    argo_b_arctic = scipy.interpolate.interp1d(argo_b_arctic_data[:,0],argo_b_arctic_data[:,1])
    argo_b_arctic_unc = scipy.interpolate.interp1d(argo_b_arctic_data[:,0],argo_b_arctic_data[:,2])

    argo_f_arctic_data = np.loadtxt('relationships/forwards_'+files[1]+'_'+str(res)+'deg.csv',delimiter=',',skiprows=1)
    argo_f_arctic = scipy.interpolate.interp1d(argo_f_arctic_data[:,0],argo_f_arctic_data[:,1])
    argo_f_arctic_unc = scipy.interpolate.interp1d(argo_f_arctic_data[:,0],argo_f_arctic_data[:,2])
    """
    Here we apply the wintertime chl-a relationships developed from the Bio-Argo profilers for the Arctic and Southern Ocean
    seperately. This uses a relative relationship between the nearest in time OC-CCI observation (either backwards or forwards) and the Bio-Argo
    profiler chl-a to construct the wintertime chl-a values.
    """
    # Here we check for pixels that flag as 0 (not dealt with yet) and 4,5,6 (for some reason gone through previous averaging - should be fixed now...), and that aren't
    # highly ice covered.
    f = np.where(((np.isnan(flag) == 1) | (flag == 5) | (flag == 6)) & (lat[np.newaxis,:,np.newaxis] < -40) & (ice <0.9))
    print(f)

    for i in range(len(f[0])):
        t = 0
        l = daystep
        if k[f[2][i]] == 1:
            while t == 0:
                t2 = 0
                while t2 == 0:
                    g = flag[f[0][i],f[1][i],f[2][i]-int(l/daystep)]
                    if (g == 2) | (g == 3):
                        #print(l)
                        if f[2][i]-l >0:
                            chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]-int(l/daystep)]) * np.abs((1-np.abs(argo_b(l)))))
                            flag[f[0][i],f[1][i],f[2][i]] = 5
                            flag_l[f[0][i],f[1][i],f[2][i]] = -l
                            unc[f[0][i],f[1][i],f[2][i]] = unc_montecarlo(chla[f[0][i],f[1][i],f[2][i]-int(l/daystep)],unc[f[0][i],f[1][i],f[2][i]-int(l/daystep)],argo_b(l),argo_b_unc(l))

                            l2 = l
                            t2 = 1
                        else:
                            l2 = 10
                            t2 = 1

                    l=l+daystep
                    if l > 9:
                        l2 = l
                        t2 =1
                t2 = 0
                l = daystep
                while t2 == 0:
                    #l2 = l
                    if f[2][i]+l <chla.shape[2]:
                        g = flag[f[0][i],f[1][i],f[2][i]+int(l/daystep)]
                        if (g == 2) | (g == 3):
                            if l < l2:
                                chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]+int(l/daystep)]) * np.abs((1-np.abs(argo_f(l)))))
                                flag[f[0][i],f[1][i],f[2][i]] = 6
                                flag_l[f[0][i],f[1][i],f[2][i]] = l
                                unc[f[0][i],f[1][i],f[2][i]] = unc_montecarlo(chla[f[0][i],f[1][i],f[2][i]+int(l/daystep)],unc[f[0][i],f[1][i],f[2][i]+int(l/daystep)],argo_f(l),argo_f_unc(l))
                                t=1
                            t2 = 1
                    l = l+daystep
                    if l > 9:
                        t2 = 1
                t = 1

    """
    Artic version
    """
    f = np.where(((np.isnan(flag) == 1) | (flag == 7) | (flag == 8)) & (lat[np.newaxis,:,np.newaxis] > 40) & (ice <0.9))
    print(f)

    for i in range(len(f[0])):
        t = 0
        l = daystep
        if k[f[2][i]] == 1:
            while t == 0:
                t2 = 0
                while t2 == 0:
                    g = flag[f[0][i],f[1][i],f[2][i]-int(l/daystep)]
                    if (g == 2) | (g == 3):
                        print(l)
                        print(int(l/daystep))
                        if f[2][i]-l >0:
                            chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]-int(l/daystep)]) * np.abs((1-np.abs(argo_b_arctic(l)))))
                            flag[f[0][i],f[1][i],f[2][i]] = 7
                            flag_l[f[0][i],f[1][i],f[2][i]] = -l
                            unc[f[0][i],f[1][i],f[2][i]] = unc_montecarlo(chla[f[0][i],f[1][i],f[2][i]-int(l/daystep)],unc[f[0][i],f[1][i],f[2][i]-int(l/daystep)],argo_b_arctic(l),argo_b_arctic_unc(l))
                            l2 = l
                            t2 = 1
                        else:
                            l2 = 10
                            t2 = 1
                    l=l+daystep
                    if l > 9:
                        t2 =1
                t2 = 0
                l = daystep
                while t2 == 0:
                    #l2 = l
                    if f[2][i]+l <chla.shape[2]:
                        g = flag[f[0][i],f[1][i],f[2][i]+int(l/daystep)]
                        if (g == 2) | (g == 3):
                            if l < l2:
                                chla[f[0][i],f[1][i],f[2][i]] = np.log10((10**chla[f[0][i],f[1][i],f[2][i]+int(l/daystep)]) * np.abs((1-np.abs(argo_f_arctic(l)))))
                                flag[f[0][i],f[1][i],f[2][i]] = 8
                                flag_l[f[0][i],f[1][i],f[2][i]] = l
                                unc[f[0][i],f[1][i],f[2][i]] = unc_montecarlo(chla[f[0][i],f[1][i],f[2][i]+int(l/daystep)],unc[f[0][i],f[1][i],f[2][i]+int(l/daystep)],argo_f_arctic(l),argo_f_arctic_unc(l))
                                t=1
                            t2 = 1
                    l = l+daystep
                    if l > 9:
                        t2 = 1
                t = 1

    c = Dataset(file,'a')
    c.variables['chl_filled'][:] = chla
    c.variables['chl_flag'][:] = flag
    c.variables['chl_filled_unc'][:] = unc

    if 'flag_l' in c.variables.keys():
        c['flag_l'][:] = flag_l
    else:
        sst_o = c.createVariable('flag_l','f4',('longitude','latitude','time'),zlib=False)
        sst_o[:] = flag_l
        sst_o.standard_name = 'Chlorophyll-a concentration flag for number of months used in Argo approach'
        sst_o.comment = 'Negative indicates backwards; positive indicates forwards'
    c.close()

    del chla, flag, unc, flag_l

    """
    """

    lon_g,lat_g = np.meshgrid(lon,lat)
    for i in range(0,chla.shape[2]):
        print('Timestep: ' + str(i) + ' of ' + str(chla.shape[2]))
        c = Dataset(file,'r')
        chl_process = c['chl_process'][i]
        if chl_process == 1:
            cci_data = np.transpose(c['chl_filled'][:,:,i])
            flag_t = np.transpose(c['chl_flag'][:,:,i])
            unc_t = np.transpose(c['chl_filled_unc'][:,:,i])
            c.close()
            cci_sum = np.zeros((ocean.shape[0]))
            for j in range(ocean.shape[0]):
                f = np.where(np.isnan(cci_data[j,:])==0)
                #print(f)
                if len(f[0])>0:
                    cci_sum[j] = len(f[0])

            total = cci_sum/ocean_sum
            print(np.nanmean(total))
            if np.nanmean(total) > 0.3:
                print(lon_g.shape)
                print(lat_g.shape)
                f = np.where((np.isnan(cci_data)==0) & (ocean > 0))
                print(len(f[0]))
                l =int(cci_data.shape[1]*0.06)
                vlon = np.ravel(lon_g[f])[0::l]
                vlat = np.ravel(lat_g[f])[0::l]
                #
                vchl = np.ravel(cci_data[f])
                print(vchl.size)
                vchl = vchl[0::l]
                print(vchl.size)
                vunc = np.ravel(unc_t[f])[0::l]


                V = skg.Variogram(np.column_stack((vlon,vlat)), vchl, n_lags=70,model='exponential', normalize=False,maxlag=100,fit_method='lm')
                if plot:
                    V.plot(show=True)
                    wait = input("Press Enter to continue.")
                ok = skg.OrdinaryKriging(V, min_points=3, max_points=6, mode='exact')
                f = np.where((np.isnan(cci_data)==1) & (ocean > 0))
                print(len(f[0]))
                field = ok.transform(lon_g[f].flatten(), lat_g[f].flatten())
                cci_data[f] = field
                flag_t[f] = 9

                ok = skg.OrdinaryKriging(V, min_points=3, max_points=6, mode='exact',values = vunc)
                field = ok.transform(lon_g[f].flatten(), lat_g[f].flatten())
                unc_t[f] =field
                c = Dataset(file,'a')
                c.variables['chl_filled'][:,:,i] = np.transpose(cci_data)
                c.variables['chl_filled_unc'][:,:,i] = np.transpose(cci_unc)
                c.variables['chl_flag'][:,:,i] = np.transpose(flag_t)
                c.variables['chl_process'][i] = 2
                #c.sync()
                c.close()
        else:
            print('Processing != 1')
            c.close()

    # c = Dataset(file,'a')
    # c['chl_filled'][:] = chla
    # c['chl_flag'][:] = flag
    # c['chl_filled_unc'][:] = unc
    # c.close()
    #
    # # c = Dataset(file,'r')
    # # chla = np.array(c['chl_filled'])
    # # flag = np.array(c['chl_flag'])
    # # flag_l = np.array(c['flag_l'])
    # # c.close()
    # # save_chla(chla,flag,flag_l,lat,lon,time,outloc,time_r)
