#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import cmocean
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from matplotlib.animation import FuncAnimation, PillowWriter,FFMpegWriter
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.colors as colors

mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Users\\df391\\OneDrive - University of Exeter\\Python\ffmpeg\\bin\\ffmpeg.exe'

# worldmap = gpd.read_file(gpd.datasets.get_path('ne_50m_land'))
# print(gpd.datasets.get_path('ne_50m_land'))
file = 'E:/SCOPE/Argo/v0-3/netcdf/Ford_et_al_OC-CCI_chlor_a_gap_filled_1997_2024_0.25deg_v0-3.nc'
logo = mpl.image.imread('C:/Users/df391/OneDrive - University of Exeter/Post_Doc_ESA_Carbon_for_climate/OC4C_icon-positive-txt.png')
logo_scope = mpl.image.imread('C:/Users/df391/OneDrive - University of Exeter/Post_Doc_ESA_SCOPE/SCOPE logo - alternative version - dark text.png')
esa_logo = mpl.image.imread('C:/Users/df391/OneDrive - University of Exeter/Post_Doc_ESA_Carbon_for_climate/ESA_logo-400x284.png')

def date_parser(date):
    return datetime.timedelta(days=int(date))+datetime.datetime(1970,1,15)
c = Dataset(file,'r')
lat = np.array(c['latitude'])
lon = np.array(c['longitude'])
flux = np.array(c['chl_filled'])
flux_unc=np.array(c['chl_filled_unc'])
time_g = np.array(c['time'])

c.close()

time_g = time_g[12:]
flux = flux[:,:,12:]
flux_unc = flux_unc[:,:,12:]

"""
"""
def animate(i,time,x_locs,y_locs,ani_data,data_add):
    print(i)
    fig.clear()
    ax = plt.axes(projection=ccrs.NearsidePerspective(x_locs[i], y_locs[i]))
    ax.set_position([0.05,0.15,0.9,0.7])
    ap = ax.pcolormesh(
        lon,
        lat,
        np.transpose(10**ani_data[:,:,i]),
        #transform=proj,
        transform=ccrs.PlateCarree(),
        cmap=cmocean.cm.algae,
        norm=colors.LogNorm(vmin=10**-2, vmax=10**1),
        zorder=1,
        )
    ax.add_feature(cfeature.GSHHSFeature(scale='l'), facecolor="grey",edgecolor=None,zorder=7)
    # ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='grey', linewidth=0.5, zorder=7)
    cax = plt.axes([0.72,0.95,0.25,0.025])
    cbar = plt.colorbar(ap,cax=cax,orientation='horizontal');cbar.set_label('Chlorophyll-a (mg m$^{-3}$)')
    loax= plt.axes([0.00,0.85,0.25,0.15])
    loax.imshow(logo)
    loax.set_axis_off()
    loax= plt.axes([0.25,0.86,0.13,0.13])
    loax.imshow(logo_scope)
    loax.set_axis_off()
    loax= plt.axes([0.0,-0.08,0.26,0.26])
    loax.imshow(esa_logo)
    loax.set_axis_off()
    ax.text(-0.15,0.97,date_parser(time[i]).strftime('%b %Y'),transform=ax.transAxes,va='top',fontsize = 18)

    pad = 3*fps
    for j in range(len(data_add)):
        if (i >= data_add[j][0]) & (i <= data_add[j][1]+pad):
            if (i>= data_add[j][0]) & (i <=data_add[j][0]+pad):
                alpha = (i - data_add[j][0]) / pad
            elif (i>= data_add[j][1]) & (i <=data_add[j][1]+pad):
                alpha = 1-((i-data_add[j][1])/pad)
            else:
                alpha = 1
            ap = ax.pcolormesh(
                lon,
                lat,
                np.transpose(data_add[j][2][:,:,i]),
                #transform=proj,
                transform=ccrs.PlateCarree(),
                cmap=data_add[j][6],
                vmin=data_add[j][3],
                vmax=data_add[j][4],
                zorder=2,
                alpha=alpha
                )
            cax = plt.axes([0.72,0.85,0.25,0.025])
            cbar = plt.colorbar(ap,cax=cax,orientation='horizontal');cbar.set_label(data_add[j][5])


# plt.show()
fps = 30
locs = [(220,-30,0),#Equatorial Pacific
    (60,-60,15),#Southern Ocean
    (0,-90,25), # More Southern Ocean
    (0,-90,40), # Amazon Plume
    (0,90,70), # North Atlantic
    (0,90,90)
    ]

f = np.where(lon < 0)
flux_unc2 = np.copy(flux_unc)
flux_unc2[f,:,:] = np.nan

data_add = [[27,35,flux_unc2,0.1,0.5,'Chlorophyll-a uncertainty\n(log$_{10}$ (mg m$^{-3}$))',cmocean.cm.thermal],
    # [65,75,flux_unc3,0.1,0.5,'Chlorophyll-a uncertainty\n(log$_{10}$ (mg m$^{-3}$))',cmocean.cm.thermal],
    # [35,45,flux_unc4,0.1,0.5,'Chlorophyll-a uncertainty\n(log$_{10}$ (mg m$^{-3}$))',cmocean.cm.thermal]
    ]

for i in range(len(data_add)):
    for j in range(2):
        data_add[i][j]= data_add[i][j]*fps

for i in range(len(locs)-1):
    x = locs[i][0]
    y = locs[i][1]
    x2 = locs[i+1][0]
    y2 = locs[i+1][1]
    t1 = locs[i][2]
    t2 = locs[i+1][2]
    if i == 0:
        x_locs = np.linspace(x,x2,(t2-t1)*fps)
        y_locs = np.linspace(y,y2,(t2-t1)*fps)
    else:
        x_locs = np.concatenate((x_locs,np.linspace(x,x2,(t2-t1)*fps)))
        y_locs = np.concatenate((y_locs,np.linspace(y,y2,(t2-t1)*fps)))
    print(x_locs.shape)

ani_data = np.zeros((flux.shape[0],flux.shape[1],x_locs.shape[0]))
time_g2 = np.linspace(time_g[0],time_g[-1],x_locs.shape[0])

for i in range(flux.shape[0]):
    for j in range(flux.shape[1]):
        ani_data[i,j,:] = np.interp(time_g2,time_g,flux[i,j,:])

for g in range(len(data_add)):
    temp = np.zeros((data_add[g][2].shape[0],data_add[g][2].shape[1],x_locs.shape[0]))
    for i in range(data_add[g][2].shape[0]):
        for j in range(data_add[g][2].shape[1]):
            temp[i,j,:] = np.interp(time_g2,time_g,data_add[g][2][i,j,:])
    data_add[g][2] = temp

fig = plt.figure(figsize=(7,7))
ani = FuncAnimation(fig, animate, interval=200, blit=False, repeat=True,frames=x_locs.shape[0],fargs=(time_g2,x_locs,y_locs,ani_data,data_add,))
ani.save('animated.mp4', dpi=300, writer=FFMpegWriter(fps=fps))
# ani.save('animated.gif', dpi=300, writer=PillowWriter(fps=10))
# fig.savefig('airseaco2flux_gcb2024_UEXPFNNU_3D.png',dpi=300)
