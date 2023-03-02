# -------------------------------------------------------------------------------
# modules
#
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import cmcrameri.cm as cmc
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import matplotlib

font = {'size': 12}
matplotlib.rc('font', **font)

def drywet(numcolors, colormap):

    colors_blue = colormap(np.linspace(0.5, 1, 5))
    colors_white = np.array([1, 1, 1, 1])
    colors_brown = [[84, 48, 5, 255],
                    [140, 81, 10, 255],
                    [191, 129, 45, 255],
                    [223, 194, 125, 255],
                    [246, 232, 195, 255]]
    rgb = []
    for i in range(len(colors_brown)):
        z = [x / 255 for x in colors_brown[i]]
        rgb.append(z)
    colors = np.vstack((rgb, colors_white, colors_blue))

    cmap = LinearSegmentedColormap.from_list(name=colormap, colors=colors, N=numcolors)

    return cmap

# -------------------------------------------------------------------------------
# read data
# %%
var_name = 'FR_SEA_ICE'

sims = ['new', 'old', 'diff']
friac = {}
labels = {'new': 'Sea Ice update', 'old': 'Sea Ice static', 'diff': 'Difference between versions'}

for s in range(len(sims)):
    sim = sims[s]
    friac[sim] = {}
    friac[sim]['label'] = labels[sim]
    data = xr.open_dataset(f'{sim}_version.nc')
    dt = data[var_name].values[0, :, :]
    friac[sim][var_name] = dt
# %%
lat = xr.open_dataset('old_version.nc')['lat'].values[:]
lon = xr.open_dataset('old_version.nc')['lon'].values[:]
lat_, lon_ = np.meshgrid(lon, lat)
print("load done")
# -------------------------------------------------------------------------------
# plot
# %%
ar = 1.0  # initial aspect ratio for first trial
wi = 12  # height in inches #15
hi = 2.5  # width in inches #10
ncol = 3  # edit here
nrow = 1
axs, cs, gl = np.empty(shape=(nrow, ncol), dtype='object'), np.empty(shape=(nrow, ncol), dtype='object'), np.empty(shape=(nrow, ncol), dtype='object')

cmap1 = cmc.davos_r
levels1 = np.linspace(0, 100, 21, endpoint=True)
norm1 = BoundaryNorm(levels1, ncolors=cmap1.N, clip=True)

cmap2 = drywet(25, cmc.vik_r)
levels2 = np.linspace(0, 40, 11, endpoint=True)
norm2 = BoundaryNorm(levels2, ncolors=cmap2.N, clip=True)

# change here the lat and lon
map_ext = [-50, 50, 40, 90]

fig = plt.figure(figsize=(wi, hi))
left, bottom, right, top = 0.07, 0.01, 0.94, 0.95
gs = gridspec.GridSpec(nrows=1, ncols=3, left=left, bottom=bottom, right=right, top=top,
                       wspace=0.1, hspace=0.15)

for i in range(3):
    sim = sims[i]
    label = friac[sim]['label']
    axs[0, i] = fig.add_subplot(gs[0, i], projection=ccrs.PlateCarree())
    axs[0, i].set_extent(map_ext, crs=ccrs.PlateCarree())
    axs[0, i].coastlines(zorder=3)
    axs[0, i].stock_img()
    gl[0, i] = axs[0, i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, linewidth=1, color='grey', alpha=0.5, linestyle='--')
    gl[0, i].right_labels = False
    gl[0, i].top_labels = False
    gl[0, i].left_labels = False
    cs[0, i] = axs[0, i].pcolormesh(lon, lat, friac[sim][var_name], cmap=cmap1, norm=norm1, shading="auto",
                                    transform=ccrs.PlateCarree())
    axs[0, i].set_title(f'{label}', fontweight='bold', pad=6, fontsize=14, loc='center')

gl[0, 0].left_labels = True

cax = fig.add_axes(
    [axs[0, 2].get_position().x1 + 0.01, axs[0, 2].get_position().y0, 0.01, axs[0, 2].get_position().height])
cbar = fig.colorbar(cs[0, 1], cax=cax, orientation='vertical',
                    ticks=np.linspace(0, 100, 6, endpoint=True))
cbar.ax.tick_params(labelsize=14)


axs[0, 0].text(-0.2, 0.5, 'Sea ice', ha='center', va='center', rotation='vertical',
               transform=axs[0, 0].transAxes, fontsize=14, fontweight='bold')
axs[0, 2].text(1.07, 1.09, '[%]', ha='center', va='center', rotation='horizontal',
               transform=axs[0, 2].transAxes, fontsize=12)

fig.show()
# plotpath = "/project/pr133/rxiang/figure/echam5/"
# fig.savefig(plotpath + 'friac' + f'{mon}.png', dpi=500)
plt.close(fig)
"""
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as mticker
import numpy as np
from netCDF4 import Dataset
import mplotutils as mpu
import xarray as xr
from matplotlib.colors import TwoSlopeNorm



def plot(filename, field, title, metric, lon, lat, max=None, min=None, 
    lat_pole=90, lon_pole=-180, coastline=True, colormap="RdBu_r", 
    centered_bar=False): 

    rotated_pole = ccrs.RotatedPole(pole_latitude=lat_pole, pole_longitude=lon_pole)
    data_crs = ccrs.PlateCarree()

    # create the plot and set the size
    plt.figure(figsize=(20,10))
    axes = plt.axes(projection= rotated_pole)


    # create country's borders and landsea mask
    #land_50m = cfeature.NaturalEarthFeature('cultural', 'admin_0_countries', '50m', edgecolor='black', facecolor='none', linewidth=0.2)
    broder_50m = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='black', facecolor='none', linewidth=0.8)
    if max is None:
        max = np.nanmax(field)
    if min is None:
        min = np.nanmin(field)
    # activate the labels and set countour of the countourf
    draw_labels = True
    reversed_cmap = False
    levels = np.arange(min, max, 0.05)
    color_map = plt.cm.get_cmap(colormap)
    if reversed_cmap:
        reversed_color_map = color_map.reversed()
    else:
        reversed_color_map = color_map
    plt.gca().set_facecolor("dimgrey")

    if centered_bar:
        norm = TwoSlopeNorm(vmin=min, vmax=max, vcenter=0.0001)
        # plot in each subplot
        h = plt.contourf(lon, lat, field[:,:], levels=levels, cmap=reversed_color_map, extend='both', norm=norm)
    else:
        h = plt.contourf(lon, lat, field[:,:], levels=levels, cmap=reversed_color_map , extend='both')
    axes.set_title(title, fontsize=25, weight="bold")

    ## add borders and landsea mask
    #axes.add_feature(land_50m)
    if coastline:
        axes.add_feature(broder_50m)
    
    gl = axes.gridlines(color='black', linestyle='--', linewidth=1., alpha=0.35, draw_labels=draw_labels, dms=True, x_inline=False, y_inline=False)
    gl.ylocator     = mticker.FixedLocator(np.arange(-60, 80, 10))
    gl.xlocator     = mticker.FixedLocator(np.arange(-100,  90 ,10))
    gl.xlabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.ylabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = True
    gl.left_labels = True

    #set colorbar
    cb = plt.colorbar(orientation="horizontal", shrink=0.4, pad=0.07, format="%.2f")
    cb.ax.tick_params(labelsize=14)
    cb.set_label(label= str(metric),fontsize=20)

    plt.tight_layout()
    plt.savefig(str(filename) + ".png")
    
# create var wherre to store others
result = np.zeros((1,224,544))


era5_old = xr.open_dataset("old_version.nc")
era5_new = xr.open_dataset("new_version.nc")
tas = xr.open_dataset("ts_delta.nc")
tos = xr.open_dataset("tos_delta.nc")


lon = tos.variables['lon'][:]
lat = tos.variables['lat'][:]
month = 0
print("loads done")
plot("old_final_temp", era5_old.variables['T_SKIN'][0,:,:].values, "Final PGW Temperature w/o sea ice update", 
     "T_SKIN [K]", lon, lat)
print("1 done")
plot("new_final_temp", era5_new.variables['T_SKIN'][0,:,:].values, "Final PGW Temperature with sea ice update", 
     "T_SKIN [K]", lon, lat)
print("2 done")
plot("old_final_sic", era5_old.variables['FR_SEA_ICE'][0,:,:].values, "ERA5 Sea Ice", 
     "Sea Ice frac [1]", lon, lat)
print("3 done")
plot("new_final_sic", era5_new.variables['FR_SEA_ICE'][0,:,:].values, "PGW Sea Ice", 
     "Sea Ice frac [1]", lon, lat)
print("4 done")

plot("diff_winter", result - tas.variables['tas'][month,:,:].values,"Differences between new and previous PGW versions for January", 
     "TAS delta [K]",lon,lat )
plot("tas_winter"+addon, tas.variables['tas'][month,:,:].values, "TAS field from previous PGW version for January", 
      "TAS delta [K]",lon,lat )

plot ("cdo_winter"+addon, cdo.variables['tos'][month,:,:].values, "SST field from bi-linear interpolation for January", 
      "SST delta [K]",lon,lat )
plot ("sst"+addon, tos.variables['sst'][month,:,:].values, "SST field from NaN-ignoring interpolation using kernel interp for January", 
      "SST delta [K]" ,lon,lat)

plot ("sst_tas_diff"+addon, tos.variables['sst'][month,:,:].values- tas.variables['tas'][month,:,:].values, "Differences between SST and TAS for January", 
      "SST delta [K]",lon,lat)
plot ("ice"+addon, era5.variables['FR_SEA_ICE'][0,:,:].values, "Sea ice fraction from ERA5 for January", 
      "Ice fraction [%]",lon,lat, colormap="Blues")
lon = christoph.variables['lon'][:]
lat = christoph.variables['lat'][:]
print(np.sum(christoph.variables['ts'][0,:,:].values))
print(np.sum(christoph.variables['ts'][0,:,:].values)/ (len(lon)*len(lat)))
plot ("heim", christoph.variables['ts'][0,:,:].values, "Differences between TS and SST for January", 
      "Temperature [K]",lon,lat)


def plot_paper(filename, field, title, metric, max=None, min=None, lat_pole=90, lon_pole=-180, 
         coastline=True, colormap="RdBu_r"): 

    rotated_pole = ccrs.RotatedPole(pole_latitude = lat_pole, pole_longitude = lon_pole)
    data_crs = ccrs.PlateCarree()
    
    # create the plot and set the size
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True , subplot_kw=dict(projection= rotated_pole), figsize = (20*3,10))
    fig.subplots_adjust( wspace=0.05, left=0.05, right=0.99, bottom=0.12, top=0.92)
    #fig.suptitle(r'$\Delta$'+ "SST comparison between native GCM data, bi-linear interpolation \n and NaN-ignoring interpolation", fontsize=30, weight="bold")

    # create country's borders and landsea mask
    #land_50m = cfeature.NaturalEarthFeature('cultural', 'admin_0_countries', '50m', edgecolor='black', facecolor='none', linewidth=0.2)
    broder_50m = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='black', facecolor='none', linewidth=0.8)
    if max == None:
        max = np.nanmax(field)
    if min == None:
        min = np.nanmin(field)
    # activate the labels and set countour of the countourf
    draw_labels = True
    reversed_cmap = False
    levels = np.arange(min, max, 0.05)
    color_map = plt.cm.get_cmap("Reds")
    if reversed_cmap == True:
        reversed_color_map = color_map.reversed()
    else:
        reversed_color_map = color_map
    #plt.gca().set_facecolor("dimgrey")
    #norm = TwoSlopeNorm(vmin=min, vmax = max, vcenter=0.0001)
    origin_dim = raw_sst.coords['longitude'].values.shape
    lon_raw = raw_sst.coords['longitude'].values.reshape(-1)
    for i in range(len(lon_raw)):
        if lon_raw[i] > 180:
            lon_raw[i] -= 360    
    # plot in each subplot
    h1 = axs[0].contourf(lon_raw.reshape(origin_dim), raw_sst.coords['latitude'].values, raw_sst.variables['tos'][month,:,:].values, levels=levels, cmap=reversed_color_map , extend='both')
    axs[0].set_title(r'$\Delta$SST on GCM ocean model grid', fontsize=25, weight="bold")
    ## add borders and landsea mask
    #axes.add_feature(land_50m)
    if coastline:
        axs[0].add_feature(broder_50m)
    
    gl = axs[0].gridlines(color='black', linestyle='--', linewidth=1., alpha=0.35, draw_labels=draw_labels, dms=True, x_inline=False, y_inline=False)
    gl.ylocator     = mticker.FixedLocator(np.arange(-60, 80, 10))
    gl.xlocator     = mticker.FixedLocator(np.arange(-100,  90 ,10))
    gl.xlabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.ylabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    
        # plot in each subplot
    h = axs[1].contourf(lon[210:-150], lat[30:-115], cdo.variables['tos'][month,:,:].values[30:-115,210:-150], levels=levels, cmap=reversed_color_map , extend='both')
    axs[1].set_title(r'$\Delta$'+ "SST using bi-linear interpolation", fontsize=25, weight="bold")
    ## add borders and landsea mask
    #axes.add_feature(land_50m)
    if coastline:
        axs[1].add_feature(broder_50m)
    
    gl = axs[1].gridlines(color='black', linestyle='--', linewidth=1., alpha=0.35, draw_labels=draw_labels, dms=True, x_inline=False, y_inline=False)
    gl.ylocator     = mticker.FixedLocator(np.arange(-60, 80, 10))
    gl.xlocator     = mticker.FixedLocator(np.arange(-100,  90 ,10))
    gl.xlabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.ylabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = True
    gl.left_labels = True

        # plot in each subplot
    h = axs[2].contourf(lon[210:-150], lat[30:-115], tos.variables['sst'][month,:,:].values[30:-115,210:-150], levels=levels, cmap=reversed_color_map , extend='both')
    axs[2].set_title(r'$\Delta$'+ "SST using NaN-ignoring interpolation", fontsize=25, weight="bold")
    ## add borders and landsea mask
    #axes.add_feature(land_50m)
    if coastline:
        axs[2].add_feature(broder_50m)
    
    gl = axs[2].gridlines(color='black', linestyle='--', linewidth=1., alpha=0.35, draw_labels=draw_labels, dms=True, x_inline=False, y_inline=False)
    gl.ylocator     = mticker.FixedLocator(np.arange(-60, 80, 10))
    gl.xlocator     = mticker.FixedLocator(np.arange(-100,  90 ,10))
    gl.xlabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.ylabel_style = {'size':12, 'rotation': 0, 'rotation_mode': 'anchor'}
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = True
    gl.left_labels = True   

    #set colorbar
    cb = mpu.colorbar(h1, axs[1], orientation = 'horizontal', pad = 0.15, aspect=50, format='%.1f') 
    cb.ax.tick_params(labelsize=14)
    cb.set_label(label="SST delta [K]",fontsize=20)

    #plt.tight_layout()
    plt.savefig(str(filename) + ".png")

#plot_paper("jonas_figure"+addon, raw_sst.variables['tos'][month,:,:].values, "Combined TAS field from NaN-ignoring interpolation for January", 
#     "TAS delta [K]")
"""
