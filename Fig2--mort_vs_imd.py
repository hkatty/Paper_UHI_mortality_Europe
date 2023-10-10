# -*- coding: utf-8 -*-
"""
Figure 2 in the main text. Plots of mortality risk as function of land 
imperviousness. Cities grouped by climate zones.

@author: Katty Huang
"""


import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from cycler import cycler
import glob
import cartopy.crs as ccrs
import matplotlib.ticker as mticker



def prep_data(datavar='imd', cityagewgt=True):
    """calculate bin averages for all cities and add group label"""
    
    def binaverage(var,data,bincount,datavar):
        binavg = data[var].groupby_bins(data[datavar],bins=databin).mean()
        # binavg = binavg.where(bincount>5)    
        return binavg
        
    #climate groups
    climate = pd.read_csv('data/city_climate_zones.csv',index_col=0)
    
    for i,group in enumerate(dict.fromkeys(climate.group)):
        for j,city in enumerate(climate.groupby('group').groups[group]):
            
            #make sure city file name uses '_' for both spaces and '-'
            cityfile = city.replace(' ','_').replace('-','_')
            
            data = xr.open_dataset('data/analysed/data_annual_seasavg_extremes_heatcold_counts_'+cityfile+'.nc')
    
            #2085.1 age group uses local age distribution. 2085.5 uses the standard age 
            age = 2085.1 if cityagewgt else 2085.5
            data = data.sel(age=age)
    
            if datavar=='imd':
                databin=np.insert(np.arange(0,105,5),0,-5)
            else:
                if data[datavar].max()<100:
                    print('Max population density difference < 100 for '+city)
                    continue
                databin=np.arange(0,data[datavar].max(),100)
            #for masking out bins containing 5 or less grids        
            bincount = data['fAD'].groupby_bins(data[datavar],bins=databin).count()
    
            binavg = binaverage('fAD',data,bincount,datavar)
            # binavg_h = binaverage('fADheat',data,bincount,datavar)
            # binavg_c = binaverage('fADcold',data,bincount,datavar)
            binavg_hex = binaverage('fAD_heat_ex',data,bincount,datavar)
            binavg_cex = binaverage('fAD_cold_ex',data,bincount,datavar)
            binavgT = binaverage('tas',data,bincount,datavar)
            binavg_hexT = binaverage('tas_heat_ex',data,bincount,datavar)
            binavg_cexT = binaverage('tas_cold_ex',data,bincount,datavar)
                
            if j == 0:
                binavg_all = binavg.assign_coords({'city':city})
                # binavg_h_all = binavg_h
                # binavg_c_all = binavg_c
                binavg_hex_all = binavg_hex.assign_coords({'city':city})
                binavg_cex_all = binavg_cex.assign_coords({'city':city})
                binavgT_all = binavgT.assign_coords({'city':city})
                binavg_hexT_all = binavg_hexT.assign_coords({'city':city})
                binavg_cexT_all = binavg_cexT.assign_coords({'city':city})
            else:
                binavg_all = xr.concat([binavg_all,binavg.assign_coords({'city':city})],dim='city')
                # binavg_h_all = xr.concat([binavg_h_all,binavg_h.assign_coords({'city':city})],dim='city')
                # binavg_c_all = xr.concat([binavg_c_all,binavg_c.assign_coords({'city':city})],dim='city')
                binavg_hex_all = xr.concat([binavg_hex_all,binavg_hex.assign_coords({'city':city})],dim='city')
                binavg_cex_all = xr.concat([binavg_cex_all,binavg_cex.assign_coords({'city':city})],dim='city')
                binavgT_all = xr.concat([binavgT_all,binavgT.assign_coords({'city':city})],dim='city')
                binavg_hexT_all = xr.concat([binavg_hexT_all,binavg_hexT.assign_coords({'city':city})],dim='city')
                binavg_cexT_all = xr.concat([binavg_cexT_all,binavg_cexT.assign_coords({'city':city})],dim='city')

        if i == 0:
            binavg_allgroups = binavg_all.assign_coords({'group':group})
            # binavg_allgroups_h = binavg_h_all.assign_coords({'group':group})
            # binavg_allgroups_c = binavg_c_all.assign_coords({'group':group})
            binavg_allgroups_hex = binavg_hex_all.assign_coords({'group':group})
            binavg_allgroups_cex = binavg_cex_all.assign_coords({'group':group})
            binavgT_allgroups = binavgT_all.assign_coords({'group':group})
            binavg_allgroups_hexT = binavg_hexT_all.assign_coords({'group':group})
            binavg_allgroups_cexT = binavg_cexT_all.assign_coords({'group':group})
        else:
            binavg_allgroups = xr.concat([binavg_allgroups,binavg_all.assign_coords({'group':group})],dim='group')
            # binavg_allgroups_h = xr.concat([binavg_allgroups_h,binavg_h_all.assign_coords({'group':group})],dim='group')
            # binavg_allgroups_c = xr.concat([binavg_allgroups_c,binavg_c_all.assign_coords({'group':group})],dim='group')
            binavg_allgroups_hex = xr.concat([binavg_allgroups_hex,binavg_hex_all.assign_coords({'group':group})],dim='group')
            binavg_allgroups_cex = xr.concat([binavg_allgroups_cex,binavg_cex_all.assign_coords({'group':group})],dim='group')
            binavgT_allgroups = xr.concat([binavgT_allgroups,binavgT_all.assign_coords({'group':group})],dim='group')
            binavg_allgroups_hexT = xr.concat([binavg_allgroups_hexT,binavg_hexT_all.assign_coords({'group':group})],dim='group')
            binavg_allgroups_cexT = xr.concat([binavg_allgroups_cexT,binavg_cexT_all.assign_coords({'group':group})],dim='group')
                    
    return xr.merge([binavg_allgroups, binavg_allgroups_hex, binavg_allgroups_cex, 
                     binavgT_allgroups, binavg_allgroups_hexT, binavg_allgroups_cexT])


def plot_imd(ax,binavg,title,iqr=True,ls='-',ylabel='$\Delta$ Mean attributable mortality\n (per 100,000 per day)',cost=True,annual=False,individual=False):
    """plot mortality vs imd subplot"""
        
    #move cities which doesn't reach deltaIMD of 80 to separate variable
    low_imd = (binavg.sel(imd_bins=binavg.imd_bins[-1]).isnull()
                & binavg.sel(imd_bins=binavg.imd_bins[-2]).isnull()
                & binavg.sel(imd_bins=binavg.imd_bins[-3]).isnull()
                & binavg.sel(imd_bins=binavg.imd_bins[-4]).isnull())
    not_null_combo = ~binavg.sel(imd_bins=binavg.imd_bins[1]).isnull()
    binavg_low_imd = binavg.where(low_imd & not_null_combo)
    #above contains 18 cities: binavg_low_imd.dropna(dim='city',how='all').city
    #1 arid, 13 cold no dry season, 1 temperatre no dry season hot summer, 
    #1 temperate no dry season warm summer, 2 temperate dry summer:
    #binavg_low_imd.dropna(dim='city',how='all').groupby('group').count(dim='city')
    binavg = binavg.where(~low_imd)
    #above contains 67 cities: binavg.dropna(dim='city',how='all').city
    #6 arid, 19 cold no dry season, 5 temperate no dry season hot summer, 
    #26 temperate no dry season warm summer, 11 temperate dry summer:
    #binavg.dropna(dim='city',how='all').groupby('group').count(dim='city')

    #mask out bins containing less than 4 cities        
    citycount = binavg.dropna(dim='city',how='all').groupby('group').count(dim='city')
    binavg = binavg.where(citycount>=4)
    
    if annual:
        binavg = binavg*365.2425
        
    y1 = binavg.groupby('group').quantile(0.25,dim='city')
    y2 = binavg.groupby('group').quantile(0.75,dim='city')  
     
    binavg.groupby('group').median(dim='city').plot(ax=ax,hue='group',ls=ls)
    # binavg_low_imd.groupby('group').median(dim='city').plot(hue='group',ls='--')
    
    if iqr:
        for g,gr in enumerate(binavg.group): 
            if individual:
                #plot all individual cities
                for city in binavg.sel(group=gr).city:
                    binavg.sel(group=gr,city=city).plot(ax=ax,color=colours[g])
            else:
                ax.fill_between(x=(np.insert(np.arange(0,105,5),0,-5)+2.5)[:-1], 
                                  y1=y1.sel(group=gr), y2=y2.sel(group=gr),
                                  alpha=0.3,label=gr.item().replace('_',' '))
    ax.set_title(title)
    ax.get_legend().remove()
    ax.set_ylabel(ylabel)
    ax.set_xlabel(None)
    #cost
    if cost:
        ax_cost = ax.secondary_yaxis('right',functions=(lambda x: x*39,lambda x: x/39))
        if annual:
            ax_cost.set_ylabel('Cost (€ per person per year)')
        else:
            ax_cost.set_ylabel('Cost (€ per person per day)')
    return binavg


def plot_climate_classification(subfig,binavg):
    """plot subfigure showing climate classification legend and map"""
    
    zones = pd.read_csv('data/city_climate_zones.csv',index_col=0)
    zones_dict = dict(zip(list(dict.fromkeys(zones.group)),range(5)))
    zones['group_number'] = zones.group.map(zones_dict)
    
    #get lat lon
    files = glob.glob('data/analysed/city_latlon_avgT_*.csv')
    latlons = pd.concat([pd.read_csv(f,index_col=0) for f in files])
    zones = zones.merge(latlons,left_index=True,right_index=True)
    
    #only select cities with deltaIMD>=80
    zones = zones[zones.index.isin(binavg.dropna(dim='city',how='all').city.values)]
    groupcount = zones.groupby('group').count().zone
    
    x=zones.lon
    y=zones.lat
    z=zones.group_number
    #change colour order to match climate groups
    cmaptab = matplotlib.cm.get_cmap('tab10')
    colours = [cmaptab(2),#green (warm summer)
               cmaptab(1),#orange (dry summer)
               cmaptab(4),#purple (Arid)
               cmaptab(0),#blue (Cold)
               cmaptab(3)]#red (hot summer)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('climategroups', colours, N=len(colours))
    
    subfig2 = subfig.subfigures(nrows=2)
    ax = subfig2[0].add_subplot(projection=ccrs.EuroPP())
    ax.set_extent([-15, 25, 30, 62],crs=ccrs.PlateCarree())#need to set to avoid issue with scatter plot coordinate
    ax.coastlines(color='grey',zorder=0)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, zorder=0)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([35, 45, 55])
    gl.xlocator = mticker.FixedLocator([-20, 0, 20])
    gl.xlabel_style = {'size': 'small', 'color': 'gray'}
    gl.ylabel_style = {'size': 'small', 'color': 'gray'}
    scatterplt = ax.scatter(x=x,y=y,c=z,s=15,alpha=1,transform=ccrs.PlateCarree(),cmap=cmap,vmin=z.min()-0.5,vmax=z.max()+0.5)
    ax.set_title('(d)')    
    ax2 = subfig2[1].add_subplot()
    ax2.axis('off')
    cbarax = subfig2[1].add_axes([ax2.get_position().x0-0.1,ax2.get_position().y0,0.05,ax2.get_position().height])
    cbar = subfig.colorbar(scatterplt,cax=cbarax,ticks=list(zones_dict.values()),shrink=0.6,aspect=10)
    cbar.ax.set_yticklabels([g+'\n('+str(groupcount[g])+' cities)' for g in zones_dict.keys()])
    cbar.ax.set_title('Köppen-Geiger climate classification',ha='left')


def output_plots(datavar='imd', cityagewgt=True):
    """
    Output two complete figures with each call. One showing IQR and another 
    showing individual cities in the mortality vs. imperviousness plots.
    All figures contain subplots of annual net, heat extreme, cold extreme, 
    and climate group classifications.
    
    """
    
    data = prep_data(datavar=datavar, cityagewgt=cityagewgt)
    
    # plot with interquartile range
    fig = plt.figure(constrained_layout=True,figsize=(9,8))
    subfigs = fig.subfigures(ncols=2,width_ratios=[2,1.3],wspace=0.07)
    axleft = subfigs[0].subplots(nrows=3,sharex=True)
    binavg=plot_imd(axleft[0],data.fAD,'(a) Annual',ylabel='$\Delta$ Mean attributable mortality\n (per 100,000 per year)',annual=True)
    plot_imd(axleft[1],data.fAD_heat_ex,'(b) Heat extreme')
    plot_imd(axleft[2],data.fAD_cold_ex,'(c) Cold extreme')
    subfigs[0].supxlabel('$\Delta$ Imperviousness (%)')
    plot_climate_classification(subfigs[1],binavg)
    agestr = 'cityage' if cityagewgt else 'stdage'    
    plt.savefig('figures/cost_vs_'+datavar+'_'+agestr+'.pdf')
    
    #plot with individual cities shown
    fig = plt.figure(constrained_layout=True,figsize=(9,8))
    subfigs = fig.subfigures(ncols=2,width_ratios=[2,1.3],wspace=0.07)
    axleft = subfigs[0].subplots(nrows=3,sharex=True)
    binavg=plot_imd(axleft[0],data.fAD,'(a) Annual',ylabel='$\Delta$ Mean attributable mortality\n (per 100,000 per year)',annual=True,individual=True)
    plot_imd(axleft[1],data.fAD_heat_ex,'(b) Heat extreme',individual=True)
    plot_imd(axleft[2],data.fAD_cold_ex,'(c) Cold extreme',individual=True)
    subfigs[0].supxlabel('$\Delta$ Imperviousness (%)')
    plot_climate_classification(subfigs[1],binavg)
    agestr = 'cityage' if cityagewgt else 'stdage'
    plt.savefig('figures/cost_vs_'+datavar+'_'+agestr+'_individual.pdf')
    

if __name__ == '__main__':

    #change colour cycle to match climate groups
    cmaptab = matplotlib.cm.get_cmap('tab10')
    colours = [cmaptab(2),#green
                cmaptab(1),#orange
                cmaptab(4),#purple
                cmaptab(0),#blue
                cmaptab(3)]#red
    custom_cycler = cycler(color=colours)
    plt.rc('axes', prop_cycle=custom_cycler)
    
    #plot figures
    output_plots(datavar='imd', cityagewgt=True)
    output_plots(datavar='imd', cityagewgt=False)
