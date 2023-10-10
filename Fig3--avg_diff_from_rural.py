# -*- coding: utf-8 -*-
"""
Figure 3 in the main text. Plots showing difference between urban and rural 
average. Boxplots and maps.

@author: Katty Huang
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

def prep_data():
    """
    Merge and cleanup data. Returns spatial-temporally aggregated values 
    for all cities, weighted by local population age structure.
    """
    
    cities_list = pd.read_csv('data/city_climate_zones.csv').city
    
    #add city lat/lon info and merge all city files
    for i,city in enumerate(cities_list):
        #make sure city file name uses '_' for both spaces and '-'
        cityfile = city.replace(' ','_').replace('-','_')
        
        citydata = xr.open_dataset('data/analysed/data_urbanruralavg_timeseries_analysis_'+cityfile+'.nc')
        
        citydesc = pd.read_csv('data/analysed/city_latlon_avgT_'+cityfile+'.csv')
        lat = citydesc.lat.item()
        lon = citydesc.lon.item()
        citydata['lon'] = lon
        citydata['lat'] = lat
        
        if i == 0:
            dataall = citydata.assign_coords({'city':city})
        else:
            dataall = xr.concat([dataall,citydata.assign_coords({'city':city})],dim='city')
    
    dataall = dataall.rename(dict([(v,v.replace('fADheat','fAD_heat')\
                                        .replace('fADcold','fAD_cold')) for v in dataall.keys()]))
        
    #local age weighting, urban-rural
    avgdiff_all_cityage = dataall.sel(age=2085.1,urbanrural='diff')
    avgdiff_all_cityage_rural = dataall.sel(age=2085.1,urbanrural='rural')
    
    #relative difference
    avgdiff_all_cityage['fADrel']=avgdiff_all_cityage.fAD/avgdiff_all_cityage_rural.fAD*100
    avgdiff_all_cityage['fADheatrel']=avgdiff_all_cityage.fAD_heat/avgdiff_all_cityage_rural.fAD_heat*100
    avgdiff_all_cityage['fADcoldrel']=avgdiff_all_cityage.fAD_cold/avgdiff_all_cityage_rural.fAD_cold*100
    avgdiff_all_cityage['fADheatexrel']=avgdiff_all_cityage.fAD_heat_ex/avgdiff_all_cityage_rural.fAD_heat_ex*100
    avgdiff_all_cityage['fADcoldexrel']=avgdiff_all_cityage.fAD_cold_ex/avgdiff_all_cityage_rural.fAD_cold_ex*100
    avgdiff_all_cityage['fADseasrel']=avgdiff_all_cityage.fAD_seasavg/avgdiff_all_cityage_rural.fAD_seasavg*100

    return avgdiff_all_cityage
    

def plot_map_scatter(x,y,z,label,title,ax=False,symmetric=False,cmap='viridis',vmax=False):
    """
    x=lons, y=lats, z=variable to plot as color
    label=colourbar label
    title=plot title
    """    
    import matplotlib.ticker as mticker
    if not ax:
        plt.figure()
        ax = plt.axes(projection=ccrs.EuroPP())
    if not vmax:
        vmax = abs(z).max()
    ax.set_extent([-15, 25, 30, 62],crs=ccrs.PlateCarree())#need to set to avoid issue with scatter plot coordinate
    ax.coastlines(color='grey',zorder=0)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, zorder=0)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([35, 45, 55])
    gl.xlocator = mticker.FixedLocator([-20, 0, 20])
    gl.xlabel_style = {'size': 'small', 'color': 'gray'}
    gl.ylabel_style = {'size': 'small', 'color': 'gray'}
    if symmetric:
        sc = ax.scatter(x=x,y=y,c=z,s=25,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmax=vmax,vmin=-vmax)
    else:
        sc = ax.scatter(x=x,y=y,c=z,s=25,transform=ccrs.PlateCarree(),cmap=cmap,vmax=vmax)
    ax.scatter(x=x,y=y,s=25,transform=ccrs.PlateCarree(),edgecolor='k',linewidth=0.7,facecolor='None')
    plt.colorbar(sc,label=label,ax=ax)
    ax.set_title(title)


def plot_box_maps(avgdiff_all_cityage): 
    """plot figure with boxplot and maps of net and extremes, local age structure"""
    
    df = avgdiff_all_cityage[['fAD','fAD_heat','fAD_cold','fAD_heat_ex','fAD_cold_ex']].to_dataframe().drop(columns=['age','crs','band','spatial_ref'])
    df_seas = avgdiff_all_cityage[['fAD_seasavg']].to_dataframe().drop(columns=['age','crs','band','spatial_ref'])
    
    df = df.merge(df_seas.unstack()['fAD_seasavg'],left_index=True,right_index=True)
    df = df.rename(columns={'fAD':'Annual \nnet','fAD_heat':'Heat','fAD_cold':'Cold',
               'fAD_heat_ex':'Heat \nextreme','fAD_cold_ex':'Cold \nextreme',
               'DJF':'Winter','MAM':'Spring','JJA':'Summer','SON':'Autumn'})
    df = df[['Annual \nnet','Heat \nextreme','Cold \nextreme','Winter','Spring','Summer','Autumn']]
    
    fig = plt.figure(figsize=(8,6))
    
    ax0 = fig.add_subplot(221)
    df.plot.box(ax=ax0,rot=70,flierprops={'marker':'.'})
    ax0.axhline(y=0,color='k',ls=':',zorder=0)
    ax0.set_title('(a) Average UHI impact')
    ax0.set_ylabel('$\Delta$ Deaths per 100,000 per day')
    ax0.text(0.65, 0.9, 'n = 85 cities', transform=ax0.transAxes)
    # #cost
    # ax_cost = ax0.secondary_yaxis('right',functions=(lambda x: x*39.1,lambda x: x/39.1))
    # ax_cost.set_ylabel('€ per person per day')
    
    ax1 = fig.add_subplot(222,projection=ccrs.EuroPP())
    plot_map_scatter(avgdiff_all_cityage.lon,avgdiff_all_cityage.lat,
                     avgdiff_all_cityage.fAD*365.24,'$\Delta$ Deaths per 100,000 per year',
                     '(b) Annual net',ax=ax1,symmetric=True)#,vmax=vmax)
    net_adverse = avgdiff_all_cityage[['lat','lon']].where(avgdiff_all_cityage.fAD>0).dropna(dim='city')
    ax1.scatter(x=net_adverse.lon,y=net_adverse.lat,s=25,transform=ccrs.PlateCarree(),edgecolor='r',linewidth=0.7,facecolor='None')
    
    
    vmax = abs(avgdiff_all_cityage[['fAD_heat_ex','fAD_cold_ex']]).max().to_array().max().item()
    
    ax2 = fig.add_subplot(223,projection=ccrs.EuroPP())
    plot_map_scatter(avgdiff_all_cityage.lon,avgdiff_all_cityage.lat,
                     avgdiff_all_cityage.fAD_heat_ex,'$\Delta$ Deaths per 100,000 per day',
                     '(c) Heat extreme days',ax=ax2,symmetric=True,vmax=vmax)
    
    ax3 = fig.add_subplot(224,projection=ccrs.EuroPP())
    plot_map_scatter(avgdiff_all_cityage.lon,avgdiff_all_cityage.lat,
                     avgdiff_all_cityage.fAD_cold_ex,'$\Delta$ Deaths per 100,000 per day',
                     '(d) Cold extreme days',ax=ax3,symmetric=True,vmax=vmax)
    
    plt.tight_layout()
    plt.savefig('figures/boxplot_netextrememap_cityage.pdf')


def plot_box_relative(avgdiff_all_cityage):
    """plot boxplot of % difference, local age structure"""
    dfrel_cityage = avgdiff_all_cityage[['fADrel','fADheatrel','fADcoldrel','fADheatexrel','fADcoldexrel']].to_dataframe().drop(columns=['age','crs','band','spatial_ref'])
    dfrel_cityage_seas = avgdiff_all_cityage[['fADseasrel']].to_dataframe().drop(columns=['age','crs','band','spatial_ref'])
    
    dfrel_cityage = dfrel_cityage.merge(dfrel_cityage_seas.unstack()['fADseasrel'],left_index=True,right_index=True)
    dfrel_cityage = dfrel_cityage.rename(columns={'fADrel':'Annual \nnet','fADheatrel':'Heat','fADcoldrel':'Cold',
               'fADheatexrel':'Heat \nextreme','fADcoldexrel':'Cold \nextreme',
               'DJF':'Winter','MAM':'Spring','JJA':'Summer','SON':'Autumn'})
    dfrel_cityage = dfrel_cityage[['Annual \nnet','Heat \nextreme','Cold \nextreme','Winter','Spring','Summer','Autumn']]
    dfrel_cityage.plot.box(rot=70,flierprops={'marker':'.'})
    plt.axhline(y=0,color='k',ls=':',zorder=0)
    plt.title('(urban - rural)/rural')
    plt.ylabel('%')
    plt.tight_layout()
    plt.savefig('figures/boxplot_percentchange_cityage.pdf')


def print_values_for_city(city,avgdiff_all_cityage):
    #temperature
    print('T')
    #rural
    print('rural: ',avgdiff_all_cityage.tas_rural_heat_ex.sel(city=city).values-273.15)
    #urban
    print('urban: ',(avgdiff_all_cityage.tas_rural_heat_ex.sel(city=city)+avgdiff_all_cityage.tas_heat_ex.sel(city=city)).values-273.15)
    #diff
    print('diff: ',avgdiff_all_cityage.tas_heat_ex.sel(city=city).values)
    #%diff
    print('% diff: ',(avgdiff_all_cityage.tas_heat_ex.sel(city=city)/avgdiff_all_cityage.tas_rural_heat_ex.sel(city=city)).values*100)
    #mortality
    print('mortality')
    #rural
    print('rural: ',avgdiff_all_cityage.fAD_rural_heat_ex.sel(city=city).values)
    #urban
    print('urban: ',(avgdiff_all_cityage.fAD_rural_heat_ex.sel(city=city)+avgdiff_all_cityage.fAD_heat_ex.sel(city=city)).values)
    #diff
    print('diff: ',avgdiff_all_cityage.fAD_heat_ex.sel(city=city).values)
    #%diff
    print('% diff: ',(avgdiff_all_cityage.fAD_heat_ex.sel(city=city)/avgdiff_all_cityage.fAD_rural_heat_ex.sel(city=city)).values*100)



if __name__ == '__main__':
    
    avgdiff_all_cityage = prep_data()
    #figure 3 in the main text
    plot_box_maps(avgdiff_all_cityage)
    #supplementary figure S10
    plot_box_relative(avgdiff_all_cityage)

    # values quoted in paper text
    print('annual mean net per year, median of all cities:')
    print(avgdiff_all_cityage.fAD.median().values*365.24)
    print('IQR:')
    print(avgdiff_all_cityage.fAD.quantile(0.25).values*365.24)
    print(avgdiff_all_cityage.fAD.quantile(0.75).values*365.24)
    
    print('number of cities with net protective effects:')
    print(avgdiff_all_cityage.fAD.where(avgdiff_all_cityage.fAD<0).count().values)
    
    print('heat extreme impact greater than cold extreme for each city:')
    print((avgdiff_all_cityage.fAD_heat_ex>avgdiff_all_cityage.fAD_cold_ex).all().values)
    
    print('daily mean during heat extreme, median of all cities:')
    print(avgdiff_all_cityage.fAD_heat_ex.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fAD_heat_ex.quantile(0.25).values)
    print(avgdiff_all_cityage.fAD_heat_ex.quantile(0.75).values)
    print('percent change, median:')
    print(avgdiff_all_cityage.fADheatexrel.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fADheatexrel.quantile(0.25).values)
    print(avgdiff_all_cityage.fADheatexrel.quantile(0.75).values)
    
    print('daily mean during cold extreme, median of all cities:')
    print(avgdiff_all_cityage.fAD_cold_ex.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fAD_cold_ex.quantile(0.25).values)
    print(avgdiff_all_cityage.fAD_cold_ex.quantile(0.75).values)
    print('percent change, median:')
    print(avgdiff_all_cityage.fADcoldexrel.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fADcoldexrel.quantile(0.25).values)
    print(avgdiff_all_cityage.fADcoldexrel.quantile(0.75).values)
    
