# -*- coding: utf-8 -*-
"""
Figure 1 in the main text. Example cities.

@author: Katty Huang
"""

import pandas as pd
import numpy as np
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
import glob
import calendar
from prep_masks import prep_eurostat, get_nuts3_codes

def extremes_mean(cityfile,vars):
    """values of vars averaged over the warmest/coldest 2% days"""

    #population density (world map) from NASA, clip to Europe
    popden = rxr.open_rasterio("data/population/NASA_SEDAC/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif",masked=True) \
            .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0)          
    popden = popden.rio.reproject('epsg:3035')

    #get temperature 
    data = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc', decode_coords="all")
    tas = data.tas
    #mask out grids with elevations more than 100m greater/less than the weighted average
    ele_mask = xr.open_dataarray('data/elevation_mask/elevation_mask_'+cityfile+'.nc')
    tas = tas.where(ele_mask)  
    #population density (per km2) 
    popden_city = popden.interp_like(tas)
    # city average daily T (population weighted & from period of consideration & masking out regions of elevation difference)
    avgT_daily = tas.weighted(popden_city.fillna(0)).mean(dim=['x','y'])
    
    #annual extreme
    #mask for warmest/coldest days over 2015-2017 time period (98th/2nd percentile threshold)
    heat_ex_mask = avgT_daily.where(lambda x: x>=x.quantile(0.98)).notnull().drop_vars('quantile')
    cold_ex_mask = avgT_daily.where(lambda x: x<=x.quantile(0.02)).notnull().drop_vars('quantile')   

    heatex = data[vars].where(heat_ex_mask).mean(dim='time')
    coldex = data[vars].where(cold_ex_mask).mean(dim='time')
    return heatex, coldex

def get_avg_mort(cityfile):
    """output average mortality per day per 100,000 for each age group, over time period 2015-2017"""

    citydesc = pd.read_csv('data/analysed/city_latlon_avgT_'+cityfile+'.csv')
    lat = citydesc.lat
    lon = citydesc.lon

    #get NUTS3 codes for city centre
    nuts3codes = get_nuts3_codes(lat, lon)

    #Eurostat annual mortality by age and sex (NUTS3) - take total of both sex
    mortEU = prep_eurostat('data/mortality/eurostat/demo_r_magec3.tsv',sex=True)
    #add age group 85+
    mortEU = mortEU.set_index(['geo\\time','age'])
    mortEU = mortEU.unstack()
    for yr in set(mortEU.columns.get_level_values(0)):
        mortEU.loc[:,(yr,'Y_GE85')] = mortEU[(yr,'Y_GE90')]+mortEU[(yr,'Y85-89')]   
    mortEU = mortEU.stack()

    #Eurostat population on 1 January by age group and sex (NUTS3)
    popEU = prep_eurostat('data/population/eurostat/demo_r_pjangrp3.tsv',sex=True)
    #add age group 85+ if not already available
    popEU = popEU.set_index(['geo\\time','age'])
    popEU = popEU.unstack()
    for yr in set(popEU.columns.get_level_values(0)):
        popEU.loc[:,(yr,'Y_GE85')] = popEU[(yr,'Y_GE85')].fillna(popEU[(yr,'Y_GE90')]+popEU[(yr,'Y85-89')])
    popEU = popEU.stack()
    
    #sum mortality and population separately across all sub-regions of city
    mortEUcity = mortEU[mortEU.index.get_level_values(0).isin(nuts3codes)].unstack().sum(axis=0)
    popEUcity = popEU[popEU.index.get_level_values(0).isin(nuts3codes)].unstack().sum(axis=0)
    #annual mortality per 100,000 (in age group)
    mort_per100000 = (mortEUcity/popEUcity*100000).replace(0,np.nan).dropna().unstack(level=0)
    
    #average mortality per day per 100,000 for each age group, over time period (years within where mort data available)
    mort_avg = (mort_per100000/[365+calendar.isleap(y) for y in mort_per100000.columns.astype(int)])[mort_per100000.columns.intersection(['2015','2016','2017'])].mean(axis=1)
    #select and rename age groups for consistency
    mort_avg['20'] = mort_avg[['Y20-24','Y25-29','Y30-34','Y35-39','Y40-44']].mean()
    mort_avg['45'] = mort_avg[['Y45-49','Y50-54','Y55-59','Y60-64']].mean()    
    mort_avg['65'] = mort_avg[['Y65-69','Y70-74']].mean() 
    mort_avg['75'] = mort_avg[['Y75-79','Y80-84']].mean() 
    mort_avg['85'] = mort_avg['Y_GE85']
    
    mort_avg_xr = xr.Dataset(data_vars=dict(mort_avg=(["age"], mort_avg[['20','45','65','75','85']])),
                              coords=dict(age=(["age"], [20,45,65,75,85])))
    return mort_avg_xr

def get_popprop(city):
    """get population age weighting"""
    #Eurostat annual population structure by age (NUTS3) -- % of each age group
    poppropEU = prep_eurostat('data/population/eurostat/demo_r_pjanind3.tsv',sex=False)
    #select only info on proportion of population in 5-year age groups (also contain, but not selected here, median age of population "MEDAGEPOP")
    indic_de_list = ["PC_Y0_4", "PC_Y5_9", "PC_Y10_14", "PC_Y15_19", 
      "PC_Y20_24", "PC_Y25_29", "PC_Y30_34", "PC_Y35_39", "PC_Y40_44", 
      "PC_Y45_49", "PC_Y50_54", "PC_Y55_59", "PC_Y60_64", "PC_Y65_69", 
      "PC_Y70_74", "PC_Y75_79", "PC_Y80_84", "PC_Y85_MAX"]
    poppropEU = poppropEU[poppropEU.indic_de.isin(indic_de_list)]
    #map to age groups of interest
    agegroup = {"PC_Y20_24":20, "PC_Y25_29":20, "PC_Y30_34":20, "PC_Y35_39":20, "PC_Y40_44":20, 
      "PC_Y45_49":45, "PC_Y50_54":45, "PC_Y55_59":45, "PC_Y60_64":45, "PC_Y65_69":65, 
      "PC_Y70_74":65, "PC_Y75_79":75, "PC_Y80_84":75, "PC_Y85_MAX":85}
    poppropEU['age']  = poppropEU.indic_de.map(agegroup)
    #drop younger age groups (<20)
    poppropEU = poppropEU.dropna(subset=['age'])
    #age group sum
    poppropEU = poppropEU.groupby(['geo\\time','age']).sum()
    #rescale to proportion of population considered (20+)
    poppropEU = poppropEU.groupby('geo\\time').apply(lambda x: x/x.sum())
    #average across years of interest
    poppropEU = poppropEU[['2015','2016','2017']].mean(axis=1)
    return poppropEU


def plot_map_with_mask(ax,city,title):
    """plot map with rural hatched"""
    heatex, coldex = extremes_mean(city,'tas')
    tasheatex = heatex - 273.15
    
    #true where rural, rounded by predominant land cover type
    isrural = xr.open_dataarray('data/copernicus_urban/urbanrural/ruralurbanmask_'+city+'_UrbClim_v1.0_500m.nc',decode_coords="all")
    isrural = isrural.where(isrural)
    
    middleT=(tasheatex.max()-tasheatex.min())/2+tasheatex.min()
    vmin=middleT-2.2
    vmax=middleT+2.2
    tasheatex.plot(ax=ax,cmap='RdBu_r',vmin=vmin,vmax=vmax,cbar_kwargs={'label':'Heat extreme T ($^\circ$C)'})
    ax.pcolor(isrural.x,isrural.y,isrural,alpha=0,hatch='//',shading='nearest')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title(title)


def plot_hist(ax,city,title):
    """plot hitogram"""
    heatex, coldex = extremes_mean(city,'fAF')
    mort_avg_xr = get_avg_mort(city)
    poppropEU = get_popprop(city)

    #true where rural, rounded by predominant land cover type
    isrural = xr.open_dataarray('data/copernicus_urban/urbanrural/ruralurbanmask_'+city+'_UrbClim_v1.0_500m.nc',decode_coords="all")
    
    #fAD per 100,000 per day = fAF * mort_avg for each age group
    fAD_heatex = heatex * mort_avg_xr.mort_avg
    
    citydesc = pd.read_csv('data/analysed/city_latlon_avgT_'+city+'.csv')
    lat = citydesc.lat
    lon = citydesc.lon
    #get NUTS3 codes for city centre
    nuts3codes = get_nuts3_codes(lat, lon)
    
    #combine age groups by local age weight
    #fill NAN with zeros to maintain weighting across age groups 
    #(if there's no risk to a certain age group, the overall risk would be lower as a result)
    age_wgt_local = poppropEU[poppropEU.index.get_level_values(0).isin(nuts3codes)].groupby('age').mean()
    fAD_heatex_local = fAD_heatex.fillna(0).weighted(xr.DataArray(age_wgt_local)).sum(dim='age').expand_dims(dim={'age':[2085.1]})
    fAD_heatex_local = fAD_heatex_local.where(fAD_heatex_local != 0)

    urban = fAD_heatex_local.where(~isrural)
    rural = fAD_heatex_local.where(isrural)

    bins = np.linspace(fAD_heatex_local.min(),fAD_heatex_local.max(),20)
    
    ##density
    # rural.plot.hist(ax=ax,bins=bins,label='rural',density=True)
    # urban.plot.hist(ax=ax,bins=bins,alpha=0.6,label='urban',density=True)
    # ax.set_ylabel('Density')
    
    #relative frequency (% of cases)
    def prep_for_hist(darray):
        #get data into shape used by np.hist(), replicating xarray wrapper (1D and w/o NaN)
        no_nan = np.ravel(darray.values)
        no_nan = no_nan[pd.notnull(no_nan)]
        return no_nan
    def rel_freq_weight(data):
        """ histogram weighting to output relative frequency """
        #weight calculation, in shape of xarray data (2D)
        darray = xr.ones_like(data).where(~data.isnull())/data.count()
        #get data into shape used by np.hist(), replicating xarray wrapper (1D and w/o NaN)
        return prep_for_hist(darray)
    ax.hist(prep_for_hist(rural),bins=bins,label='rural',weights=rel_freq_weight(rural)*100)
    ax.hist(prep_for_hist(urban),bins=bins,alpha=0.6,label='urban',weights=rel_freq_weight(urban)*100)
    ax.set_ylabel('Relative frequency (%)')
    
    ax.axvline(x=rural.mean(),c='k',ls=':')
    ax.axvline(x=urban.mean(),c='k',ls='--')
    ax.set_title(title)


def set_share_axes(axs, target=None, sharex=False, sharey=False):
    """set axes for subplots"""
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_x_axes.join(target, ax)
        if sharey:
            target._shared_y_axes.join(target, ax)
    # Turn off x tick labels and offset text for all but the bottom row
    if sharex and axs.ndim > 1:
        for ax in axs[:-1,:].flat:
            ax.xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            ax.xaxis.offsetText.set_visible(False)
    # Turn off y tick labels and offset text for all but the left most column
    if sharey and axs.ndim > 1:
        for ax in axs[:,1:].flat:
            ax.yaxis.set_tick_params(which='both', labelleft=False, labelright=False)
            ax.yaxis.offsetText.set_visible(False)


def ERF_Tlines(ax,city,label):
    """plot exposure-response relationships"""
    
    #RR vs T relationship for all age groups into one dataframe
    RRfiles = glob.glob('data/from_Pierre/RRfit_zenodo/RRfit_urbclimT_'+city+'_[!p]*.csv')
    RRall = pd.DataFrame()
    for RRfile in RRfiles:
        agegroup = RRfile.split(city+'_')[1].split('.csv')[0]
        RR = pd.read_csv(RRfile,index_col=0)
        RR.columns = pd.MultiIndex.from_product([RR.columns,[agegroup]])
        RRall = pd.concat([RRall,RR],axis=1)
    
    #average T over warmest/coldest 2% days
    tasheatex, tascoldex = extremes_mean(city,'tas')
    dataT = xr.merge([tasheatex.rename('heat_ex'),tascoldex.rename('cold_ex')]) - 273.15

    #true where rural, rounded by predominant land cover type
    isrural = xr.open_dataarray('data/copernicus_urban/urbanrural/ruralurbanmask_'+city+'_UrbClim_v1.0_500m.nc',decode_coords="all")
    
    urbanT = dataT.where(~isrural).mean(dim=['x','y'])
    ruralT = dataT.where(isrural).mean(dim=['x','y'])
    
    #exposure-response plot
    RRall.allfit.plot(ax=ax,legend=False)
    ax.axvline(x=urbanT.heat_ex,color='k',ls='--',label='urban')
    ax.axvline(x=urbanT.cold_ex,c='grey',ls='--')
    ax.axvline(x=ruralT.heat_ex,color='k',ls=':',label='rural')
    ax.axvline(x=ruralT.cold_ex,c='grey',ls=':')
    ax.set_ylabel('Relative risk')
    ax.set_title(label+' '+city)
    ax.set_ylim([0.7,3.5])


if __name__ == '__main__':

    fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(9,7))
    set_share_axes(ax[:,0], sharex=True)
    set_share_axes(ax[:,2], sharex=True)
    
    ERF_Tlines(ax[0,0],'Milan','(a)')
    ERF_Tlines(ax[1,0],'London','(b)')
    ERF_Tlines(ax[2,0],'Stockholm','(c)')
    
    ax[0,0].legend(fontsize='small',ncol=2,loc='upper left')
    ax[2,0].set_xlabel('Daily mean temperature ($^\circ$C)')
    plot_map_with_mask(ax[0,1],'Milan','(d)')
    plot_map_with_mask(ax[1,1],'London','(e)')
    plot_map_with_mask(ax[2,1],'Stockholm','(f)')
    plot_hist(ax[0,2],'Milan','(g)')
    plot_hist(ax[1,2],'London','(h)')
    plot_hist(ax[2,2],'Stockholm','(i)')
    ax[0,2].legend(fontsize='small')
    ax[2,2].set_xlabel('Deaths per day per 100,000')
    plt.tight_layout()
    plt.savefig('figures/example_cities.pdf')
