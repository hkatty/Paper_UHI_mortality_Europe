# -*- coding: utf-8 -*-
"""
Temporal analysis of urban-rural daily mean timeseries data

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import rioxarray as rxr


#population density (world map) from NASA, clip to Europe
popden = rxr.open_rasterio("data/population/NASA_SEDAC/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif",masked=True) \
            .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0)          
popden = popden.rio.reproject('epsg:3035')

#cities list
cities_list = pd.read_csv('data/city_climate_zones.csv').city

for i,cityfile in enumerate(cities_list):
    
    #make sure city file name uses '_' for both spaces and '-'
    cityfile = cityfile.replace(' ','_').replace('-','_')
    #city name in Copernicus dataset (use spaces rather than '_' used in filename)
    city = (cityfile.replace('_','-') if cityfile=='Cluj_Napoca' else cityfile.replace('_',' '))    
    
    #timeseries of daily difference in urban-rural averages
    timeseries = xr.open_dataset('data/analysed/data_urbanruralavg_timeseries_'+cityfile+'.nc')
    
    #city domain mean temperature
    tasm = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_fldmean.nc').tas.squeeze(['lat','lon'])
          
    #minimum mortality temperature
    MMT = pd.read_csv('data/RRfit/MMT_urbclimT_'+cityfile+'.csv',index_col=0)
    MMT.index = [int(age.split('-')[0].split('+')[0]) for age in MMT.Age]
    MMT.index.rename('age',inplace=True)
    MMT_xr = xr.Dataset.from_dataframe(MMT) 
    
    #list of datasets to be merged for output
    outlist = []
    
    #annual average over whole time period
    timeseries_annual = timeseries.mean(dim='time')
    outlist.append(timeseries_annual)
    
    #seasonal average
    timeseries_seasavg = timeseries.groupby('time.season').mean()
    timeseries_seasavg = timeseries_seasavg.rename(dict([(v,(v+'_seasavg')) for v in timeseries_seasavg.keys()]))
    outlist.append(timeseries_seasavg)

    #count of total days warm/cold for each age group
    annual_counts = (~timeseries[['fADheat','fADcold']].isnull()).sum(dim='time').rename({'fADheat':'heat_count','fADcold':'cold_count'})
    outlist.append(annual_counts)    
    
    #count of days warm/cold per season per age group
    seas_counts = (~timeseries[['fADheat','fADcold']].isnull()).groupby('time.season').sum(dim='time').rename({'fADheat':'heat_seas_count','fADcold':'cold_seas_count'})
    outlist.append(seas_counts)    
    
    #temperature data
    tas = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc', decode_coords="all")
    tas = tas.tas
    #mask out grids with elevations more than 100m greater/less than the weighted average
    ele_mask = xr.open_dataarray('data/elevation_mask/elevation_mask_'+cityfile+'.nc')
    tas = tas.where(ele_mask)  
    #population density (per km2) 
    popden_city = popden.interp_like(tas)
    #city average daily T (population weighted & from period of consideration & masking out regions of elevation difference)
    avgT_daily = tas.weighted(popden_city.fillna(0)).mean(dim=['x','y'])
        
    #seasonal extreme
    #mask for seasonal warm/cold days (98th/2nd percentile threshold over three year period)
    heat_ex_mask = avgT_daily.groupby('time.season').map(lambda x: x>=x.quantile(0.98)).drop_vars('quantile')
    cold_ex_mask = avgT_daily.groupby('time.season').map(lambda x: x<=x.quantile(0.02)).drop_vars('quantile')
    #average of seasonal extreme days
    vars_ex = ['fAD']
    fADtimeseries_heat_ex_seasavg = timeseries[vars_ex].where(heat_ex_mask).groupby('time.season').mean()
    fADtimeseries_heat_ex_seasavg = fADtimeseries_heat_ex_seasavg.rename(dict([(v,(v+'_seas_heat_ex')) for v in fADtimeseries_heat_ex_seasavg.keys()]))
    fADtimeseries_cold_ex_seasavg = timeseries[vars_ex].where(cold_ex_mask).groupby('time.season').mean()
    fADtimeseries_cold_ex_seasavg = fADtimeseries_cold_ex_seasavg.rename(dict([(v,(v+'_seas_cold_ex')) for v in fADtimeseries_cold_ex_seasavg.keys()]))
    outlist.extend([fADtimeseries_heat_ex_seasavg,fADtimeseries_cold_ex_seasavg])
    #annual extreme
    #mask for warmest/coldest days over 2015-2017 time period (98th/2nd percentile threshold)
    heat_ex_mask = avgT_daily.where(lambda x: x>=x.quantile(0.98)).notnull().drop_vars('quantile')
    cold_ex_mask = avgT_daily.where(lambda x: x<=x.quantile(0.02)).notnull().drop_vars('quantile')
    #average of extreme days over time period
    vars_ex = ['fAD'] 
    fADtimeseries_heat_ex = timeseries[vars_ex].where(heat_ex_mask).mean(dim='time')
    fADtimeseries_heat_ex = fADtimeseries_heat_ex.rename(dict([(v,(v+'_heat_ex')) for v in fADtimeseries_heat_ex.keys()]))
    fADtimeseries_cold_ex = timeseries[vars_ex].where(cold_ex_mask).mean(dim='time')
    fADtimeseries_cold_ex = fADtimeseries_cold_ex.rename(dict([(v,(v+'_cold_ex')) for v in fADtimeseries_cold_ex.keys()]))
    outlist.extend([fADtimeseries_heat_ex,fADtimeseries_cold_ex])
        
    out = xr.merge(outlist)
    out.to_netcdf(path='data/analysed/data_urbanruralavg_timeseries_analysis_'+cityfile+'.nc')
                                     
    citydesc = pd.read_csv('data/analysed/city_latlon_avgT_'+cityfile+'.csv',index_col=0)
    lat = citydesc.lat.item()
    lon = citydesc.lon.item()
    out['lon'] = lon
    out['lat'] = lat
    
    if i == 0:
        dataall = out.assign_coords({'city':city})
    else:
        dataall = xr.concat([dataall,out.assign_coords({'city':city})],dim='city')

dataall.to_netcdf('data/analysed/avg_diff_from_rural_urbanrural.nc')