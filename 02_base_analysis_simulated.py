# -*- coding: utf-8 -*-
"""
Read in timeseries of daily urban-rural difference in attributable fraction,
multiply by mean mortality to get attributable number, and perform time 
selection and averaging (annual, seasonal, heat/cold, extremes)

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import rioxarray as rxr
import numpy as np
import calendar
import glob
import os
from prep_masks import prep_eurostat, get_nuts3_codes


#options:
#output annual averages
annual_out = True
#output seasonal averages
seas_out = True
#output extremes (annual/seasonal extremes, depending on annual_out/seas_out)
extremes_out = True
#output heat/cold variables separated by optimal T
heatcold_out = True
#output counts of heat/cold days (only relevant if heatcold_out==True)
counts_out = True 
if not heatcold_out:
    counts_out = False
    

#%%### Prep data ###########################

#Eurostat annual mortality by age and sex (NUTS3) - take total of both sex
mortEU = prep_eurostat('data/mortality/eurostat/demo_r_magec3.tsv',sex=True)
#add age group 85+
mortEU = mortEU.set_index(['geo\\time','age'])
mortEU = mortEU.unstack()
for yr in set(mortEU.columns.get_level_values(0)):
    mortEU.loc[:,(yr,'Y_GE85')]=mortEU[(yr,'Y_GE90')]+mortEU[(yr,'Y85-89')]   
mortEU = mortEU.stack()

#Eurostat population on 1 January by age group and sex (NUTS3)
popEU = prep_eurostat('data/population/eurostat/demo_r_pjangrp3.tsv',sex=True)
#add age group 85+ if not already available
popEU = popEU.set_index(['geo\\time','age'])
popEU = popEU.unstack()
for yr in set(popEU.columns.get_level_values(0)):
    popEU.loc[:,(yr,'Y_GE85')] = popEU[(yr,'Y_GE85')].fillna(popEU[(yr,'Y_GE90')]+popEU[(yr,'Y85-89')])
popEU = popEU.stack()

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

#population density (world map) from NASA, clip to Europe
popden = rxr.open_rasterio("data/population/NASA_SEDAC/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif",masked=True) \
            .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0)          
popden = popden.rio.reproject('epsg:3035')

#city information from Pierre (to extract lat-lon)
citydesc = pd.read_csv('data/citydesc_v3.csv',index_col=0,encoding='latin1')

#names matches between Copernicus and Pierre's data
epinames = pd.read_csv('data/copernicus_epiv3_cities_match.csv',index_col=0,encoding='latin1',keep_default_na=False)

#table of different geographical classification and codes
geotable = pd.read_excel('data/eurostat--EU-28-LAU-2019-NUTS-2016.xlsx',sheet_name='Combined')

#age group weighting by European standard population 2013
#for age groups considered here (does not add up to 100% because missing younger age groups)
ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
age_wgt_std = ESP13/ESP13.sum()
age_wgt_std.index.rename('age',inplace=True)


#%%### city loop ####################

#output netcdf file name label based on options at the top
options = [annual_out,seas_out,extremes_out,heatcold_out,counts_out]
out_labels = np.array(['_annual','_seasavg','_extremes','_heatcold','_counts'])
out_labels = ''.join(out_labels[options])

#cities from file names (use '_' instead of spaces)
cities_list = glob.glob(os.path.join('data','copernicus_urban','tas_*_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc'))
cities_list = [x.split(os.path.join('data','copernicus_urban','tas_'))[1].split('_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc')[0] for x in cities_list]

for cityfile in cities_list:
    
    print(cityfile)
    
    #make sure city file name uses '_' for both spaces and '-'
    cityfile = cityfile.replace(' ','_').replace('-','_')
    #city name in Copernicus dataset (use spaces rather than '_' used in filename)
    city = (cityfile.replace('_','-') if cityfile=='Cluj_Napoca' else cityfile.replace('_',' '))
    #city name used in epidemiology model
    city2 = epinames[epinames.copernicus==city].epi.values[0]
    if city2=='':
        print('epi model not found for '+city)
        continue

    lat = citydesc[citydesc.LABEL==city2]['lat'].iloc[0]
    lon = citydesc[citydesc.LABEL==city2]['lon'].iloc[0]

    # calculate average mortality for each age group ####################
    
    #get NUTS3 codes for city centre
    nuts3codes = get_nuts3_codes(lat, lon)

    #sum mortality and population separately across all sub-regions of city
    mortEUcity=mortEU[mortEU.index.get_level_values(0).isin(nuts3codes)].unstack().sum(axis=0)
    popEUcity=popEU[popEU.index.get_level_values(0).isin(nuts3codes)].unstack().sum(axis=0)
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
    #######################################################################
    
    #city domain mean temperature
    tasm = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_fldmean.nc').tas.squeeze(['lat','lon'])
    #urban-rural difference in attributable fraction
    fAFdiff = xr.open_mfdataset('data/simulated/simulated_fAFdiff_'+cityfile+'*.nc')#, decode_coords="all")
    if fAFdiff.simulation.size != 1000:
        print('simulation count '+str(fAFdiff.simulation.size)+' for '+cityfile)
          
    #minimum mortality temperature
    MMT = pd.read_csv('data/simulated/MMT/simulated_MMT_urbclimT_'+cityfile+'.csv',index_col=1).iloc[:,1:]
    MMT2 = pd.read_csv('data/simulated/MMT/simulated_MMT_urbclimT_'+cityfile+'_s501to1000.csv',index_col=1).iloc[:,1:]
    MMT2.index = MMT2.index + 500
    MMT = pd.concat([MMT,MMT2])
    MMT.columns = [int(age.split('-')[0].split('+')[0]) for age in MMT.columns]
    MMT.columns.rename('age',inplace=True)
    MMT_xr = MMT.stack().to_xarray().rename('MMT')
    
    if heatcold_out:
        #separate by heat and cold days
        heat = fAFdiff.where(tasm-273.15>=MMT_xr).rename({'fAFdiff':'fAFdiff_heat'})
        cold = fAFdiff.where(tasm-273.15<MMT_xr).rename({'fAFdiff':'fAFdiff_cold'})
        fAFdiff = xr.merge([fAFdiff,heat,cold])
    
    #weight age groups by local (NUTS3) population proportion
    #average aross all NUT3 codes corresponding to city
    age_wgt_local = poppropEU[poppropEU.index.get_level_values(0).isin(nuts3codes)].groupby('age').mean()
    
    #fAD per 100,000 per day
    #fAF*mort by age group
    fADdiff = fAFdiff*mort_avg_xr.mort_avg
    
    #combine age groups by population ratio
    #fill NAN with zeros to maintain weighting across age groups 
    #(if there's no risk to a certain age group, the overall risk would be lower as a result)
    fADdiffall_local = fADdiff.fillna(0).weighted(xr.DataArray(age_wgt_local)).sum(dim='age').expand_dims(dim={'age':[2085.1]})
    fADdiffall_local = fADdiffall_local.where(fADdiffall_local != 0)
    fADdiffall_std = fADdiff.fillna(0).weighted(xr.DataArray(age_wgt_std)).sum(dim='age').expand_dims(dim={'age':[2085.5]})
    fADdiffall_std = fADdiffall_std.where(fADdiffall_std != 0)
    fADdiff = xr.merge([fADdiffall_local,fADdiffall_std,fADdiff])    
    fADdiff = fADdiff.rename(dict([(v,v.replace('AF','AD')) for v in fADdiff.keys()]))

    
    #list of datasets to be merged for output
    outlist = []
    
    if annual_out:
        #annual average over whole time period
        fADdiff_annual = fADdiff.mean(dim='time')
        outlist.append(fADdiff_annual)
    
    if seas_out:
        #seasonal average
        fADdiff_seasavg = fADdiff.groupby('time.season').mean()
        fADdiff_seasavg = fADdiff_seasavg.rename(dict([(v,(v+'_seasavg')) for v in fADdiff_seasavg.keys()]))
        outlist.append(fADdiff_seasavg)

    if heatcold_out and counts_out:
        if annual_out:
            #count of total days warm/cold for each age group
            annual_counts = (~fADdiff[['fADdiff_heat','fADdiff_cold']].isnull()).sum(dim='time').rename({'fADdiff_heat':'heat_count','fADdiff_cold':'cold_count'})
            outlist.append(annual_counts)    
        if seas_out:
            #count of days warm/cold per season per age group
            seas_counts = (~fADdiff[['fADdiff_heat','fADdiff_cold']].isnull()).groupby('time.season').sum(dim='time').rename({'fADdiff_heat':'heat_seas_count','fADdiff_cold':'cold_seas_count'})
            outlist.append(seas_counts)    
    
        
    tas = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc', decode_coords="all")
    tas = tas.tas
    
    #population density (per km2) 
    popden_city=popden.interp_like(tas)
        
    #mask out grids with elevations more than 100m greater/less than the weighted average
    ele_mask = xr.open_dataarray('data/elevation_mask/elevation_mask_'+cityfile+'.nc')
    tas = tas.where(ele_mask)  
    # city average daily T (population weighted & from period of consideration & masking out regions of elevation difference)
    avgT_daily = tas.weighted(popden_city.fillna(0)).mean(dim=['x','y'])
        
    if extremes_out:
        if seas_out:
            #seasonal extreme
            #mask for seasonal warm/cold days (98th/2nd percentile threshold over three year period)
            heat_ex_mask = avgT_daily.groupby('time.season').map(lambda x: x>=x.quantile(0.98)).drop_vars('quantile')
            cold_ex_mask = avgT_daily.groupby('time.season').map(lambda x: x<=x.quantile(0.02)).drop_vars('quantile')
            #average of seasonal extreme days
            vars_ex = ['fADdiff']
            fADdiff_heat_ex_seasavg = fADdiff[vars_ex].where(heat_ex_mask).groupby('time.season').mean()
            fADdiff_heat_ex_seasavg = fADdiff_heat_ex_seasavg.rename(dict([(v,(v+'_seas_heat_ex')) for v in fADdiff_heat_ex_seasavg.keys()]))
            fADdiff_cold_ex_seasavg = fADdiff[vars_ex].where(cold_ex_mask).groupby('time.season').mean()
            fADdiff_cold_ex_seasavg = fADdiff_cold_ex_seasavg.rename(dict([(v,(v+'_seas_cold_ex')) for v in fADdiff_cold_ex_seasavg.keys()]))
            outlist.extend([fADdiff_heat_ex_seasavg,fADdiff_cold_ex_seasavg])
        if annual_out:
            #annual extreme
            #mask for warmest/coldest days over 2015-2017 time period (98th/2nd percentile threshold)
            heat_ex_mask = avgT_daily.where(lambda x: x>=x.quantile(0.98)).notnull().drop_vars('quantile')
            cold_ex_mask = avgT_daily.where(lambda x: x<=x.quantile(0.02)).notnull().drop_vars('quantile')
            #average of extreme days over time period
            vars_ex = ['fADdiff'] 
            fADdiff_heat_ex = fADdiff[vars_ex].where(heat_ex_mask).mean(dim='time')
            fADdiff_heat_ex = fADdiff_heat_ex.rename(dict([(v,(v+'_heat_ex')) for v in fADdiff_heat_ex.keys()]))
            fADdiff_cold_ex = fADdiff[vars_ex].where(cold_ex_mask).mean(dim='time')
            fADdiff_cold_ex = fADdiff_cold_ex.rename(dict([(v,(v+'_cold_ex')) for v in fADdiff_cold_ex.keys()]))
            outlist.extend([fADdiff_heat_ex,fADdiff_cold_ex])
        
    out = xr.merge(outlist)
    out.to_netcdf(path='data/analysed/simulated_data'+out_labels+'_'+cityfile+'.nc')
                                 
