# -*- coding: utf-8 -*-
"""
Output time series of urban and rural daily averages, as well as their difference

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

#city information from Pierre (to extract lat-lon)
citydesc = pd.read_csv('data/citydesc_v3.csv',index_col=0,encoding='latin1')

#names matches between Copernicus and Pierre's data
epinames = pd.read_csv('data/copernicus_epiv3_cities_match.csv',index_col=0,encoding='latin1',keep_default_na=False)

#table of different geographical classification and codes
geotable = pd.read_excel('data/eurostat--EU-28-LAU-2019-NUTS-2016.xlsx',sheet_name='Combined')

#population density (world map) from NASA, clip to Europe
popden = rxr.open_rasterio("data/population/NASA_SEDAC/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif",masked=True) \
            .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0)          
popden = popden.rio.reproject('epsg:3035')

#age group weighting by European standard population 2013
#for age groups considered here (does not add up to 100% because missing younger age groups)
ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
age_wgt_std = ESP13/ESP13.sum()
age_wgt_std.index.rename('age',inplace=True)


#%%### city loop #####################

#cities from file names (use '_' instead of spaces)
cities_list = glob.glob(os.path.join('data','copernicus_urban','tas_*_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc'))
cities_list = [x.split(os.path.join('data','copernicus_urban','tas_'))[1].split('_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc')[0] for x in cities_list]

i=-1

for cityfile in cities_list:
    
    #make sure city file name uses '_' for both spaces and '-'
    cityfile = cityfile.replace(' ','_').replace('-','_')
    #city name in Copernicus dataset (use spaces rather than '_' used in filename)
    city = (cityfile.replace('_','-') if cityfile=='Cluj_Napoca' else cityfile.replace('_',' '))
    #city name used by Pierre
    city2 = epinames[epinames.copernicus==city].epi.values[0]
    if city2=='':
        print('epi model not found for '+city)
        continue

    lat = citydesc[citydesc.LABEL==city2]['lat'].iloc[0]
    lon = citydesc[citydesc.LABEL==city2]['lon'].iloc[0]
    i+=1
    
    tas = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc', decode_coords="all")
    #true where rural
    isrural = xr.open_dataset('data/copernicus_urban/urbanrural/ruralurbanmask_'+cityfile+'_UrbClim_v1.0.nc',decode_coords="all")
    #round by predominant land cover type
    isrural = isrural.fillna(0).interp_like(tas).ruralurbanmask.round()
    isrural = (isrural==1)
    #1 where land, nan otherwise
    landseamask = xr.open_dataset('data/copernicus_urban/landseamask/landseamask_'+cityfile+'_UrbClim_v1.0.nc',decode_coords="all")
    landseamask = landseamask.interp_like(tas).landseamask
    #minimum mortality temperature
    MMT = pd.read_csv('data/RRfit/MMT_urbclimT_'+cityfile+'.csv',index_col=0)
    MMT.index = [int(age.split('-')[0].split('+')[0]) for age in MMT.Age]
    MMT.index.rename('age',inplace=True)

    #get NUTS3 codes for city centre
    nuts3codes = get_nuts3_codes(lat, lon)
            
    #population density (per km2) 
    popden_city=popden.interp_like(tas)
        
    #mask out grids with elevations more than 100m greater/less than the weighted average
    ele_mask = xr.open_dataarray('data/elevation_mask/elevation_mask_'+cityfile+'.nc')
    tas = tas.where(ele_mask)
              
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
    MMT_xr = xr.Dataset.from_dataframe(MMT)
    
    #separate by heat and cold
    heat = tas[['tas','fAF']].where(tas.tas.mean(dim=['x','y'])-273.15>=MMT_xr.MMT).rename({'tas':'tasheat','fAF':'fAFheat'})
    cold = tas[['tas','fAF']].where(tas.tas.mean(dim=['x','y'])-273.15<MMT_xr.MMT).rename({'tas':'tascold','fAF':'fAFcold'})
    tas = xr.merge([tas,heat,cold])
    
    #weight age groups by local (NUTS3) population proportion
    #average aross all NUT3 codes corresponding to city
    age_wgt_local = poppropEU[poppropEU.index.get_level_values(0).isin(nuts3codes)].groupby('age').mean()
    
    
    #fAD per 100,000 per day
    #fAF*mort by age group
    fAD = tas[[x for x in tas.keys() if x.startswith('fAF')]]*mort_avg_xr.mort_avg
    fAD = fAD.rename(dict([(v,v.replace('AF','AD')) for v in fAD.keys()]))
    
    #combine all data (spatial maps)
    data = xr.merge([tas,fAD])
    data = data.where(~np.isnan(tas.tas.isel(time=0)))
    
    #population weighted mean across age groups
    #fill NAN with zeros to maintain weighting across age groups 
    #(if there's no risk to a certain age group, the overall risk would be lower as a result)
    #then set zeros back to NAN
    #weighted sum for attributable number
    tasall_local = data[[x for x in data.keys() if x.startswith('fAD')]].fillna(0).weighted(xr.DataArray(age_wgt_local)).sum(dim='age').expand_dims(dim={'age':[2085.1]})
    tasall_local = tasall_local.where(tasall_local != 0)
    tasall_std = data[[x for x in data.keys() if x.startswith('fAD')]].fillna(0).weighted(xr.DataArray(age_wgt_std)).sum(dim='age').expand_dims(dim={'age':[2085.5]})
    tasall_std = tasall_std.where(tasall_std != 0)
    #weighted mean for attributable fraction
    tasall_local_AF = data[[x for x in data.keys() if x.startswith('fAF')]].fillna(0).weighted(xr.DataArray(age_wgt_local)).mean(dim='age').expand_dims(dim={'age':[2085.1]})
    tasall_local_AF = tasall_local_AF.where(tasall_local_AF != 0)
    tasall_std_AF = data[[x for x in data.keys() if x.startswith('fAF')]].fillna(0).weighted(xr.DataArray(age_wgt_std)).mean(dim='age').expand_dims(dim={'age':[2085.5]})
    tasall_std_AF = tasall_std_AF.where(tasall_std_AF != 0)    
    data = xr.merge([tasall_local_AF,tasall_std_AF,tasall_local,tasall_std,data])


    #urban rural averages at each time step
    rural = data.where(isrural).mean(dim=['x','y'])#.rename(dict([(v,(v+'_rural')) for v in data.keys()]))
    urban = data.where(~isrural).mean(dim=['x','y'])#.rename(dict([(v,(v+'_urban')) for v in data.keys()]))
    diff = urban - rural
    out = xr.concat([rural.expand_dims(dim={'urbanrural':['rural']}),
                     urban.expand_dims(dim={'urbanrural':['urban']}),
                     diff.expand_dims(dim={'urbanrural':['diff']})],
                    dim='urbanrural')
        
    #combine all cities outputs
    if i == 0:
        data_all=out.assign_coords({'city':[city]})
    else:
        data_all=xr.concat([data_all,out.assign_coords({'city':[city]})],
                            dim='city')
        
data_all.to_netcdf(path='data/analysed/data_urbanruralavg_timeseries.nc')
                                 
