# -*- coding: utf-8 -*-
"""
Classify cities into Koppen-Geiger climate classes, with some similar classes grouped 
together where too few cities are classified in a class

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


#climate classification (Koppen-Geiger classes, version by Beck et al. 2018), clip to Europe
climate = rxr.open_rasterio("data/climate_zones/Beck_KG_V1_present_0p0083.tif",decode_coords="all") \
            .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0).astype('float').rio.write_nodata(np.nan)
climate = climate.where(climate!=0).rio.reproject('epsg:3035')

#cities from file names (use '_' instead of spaces)
cities = glob.glob(os.path.join('data','copernicus_urban','tas_*_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc'))
cities = [x.split(os.path.join('data','copernicus_urban','tas_'))[1].split('_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc')[0] for x in cities]

# #cities where the climate classification spatial resolution has an effect on the zone
# cities = ['Helsinki','Lisbon','Murcia','Nice','Zagreb']

out = []

for cityfile in cities:
    
    #city name in Copernicus dataset (use spaces rather than '_' used in filename)
    city = (cityfile.replace('_','-') if cityfile=='Cluj_Napoca' else cityfile.replace('_',' '))

    #true where urban
    isurban = xr.open_dataset('data/copernicus_urban/urbanrural/ruralurbanmask_'+cityfile+'_UrbClim_v1.0.nc',decode_coords="all")
    isurban = np.isnan(isurban.ruralurbanmask)
    
    #climate zone for the urban subset
    climate_city = climate.interp_like(isurban,method='nearest').where(isurban)
    
    # plt.figure()
    # climate_city.plot()
    
    #find the most common climate classification for city
    count = climate_city.groupby(climate_city).count()
    zone = count.sortby(count,ascending=False).group[0].item()

    out.append({'city':city,'zone':zone}) 
    
climate_zones = pd.DataFrame(out).set_index('city')

# #manually check classifications
# for z in set(climate_zones.zone):
#     print('zone: '+str(z))
#     print(climate_zones.index[climate_zones.zone==z].values)
    
#group similar climate zones
regroup = {4:'Arid',7:'Arid',8:'Temperate dry summer',9:'Temperate dry summer', 
           14:'Temperate no dry season hot summer',15:'Temperate no dry season warm summer', 
           25:'Cold no dry season',26:'Cold no dry season'}
climate_zones['group'] = climate_zones.zone.map(regroup)

#total number of cities within each group
climate_zones['group_total'] = climate_zones.groupby('group').group.transform('count')
#member number for each city within group
climate_zones['city_number'] = climate_zones.groupby('group').cumcount()

#split group into A and B subgroups if group contains more than 20 cities
climate_zones['subgroup'] = climate_zones['group']
climate_zones.loc[(climate_zones.group_total>20)&
                  (climate_zones.city_number<climate_zones.group_total/2),
                  'subgroup'] = climate_zones.subgroup+' A'
climate_zones.loc[(climate_zones.group_total>20)&
                  (climate_zones.city_number>=climate_zones.group_total/2),
                  'subgroup'] = climate_zones.subgroup+' B'

climate_zones.drop(columns=['group_total', 'city_number'],inplace=True)
climate_zones.to_csv('data/city_climate_zones.csv')


###### classification differs with spatial resolution for following cities

#high resolution (Beck_KG_V1_present_0p0083.tif)
# climate_zones.where(climate_zones.zone!=climate_zones_coarse.zone).dropna(how='all')
# Out[335]: 
#         city  zone
# 33  Helsinki  26.0
# 41    Lisbon   8.0
# 53    Murcia   4.0
# 56      Nice   8.0
# 85    Zagreb  14.0

#coarse: 0.5 degree resolution (Beck_KG_V1_present_0p5.tif)
# climate_zones_coarse.where(climate_zones.zone!=climate_zones_coarse.zone).dropna(how='all')
# Out[336]: 
#         city  zone
# 33  Helsinki  27.0
# 41    Lisbon   9.0
# 53    Murcia   7.0
# 56      Nice   9.0
# 85    Zagreb  26.0