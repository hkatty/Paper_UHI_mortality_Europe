# -*- coding: utf-8 -*-
"""
Regrid urban-rural masks to 500m and create elevation masks that 
mask out grids with elevations more than 100m greater/less than 
the weighted average for each city. Output to file.

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import rioxarray as rxr
from rioxarray.merge import merge_arrays
import numpy as np
from nuts_finder import NutsFinder  
import glob
import os


def prep_eurostat(file,sex=False):
    """
    read in and cleanup Eurostat dataset with file path as input
    outputs data as dataframe
    set sex to True if data is sex-separated to output only the total count across both sex
    """
    #read in data
    x = pd.read_csv(file,sep='[\t,]',engine='python',na_values=(':',': '))
    #strip white spaces at beginning/end of column names
    x.columns = x.columns.str.strip()
    #remove all notations for estimated/provisional values/values with breaks
    x.replace({' e':'',' b':'',' p':'',' bep':'',' ep':''},regex=True,inplace=True)
    #make sure all year columns are floats and not strings/objects
    years = np.arange(1990,2021).astype(str)
    for y in years:
        if y in x.columns:
            # x[y] = x[y].astype(float)     
            #convert to numeric and coerce errors to NaN
            x[y] = pd.to_numeric(x[y], errors='coerce')
    #return the total number across both sex only if data includes sex-separated counts
    return x[x.sex=='T'].drop(columns=['unit','sex']) if sex else x.drop(columns='unit')


def get_nuts3_codes(lat,lon): 
    """find NUTS 3 code for the city centre using nuts-finder package (uses point-in-shape method)"""
    try:
        nuts3code = NutsFinder(year=2016).find(lat,lon)[3]['NUTS_ID']
    except Exception:
        nuts3code = NutsFinder().find(lat,lon)[3]['NUTS_ID']
    
    #table of different geographical classification and codes
    geotable = pd.read_excel('data/eurostat--EU-28-LAU-2019-NUTS-2016.xlsx',sheet_name='Combined')

    #find all NUTS 3 codes that correspond to the same city/greater city
    nuts_cities = geotable[geotable['NUTS 3 CODE']==nuts3code].CITY_NAME.drop_duplicates().dropna()
    nuts_greater_cities = geotable[geotable['NUTS 3 CODE']==nuts3code].GREATER_CITY_NAME.drop_duplicates().dropna()
    nuts3codes = geotable[(geotable.CITY_NAME.isin(nuts_cities))|(geotable.GREATER_CITY_NAME.isin(nuts_greater_cities))]['NUTS 3 CODE'].drop_duplicates()
    #if not found in excel table, use the one from nuts-finder (may miss other codes that also cover the greater city area)
    if nuts3codes.empty:
        print('NUTS 3 code '+str(nuts3code)+' not found in excel table')
        nuts3codes = pd.Series(nuts3code)    
    return nuts3codes


def get_ele_grids(lat,lon):
    """
    Get elevation grid containing city coordinate and three neighbouring grids (north/south, west/east, diagonal)
    """
    #raster grid containing city coordinate and three neighbouring grids (north/south, west/east, and diagonal)
    latgrid = 5*round(lat/5) #closest, doesn't necessarily contain city coordinate
    latgrid2 = latgrid-5 if lat<latgrid else latgrid+5
    longrid = 5*round(lon/5)
    longrid2 = longrid-5 if lon<longrid else longrid+5
    elegrid = 'n'+str(int(latgrid))+('w' if longrid<0 else 'e')+f'{abs(longrid):03.0f}'
    elegrid2 = 'n'+str(int(latgrid2))+('w' if longrid<0 else 'e')+f'{abs(longrid):03.0f}'
    elegrid3 = 'n'+str(int(latgrid))+('w' if longrid2<0 else 'e')+f'{abs(longrid2):03.0f}'
    elegrid4 = 'n'+str(int(latgrid2))+('w' if longrid2<0 else 'e')+f'{abs(longrid2):03.0f}'
    #grid may not exist in some cases (full ocean grid)
    elelist = []
    for grid in [elegrid,elegrid2,elegrid3,elegrid4]:
        try:
            ele1 = rxr.open_rasterio('data/elevation/MERIT_DEM/'+grid+'_dem.tif').isel(band=0)
            elelist.append(ele1)
        except Exception:
            pass
    #if city central coordinate falls within 1.5 degrees of the border, also include grid on the other side
    if abs(lat-latgrid)<1.5:
        latgrid3 = latgrid-5 if lat>latgrid else latgrid+5
        elegrid5 = 'n'+str(int(latgrid3))+('w' if longrid<0 else 'e')+f'{abs(longrid):03.0f}'
        elegrid6 = 'n'+str(int(latgrid3))+('w' if longrid2<0 else 'e')+f'{abs(longrid2):03.0f}'
        for grid in [elegrid5,elegrid6]:
            try:
                ele1 = rxr.open_rasterio('data/elevation/MERIT_DEM/'+grid+'_dem.tif').isel(band=0)
                elelist.append(ele1)
            except Exception:
                pass
    if abs(lon-longrid)<1.5:
        longrid3 = longrid-5 if lon>longrid else longrid+5
        elegrid7 = 'n'+str(int(latgrid))+('w' if longrid3<0 else 'e')+f'{abs(longrid3):03.0f}'
        elegrid8 = 'n'+str(int(latgrid2))+('w' if longrid3<0 else 'e')+f'{abs(longrid3):03.0f}'
        for grid in [elegrid7,elegrid8]:
            try:
                ele1 = rxr.open_rasterio('data/elevation/MERIT_DEM/'+grid+'_dem.tif').isel(band=0)
                elelist.append(ele1)
            except Exception:
                pass
    if (abs(lon-longrid)<1.5)&(abs(lat-latgrid)<1.5):
        elegrid9 = 'n'+str(int(latgrid3))+('w' if longrid3<0 else 'e')+f'{abs(longrid3):03.0f}'
        elegrid10 = 'n'+str(int(latgrid3))+('w' if longrid3<0 else 'e')+f'{abs(longrid3):03.0f}'
        for grid in [elegrid9,elegrid10]:
            try:
                ele1 = rxr.open_rasterio('data/elevation/MERIT_DEM/'+grid+'_dem.tif').isel(band=0)
                elelist.append(ele1)
            except Exception:
                pass
    ele = merge_arrays(elelist)
    return ele
    
    
def get_ele_mask(cityfile,lat,lon,tas,popden_city):
    """
    Calculate and output a mask which excludes grids with elevations more than 100m
    greater/less than the population weighted average for the city 

    cityfile = name of city in Copernicus data files
    lat = latitude of city centre
    lon = longitude of city centre
    tas = temperature file, for interpolating to the same grid
    popden_city = population density map for the city
    """
    
    #1 where land, nan otherwise
    landseamask = xr.open_dataset('data/copernicus_urban/landseamask/landseamask_'+cityfile+'_UrbClim_v1.0.nc',decode_coords="all")
    #any grid that contains some water body is set to nan
    landseamask = landseamask.interp_like(tas).landseamask

    #elevation raster grids surrounding the city
    ele = get_ele_grids(lat,lon)
    
    #elevation
    ele_re = ele.rio.reproject('epsg:3035')
    ele_city = ele_re.interp_like(tas)*landseamask
    #population weighted mean, excluding NAN (e.g. non-land) grids
    ele_city_wgtd_mean=ele_city.weighted(popden_city.fillna(0)).mean()
    #mask to mask out grids with elevations more than 100m greater/less than the weighted average
    ele_mask = (abs(ele_city-ele_city_wgtd_mean)<=100)

    return ele_mask



if __name__ == '__main__':
    
    #city information from Pierre (to extract lat-lon)
    citydesc = pd.read_csv('data/citydesc_v3.csv',index_col=0,encoding='latin1')
    
    #names matches between Copernicus and Pierre's data
    epinames = pd.read_csv('data/copernicus_epiv3_cities_match.csv',index_col=0,encoding='latin1',keep_default_na=False)
    
    #population density (world map) from NASA, clip to Europe
    popden = rxr.open_rasterio("data/population/NASA_SEDAC/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif",masked=True) \
                .rio.clip_box(miny=35,maxy=72,minx=-25,maxx=65).isel(band=0)          
    popden = popden.rio.reproject('epsg:3035')
        
    #cities from file names (use '_' instead of spaces)
    cities_list = glob.glob(os.path.join('data','copernicus_urban','tas_*_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc'))
    cities_list = [x.split(os.path.join('data','copernicus_urban','tas_'))[1].split('_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc')[0] for x in cities_list]

    
    for cityfile in cities_list:
        
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
       
        #temperature at 500m resolution
        tas = xr.open_dataset('data/copernicus_urban/tas_'+cityfile+'_mergetime_daymean_mask_500m_2015to2017_fAF.nc', decode_coords="all")
        
        #population density (per km2) 
        popden_city=popden.interp_like(tas)
            
        #mask out grids with elevations more than 100m greater/less than the weighted average
        ele_mask = get_ele_mask(cityfile,lat,lon,tas,popden_city)
        #write elevation mask to file
        ele_mask.rename('elevation_and_land_mask').to_netcdf('data/elevation_mask/elevation_mask_'+cityfile+'.nc')
    
        #true where rural
        isrural = xr.open_dataset('data/copernicus_urban/urbanrural/ruralurbanmask_'+cityfile+'_UrbClim_v1.0.nc',decode_coords="all")
        #round by predominant land cover type
        isrural = isrural.fillna(0).interp_like(tas).ruralurbanmask.round()
        isrural = (isrural==1)
        #write rural mask to file
        isrural.to_netcdf('data/copernicus_urban/urbanrural/ruralurbanmask_'+cityfile+'_UrbClim_v1.0_500m.nc')    
