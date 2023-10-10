#!/usr/bin/env bash

cd data/copernicus_urban/

#echo "City:"
#read city

#loop through all cities in raw
cities=$(ls raw/tas_*_UrbClim_2008_01_v1.0.nc | sed 's/^.*tas_\(.*\)_UrbClim.*/\1/')
##specify list of cities
#cities=(Athens Antwerp)

for city in $cities

do

#merge all time steps, calculate daily mean, mask out non-land
cdo mul landseamask/landseamask_${city^}_UrbClim_v1.0.nc -daymean -mergetime raw/tas_${city^}_UrbClim_20*v1.0.nc tas_${city^}_mergetime_daymean_mask.nc

#calculate daily city domain mean for constructing model basis
cdo fldmean -setctomiss,nan tas_${city^}_mergetime_daymean_mask.nc tas_${city^}_mergetime_daymean_mask_fldmean.nc

#get original grid description
xinc=$(cdo griddes tas_${city^}_mergetime_daymean_mask.nc | grep -Po 'xinc\s*= \K.*')
xsize=$(cdo griddes tas_${city^}_mergetime_daymean_mask.nc | grep -Pom 1 'xsize\s*= \K.*')
#make copy of original grid description, remove unnecessary lines, and modify xinc/yinc and xsize/ysize
cdo griddes tas_${city^}_mergetime_daymean_mask.nc | sed -n '/gridtype  = projection/,$p' | sed '/gridsize/d;/datatype/d;/scanningMode/d;/grid_mapping_name/d;/name/d' | sed "s/inc.*=$PARTITION_COLUMN.*/inc = $((${xinc}*5))/" | sed "s/size.*=$PARTITION_COLUMN.*/size = $((${xsize}/5))/" > ${city,,}grid_500m.txt
#add additional information
echo "grid_mapping = crs
grid_mapping_name = \"EPSG:3035 (ETRS89 / LAEA Europe)\" 
proj_params = \"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs\"" >> ${city,,}grid_500m.txt

#for attribution
#regrid to 500m
cdo remapbil,${city,,}grid_500m.txt tas_${city^}_mergetime_daymean_mask.nc tas_${city^}_mergetime_daymean_mask_500m.nc
#select last 3 years (2015-2017)
cdo selyear,2015/2017 tas_${city^}_mergetime_daymean_mask_500m.nc tas_${city^}_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc

done #loop through cities
