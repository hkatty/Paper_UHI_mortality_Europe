# -*- coding: utf-8 -*-
"""
Output csv tables of mortality risk and cost 

@author: Katty Huang
"""

import xarray as xr
import pandas as pd

#age weighting option:
cityage = True

#urban-rural difference, average per 100,000 per day
uhi = xr.open_dataset('data/analysed/avg_diff_from_rural_urbanrural.nc').sel(urbanrural='diff')

if cityage:
    #local age distribution
    age_wgt = pd.read_csv('data/age_weights.csv',index_col=0)[1:].T
    age_wgt.columns.name = 'city'
    age_wgt.index = [float(age) for age in age_wgt.index]
    age_wgt.index.name = 'age'
    age_wgt_xr = age_wgt.stack().swaplevel().to_xarray()
else:
    #age group weighting by European standard population 2013
    #for age groups considered here (does not add up to 100% because missing younger age groups)
    ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
    age_wgt_std = ESP13/ESP13.sum()
    age_wgt_std.index.rename('age',inplace=True)
    age_wgt_xr = age_wgt_std.to_xarray()
uhi = uhi.assign(age_wgt=age_wgt_xr)

#per year for net, heat, cold; per day for extremes
uhi['fAD'] = (uhi.fAD).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')*365.24
uhi['fADheat'] = (uhi.fADheat*uhi.heat_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
uhi['fADcold'] = (uhi.fADcold*uhi.cold_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
uhi['fAD_heat_ex'] = uhi.fAD_heat_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
uhi['fAD_cold_ex'] = uhi.fAD_cold_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')

outlist = ['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']
outdict = {'fAD':'Net (per year)','fADheat':'Heat (per year)','fADcold':'Cold (per year)','fAD_heat_ex':'Heat extreme (per day)','fAD_cold_ex':'Cold extreme (per day)'}
#sort by net, descending
mort_table = uhi[outlist].sortby('fAD',ascending=False)
mort_table = mort_table.rename(outdict)
mort_table = mort_table.to_dataframe()[[outdict[v] for v in outlist]]
mort_table.to_csv('data/analysed/uhi_mortality_table_cityage.csv' if cityage else
                  'data/analysed/uhi_mortality_table_stdage.csv')

# output net, heat, cold cost to csv (VSL, local population)   
cost_table = mort_table / 100000 * 3910000
cost_table.to_csv('data/analysed/uhi_mortality_cost_table_cityage.csv' if cityage else
                  'data/analysed/uhi_mortality_cost_table_stdage.csv')
