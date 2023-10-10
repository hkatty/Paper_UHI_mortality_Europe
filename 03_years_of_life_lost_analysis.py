# -*- coding: utf-8 -*-
"""
Years of life lost computations

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import numpy as np
from nuts_finder import NutsFinder
from prep_masks import prep_eurostat


def avg_yrs_lost(avgdiff_all):
    """output average years lost by age group and city to csv file"""
    
    #Eurostat life expectancy by age and sex (NUTS2), life expectancy is the number of years remaining
    lexpEU = prep_eurostat('data/life_expectancy/demo_r_mlifexp.tsv',sex=True)
    #define age groups
    agegroup = dict([('Y'+str(n),20) for n in np.arange(20,45)]+\
               [('Y'+str(n),45) for n in np.arange(45,65)]+\
               [('Y'+str(n),65) for n in np.arange(65,75)]+\
               [('Y'+str(n),75) for n in np.arange(75,85)]+\
               [('Y_GE85',85)])
    #group age years into age groups
    lexpEU['agegroup'] = lexpEU.age.map(agegroup)   
    #drop younger ages (<20) 
    lexpEU = lexpEU.dropna(subset=['agegroup'])
    #age group average
    lexpEU = lexpEU.groupby(['geo\\time','agegroup']).mean()
    #average across years of interest
    lexpEU = lexpEU[['2015','2016','2017']].mean(axis=1)
    
    yrs_lost = pd.DataFrame()
    for lat,lon,city in zip(avgdiff_all.lat,avgdiff_all.lon,avgdiff_all.city):
    
        #find NUTS 2 code for the city centre using nuts-finder package (uses point-in-shape method)
        try:
            nuts2code = NutsFinder(year=2016).find(lat,lon)[2]['NUTS_ID']
        except Exception:
            nuts2code = NutsFinder().find(lat,lon)[2]['NUTS_ID']
                    
        #extract city value by NUTS2 code
        yrs_lost[city.item()]=lexpEU[lexpEU.index.get_level_values(0)==nuts2code].droplevel(level=0,axis=0)
    
    yrs_lost.to_csv('data/avg_yrs_lost_by_age_group.csv')


def yll_tables(uhi, cityagewgt=True):
    """
    calculate years of life lost based on age-dependent UHI mortality impact
    output YLL to tables
    option to weight by local age distribution or standard population distribution
    input: uhi = urban-rural difference, average per 100,000 per day
    """   
        
    if cityagewgt:
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
    
    #add years of life lost to dataset
    yrs_lost=pd.read_csv('data/avg_yrs_lost_by_age_group.csv',index_col=0)
    yrs_lost.columns.name = 'city'
    yrs_lost.index.name = 'age'
    uhi = uhi.assign(yrs_lost=yrs_lost)
    #multiply each death by years of life lost
    #for fAD, fADheat, fADcold, fAD_heat_ex, and fAD_cold_ex only
    uhi[['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']] = uhi[['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']]*uhi.yrs_lost
    
    #per year for net, heat, cold; per day for extremes
    uhi['fAD'] = (uhi.fAD).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')*365.24
    uhi['fADheat'] = (uhi.fADheat*uhi.heat_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fADcold'] = (uhi.fADcold*uhi.cold_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_heat_ex'] = uhi.fAD_heat_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_cold_ex'] = uhi.fAD_cold_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    
    #sort by net, descending
    out2 = uhi[['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']].sortby('fAD',ascending=False)
    out2 = out2.rename({'fAD':'Net (per year)','fADheat':'Heat (per year)','fADcold':'Cold (per year)','fAD_heat_ex':'Heat extreme (per day)','fAD_cold_ex':'Cold extreme (per day)'})
    out2.to_dataframe()[['Net (per year)','Heat (per year)','Cold (per year)','Heat extreme (per day)','Cold extreme (per day)']]\
        .to_csv('data/analysed/uhi_YLL_table_cityage.csv' if cityagewgt else
                'data/analysed/uhi_YLL_table_stdage.csv')


def mort_to_cost(filein):
    """convert YLL table to cost by value of life year"""
    mort_table = pd.read_csv(filein,index_col=0)
    cost_table = mort_table / 100000 * 46000
    cost_table.to_csv(filein.replace('YLL_table','YLL_cost_table'))


def print_yll_vs_vsl():
    """
    print out some comparisons of assessments by years of life lost vs
    by value of statistical life
    """
    yll = pd.read_csv('data/analysed/uhi_YLL_table_cityage.csv',index_col=0)
    yllcost = pd.read_csv('data/analysed/uhi_YLL_cost_table_cityage.csv',index_col=0)
    yll_std = pd.read_csv('data/analysed/uhi_YLL_table_stdage.csv',index_col=0)
    yllcost_std = pd.read_csv('data/analysed/uhi_YLL_cost_table_stdage.csv',index_col=0)
    yllall = yll.join(yll_std,rsuffix='_std').join(yllcost,rsuffix='_cost').join(yllcost_std,rsuffix='_stdcost')
    
    vsl = pd.read_csv('data/analysed/uhi_mortality_table_cityage.csv',index_col=0)
    vslcost = pd.read_csv('data/analysed/uhi_mortality_cost_table_cityage.csv',index_col=0)
    vsl_std = pd.read_csv('data/analysed/uhi_mortality_table_stdage.csv',index_col=0)
    vslcost_std = pd.read_csv('data/analysed/uhi_mortality_cost_table_stdage.csv',index_col=0)
    vslall = vsl.join(vsl_std,rsuffix='_std').join(vslcost,rsuffix='_cost').join(vslcost_std,rsuffix='_stdcost')
    
    #upper estimate of VOLY for YLL
    yllcost2 = yllcost/46000*116000
    
    #ratio of YLL/VSL costs for both definitions of VOLY
    yllvslcost = pd.concat([yllcost/vslcost,yllcost2/vslcost],keys=['best','upper'],axis=1)
    
    # values to quote in text
    print('YLLcost/VSLcost, median:')
    print((yllcost/vslcost).median())
    print('IQR:')
    print((yllcost/vslcost).quantile(0.25))
    print((yllcost/vslcost).quantile(0.75))
    print('median of net, heat, cold medians:')
    print((yllcost/vslcost).median()[['Net (per year)','Heat (per year)','Cold (per year)']].median())
    
    print('YLLcost higher estimate/VSLcost, median:')
    print((yllcost2/vslcost).median())
    print('IQR:')
    print((yllcost2/vslcost).quantile(0.25))
    print((yllcost2/vslcost).quantile(0.75))
    print('median of net, heat, cold medians:')
    print((yllcost2/vslcost).median()[['Net (per year)','Heat (per year)','Cold (per year)']].median())
    
    print('number of cities with adverse net effect, YLL:')
    print(yll[(yll['Net (per year)']>0)].count()[0])
    print('number of cities with adverse net effect, VSL:')
    print(vsl[(vsl['Net (per year)']>0)].count()[0])
    

if __name__ == '__main__':
    
    #urban-rural difference, average per 100,000 per day
    uhi = xr.open_dataset('data/analysed/avg_diff_from_rural_urbanrural.nc').sel(urbanrural='diff')
    
    #output average life years lost
    avg_yrs_lost(uhi)
    #output years of life lost tables
    yll_tables(uhi, cityagewgt=True)
    yll_tables(uhi, cityagewgt=False)
    #output cost tables by value of life year
    mort_to_cost('data/analysed/uhi_YLL_table_cityage.csv')
    mort_to_cost('data/analysed/uhi_YLL_table_stdage.csv')
    #print out comparisons between VSL and YLL(VOLY) 
    print_yll_vs_vsl()
    