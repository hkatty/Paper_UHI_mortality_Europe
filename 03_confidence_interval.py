# -*- coding: utf-8 -*-
"""
calculate confidence interval for each city estimate in supplement table
and for estimate of number of cities with adverse annual net impact
based on uncertainty in exposure-response relationship, represented by 
1000 Monte Carlo simulations

@author: Katty Huang
"""

import xarray as xr
import pandas as pd


def table_calc(uhi):
    #per year for net, heat, cold; per day for extremes
    uhi['fAD'] = (uhi.fADdiff).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')*365.24
    uhi['fADheat'] = (uhi.fADdiff_heat*uhi.heat_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fADcold'] = (uhi.fADdiff_cold*uhi.cold_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_heat_ex'] = uhi.fADdiff_heat_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_cold_ex'] = uhi.fADdiff_cold_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi = uhi[['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']]    
    return uhi
    
def table_calc2(uhi):
    #per year for net, heat, cold; per day for extremes
    uhi['fAD'] = (uhi.fAD).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')*365.24
    uhi['fADheat'] = (uhi.fADheat*uhi.heat_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fADcold'] = (uhi.fADcold*uhi.cold_count/3).weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_heat_ex'] = uhi.fAD_heat_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi['fAD_cold_ex'] = uhi.fAD_cold_ex.weighted(uhi.age_wgt.fillna(0)).sum(dim='age')
    uhi = uhi[['fAD','fADheat','fADcold','fAD_heat_ex','fAD_cold_ex']]    
    return uhi

def get_yll(uhi,city):
    #add years of life lost to dataset
    yrs_lost=pd.read_csv('data/avg_yrs_lost_by_age_group.csv',index_col=0)
    yrs_lost.columns.name = 'city'
    yrs_lost.index.name = 'age'
    uhi = uhi.assign(yrs_lost=yrs_lost[city])
    #multiply each death by years of life lost
    #for fAD, fADheat, fADcold, fAD_heat_ex, and fAD_cold_ex only
    keys = [k for k in uhi.keys() if k.startswith('fAD')]
    uhi[keys] = uhi[keys]*uhi.yrs_lost
    return uhi


def do_analysis(cityage=True,yll=False):
    
    cities_list = pd.read_csv('data/city_climate_zones.csv').city
    
    for i,cityfile in enumerate(cities_list):
    
        cityfile = cityfile.replace(' ','_').replace('-','_')
        #city name in Copernicus dataset (use spaces rather than '_' used in filename)
        city = (cityfile.replace('_','-') if cityfile=='Cluj_Napoca' else cityfile.replace('_',' '))
            
        if cityage:
            #local age distribution
            age_wgt = pd.read_csv('data/age_weights.csv',index_col=0)[1:].T
            age_wgt.columns.name = 'city'
            age_wgt.index = [float(age) for age in age_wgt.index]
            age_wgt.index.name = 'age'
            age_wgt_xr = age_wgt.stack().swaplevel().to_xarray().sel(city=city)
        else:
            #age group weighting by European standard population 2013
            #for age groups considered here (does not add up to 100% because missing younger age groups)
            ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
            age_wgt_std = ESP13/ESP13.sum()
            age_wgt_std.index.rename('age',inplace=True)
            age_wgt_xr = age_wgt_std.to_xarray()
        
        #urban-rural difference, average per 100,000 per day
        uhi = xr.open_dataset('data/analysed/simulated_data_annual_seasavg_extremes_heatcold_counts_'+cityfile+'.nc')
        uhi = uhi.assign(age_wgt=age_wgt_xr)
        
        timeseriesuhi = xr.open_dataset('data/analysed/data_urbanruralavg_timeseries_analysis_'+cityfile+'.nc').sel(urbanrural='diff')
        timeseriesuhi = timeseriesuhi.assign(age_wgt=age_wgt_xr)

        if yll:
            #multiply by avg years lost to get YLL
            uhi = get_yll(uhi,city)
            timeseriesuhi = get_yll(timeseriesuhi, city)
                
        uhi = table_calc(uhi)
        qs = uhi.quantile([0.05,0.5,0.95],dim='simulation')

        timeseriesuhi = table_calc2(timeseriesuhi)    
        
        qs = xr.concat([qs,timeseriesuhi.expand_dims(dim={'quantile':['timeseries']})],dim='quantile')
        
        #combine all cities outputs for later analysis & plotting
        if i == 0:
            data_all = qs.expand_dims(dim={'city':[city]})
            fAD_all = uhi['fAD'].expand_dims(dim={'city':[city]})
            timeseriesuhi_all = timeseriesuhi.expand_dims(dim={'city':[city]})
            uhi_all = uhi.expand_dims(dim={'city':[city]})
        else:
            data_all = xr.concat([data_all,qs.expand_dims(dim={'city':[city]})],
                               dim='city')
            fAD_all = xr.concat([fAD_all,uhi['fAD'].expand_dims(dim={'city':[city]})],
                               dim='city')
            timeseriesuhi_all = xr.concat([timeseriesuhi_all,timeseriesuhi.expand_dims(dim={'city':[city]})],
                                dim='city')
            uhi_all = xr.concat([uhi_all,uhi.expand_dims(dim={'city':[city]})],
                               dim='city')
        
    return data_all, fAD_all, timeseriesuhi_all, uhi_all


def quantiles_data(data_all):
    #cities sorted by median net, descending
    data_all = data_all.sortby(data_all.fAD.sel(quantile=0.5),ascending=False)
    
    name_dict = {'fAD':'Net [per year]','fADheat':'Heat [per year]','fADcold':'Cold [per year]','fAD_heat_ex':'Heat extreme [per day]','fAD_cold_ex':'Cold extreme [per day]'}

    return data_all.to_dataframe()[name_dict.keys()].rename(name_dict)


def to_outstr(data_all,cost=''):
    
    if cost.casefold() == ('vsl').casefold():
        #cost by VSL (datain needs to be mortality counts)
        data_all = data_all / 100000 * 3910000
    elif cost.casefold() == ('yll').casefold():
        #cost by YLL (datain needs to be years of life lost)
        data_all = data_all / 100000 * 46000
        
    #specify number of decimal places to round to
    round_dict = {'fAD':1,'fADheat':1,'fADcold':1,'fAD_heat_ex':2,'fAD_cold_ex':2} 
    name_dict = {'fAD':'Net [per year]','fADheat':'Heat [per year]','fADcold':'Cold [per year]','fAD_heat_ex':'Heat extreme [per day]','fAD_cold_ex':'Cold extreme [per day]'}
        
    #for supplementary table, show confidence interval in brackets ###################
    upper = data_all.sel(quantile=0.95).to_dataframe()[list(name_dict.keys())]
    lower = data_all.sel(quantile=0.05).to_dataframe()[list(name_dict.keys())]
    timeseries = data_all.sel(quantile='timeseries').to_dataframe()[list(name_dict.keys())]  

    quantiles = '(' + lower.round(round_dict).astype('str') + ', ' + upper.round(round_dict).astype('str') + ')'
    outstr = timeseries.round(round_dict).astype('str') + ' ' + quantiles 
    outstr = outstr.rename(columns=name_dict)
    
    # sort by best estimate net, descending
    cityorder = timeseries.sort_values(by='fAD',ascending=False).index
    outstr = outstr.reindex(cityorder)
    
    return outstr



if __name__ == '__main__':

    # tables of risks and economic assessments in supplementary materials
    for cityage in [True,False]:
        data_all = do_analysis(cityage=cityage,yll=False)[0]
        outstr = to_outstr(data_all,cost='')
        if cityage:
            outstr.to_csv('data/analysed/uhi_mortality_table_ci_cityage.csv')
            # quantiles_data(data_all).to_csv('data/analysed/simulated_quantiles_table_cityage.csv')
        else:
            outstr.to_csv('data/analysed/uhi_mortality_table_ci_stdage.csv')
            # quantiles_data(data_all).to_csv('data/analysed/simulated_quantiles_table_stdage.csv')
            
        outstr = to_outstr(data_all,cost='vsl')
        if cityage:
            outstr.to_csv('data/analysed/uhi_mortality_cost_table_ci_cityage.csv')
        else:
            outstr.to_csv('data/analysed/uhi_mortality_cost_table_ci_stdage.csv')
            
        data_all = do_analysis(cityage=cityage,yll=True)[0]
        outstr = to_outstr(data_all,cost='')
        if cityage:
            outstr.to_csv('data/analysed/uhi_yll_table_ci_cityage.csv')
        else:
            outstr.to_csv('data/analysed/uhi_yll_table_ci_stdage.csv')
            
        outstr = to_outstr(data_all,cost='yll')
        if cityage:
            outstr.to_csv('data/analysed/uhi_yll_cost_table_ci_cityage.csv')
        else:
            outstr.to_csv('data/analysed/uhi_yll_cost_table_ci_stdage.csv')
    
    
    # values quoted in paper text
    
    #confidence interval of number of cities with adverse annual net impact #######
    data_all, fAD_all, best_all, timeseriesuhi_all, uhi_all = do_analysis(cityage=True,yll=False)
    adversenum_ci = (fAD_all>0).sum(dim='city')
    print(adversenum_ci.quantile([0.05,0.5,0.95],dim='simulation'))
    adversenum = (best_all.fAD>0).sum(dim='city')
    print('best fit:', adversenum.item())
    adversenum_timeseries = (timeseriesuhi_all.fAD>0).sum(dim='city')
    print('timeseries best fit:', adversenum_timeseries.item())
    ###############################################################################
    
    # median values quoted in text
    avgdiff_all_cityage = timeseriesuhi_all #simulation to use in text for best estimate
    print('annual mean net per year, median of all cities:')
    print(avgdiff_all_cityage.fAD.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fAD.quantile(0.25).values)
    print(avgdiff_all_cityage.fAD.quantile(0.75).values)
    print('CI of cities\' median:')
    print(uhi_all.fAD.median(dim='city').quantile([0.05,0.5,0.95],dim='simulation'))
    
    print('heat extreme impact greater than cold extreme for each city:')
    print((avgdiff_all_cityage.fAD_heat_ex>avgdiff_all_cityage.fAD_cold_ex).all().values)
    
    print('daily mean during heat extreme, median of all cities:')
    print(avgdiff_all_cityage.fAD_heat_ex.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fAD_heat_ex.quantile(0.25).values)
    print(avgdiff_all_cityage.fAD_heat_ex.quantile(0.75).values)
    print('CI of cities\' median:')
    print(uhi_all.fAD_heat_ex.median(dim='city').quantile([0.05,0.5,0.95],dim='simulation'))
    
    print('daily mean during cold extreme, median of all cities:')
    print(avgdiff_all_cityage.fAD_cold_ex.median().values)
    print('IQR:')
    print(avgdiff_all_cityage.fAD_cold_ex.quantile(0.25).values)
    print(avgdiff_all_cityage.fAD_cold_ex.quantile(0.75).values)
    print('CI of cities\' median:')
    print(uhi_all.fAD_cold_ex.median(dim='city').quantile([0.05,0.5,0.95],dim='simulation'))

