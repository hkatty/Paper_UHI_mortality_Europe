# -*- coding: utf-8 -*-
"""
Figure 4 in the main text
Calculate correlation between attributable mortality measures and city properties

@author: Katty Huang
"""


import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from heatmapz import heatmap
import seaborn as sns


def cal_corr(age,agerisk=False,annualrisk=False,cityagewgt=True):
    """
    age: age group used for factors such as MMT, heat/cold days
    agerisk: True to show correlation with fAD of a specific age group only
    annualrisk: True to show also the annual risk (in addition to risk per day)
    cityagewgt: True to use local age distribution weighting
    posineg: True to use data where heat/cold are defined as when urbanfAF>/<ruralfAF
    """

    avgdiff_all_ages = xr.open_dataset('data/analysed/avg_diff_from_rural_urbanrural.nc')
    avgdiff_all = avgdiff_all_ages.sel(age=2085.1 if cityagewgt else 2085.5,urbanrural='diff')
    avgdiff_all_age = avgdiff_all_ages.sel(age=float(age.split(sep='-')[0].split('+')[0]))
    citydesc = pd.read_csv('data/analysed/city_latlon_avgT.csv',index_col=0)
    #heat/cold days count
    counts = avgdiff_all_age[['heat_count','cold_count']].sel(urbanrural='diff').to_dataframe()[['heat_count','cold_count']]
        
    #add age weight to xarray
    if cityagewgt:
        #local age distribution
        age_wgt = pd.read_csv('data/age_weights.csv',index_col=0).drop(index='std')
        age_wgt.index.name = 'city'
        age_wgt.columns.name = 'age'
        age_wgt.columns = age_wgt.columns.astype('float')
        age_wgt_xr = age_wgt.stack().to_xarray()
        age85to20 = (age_wgt[85]/age_wgt[20]).rename('Age ratio')
    else:
        #age group weighting by European standard population 2013
        #for age groups considered here (does not add up to 100% because missing younger age groups)
        ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
        age_wgt_std = ESP13/ESP13.sum()
        age_wgt_std.index.rename('age',inplace=True)
        age_wgt_xr = age_wgt_std.to_xarray()
        
    avgdiff_all_ages = avgdiff_all_ages.assign(age_wgt=age_wgt_xr)
    # #weighted sum across age groups
    # avgdiff_all = avgdiff_all_ages.weighted(avgdiff_all_ages.age_wgt.fillna(0)).sum(dim='age')
    
    if agerisk:
        avgdiff_all=avgdiff_all_age
    
    #fAD measures
    fAD = avgdiff_all[['fAD','fAD_heat_ex','fAD_cold_ex']].to_dataframe()[['fAD','fAD_heat_ex','fAD_cold_ex']]
    fAD_seas = avgdiff_all[['fAD_seasavg']].to_dataframe()[['fAD_seasavg']]
    fAD = fAD.merge(fAD_seas.unstack()['fAD_seasavg'],left_index=True,right_index=True)
    fAD = fAD[['fAD','fAD_heat_ex','fAD_cold_ex','DJF','MAM','JJA','SON']]
    fAD = fAD.rename(columns={'fAD':'annual net',#'fADheat':'heat','fADcold':'cold',
               'fAD_heat_ex':'heat extreme','fAD_cold_ex':'cold extreme',
               'DJF':'winter','MAM':'spring','JJA':'summer','SON':'autumn'})   
    
    if annualrisk:
        
        if agerisk:
            fADheat_annual = (avgdiff_all_age.fADheat*avgdiff_all_age.heat_count/3).rename('annual heat').to_series()
            fADcold_annual = (avgdiff_all_age.fADcold*avgdiff_all_age.cold_count/3).rename('annual cold').to_series()
            # fADnet_annual = (avgdiff_all_age.fAD*365.25).rename('annual net').to_series()
        else:
        
            fADheat_annual = (avgdiff_all_ages.fADheat*avgdiff_all_ages.heat_count/3).weighted(avgdiff_all_ages.age_wgt.fillna(0)).sum(dim='age')
            fADcold_annual = (avgdiff_all_ages.fADcold*avgdiff_all_ages.cold_count/3).weighted(avgdiff_all_ages.age_wgt.fillna(0)).sum(dim='age')
            # fADnet_annual = (avgdiff_all_ages.fAD*365.25).weighted(avgdiff_all_ages.age_wgt.fillna(0)).sum(dim='age')
            fADheat_annual = fADheat_annual.rename('annual heat').to_series()
            fADcold_annual = fADcold_annual.rename('annual cold').to_series()
            # fADnet_annual = fADnet_annual.rename('annual net').to_series()   
        
        fAD = pd.concat([fAD,fADheat_annual,fADcold_annual],axis=1)
    
    ##%%possible controlling factors
    #RR
    RRtmin=[]
    RRtmax=[]
    tmin=[]
    tmax=[]
    MMT=[]
    for city in fAD.index:
        RRfit = pd.read_csv('data/RRfit/RRfit_urbclimT_'+city.replace(' ','_').replace('-','_')+'_'+age+'.csv',index_col=0).allfit
        RRtmin.append((city,RRfit.iloc[0]))
        RRtmax.append((city,RRfit.iloc[-1]))
        tmin.append((city,RRfit.index[0]))
        tmax.append((city,RRfit.index[-1]))
        MMTcity = pd.read_csv('data/RRfit/MMT_urbclimT_'+city.replace(' ','_').replace('-','_')+'.csv',index_col=0)
        MMT.append((city,MMTcity[MMTcity.Age==age].MMT.values[0]))
        
    RRtmin = pd.Series(dict(RRtmin)).rename('RR$_\mathregular{Tmin}$')
    RRtmax = pd.Series(dict(RRtmax)).rename('RR$_\mathregular{Tmax}$')
    # RRratio = (RRtmax/RRtmin).rename('RR_ratio')
    tmin = pd.Series(dict(tmin)).rename('Tmin')
    tmax = pd.Series(dict(tmax)).rename('Tmax')
    # trange = (tmax-tmin).rename('T_range')
    MMT = pd.Series(dict(MMT)).rename('MMT')
    # MMT_percentile = ((MMT-tmin)/(tmax-tmin)).rename('MMT$_\mathregular{percentile}$')
    
    tas = citydesc.avgT.rename('T$_\mathregular{avg}$')
    # summerT = citydesc.JJA.rename('T$_\mathregular{summer}$')
    # winterT = citydesc.DJF.rename('T$_\mathregular{winter}$')
    uhi = (avgdiff_all.tas-273.15).rename('UHI').to_dataframe().UHI
    # uhisummer = (avgdiff_all.tas_seasavg.sel(season='JJA')-273.15).rename('UHI_summer').to_dataframe().UHI_summer.rename('UHI$_\mathregular{summer}$')
    # uhiwinter = (avgdiff_all.tas_seasavg.sel(season='DJF')-273.15).rename('UHI_winter').to_dataframe().UHI_winter.rename('UHI$_\mathregular{winter}$')
    lon = avgdiff_all.lon.to_dataframe().lon.rename('Longitude')
    lat = avgdiff_all.lat.to_dataframe().lat.rename('Latitude')
    # diag = np.sqrt((lon.max()-lon)**2+lat**2).rename('Northwest')
    
    heatdays = counts.heat_count.rename('Warm days')
    # colddays = counts.cold_count.rename('Cold days')
    
    if cityagewgt:
        factors = pd.concat([RRtmin,RRtmax,heatdays,uhi,age85to20,tas,lon,lat],axis=1)    
    else:
        # factors = pd.concat([RRtmin,RRtmax,RRratio,tmin,tmax,trange,MMT,tas,summerT,winterT,uhi,uhisummer,uhiwinter,lon,lat,diag],axis=1)
        factors = pd.concat([RRtmin,RRtmax,heatdays,uhi,tas,lon,lat],axis=1)        
    
    ##%%correlations
    spearman = fAD.apply(lambda s: factors.corrwith(s,method='spearman'))
    # pearson = fAD.apply(lambda s: factors.corrwith(s,method='pearson'))
    
    ##%%include p-value considerations
    coef=np.empty([fAD.shape[1],factors.shape[1]])
    coef[:]=np.nan
    pval=np.empty([fAD.shape[1],factors.shape[1]])
    pval[:]=np.nan
    for i in range(fAD.shape[1]):
        for j in range(factors.shape[1]):
            coef[i,j],pval[i,j]=scipy.stats.spearmanr(fAD.iloc[:,i],factors.iloc[:,j])
    
    spearman = pd.DataFrame(data=coef,index=fAD.columns,columns=factors.columns).T 
    p = pd.DataFrame(data=pval,index=fAD.columns,columns=factors.columns).T    
    return spearman,p


def plot_corr(spearman,p,fig,title,vmax=0):
    """plot correlations as heatmap"""
    
    spearman = spearman.where(p<0.01).fillna(0)
    
    vmax = abs(spearman.unstack()).max()
    if vmax != 0:
        vmax = vmax
    
    heatmap(x=spearman.unstack().index.get_level_values(0).str.capitalize(),
            y=spearman.unstack().index.get_level_values(1),
            x_order=spearman.unstack().index.get_level_values(0).unique().str.capitalize(),
            y_order=spearman.unstack().index.get_level_values(1).unique()[::-1],
            color=spearman.unstack(),
            size=abs(spearman.unstack()),
            color_range=(-vmax,vmax),
            palette=sns.diverging_palette(220, 20, n=200)[:],
            size_scale=400,
            fig=fig,
            clabel='Spearman\'s correlation',
            title=title)


def cal_corr_factors(age,cityagewgt=False):
    """
    correlation between factors
    age: age group used for factors such as MMT, heat/cold days
    cityagewgt: True to use local age distribution weighting
    """
    
    avgdiff_all_ages = xr.open_dataset('data/analysed/avg_diff_from_rural_urbanrural.nc')
    avgdiff_all = avgdiff_all_ages.sel(age=2085.1 if cityagewgt else 2085.5,urbanrural='diff')
    avgdiff_all_age = avgdiff_all_ages.sel(age=float(age.split(sep='-')[0].split('+')[0]))
    avgT = pd.read_csv('data/analysed/city_latlon_avgT.csv',index_col=0).avgT
    #heat/cold days count
    counts = avgdiff_all_age[['heat_count','cold_count']].sel(urbanrural='diff').to_dataframe()[['heat_count','cold_count']]
    
    #age weight
    if cityagewgt:
        #local age distribution
        age_wgt = pd.read_csv('data/age_weights.csv',index_col=0).drop(index='std')
        age_wgt.index.name = 'city'
        age_wgt.columns.name = 'age'
        age_wgt.columns = age_wgt.columns.astype('float')
        age85to20 = (age_wgt[85]/age_wgt[20]).rename('Age ratio')
    else:
        #age group weighting by European standard population 2013
        #for age groups considered here (does not add up to 100% because missing younger age groups)
        ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
        age_wgt_std = ESP13/ESP13.sum()
        age_wgt_std.index.rename('age',inplace=True)
   
    #possible controlling factors
    #RR
    RRtmin=[]
    RRtmax=[]
    tmin=[]
    tmax=[]
    MMT=[]
    for city in avgdiff_all.city.values:
        RRfit = pd.read_csv('data/RRfit/RRfit_urbclimT_'+city.replace(' ','_').replace('-','_')+'_'+age+'.csv',index_col=0).allfit
        RRtmin.append((city,RRfit.iloc[0]))
        RRtmax.append((city,RRfit.iloc[-1]))
        tmin.append((city,RRfit.index[0]))
        tmax.append((city,RRfit.index[-1]))
        MMTcity = pd.read_csv('data/RRfit/MMT_urbclimT_'+city.replace(' ','_').replace('-','_')+'.csv',index_col=0)
        MMT.append((city,MMTcity[MMTcity.Age==age].MMT.values[0]))
        
    RRtmin = pd.Series(dict(RRtmin)).rename('RR$_\mathregular{Tmin}$')
    RRtmax = pd.Series(dict(RRtmax)).rename('RR$_\mathregular{Tmax}$')
    # RRratio = (RRtmax/RRtmin).rename('RR_ratio')
    tmin = pd.Series(dict(tmin)).rename('Tmin')
    tmax = pd.Series(dict(tmax)).rename('Tmax')
    # trange = (tmax-tmin).rename('T_range')
    MMT = pd.Series(dict(MMT)).rename('MMT')
    # MMT_percentile = ((MMT-tmin)/(tmax-tmin)).rename('MMT$_\mathregular{percentile}$')
    
    tas = avgT.rename('T$_\mathregular{avg}$')
    # summerT = (avgdiff_all.tas_seasavg+avgdiff_all.tas_rural_seasavg-273.15).sel(season='JJA').rename('T_summer').to_dataframe().T_summer.rename('T$_\mathregular{summer}$')
    # winterT = (avgdiff_all.tas_seasavg+avgdiff_all.tas_rural_seasavg-273.15).sel(season='DJF').rename('T_winter').to_dataframe().T_winter.rename('T$_\mathregular{winter}$')
    uhi = (avgdiff_all.tas-273.15).rename('UHI').to_dataframe().UHI
    # uhisummer = (avgdiff_all.tas_seasavg.sel(season='JJA')-273.15).rename('UHI_summer').to_dataframe().UHI_summer.rename('UHI$_\mathregular{summer}$')
    # uhiwinter = (avgdiff_all.tas_seasavg.sel(season='DJF')-273.15).rename('UHI_winter').to_dataframe().UHI_winter.rename('UHI$_\mathregular{winter}$')
    lon = avgdiff_all.lon.to_dataframe().lon.rename('Longitude')
    lat = avgdiff_all.lat.to_dataframe().lat.rename('Latitude')
    # diag = np.sqrt((lon.max()-lon)**2+lat**2).rename('Northwest')
    
    heatdays = counts.heat_count.rename('Warm days')
    # colddays = counts.cold_count.rename('Cold days')
    
    if cityagewgt:
        factors = pd.concat([RRtmin,RRtmax,heatdays,uhi,age85to20,tas,lon,lat],axis=1)    
    else:
        # factors = pd.concat([RRtmin,RRtmax,RRratio,tmin,tmax,trange,MMT,tas,summerT,winterT,uhi,uhisummer,uhiwinter,lon,lat,diag],axis=1)
        factors = pd.concat([RRtmin,RRtmax,heatdays,uhi,tas,lon,lat],axis=1)        
       
    #correlations including p-value considerations
    coef=np.empty([factors.shape[1],factors.shape[1]])
    coef[:]=np.nan
    pval=np.empty([factors.shape[1],factors.shape[1]])
    pval[:]=np.nan
    for i in range(factors.shape[1]):
        for j in range(factors.shape[1]):
            coef[i,j],pval[i,j]=scipy.stats.spearmanr(factors.iloc[:,i],factors.iloc[:,j])
    
    spearman = pd.DataFrame(data=coef,index=factors.columns,columns=factors.columns).T 
    p = pd.DataFrame(data=pval,index=factors.columns,columns=factors.columns).T
    
    return spearman,p


def plot_corr_factors(spearman,p,fig,title,vmax=0):
    """plot heatmap showing correlation between factors"""
    
    spearman = spearman.where(p<0.01).fillna(0)
    
    vmax = abs(spearman.unstack()).max()
    if vmax != 0:
        vmax = vmax
    
    heatmap(x=spearman.unstack().index.get_level_values(0),
            y=spearman.unstack().index.get_level_values(1),
            x_order=spearman.unstack().index.get_level_values(0).unique(),
            y_order=spearman.unstack().index.get_level_values(1).unique()[::-1],
            color=spearman.unstack(),
            size=abs(spearman.unstack()),
            color_range=(-vmax,vmax),
            palette=sns.diverging_palette(220, 20, n=200)[:],
            size_scale=400,
            fig=fig,
            clabel='Spearman\'s correlation',
            title=title)
    


if __name__ == '__main__':
        
    #correlation between mortality risks and factors, local age structure (Fig. 4)
    fig=plt.figure()
    spearman,p = cal_corr(age='65-75',agerisk=False,annualrisk=False,cityagewgt=True)
    plot_corr(spearman,p,fig,'')
    plt.savefig('figures/correlation_cityage.pdf')
    
    #correlation between factors (Fig. S14)
    fig=plt.figure()
    spearman,p = cal_corr_factors(age='65-75',cityagewgt=True)
    plot_corr_factors(spearman,p,fig,'')
    plt.savefig('figures/correlation_factors.pdf')
