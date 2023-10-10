# -*- coding: utf-8 -*-
"""
Figure 5 in the main text. Economic assessments and comparison to other costs.

@author: Katty Huang
"""

import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from scipy.stats import linregress


def uhi_cost():
    """mortality cost associated with UHI based on VSL"""
    #urban-rural difference, average per 100,000 per day
    avgdiff_all = xr.open_dataset('data/analysed/avg_diff_from_rural_urbanrural.nc').sel(urbanrural='diff')
    
    for i,cityagewgt in enumerate([True,False]):#True to use local age distribution weighting
        #add age weight to xarray
        if cityagewgt:
            #local age distribution
            age_wgt = pd.read_csv('data/age_weights.csv',index_col=0)[1:].T
            age_wgt.columns.name = 'city'
            age_wgt.index = [float(age) for age in age_wgt.index]
            age_wgt.index.name = 'age'
            age_wgt_xr = age_wgt.stack().swaplevel().to_xarray()        
            avgdiff_all = avgdiff_all.assign(age_wgt=age_wgt_xr)
        else:
            #age group weighting by European standard population 2013
            #for age groups considered here (does not add up to 100% because missing younger age groups)
            ESP13 = pd.Series(data=[32.5,26.5,10.5,6.5,2.5],index=[20,45,65,75,85]) 
            age_wgt_std = ESP13/ESP13.sum()
            age_wgt_std.index.rename('age',inplace=True)
            avgdiff_all = avgdiff_all.assign(age_wgt=age_wgt_std.to_xarray())
        
        #### by value of statistical life ####################
        #cost per capita per year associated with heat/cold, calculated separately for each age group (because heat/cold threshold differs) and combined by weighting
        #cost by value of statisticl life
        #$3.6 million 2005 USD -> 3.0708 million 2005 EUR -> 3.91 million 2021 EUR
        #avg fADheat/cold per 100,000 per day /100,000 * 3,910,000 EUR per person * number of heat/cold days over whole period (3 years) / 3 years = cost per year
        fADheat_cost = (avgdiff_all.fADheat*39.1*avgdiff_all.heat_count/3).weighted(avgdiff_all.age_wgt.fillna(0)).sum(dim='age').expand_dims({'age_structure':['local' if cityagewgt else 'std']})
        fADcold_cost = (avgdiff_all.fADcold*39.1*avgdiff_all.cold_count/3).weighted(avgdiff_all.age_wgt.fillna(0)).sum(dim='age').expand_dims({'age_structure':['local' if cityagewgt else 'std']})
        
        #net cost per year
        fADnet_cost = (avgdiff_all.fAD*39.1*365.25).weighted(avgdiff_all.age_wgt.fillna(0)).sum(dim='age').expand_dims({'age_structure':['local' if cityagewgt else 'std']})
        #####################################################
        
        if i==0:
            fADheat_cost_all = fADheat_cost
            fADcold_cost_all = fADcold_cost
            fADnet_cost_all = fADnet_cost
        else:
            fADheat_cost_all = xr.concat([fADheat_cost_all,fADheat_cost],dim='age_structure')
            fADcold_cost_all = xr.concat([fADcold_cost_all,fADcold_cost],dim='age_structure')
            fADnet_cost_all = xr.concat([fADnet_cost_all,fADnet_cost],dim='age_structure')
    
    
    #match air pollution and heat/cold cities, standard age for all unless specified
    tcost = xr.merge([fADheat_cost_all.sel(age_structure='std').rename('fADheat'),
                      fADcold_cost_all.sel(age_structure='std').rename('fADcold'),
                      fADnet_cost_all.sel(age_structure='std').rename('fADnet'),
                      fADheat_cost_all.sel(age_structure='local').rename('fADheat_city'),
                      fADcold_cost_all.sel(age_structure='local').rename('fADcold_city'),
                      fADnet_cost_all.sel(age_structure='local').rename('fADnet_city'),
                      avgdiff_all.lon,avgdiff_all.lat],compat='override')
    return tcost,avgdiff_all


def lancet_air_pollution(pollutant):
    """Air pollution costs from Khomenko et al 2021"""
    #city information from Pierre (to extract lat-lon)
    citydesc = pd.read_csv('data/citydesc.csv',index_col=0,encoding='latin1')
    #list of copernicus and matching epi city names
    citynames = pd.read_csv('data/copernicus_epiv3_cities_match.csv',index_col=0,encoding='latin1',keep_default_na=False)
    
    #city and greater city codes 
    citycode_subset = citydesc[citydesc.LABEL.isin(citynames.epi)].URAU_CODE
    
    if pollutant == 'pm2p5':
        files = ['pm2p5_cities.csv']
    elif pollutant == 'no2':
        files = ['no2_cities.csv']
    elif pollutant == 'both':
        files = ['pm2p5_cities.csv','no2_cities.csv']
    groupcode = 'City code'
    var = 'Deaths/100,000 population (95% CI)'
    
    dataall = pd.DataFrame()
    for file in files:
        data=pd.read_csv('data/air_pollution/'+file,index_col=groupcode)
        datasum=data.groupby(data.index)[var].sum()
        datasum_subset = datasum[datasum.index.isin(citycode_subset)]
        lonlats = citydesc[citydesc.URAU_CODE.isin(datasum_subset.index)][['URAU_CODE','lon','lat']]
        lonlats = lonlats.groupby('URAU_CODE').mean()
        name = citynames.merge(citydesc,left_on='epi',right_on='LABEL')[['copernicus','URAU_CODE']].drop_duplicates()
        # datacities = pd.merge(datasum_subset,lonlats,left_index=True,right_index=True)
        datacities = pd.merge(pd.merge(datasum_subset,lonlats,left_index=True,right_index=True),name,left_index=True,right_on='URAU_CODE').set_index('URAU_CODE')
        datacities['Euro per capita per year'] = datacities[var]*39
        dataall = pd.concat([dataall,datacities])
    
    if pollutant == 'both':
        pm25no2=dataall[['Deaths/100,000 population (95% CI)','Euro per capita per year']].groupby(dataall.index).sum()
        pm25no2=pm25no2.merge(dataall[['lon','lat']].groupby(dataall.index).mean(),left_index=True,right_index=True)
        pm25no2=pm25no2.merge(dataall['copernicus'],left_index=True,right_index=True).drop_duplicates()
    
    if pollutant == 'pm2p5':
        airpollution = dataall
        title = 'PM2.5'
    elif pollutant == 'no2':
        airpollution = dataall
        title = 'NO2'
    else:
        airpollution = pm25no2
        title = 'PM2.5+NO2'
    return airpollution,title


def get_citycodes():
    """list of city codes"""
    #city information from Pierre (to extract lat-lon)
    citydesc = pd.read_csv('data/citydesc_v3.csv',index_col=0,encoding='latin1')
    #list of copernicus and matching epi city names
    citynames = pd.read_csv('data/copernicus_epiv3_cities_match.csv',index_col=0,encoding='latin1',keep_default_na=False)
    #city and greater city codes 
    # citycodes = citydesc[citydesc.LABEL.isin(citynames.epi)].URAU_CODE.drop_duplicates()
    #list of matching copernicus & epi city names with city code as index 
    citycodesnames = citydesc.merge(citynames,left_on='LABEL',right_on='epi')[['URAU_CODE','copernicus','epi']].drop_duplicates().set_index('URAU_CODE')
    #there seems to be a bug in the lookup table, ES552K1 mislabelled as Barcelona, which should be ES002K2 and is already accounted for
    return citycodesnames.drop('ES552K1')

    
def add_air_pollution_costs(tcost_df):
    """
    O3 (and PM2.5 but not used) air pollution cost from EEA 
    and PM2.5 cost from Khomenko et al. 2021
    """
    
    #air pollution Khomenko et al 2021 (for PM2.5)
    airpollution,pollutant = lancet_air_pollution(pollutant='pm2p5')
    bothcosts = pd.merge(tcost_df,airpollution,left_on='city',right_on='copernicus')

    #European Environment Agency estimate of O3 and PM2.5 cost
    eea = pd.read_csv('data/air_pollution/eea_air_quality_health_risk_assessment2015to2017.csv')
    eea['Cost per person'] = eea['Premature Deaths']/eea['Total City Population *']*3.91*10**6
    eea_mcost = eea.groupby(['City Code','Air Pollutant'])['Cost per person'].mean().unstack()
    
    subset_cities = get_citycodes()
    eea_mcost_match = eea_mcost[['O3','PM2.5']].merge(subset_cities['copernicus'],left_index=True,right_index=True)
    return bothcosts.merge(eea_mcost_match,left_on='copernicus',right_on='copernicus')


def prep_eurostat(file,sex=False,unit=True):
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
    return x[x.sex=='T'].drop(columns=np.array(['unit','sex'])[[unit,sex]].tolist()) if sex else x.drop(columns=np.array(['unit'])[[unit]].tolist())


def rent_and_transit():
    """rent and public transit costs"""
    #rent calculated by multiplying area per person with rent per area
    livcon = prep_eurostat('data/rent_transport/urb_clivcon.tsv',unit=False)
    #SA1049V average annual rent for housing per m2 (EUR)
    #SA1022V average area of living accommodation (m2/person)
    rentperm2 = livcon[livcon.indic_ur=='SA1049V'].set_index('cities\\time').drop(columns=['indic_ur'])
    areaperperson = livcon[livcon.indic_ur=='SA1022V'].set_index('cities\\time').drop(columns=['indic_ur'])
    rent = (rentperm2*areaperperson).dropna(axis=0,how='all')
    
    #public transit cost
    tran = prep_eurostat('data/rent_transport/urb_ctran.tsv',unit=False)
    #TT1080V Cost of a combined monthly ticket (all modes of public transport) for 5-10 km in the central zone - EUR
    transit = tran[tran.indic_ur=='TT1080V'].set_index('cities\\time').drop(columns=['indic_ur'])

    #cost for all cities with data (exclude country level data)
    rent_all_cities = rent[rent.index.str.endswith('C')][['2021','2020','2019','2018','2017','2016','2015','2014','2013','2012']].dropna(how='all').mean(axis=1)/12
    transit_all_cities = transit[['2021','2020','2019','2018','2017','2016','2015','2014','2013','2012']].dropna(how='all').mean(axis=1)
    
    return pd.merge(rent_all_cities.rename('Rent\n(n=149)'),transit_all_cities.rename('Transit\n(n=357)'),how='outer',left_index=True,right_index=True)


def plot_map_scatter_cmap(x,y,z,label,title,cmap,symmetric=False,ax=False):
    """
    x=lons, y=lats, z=variable to plot as color
    label=colourbar label
    title=plot title
    cmap=colourmap
    symmetric=colourbar symmetric around zero if True
    create figure if ax argument not passed
    """
    
    if not ax:
        # plt.figure()
        # ax = plt.axes(projection=ccrs.EuroPP())
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.EuroPP()})
    ax.set_extent([-15, 25, 30, 62],crs=ccrs.PlateCarree())#need to set to avoid issue with scatter plot coordinate
    ax.coastlines(color='grey',zorder=0)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, zorder=0)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([35, 45, 55])
    gl.xlocator = mticker.FixedLocator([-20, 0, 20])
    gl.xlabel_style = {'size': 'small', 'color': 'gray'}
    gl.ylabel_style = {'size': 'small', 'color': 'gray'}
    if symmetric:
        vmax = max(abs(z.max()),abs(z.min()))
        im = ax.scatter(x=x,y=y,c=z,s=25,transform=ccrs.PlateCarree(),cmap=cmap,vmax=vmax,vmin=-vmax)
    else:
        im = ax.scatter(x=x,y=y,c=z,s=25,transform=ccrs.PlateCarree(),cmap=cmap)
    ax.scatter(x=x,y=y,s=25,transform=ccrs.PlateCarree(),edgecolor='k',linewidth=0.6,facecolor='None')
    plt.colorbar(im,label=label,ax=ax)
    ax.set_title(title)
    
def plotlinearfit(ax,title,x,y,xlabel,ylabel,dotlabel='',ls='-'):
    """plot scatter with linear regression"""
    slope, intercept, r, p, std_err = linregress(x,y)
    mn=np.min(x)
    mx=np.max(x)
    x1=np.linspace(mn,mx,500)
    y1=slope*x1+intercept
    ax.plot(x,y,'o',label=dotlabel)
    ax.plot(x1,y1,color='r',ls=ls,label='r = '+('%.2g' % r)+', p = '+('%.1g' % p))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.set_title(title)


def plot_cost_comparison(avgdiff_all,tcost):
    """
    Figure 5 in main text. Cost comparisons 
    Heat/cold/net maps (local age structure), heat vs air pollution
    """
    fig = plt.figure(figsize=(12,6.5))
    ax1 = fig.add_subplot(231,projection=ccrs.EuroPP())
    ax2 = fig.add_subplot(232,projection=ccrs.EuroPP())
    ax3 = fig.add_subplot(233,projection=ccrs.EuroPP())
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    
    #statistical life
    #heat, local age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADheat_city,'€ per capita per year','(a) Heat-related deaths',cmap='RdPu',ax=ax1)
    
    #cold, local age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADcold_city,'€ per capita per year','(b) Cold-related deaths',cmap='BuPu_r',ax=ax2)
   
    #annual net, local age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADnet_city,'€ per capita per year','(c) Annual T-related deaths',cmap='seismic',symmetric=True,ax=ax3)
    net_adverse = avgdiff_all[['lat','lon']].where(tcost.fADnet_city>0).dropna(dim='city')
    ax3.scatter(x=net_adverse.lon,y=net_adverse.lat,s=25,transform=ccrs.PlateCarree(),edgecolor='r',linewidth=0.7,facecolor='None')
    
    #cold vs heat, local age
    plotlinearfit(ax=ax4,title='(d) Heat vs. cold',x=tcost.fADcold_city.to_series(),y=tcost.fADheat_city.to_series(),xlabel='Cold-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Heat-related deaths (€ capita$^{-1}$ y$^{-1}$)')
    lims = [np.min(np.abs([ax4.get_xlim(), ax4.get_ylim()])),np.max(np.abs([ax4.get_xlim(), ax4.get_ylim()]))]
    ax4.plot([-x for x in lims],lims,'k:',label='one-to-one line')
    ax4.legend()
    
    # #compare to EEA air pollution, local age
    # plotlinearfit(ax=ax5,title='(e) Heat vs. O3',x=bothcosts.O3,y=bothcosts.fADheat_city,xlabel='O3-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Heat-related deaths (€ capita$^{-1}$ y$^{-1}$)')
    # #compare to Khomenko air pollution, standard age
    # plotlinearfit(ax=ax6,title='(f) Heat vs. PM2.5\nage standardised',x=bothcosts['Euro per capita per year'],y=bothcosts.fADheat,xlabel=pollutant+'-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Heat-related deaths (€ capita$^{-1}$ y$^{-1}$)')
    
    #compare to various air pollution, one plot
    #compare to EEA air pollution, local age
    plotlinearfit(ax=ax5,title='',x=bothcosts.O3,y=bothcosts.fADheat_city,xlabel='',ylabel='',dotlabel='ozone',ls='--')
    # plotlinearfit(ax=ax5,title='',x=bothcosts['PM2.5'],y=bothcosts.fADheat_city,xlabel='',ylabel='',dotlabel='PM2.5 all',ls='-.')
    #compare to Khomenko air pollution, standard age
    plotlinearfit(ax=ax5,title='(e) Heat vs. air pollution',x=bothcosts['Euro per capita per year'],y=bothcosts.fADheat,
                  xlabel='Air pollution-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Heat-related deaths (€ capita$^{-1}$ y$^{-1}$)',
                  dotlabel='PM2.5 (std age)')
    
    #cost comparison
    all_costs.plot.box(ax=ax6)
    ax6.set_title('(f) Cost comparison')
    ax6.set_ylabel('Average cost (€ capita$^{-1}$ month$^{-1}$)')
    
    plt.tight_layout()
    plt.savefig('figures/cost_compare.pdf')


def plot_cost_comparison_std(avgdiff_all,tcost):
    """figure S16 for supplement (standard age heat/cold/net, & cold vs pollution)"""
    fig = plt.figure(figsize=(12,6.5))
    ax1 = fig.add_subplot(231,projection=ccrs.EuroPP())
    ax2 = fig.add_subplot(232,projection=ccrs.EuroPP())
    ax3 = fig.add_subplot(233,projection=ccrs.EuroPP())
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    
    #statistical life
    #heat, std age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADheat,'€ per capita per year','(a) Heat-related deaths\nage standardised',cmap='RdPu',ax=ax1)
    
    #cold, std age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADcold,'€ per capita per year','(b) Cold-related deaths\nage standardised',cmap='BuPu_r',ax=ax2)
    
    #annual net, std age
    plot_map_scatter_cmap(avgdiff_all.lon,avgdiff_all.lat,tcost.fADnet,'€ per capita per year','(c) Annual T-related deaths\nage standardised',cmap='seismic',symmetric=True,ax=ax3)
    net_adverse = avgdiff_all[['lat','lon']].where(tcost.fADnet>0).dropna(dim='city')
    ax3.scatter(x=net_adverse.lon,y=net_adverse.lat,s=25,transform=ccrs.PlateCarree(),edgecolor='r',linewidth=0.7,facecolor='None')
    
    #cold vs heat, std age
    plotlinearfit(ax=ax4,title='(d) Heat vs. cold\nage standardised',x=tcost.fADcold.to_series(),y=tcost.fADheat.to_series(),xlabel='Cold-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Heat-related deaths (€ capita$^{-1}$ y$^{-1}$)')
    lims = [np.min(np.abs([ax4.get_xlim(), ax4.get_ylim()])),np.max(np.abs([ax4.get_xlim(), ax4.get_ylim()]))]
    ax4.plot([-x for x in lims],lims,'k:',label='one-to-one line')
    ax4.legend()
    
    #compare to various air pollution, one plot
    #compare to EEA air pollution, local age
    plotlinearfit(ax=ax5,title='',x=bothcosts.O3,y=bothcosts.fADcold_city,xlabel='',ylabel='',dotlabel='ozone',ls='--')
    #compare to Khomenko air pollution, standard age
    plotlinearfit(ax=ax5,title='(e) Cold vs. air pollution',x=bothcosts['Euro per capita per year'],y=bothcosts.fADcold,
                  xlabel='Air pollution-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='Cold-related deaths (€ capita$^{-1}$ y$^{-1}$)',
                  dotlabel='PM2.5 (std age)')
    
    #compare O3 and PM2.5
    plotlinearfit(ax=ax6,title='(f) PM2.5 vs. ozone pollution',x=bothcosts.O3,y=bothcosts['Euro per capita per year'],
                  xlabel='Ozone-related deaths (€ capita$^{-1}$ y$^{-1}$)',ylabel='PM2.5-related deaths (€ capita$^{-1}$ y$^{-1}$)',
                  dotlabel='')
    
    plt.tight_layout()
    plt.savefig('figures/cost_heatcoldannual_std_age.pdf')
    


if __name__ == '__main__':
    
    #UHI costs
    tcost, avgdiff_all = uhi_cost()
    tcost_df = tcost.to_dataframe().drop(columns=['crs','band','spatial_ref','age_structure'])
    #air pollution costs
    bothcosts = add_air_pollution_costs(tcost_df)
    #rent and transit costs
    all_costs = rent_and_transit()
    #combine all
    all_costs = all_costs.merge(bothcosts[['fADheat_city','Euro per capita per year','O3','copernicus']].set_index('copernicus')/12,how='outer',left_index=True,right_index=True)
    all_costs.rename(columns={'fADheat_city':'Heat\n(n=70)','Euro per capita per year':'PM2.5\nstd age\n(n=70)','O3':'Ozone\nall ages\n(n=70)'}, inplace=True)
    
    #Fig. 5
    plot_cost_comparison(avgdiff_all,tcost)
    #Fig. S16
    plot_cost_comparison_std(avgdiff_all,tcost)
    
    # values quoted in paper text
    print('cost of heat, per year:')
    print('mean:')
    print(tcost.fADheat_city.mean().values)
    print('median:')
    print(tcost.fADheat_city.median().values)
    print('IQR:')
    print(tcost.fADheat_city.quantile(0.25).values)
    print(tcost.fADheat_city.quantile(0.75).values)
    
    print('cost of cold, per year:')
    print('mean:')
    print(tcost.fADcold_city.mean().values)
    print('median:')
    print(tcost.fADcold_city.median().values)
    print('IQR:')
    print(tcost.fADcold_city.quantile(0.25).values)
    print(tcost.fADcold_city.quantile(0.75).values)
    
    print('cities with low heat risk but high cold protection:')
    print(tcost_df[(tcost_df.fADheat_city<200)&(tcost_df.fADcold_city<-500)].dropna()[['fADheat_city','fADcold_city']] )
    
    print('number of countries with net adverse impact:')
    print((tcost.fADnet_city>0).sum().values)
    
    print('ratio of heat vs pm2.5:')
    print(all_costs['Heat\n(n=70)'].median()/all_costs['PM2.5\nstd age\n(n=70)'].median())
    print('ratio of heat to o3:')
    print(all_costs['Heat\n(n=70)'].median()/all_costs['Ozone\nall ages\n(n=70)'].median())