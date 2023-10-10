# Apply exposure-response relationships to UrbClim temperature
# spatial time series data to get spatial time series of
# forward attributable fraction (fAF, i.e. lag-cumulative)

# Attribute fAF for spatial time series
# by reordering order of matrix dimensions and truncating 
# lon-lat-time to a single vector (array) for attribution

# use age-separated risk relationships from Masselot et al. 2023
# output NetCDF with additional age dimension


library(ncdf4)
library(abind)
library(dlnm)
library(mixmeta)
source("attrdl_onebasis.R")

#cities with both urbclim and epi data
cities <- c('Amsterdam','Antwerp','Athens','Barcelona','Bari','Basel','Berlin',
            'Bilbao','Bologna','Bordeaux','Brasov','Bratislava','Brussels','Bucharest',
            'Budapest','Charleroi','Cluj-Napoca','Cologne','Copenhagen','Debrecen',
            'Dublin','Dusseldorf','Edinburgh','Frankfurt am Main','Gdansk','Geneva',
            'Genoa','Ghent','Glasgow','Graz','Gyor','Hamburg','Helsinki',
            'Klaipeda','Kosice','Krakow','Leeds','Leipzig','Liege','Lille','Lisbon',
            'Ljubljana','London','Luxembourg','Lyon','Madrid','Malaga','Marseille',
            'Milan','Miskolc','Montpellier','Munich','Murcia','Nantes','Naples',
            'Nice','Oslo','Palermo','Paris','Pecs','Porto','Prague','Riga','Rome',
            'Rotterdam','Sevilla','Sofia','Split','Stockholm','Strasbourg','Szeged',
            'Tallinn','Thessaloniki','Toulouse','Trieste','Turin','Utrecht',
            'Valencia','Varna','Vienna','Vilnius','Warsaw','Wroclaw','Zagreb','Zurich')

# model coefficients (1 per city-age)
coefs <- read.csv("../data/Masselot2023coefs.csv")
# model variance-covariance matrix (1 per city-age)
vcovs <- read.csv("../data/Masselot2023vcov.csv")
# city descriptions (contains URAU_CODES matched with city names)
citydesc <- read.csv("../data/citydesc_v3.csv",row.names=1)
citydesc2 <- read.csv('../data/Masselot2023metadata.csv', encoding="UTF-8") 

#age groups to loop through
ages <- c('_20-45','_45-65','_65-75','_75-85','_85+')
#age dimension to be added to NetCDF file
agedim<-ncdim_def(name="age", units="age_group", vals=c(20,45,65,75,85))

# range of percentile values allowed to be the MMT
predper <- c(seq(0.1,1,0.1), 2:98, seq(99,99.9,0.1))
mmprange <- c(1, 99)
inrange <- predper >= mmprange[1] & predper <= mmprange[2]


for (city in cities){
  
  # city name used by Pierre
  city_names <- read.csv("../data/copernicus_epiv3_cities_match.csv")#, encoding="UTF-8")
  city_epi <- city_names$epi[city_names$copernicus==city]
  if (identical(city_epi,"")){
    stop(paste0("epi model not available for ",city))
  }
  # city URAU code
  city_code <- citydesc[citydesc$LABEL==city_epi,][1,1]
  # first 5 characters of the code, then ending with 'C' 
  city_code <- paste0(substr(city_code, 1, 5),'C')
  # special case for Klaipeda
  if (city_epi == 'Klaipeda'){
    city_code <- citydesc2[startsWith(citydesc2$LABEL,'Klaip'),][1,1]
  }
  # if city_code from above is not in Pierre's list, try to get it from Pierre's metadata using the city name
  if (!(city_code %in% vcovs$URAU_CODE)){
    city_code <- citydesc2[citydesc2$LABEL==city_epi,][1,1]
    if (!(city_code %in% vcovs$URAU_CODE)){
      city_code <- citydesc2[(citydesc2$cityname==city_epi)&(!is.na(citydesc2$cityname)),][1,1]
    }
  }
  
  # NetCDF data containing daily mean temperature to be modified with
  # additional variables containing fAF
  file <- paste0("../data/copernicus_urban/tas_",gsub('-','_',gsub(' ','_',city)),"_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc")
  filem <- paste0("../data/copernicus_urban/tas_",gsub('-','_',gsub(' ','_',city)),"_mergetime_daymean_mask_fldmean.nc")
  
  # get temperature data
  nc<-nc_open(file,write=TRUE)
  tas<-ncvar_get(nc,'tas')
  tas<-tas-273.15  
  # daily field average temperature
  nc2<-nc_open(filem,write=FALSE)
  tasm<-ncvar_get(nc2,'tas')
  tasm<-tasm-273.15
  
  # define new variable for NetCDF file
  attfrac_nc<-ncvar_def( name=paste0("fAF"), units="", 
                         dim=list(nc$dim[['x']],nc$dim[['y']],nc$dim[['time']],agedim), 
                         longname="forward attributed fraction", prec="float",
                         verbose=TRUE,chunksizes=c(nc$dim[['x']]$len,nc$dim[['y']]$len,1,1))
  nc<-ncvar_add( nc, attfrac_nc, verbose=TRUE)
  
  #counter for age group number
  iage<-1
  
  #initialise variable to store MMT for each age group
  MMT <- c()
  
  for (age in ages){
    
    # DEFINE THE CROSS/ONE-BASES
    ob <- onebasis(c(tasm), fun = "bs", degree = 2, knots = quantile(c(tasm), c(.1, .75, .9), na.rm = T))
    attr(ob,"argvar") <- list(fun = "bs", degree = 2, knots = quantile(c(tasm), c(.1, .75, .9), na.rm = T))
    
    # variance-covariance for city-age group
    vcov <- xpndMat(vcovs[vcovs$URAU_CODE==city_code & vcovs$age==strsplit(age,'_')[[1]][2],-(1:2)])
    
    # model coefficients for city-age-simulation, as matrix
    coeff <- t(coefs[coefs$URAU_CODE==city_code & coefs$age==strsplit(age,'_')[[1]][2],-(1:2)])
    
    ##############################  
    # run first prediction to get MMT for later rescaling
    cp <- crosspred(ob,coef=coeff,vcov=vcov,at=quantile(tasm,predper/100),cen=19)
    cen <- cp$predvar[inrange][which.min(cp$allfit[inrange])] 
    MMT <- c(MMT,cen)
    
    # re-center % get predictions for each model centering on the MMT
    bound <- range(tasm) #range of domain daily mean T
    crall <- crosspred(ob,coef=coeff,vcov=vcov,from=bound[1],to=bound[2],by=0.1,cen=cen)
    ##################################
    
    # flatten to 1D, attribute, then reshape back to 3D 
    attfrac1d <- attrdl(c(tas),ob,coef=coeff,vcov=vcov,tot=F,cen=cen,dir='forw')
    AF <- replace(tas,TRUE,attfrac1d)
    
    # write to NetCDF file
    ncvar_put(nc, attfrac_nc, AF, start=c(1,1,1,iage), count=c(-1,-1,-1,1), verbose=TRUE)
    
    #save cumulative RR fit to file
    write.csv(lapply(crall[c("allfit","allhigh","alllow")],exp),file=paste0('../data/RRfit/RRfit_urbclimT_',gsub('-','_',gsub(' ','_',city)),'_',strsplit(age,'_')[[1]][2],'.csv'))
    
    #advance age group counter
    iage<-iage+1
    
  } #loop through age groups
  
  
  #save optimal T to file
  write.csv(setNames(data.frame(gsub("_","",ages),MMT),c("Age","MMT")),file=paste0('../data/RRfit/MMT_urbclimT_',gsub('-','_',gsub(' ','_',city)),'.csv'))
  
  nc_close(nc)
  nc_close(nc2)
  
} #loop through cities
