# Output an ensemble timeseries of urban-rural difference in attributable 
# fraction, based on 1000 Monte Carlo simulated members of exposure-response 
# relationships for each city from Masselot et al. 2023. Used to capture 
# uncertainties in the exposure-response relationships.


attr_city <- function(city){
  
  library(ncdf4)
  library(abind)
  library(dlnm)
  library(mixmeta)
  source("attrdl_onebasis.R")

  # number of simulations to use
  numsims <- 1000
  # number of simulations to skip (i.e. start from simstart+1)
  simstart <- 0

  # model coefficients (1000 simulated per city-age)
  coefs <- read.csv("../data/Masselot2023coef_simu.csv")
  # model variance-covariance matrix (1 per city-age)
  vcovs <- read.csv("../data/Masselot2023vcov.csv")
  # city descriptions (contains URAU_CODES matched with city names)
  citydesc <- read.csv("../data/citydesc_v3.csv",row.names=1)
  citydesc2 <- read.csv('../data/Masselot2023metadata.csv', encoding="UTF-8")
  
  # age groups to loop through
  ages <- c('_20-45','_45-65','_65-75','_75-85','_85+')
  
  # age dimension to be added to NetCDF file
  agedim <- ncdim_def(name="age", units="age_group", vals=c(20,45,65,75,85))
  # simulation dimension to be added to NetCDF file
  simdim <- ncdim_def(name="simulation", units="", vals=(1:numsims)+simstart)
  
  # range of percentile values allowed to be the MMT
  predper <- c(seq(0.1,1,0.1), 2:98, seq(99,99.9,0.1))
  mmprange <- c(1, 99)
  inrange <- predper >= mmprange[1] & predper <= mmprange[2]
  
  # city name used by Pierre
  city_names <- read.csv("../data/copernicus_epiv3_cities_match.csv")#, encoding="UTF-8")
  city_epi <- city_names$epi[city_names$copernicus==city]
  # if (identical(city_epi,character(0))){
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
  
  # NetCDF file containing daily mean temperature
  file <- paste0("../data/copernicus_urban/tas_",gsub('-','_',gsub(' ','_',city)),"_mergetime_daymean_mask_500m_2015to2017_fAFzenodo.nc")
  # daily domain-mean temperature
  filem <- paste0("../data/copernicus_urban/tas_",gsub('-','_',gsub(' ','_',city)),"_mergetime_daymean_mask_fldmean.nc")
  # urban-rural mask (1 where rural, 0 where urban)
  filemask <- paste0("../data/copernicus_urban/urbanrural/ruralurbanmask_",gsub('-','_',gsub(' ','_',city)),"_UrbClim_v1.0_500m.nc")
  # elevation and land mask (1 where land and elevation is within 100m of the population weighted average, 0 otherwise)
  fileelemask <- paste0("../data/elevation_mask/elevation_mask_",gsub('-','_',gsub(' ','_',city)),".nc")
  
  # get temperature data
  nc <- nc_open(file,write=FALSE)
  tas <- ncvar_get(nc,'tas')
  tas <- tas-273.15  
  # daily domain average temperature
  nc2 <- nc_open(filem,write=FALSE)
  tasm <- ncvar_get(nc2,'tas')
  tasm <- tasm-273.15
  # urban-rural mask (1 where rural, 0 where urban)
  nc3 <- nc_open(filemask,write=FALSE)
  isrural <- ncvar_get(nc3,'ruralurbanmask')
  # elevation and land mask (1 where land and elevation is within 100m of the population weighted average, 0 otherwise)
  nc4 <- nc_open(fileelemask,write=FALSE)
  elemask <- ncvar_get(nc4,'elevation_and_land_mask')
  
  # define variable for output NetCDF file
  AFdiff_nc <- ncvar_def( name="fAFdiff", units="", 
                          dim=list(agedim,simdim,nc$dim[['time']]), 
                          longname="urban-rural difference in forward attributed fraction", prec="float",
                          verbose=TRUE,chunksizes=c(agedim$len,simdim$len,1))
  # define output NetCDF file
  ncout <- nc_create(paste0('../data/simulated/simulated_fAFdiff_',gsub('-','_',gsub(' ','_',city)),'_s',simstart+1,'to',simstart+numsims,'.nc'), AFdiff_nc)
  
  
  # counter for age group number
  iage <- 1

  # initialise variable to store MMT for each age group and all simulations
  MMT <- data.frame(simulation=1:numsims)
  
  for (age in ages){
    
    # define the cross/one-basis
    ob <- onebasis(c(tasm), fun = "bs", degree = 2, knots = quantile(c(tasm), c(.1, .75, .9), na.rm = T))
    attr(ob,"argvar") <- list(fun = "bs", degree = 2, knots = quantile(c(tasm), c(.1, .75, .9), na.rm = T))
    
    # variance-covariance for city-age group
    vcov <- xpndMat(vcovs[vcovs$URAU_CODE==city_code & vcovs$age==strsplit(age,'_')[[1]][2],-(1:2)])

    # to store all simulation for same age group
    MMTsims <- c()
    
    for (sim in (1:numsims)+simstart){
      # model coefficients for city-age-simulation, as matrix
      coeff <- t(coefs[coefs$URAU_CODE==city_code & coefs$age==strsplit(age,'_')[[1]][2] & coefs$sim==sim,-(1:3)])
      
      # run first prediction to get MMT for later rescaling
      cp <- crosspred(ob,coef=coeff,vcov=vcov,at=quantile(tasm,predper/100),cen=19)
      cen <- cp$predvar[inrange][which.min(cp$allfit[inrange])] 

      # collect all MMT simulated members for same age group into list
      MMTsims <- c(MMTsims,cen)
      
      # re-center % get predictions for each model centering on the MMT
      bound <- range(tasm) #range of domain daily mean T
      crall <- crosspred(ob,coef=coeff,vcov=vcov,from=bound[1],to=bound[2],by=0.1,cen=cen)      
      
      # flatten to 1D, attribute, then reshape back to 3D 
      attfrac1d <- attrdl(c(tas),ob,coef=coeff,vcov=vcov,tot=F,cen=cen,dir='forw') 
      AF <- replace(tas,TRUE,attfrac1d)
      
      # mask out grids with any water bodies and with elevation greater/less 
      # than 100m from the population-weighted average
      ## duplicate mask across time steps so it's the same size as AF
      elemask3d <- array(elemask,dim=dim(AF))
      ## remove grids that don't satisfy the above criteria from AF array
      AF[elemask3d==0] <- NA
      
      # mask by urban/rural grids
      ## duplicate urban-rural across time steps so that it's the same size as AF
      isrural3d <- array(isrural,dim=dim(AF))
      ## make a copy of AF and overwrite urban grids with NA (get rural grids only)
      AFrural <- AF
      AFrural[isrural3d==0] <- NA
      ## same as above but remove rural grids (get urban grids only)
      AFurban <- AF
      AFurban[isrural3d==1] <- NA
      
      # reshape arrays by flattening x-y dimensions into one spatial dimension
      AFrural <- array(AFrural,dim=c(nc$dim$y$len*nc$dim$x$len,nc$dim[['time']]$len))
      AFurban <- array(AFurban,dim=c(nc$dim$y$len*nc$dim$x$len,nc$dim[['time']]$len))
      
      # urban-rural difference (urban mean minus rural mean at each time step)
      AFdiff <- colMeans(AFurban,dims=1,na.rm=TRUE)-colMeans(AFrural,dims=1,na.rm=TRUE)
      
      # write to NetCDF file
      ncvar_put(ncout, AFdiff_nc, AFdiff, start=c(iage,sim-simstart,1), count=c(1,1,-1), verbose=FALSE)
      
    } # loop through simulations

    # write all simulations for age group to dataframe
    MMT[strsplit(age,'_')[[1]][2]] <- MMTsims
    
    # advance age group counter
    iage <- iage+1
    
  } # loop through age groups

  #save MMT to file
  write.csv(MMT,
            file=paste0('../data/simulated/MMT/simulated_MMT_urbclimT_',gsub('-','_',gsub(' ','_',city)),'_s',simstart+1,'to',simstart+numsims,'.csv'))
  
  nc_close(nc)
  nc_close(nc2)
  nc_close(nc3)
  nc_close(nc4)
  nc_close(ncout)
    
}


library(parallel)
library(doParallel)

# cities with both urbclim and epi data
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

# number of cores available for parallel computing
numCores <- detectCores()

# make cluster on all available cores (or minus one to leave some processing space for other tasks)
cl <- makeCluster(numCores-1)
registerDoParallel(cl)

# call function defined at top to compute in parallel
parLapply(cl,cities,fun=attr_city)

# turn parallel processing off and run sequentially again:
registerDoSEQ()


