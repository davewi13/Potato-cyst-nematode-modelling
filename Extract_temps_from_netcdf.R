library(ncdf4)
library(fields)
library(sp)
library(rgdal)
library(abind)
library(stringr)
library(SpatialEpi)

setwd("C:/Users/dewing/PCN_climate_data/")

extract_temps <- function(lat, lon, start_time, end_time, farm_name){
  # Take desired latitude and longitude and start and end times and construct temperature time series by stitching together the minimum and maximum temperatures from the Hadley Centre model output for the relevant 1km square.
  # Inputs:   lat, lon - latitude and longitude
  #           start_time, end_time - format "yyyymmdd" for consistency with Hadley input files.  Code won't work if these span multiple years or include february (naff, I know but not relevant for what I'm doing)
  # Outputs:  .csv file with temperature time series
  
  # Create a space to store the temperature outputs later
  fulltemps <- NULL
  
  # Create vector with filenames for required files
  base_string_min <- "tasmin_hadukgrid_uk_1km_day_"
  base_string_max <- "tasmax_hadukgrid_uk_1km_day_"
  year <- substr(start_time, 1, 4)
  start_month <- substr(start_time, 5, 6)
  end_month <- substr(end_time, 5, 6)
  start_day <- substr(start_time, 7, 8)
  end_day <- substr(end_time, 7, 8)
  
  monthvec <- str_pad(as.character(seq(as.numeric(start_month), as.numeric(end_month))), 2, pad="0")
  
  for(month in monthvec){
    if(month %in% c("04", "06", "09", "11")){
      daysinmonth <- "30"
    } else {
      daysinmonth <- "31"
    }
    filename_min <- paste0(base_string_min, year, month, "01-", year, month, daysinmonth, ".nc")
    filename_max <- paste0(base_string_max, year, month, "01-", year, month, daysinmonth, ".nc")
    # Finally open the file
    current_min <- nc_open(filename_min)
    current_max <- nc_open(filename_max)
    
    # Find correct grid square
    farm.coords.latlon <- SpatialPoints(matrix(c(lon, lat), ncol=2), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    farm.coords.tmerc <- spTransform(farm.coords.latlon, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000"))

    lonvec <- current_min$var$longitude$dim[[1]]$vals
    latvec <- current_min$var$longitude$dim[[2]]$vals
    LonIdx <- which.min(abs(lonvec - farm.coords.tmerc@coords[1]))
    LatIdx <- which.min(abs(latvec - farm.coords.tmerc@coords[2]))
    
    # Extract temperatures from correct grid square
    mintemps <- ncvar_get(current_min, "tasmin")[LonIdx, LatIdx,]
    maxtemps <- ncvar_get(current_max, "tasmax")[LonIdx, LatIdx,]
    
    # If it is the first or last month then remove days before start or after end
    if(month == monthvec[1]){
      start_day <- as.numeric(start_day)
      TimeIdx <- seq(start_day, length(mintemps))
      mintemps <- mintemps[TimeIdx]
      maxtemps <- maxtemps[TimeIdx]
    } else if(month == tail(monthvec,1)) {
      end_day <- as.numeric(end_day)
      TimeIdx <- seq(1, end_day)
      mintemps <- mintemps[TimeIdx]
      maxtemps <- maxtemps[TimeIdx]
    }
    
    # Combine min and max temperatures into one series
    temps <- numeric(2*length(mintemps))
    for(i in 1:length(mintemps)){
      temps[(2*i)-1] <- mintemps[i]
      temps[2*i] <- maxtemps[i]
    }
    
    # The above will be overwritten each month.  Append this to previous months
    fulltemps <- c(fulltemps, temps)
  }
  temps_df <- data.frame(Time = seq(0, length(fulltemps)/2-0.5, 0.5), Temperature=fulltemps)
  write.csv(temps_df, file = paste0("//212.219.57.223/HomeDirVol/dewing/My_files/Research/PCN/Model_code/Data/",farm_name, year, ".csv"))
}

# Roseisle, Moray, 2012
extract_temps(lat=57.684010, lon=-3.439360, start_time = "20120517", end_time = "20120912", farm_name="Roseisle")
# Roseisle, Moray, 2015
extract_temps(lat=57.684010, lon=-3.439360, start_time = "20150527", end_time = "20151020", farm_name="Roseisle")
# Rothes, Moray
extract_temps(lat=57.53927778, lon=-3.16777778, start_time = "20130611", end_time = "20130905", farm_name="Rothes")
# Duffus, Moray
extract_temps(lat=57.69837222, lon=-3.36305556, start_time = "20160620", end_time = "20161117", farm_name="Duffus")

# Derachie Farm, Angus, 2011
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20110519", end_time = "20111031", farm_name="Derachie")
# Derachie Farm, Angus, 2012
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20120518", end_time = "20120921", farm_name="Derachie")
# Derachie Farm, Angus, 2013
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20130521", end_time = "20130916", farm_name="Derachie")
# Derachie Farm, Angus, 2014
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20140430", end_time = "20140919", farm_name="Derachie")
# Derachie Farm, Angus, 2015
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20150517", end_time = "20150829", farm_name="Derachie")
# Derachie Farm, Angus, 2016
extract_temps(lat=56.747541, lon=-2.903272, start_time = "20160519", end_time = "20160829", farm_name="Derachie")

# Balruddery Farm, Angus, 2011
extract_temps(lat=56.4723, lon=-3.1148, start_time = "20110518", end_time = "20111031", farm_name="Balruddery")
# Balruddery Farm, Angus, 2012
extract_temps(lat=56.4723, lon=-3.1148, start_time = "20120517", end_time = "20120923", farm_name="Balruddery")
# Balruddery Farm, Angus, 2015
extract_temps(lat=56.4723, lon=-3.1148, start_time = "20150508", end_time = "20151001", farm_name="Balruddery")

# Luffness Mains, East Lothian, 2010
extract_temps(lat=56.00836, lon=-2.832121, start_time = "20100527", end_time = "20101001", farm_name="Luffness")
# Luffness Mains, East Lothian, 2011
extract_temps(lat=56.00836, lon=-2.832121, start_time = "20110504", end_time = "20110930", farm_name="Luffness")

# Tranent, East Lothian, 2014
extract_temps(lat=55.944569, lon=-2.954030, start_time = "20140428", end_time = "20141008", farm_name="Tranent")

# Dalrymple, Ayrshire, 2013
extract_temps(lat=55.397310, lon=-4.592270, start_time = "20130604", end_time = "20131010", farm_name="Dalrymple")
# Dalrymple, Ayrshire, 2014
extract_temps(lat=55.397310, lon=-4.592270, start_time = "20140612", end_time = "20141013", farm_name="Dalrymple")
# Dalrymple, Ayrshire, 2015
extract_temps(lat=55.397310, lon=-4.592270, start_time = "20150610", end_time = "20151005", farm_name="Dalrymple")
# Dalrymple, Ayrshire, 2016
extract_temps(lat=55.397310, lon=-4.592270, start_time = "20160524", end_time = "20160928", farm_name="Dalrymple")

# Bold Farm, Lancashire, 2011
extract_temps(lat=53.579477, lon=-2.914892, start_time = "20110531", end_time = "20110926", farm_name="Bold Farm")

# Howell, Lincolnshire, 2013
extract_temps(lat=53.000839, lon=-0.289299, start_time = "20130412", end_time = "20130920", farm_name="Howell")

# Spalding, Lincolnshire, 2012
extract_temps(lat=52.818844, lon=-0.107189, start_time = "20120503", end_time = "20120813", farm_name="Spalding")

# Holbeach, Lincolnshire, 2014
extract_temps(lat=52.801887, lon=0.0717647, start_time = "20140411", end_time = "20140905", farm_name="Holbeach")

# Harper Adams, Shropshire, 2011
extract_temps(lat=52.7797, lon=-2.4275, start_time = "20110421", end_time = "20110921", farm_name="Harper Adams")
# Harper Adams, Shropshire, 2012
extract_temps(lat=52.7797, lon=-2.4275, start_time = "20120402", end_time = "20121004", farm_name="Harper Adams")
# Harper Adams, Shropshire, 2014
extract_temps(lat=52.7797, lon=-2.4275, start_time = "20140430", end_time = "20141121", farm_name="Harper Adams")

# Aylsham, Norfolk, 2014
extract_temps(lat=52.7770, lon=-1.3103, start_time = "20140331", end_time = "20140922", farm_name="Aylsham")

# Leominster, Herefordshire, 2011
extract_temps(lat=52.305340, lon=-2.818770, start_time = "20110419", end_time = "20110808", farm_name="Leominster")

# NIAB, Cambridgeshire, 2014
extract_temps(lat=52.2203, lon=-0.0950, start_time = "20140416", end_time = "20140908", farm_name="NIAB")

# Kings Caple, Herefordshire, 2011
extract_temps(lat=51.959043, lon=-2.634648, start_time = "20110414", end_time = "20110823", farm_name="Kings Caple")

# Aythorpe, Essex, 2014
extract_temps(lat=51.808528, lon=-0.289972, start_time = "20140404", end_time = "20140902", farm_name="Aythorpe")

# Tetbury, Gloucestershire, 2011
extract_temps(lat=51.639190, lon=-2.161870, start_time = "20110414", end_time = "20110802", farm_name="Tetbury")

# Bold Farm, Lancashire, 2011
extract_temps(lat=53.579477, lon=-2.914892, start_time = "20110531", end_time = "20110926", farm_name="Bold Farm")

# NIAB, Cambridgeshire, 2014
extract_temps(lat=52.2203, lon=-0.0950, start_time = "20140416", end_time = "20140908", farm_name="NIAB")



# Load data files.  Can be found at address below.  In this case I have loaded 2070-2080 for climate run 13.
# Separate netcdf files for minimum and maximum daily temperature which I want to combine into a series of
# .csv files for each grid square because that's what the WNV model code is set up to work with.
#Citable as:  Met Office Hadley Centre (2018): UKCP18 Regional Projections on a 12km grid over the UK for 1980-2080.
#Centre for Environmental Data Analysis, 15/01/2021. 
#https://catalogue.ceda.ac.uk/uuid/589211abeb844070a95d061c8cc7f604
min_file <- nc_open("tasmin_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")
max_file <- nc_open("tasmax_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")

# Extract temperature variables
temp_min <- ncvar_get(min_file,"tasmin")
temp_max <- ncvar_get(max_file,"tasmax")

# Extract longitude and latitude.  There is some inconsistency in the way different versions of the climate files
# are presented so the commented out code may be required.
# #lon <- min_file$var$grid_longitude$dim[[1]]$vals
# #lat <- min_file$var$grid_longitude$dim[[2]]$vals
lon <- min_file$var$longitude$dim[[1]]$vals
lat <- min_file$var$longitude$dim[[2]]$vals