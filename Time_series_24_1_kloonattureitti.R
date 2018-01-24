library(rmatio)
library(raster)
library(lubridate)
library(ggplot2)
library(hms)
library(scales)
library(grid)
library(gridExtra)
library(quantreg)
#Clear objects from workspace
rm(list = ls())

#Script to extract chosen wave and wind parameters for chosen ship route and plotting timeseries
#Errors when running script most likely are result of chosen ship not having any locations in AIS data in the chosen timeframe
#c. Sakari Äärilä

####################################################################################################################

#INPUTS 

mmsi = 265611110  # Ship MMSI (Maritime Mobile Service Identity)

# https://www.marinetraffic.com/en/ais/index/ships/all   for MMSI search
#Examples
#Silja Serenade 230184000 Helsinki-Mariehamn-Stockholm
#Silja Europa 276807000 Helsinki-Tallinn 


startdate = 20121130# Date (YYYYMMDD)
starthour = 15       # Hour of day (0-23)
duration =  3      # Duration in hours (1-72)

#Set parameters: "VHM0" = Hs, "VTPK" = Tp, "VPED" = pdr, "VTM01" = Tm01, "VTM10" = Tm10, "TpI" = Tp .
parameters = c("VHM0", "VTPK", "VPED", "VTM01", "VTM10", "TpI")   #Give inputs like this: c("X", "Y", "Z")   

#Set paths to data and output directories
route_file_path = "F:/output/" # AIS data directory.  #Should be .mat files
ncdfpath = "F:/2012_waves/"    # NetCDF wavemodel data directory  #The data should be stored one month per file with a band for every hour.
output_dir =  "F:/output/"    # Output directory

#####################################################################################################################

#FUNCTIONS


#Function to choose ship route based on MMSI, date, time and duration
chooseroute = function(startdate, starthour, duration, mmsi) {
  
  #Create list of files in the AIS data folder
  aislist = list.files(path=route_file_path, pattern = '\\.mat$')
  
  #Create empty dataframe for the routes
  DFroute = data.frame()
  
  #Loop through the route files and add the points specified in inputs to the empty dataframe
  for (i in aislist) {
    file_name = i
    print(i)
    file_date <- strsplit(file_name,"_")[[1]][2]
    
    #Check if filename matches the chosen input date, if so create dataframe from the file
    if (as.numeric(file_date) == as.numeric(startdate))  {
      read = read.mat(paste(route_file_path, i, sep=""))
      print(paste("read .mat", file_name))
      DF = do.call(rbind.data.frame, read )
      DF = unique(DF) #Remove duplicates
      print("Transformed into DF")
      
      # Pick all rows with input MMSI into new dataframe
      DFpoint = DF[DF$V1 %in% mmsi,] 
      #if (nrow(DFpoint) < 1) {next}   #??????
      DFpoint$date = file_date # Add date column 
      
      #Calculate duration from starthour and time column
      DFpoint$travel_dur = round(DFpoint$V2 - starthour, 4)    
      DFpoint = DFpoint[DFpoint$travel_dur >= 0, ]
      if (nrow(DFpoint) < 1) {next} #If no matching MMSI. Skip to next file
      
      #Loop over the dataframe and bind rows to DFroute if they match the chosen time frame
      for (z in 1:nrow(DFpoint)) {
        if (DFpoint[z, 2] >= starthour && DFpoint[z, 2] <= starthour + duration ) {
          DFroute = rbind(DFroute, DFpoint[z,])}} 
      
      
      
      #If the route falls on 2 separate days, loop over the next day file. Same procedure as for 1 day route. 
    } else if (starthour + duration > 24 && as.numeric(file_date) == as.numeric(startdate) + 1) {
      read = read.mat(paste(route_file_path, i, sep=""))
      print("read .mat")
      DF = do.call(rbind.data.frame, read )
      print("Transform into DF")
      DFpoint = DF[DF$V1 %in% mmsi,] 
      
      DFpoint$date = file_date
      DFpoint$travel_dur = round(24 - starthour + DFpoint$V2, 4)
      DFpoint = DFpoint[DFpoint$travel_dur < duration, ]
      if (nrow(DFpoint) < 1) {next}
      for (z in 1:nrow(DFpoint)) {
        if (DFpoint[z, 2] <= starthour+duration - 24 ) {
          DFroute = rbind(DFroute, DFpoint[z,])}}
      
      
      
      #If the route falls on 3 separate days, loop over the next day file. Same procedure as for 1 and 2 day route.
    } else if (starthour + duration > 48 && as.numeric(file_date) == as.numeric(startdate) + 2) {
      read = read.mat(paste(route_file_path, i, sep=""))
      print("read .mat")
      DF = do.call(rbind.data.frame, read )
      
      print("Transform into DF")
      DFpoint = DF[DF$V1 %in% mmsi,] 
      
      DFpoint$date = file_date
      DFpoint$travel_dur = round(48 - starthour + DFpoint$V2, 4)
      DFpoint = DFpoint[DFpoint$travel_dur < duration, ]
      if (nrow(DFpoint) < 1) {next}
      for (z in 1:nrow(DFpoint)) {
        if (DFpoint[z, 2] <= starthour+duration - 48 ) {
          DFroute = rbind(DFroute, DFpoint[z,])}}
      
    }}
  if (nrow(DFroute) > 0) {
  #Name the columns  #KATO ET TOIMII UUSIEN NIMIEN KANSSA            
  DFroute = DFroute[ -c(6:7) ]
  colnames(DFroute) = c("mmsi", "time", "velocity", "longitude", "latitude", "VHM0", "VTPK", "VPED", "VTM01", "VTM10", "TpI", "date", "duration")
  
  
  #Create new columns for raster name, band number and dummytimedate(used for timeseries)
  DFroute["raster"] = NA
  DFroute["band"] = NA
  DFroute["dummytimedate"] = ymd_hms(paste(DFroute$date, hms(hours=DFroute$time)))
  
  #Remove row if parameter value is NA
  DFroute = DFroute[complete.cases(DFroute[,6:11]),]
}
  return(DFroute)
}

# Function to plot the chosen route. # Its ugly and wrong CRS but main point is to see if the desired cities are connected by the route before continuing
plot_route_func = function(DFroute) {
  
  #Coordinates for Helsinki, Tallinn, Mariehamn and Stockholm
  dfcity <- data.frame(longitude=c(24.9598, 24.7726, 19.9269, 18.1970),
                       latitude=c(60.1613, 59.4450, 60.0921, 59.3347))
  #Name the rows
  row.names(dfcity) = c("Helsinki", "Tallinn", "Mariehamn", "Stockholm" )
  
  #Pick the first and last point of the route
  DurationMin = DFroute[which.min(DFroute$duration),]
  DurationMax = DFroute[which.max(DFroute$duration),]
  
  #Create plot of the route 
  p = ggplot(DFroute, aes(longitude, latitude)) +geom_point()+
    geom_point(data = dfcity)+ geom_text(data=dfcity, label=row.names(dfcity))+
    geom_point(data = DurationMin, color = "red")+
    geom_point(data = DurationMax, color = "blue")+
    ggtitle(paste("MMSI:", mmsi, "Date:", startdate, "\n Duration:", duration, "hours | (Plot not in right CRS)", "\n Red = Start | Blue = End"))
  graphics.off()
  #Return plot
  return(p) 
}

#Function to create identical route for everyday of month and Create time/date stamps for each route.
timeseries_month_func = function(startdate, DFroute  ) {
  
  #the function adds the route for each day of month to DFroute with corresponding timestamps
  #The route itself is the same, only time and date changes
  
  #Set startdate to first of the month and stopdate to last of the month
  startdateM = format(floor_date(ymd(startdate), "month"), "%Y%m%d")
  stopdate = format(ceiling_date(ymd(startdate), "month")-days(1), "%Y%m%d")
  dates <- format(seq(ymd(startdateM), ymd(stopdate), by=1), "%Y%m%d")
  DFmonth = data.frame()
  DFroute <- subset(DFroute, select = -c(VHM0, VTPK, VPED, VTM01, VTM10, TpI) )
  #Start looping the dates 
  for (a in 1:length(dates)) {  
    
    d = dates[a]
    DFday = DFroute
    #Give new time/date stamps for points
    DFday$date = d
    DFday$time = starthour + DFday$duration
    DFday$routeno = a
    DFmonth = rbind(DFmonth, DFday)
    
    print(paste0("date ", d, " / ", stopdate, " added to DFroute" ))
  }
  # Check that change of day (time > 24, 48, 72) adds  day to $date and sets $time back to 00
  DFmonth$date = ifelse(DFmonth$time>=72, format(ymd(DFmonth$date)+days(3), "%Y%m%d") , DFmonth$date)
  DFmonth$time = ifelse(DFmonth$time >=72, starthour-72+DFmonth$duration, DFmonth$time)
  
  DFmonth$date = ifelse(DFmonth$time>=48, format(ymd(DFmonth$date)+days(2), "%Y%m%d") , DFmonth$date)
  DFmonth$time = ifelse(DFmonth$time >=48, starthour-48+DFmonth$duration, DFmonth$time)
  
  DFmonth$date = ifelse(DFmonth$time>=24, format(ymd(DFmonth$date)+days(1), "%Y%m%d") , DFmonth$date)
  DFmonth$time = ifelse(DFmonth$time >=24, starthour-24+DFmonth$duration, DFmonth$time)
  
  #Give timedate values from $date and $time
  # DFmonth["timedate"] = ymd_hms(paste(DFmonth$date, hms(hours=DFmonth$time)))  
  
  #Now we have a dataframe with the same route 30 times with timestamps for each day. 
  
  #Lets fill raster and band columns
  #Check the raster for startdate and see how many bands(layers) it has
  
  nband = nbands(raster(paste0(ncdfpath,"WAVE", substr(startdate, 1, 6), "0100.nc") ))
  
  #Fill the raster column according to date
  DFmonth$raster = paste0("WAVE", substr(DFmonth$date, 1, 6), "0100.nc")
  
  #Fill the band column according to date and time
  DFmonth$band = 24 * (day(ymd(DFmonth$date))-1) + round(as.numeric(DFmonth$time))
  DFmonth$band<-ifelse(DFmonth$band<1,1, DFmonth$band)
  DFmonth$band<- ifelse(DFmonth$band > nband, nband, DFmonth$band)
  #Jos band > nbands -> Käytetään viimeistä bandia. Tämä johtaa siihen, että kuun viimeinen puoli tuntia poimitaan edellisen tasatunnin kohdalta
  print(paste0("Raster/band columns filled"))  
  
  
  #Return is a dataframe with the route multiplied by days in the month with raster and band(layer) columns filled to know which raster and layer to use for extract values  
  return(DFmonth)
}    

#Function to extract the wave height from Netcdf. 
extractparametersMONTH = function(DFpoint, parameters) {
  
  
  #Choose right raster/netcdf taking the first row's raster column. Works since all will be using same raster file
  r = DFpoint$raster[1]
  
  #Loop through parameters
  for (para in parameters) {
    
    #Create rasterstack from netcdf file  
    STACK <- stack(paste0(ncdfpath, r), varname = para)
    
    #Create matrix from rasterstack # This takes some time
    print("Creating matrix. Please Wait.")
    mat = as.matrix(STACK)
    
    #Create empty dataframe 
    DFfilled= data.frame()
    
    # Loop through band numbers in data and extract values to df
    for (layer in unique(DFpoint$band)){
      
      #Select rows with certain band number     
      DFlayer = DFpoint[DFpoint$band == layer,]
      
      #Get coordinates from those rows
      dfcoords = data.frame( "lon" = DFlayer$longitude, "lat" = DFlayer$latitude)
      #Create SpatialPoints object from coordinates of the points
      points = SpatialPoints(dfcoords)
      
      #Create ID raster to be used for extracting
      ID_Raster <- raster(STACK[[1]])
      #Fill it with integers from 1 to the maximum number of cells.
      ID_Raster[] <- 1:ncell(STACK[[1]])
      
      #Use the extract on this raster to identify the correct cell and then extract the corresponding values from the matrix
      ext_ID <- extract(ID_Raster, points)
      DFlayer[para] <- mat[as.numeric(ext_ID), layer]
      
      #Bind rows to dataframe
      DFfilled = rbind(DFfilled, DFlayer)
      print(paste("Band", layer, "done"))
      
    }
    #Bind the parameter column to original dataframe 
    DFpoint = cbind(DFpoint, DFfilled[para])
    print(paste("Values extracted for", para, "| Next parameter/Day..."))
    
    #Give TpI the right value: 1/Tpi
    if (para == "TpI") {DFpoint["TpI"] = 1/DFpoint["TpI"]}
  }
  
  return(DFpoint)
}

# Function to create output CSV and PDF for timeserie data and plots. Pitäskö lisätä dummytimedate?
outputDF_func = function(parameters, startdate, mmsi, DFfinal_month, DF1){
  
  #Function to create multiplot with shared legend. 
  grid_arrange_shared_legend = function(...) {
    
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
      do.call(arrangeGrob, lapply(plots, function(x)
        x + theme(legend.position="none"))),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight))
  } 
  # List for plots to be added to grid_arrange multiplot
  plist = list() 
  
  #List dates in the month of startdate by first finding first and last day of the month, then listing dates between them
  startdateM = format(floor_date(ymd(startdate), "month"), "%Y%m%d")
  stopdate = format(ceiling_date(ymd(startdate), "month")-days(1), "%Y%m%d")
  #Create object including the dates
  dates <- format(seq(ymd(startdateM), ymd(stopdate), by=1), "%Y%m%d")
  
  
  #Loop through chosen parameters 
  for (i in 1:length(parameters)) 
    local({
      i = i 
      para =  parameters[i]
      
      #Clean DFfinal_month and choose only points of one day (unique) 
      DFts = unique(subset(DFfinal_month, select = c(mmsi, time, velocity, longitude, latitude, dummytimedate) )) 
      
      #Loop over dates and add a column of values for each day in month
      for (d in 1:length(dates)){
        
        #Temporary dataframe for day data
        DFday = DFfinal_month[DFfinal_month$routeno %in% d,]   
        
        #And bind parameter column to the final dataframe
        DFts = cbind(DFts, DFday[,para])
        
        #Name the column by the day of the month
        x = paste0("startdate_", d)
        names(DFts)[[6+d]] = x
        
      }
      
      print("Calculating quantiles... Please wait.")
      #Loop over rows and calculate 0.90 and 0.60 quantiles for parameter
      for (j in 1:nrow(DFts)){  
        
        DFts$q60[j] = quantile(DFts[j, c(7:ncol(DFts))], probs = 0.60, na.rm = T)
        DFts$q90[j] = quantile(DFts[j, c(7:ncol(DFts))], probs = 0.90, na.rm = T)
      }
      
      #Write CSV to output folder 
      write.csv2(DFts, paste0(output_dir, "fullmonth_", mmsi,"_",startdate, "_" ,para, ".csv"), row.names = F)
      
      #Next the Time-Series plot
      
      #The plot itself. Includes The chosen day, month min/max and month average and quantiles 0.90 and 0.60
      p = ggplot(DFfinal_month, aes(dummytimedate, DFfinal_month[para]))+
        geom_line(aes(color = "Month Min/Max"), alpha = 0.4)+
        stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), aes(color="SD +-1"), alpha=0.05)+
        geom_line(data=DFts,linetype = 2, aes(color = "q0.90", y = DFts$q90))+
        geom_line(data=DFts, linetype = 2, aes(color = "q0.60", y = DFts$q60))+
        geom_line(data=DF1 ,aes(color= "Day", y = DF1[para]))+
        stat_summary(fun.y="mean", geom = "line", size = 1, aes(color = "Month average"))+
        
        labs(colour="Legend",x="Time",y=para)+
        scale_x_datetime(name = "Time", labels = date_format("%H"), breaks = date_breaks("1 hours") )+
        ggtitle(paste(para, "MMSI:", mmsi, "\n Date:", startdate))+
        theme(legend.position = c(1,1),legend.justification = c(1, 1))+
        scale_color_manual(values = c("red","gray40", "gray85","royalblue1", "royalblue4","gray70"))+
        theme_minimal()
      
      
      
      # Using <<- Instead of local assignment to not overwrite plot every time.
      plist[[i]] <<- p
      
      
    })
  
  #Call the shared_legend plot function, insert plist as arguments
  pAll = do.call(grid_arrange_shared_legend, plist)
  
  #Save_the_plot
  ggsave(filename= paste0(output_dir, "timeseries_", mmsi, "_", startdate, ".pdf" ), plot = pAll, device = "pdf",
         scale = 1, width = 7, height = 7, units = "in",
         dpi = 300, limitsize = TRUE)
  
}



#####################################################################################################################
#RUN FUNCTIONS
#Silence warnings
options(warn = -1)


# Choose the route from .mat files and create dataframe 
DF1 = chooseroute(startdate, starthour, duration, mmsi)  #If prints error. No such MMSI in the chosen timeframe
print("Function chooseroute done!")

#Delete every second row to make things faster. Good for long routes. Mark lines with # to not do this 
  #toDelete <- seq(1, nrow(DF1), 2)
  #DF1  = DF1[ toDelete ,]

#Delete every second row to make things faster. Good for long routes. Mark lines with # to not do this
  #toDelete <- seq(1, nrow(DF1), 2)
  #DF1 = DF1[ toDelete ,]


#Plot the desired route  #Check that the route is ok before proceeding to next!    
# Route plots rarely go up until Stockholm and also there are no raster values close to Stockholm so the points get deleted as NAs.
plot_route_func(DF1)  

# And same procedure for the whole month and add info about the right raster and band to use for extract for each row
DFmonth = timeseries_month_func(startdate, DF1)

#Extract chosen parameters 
print("Extracting values from rasters")
DFfinal_month = extractparametersMONTH(DFmonth, parameters)


DFfinal_month  =DFfinal_month[complete.cases(DFfinal_month), ]

#Create own .CSV for each parameter and create the time-series plots 
print("Now creating outputs...")
outputDF_func(parameters, startdate, mmsi, DFfinal_month, DF1)


options(warn = 0)   # Warnings back on
##########################################################################################
#TESTAILUA

