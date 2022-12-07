library(Rfa)
library(belgium )
library(ncdf4,lib.loc="/mnt/netapp/group/mealadin/NORI/software/R/3.6.2/lib64/R/library")
library(RNetCDF)
source("functions.R")
# - easy leading zeros
lz=function(x,n=2) formatC(x,flag="0",format="d",width=n)


ref ="CONTROL_SUM"
model="AR13"
year=2021
month=7
day=14
starth=0
basedir="/home/idehmous/Desktop/emaddc_stat/exps"
id=""

files   =paste0(basedir,"/",ref,"/",year ,lz(month),lz(day),lz(starth),"/","PF",model,"be13b_l+",lz(seq(1,12,1),4))

# GET THE DATES ( UNIX FORMAT )
validdates=lapply(files ,function(fa) attr(FAopen(fa),"time")$validdate)
dates=c()
for (i in seq_along(validdates)){
epoch <- as.numeric(validdates[[i]])
dates =c(dates , epoch)  
}

# GET PCP (1st FILE TO GET COORDINATES )  USE GEOFIELD  
pcp=lapply(files[1], function(fa) FAdec(fa,"SURFACCPLUIE") )


# GET COORDINATES & RESHAPE 
p     =as.geofield( pcp  , domain=domain_list$be13b_l )
latlon=DomainPoints(p  ,type='lalo')
lons=as.array(latlon[1] )
lats=as.array(latlon[2] )
nx=564
ny=564
y  = as.array( as.numeric(lats[[1]]) )
x  = as.array( as.numeric(lons[[1]]) )
yy = matrix( as.double(y) , nrow=ny , ncol=nx, byrow = T)
xx = matrix( as.double(x) , nrow=ny , ncol=nx, byrow = T)


yy.v <- apply(yy , 2, rev)
yy= yy.v

#-----------------------------------------------------------------



# GET THE PCP FIELDS 
fieldname="SURFACCPLUIE"
pcp_all=lapply(files, function(fa) FAdec(fa,fieldname ) )

pcp_1h=list()
first=TRUE
for(i in seq_along(pcp_all)){
  print(i)
  if (first){
    pcp_1h[[i]]=pcp_all[[i]]
    first=FALSE
  } else {
    pcp_1h[[i]]=pcp_all[[i]]-pcp_all[[i-1]]
    pcp_1h[[i]][pcp_1h[[i]]<0]=0
  }
}

#RADAR (OPTIONAL ) TO USE ONLY WE NEED TO REGRID TO THE RADAR QPE GRID  THE DATA BEFORE WRITTING TO NETCDF 
library(rhdf5)
library(belgium)
radarbase="/mnt/HDS_RADAR_EDP/realtime"
validdates=lapply(files[1]  ,function(fa) attr(FAopen(fa),"time")$validdate)

radar_1h=list()

for (i in seq_along(validdates)){
  date=validdates[[i]]
#  print(date)
  file=paste0(
    radarbase,"/",
    format.Date(date,"%Y"),"/",
    format.Date(date,"%m"),"/",
    format.Date(date,"%d"),"/",
    "bhbjbwdnfa/comp/acrr/qpe2_1h/hdf/",
    format.Date(date,"%Y%m%d%H%M%S"),
    ".rad.bhbjbwdnfa.comp.acrr.qpe2_1h.hdf"
  )
 
  pcp=h5read(file,name='dataset1/data1/data')
  # hdf5 files have a different row/column dominance, so we need to transform the matrix given
  # by h5read
  # at the sametime we create the geofield by providing the qpe domain, which is included in the 
  # belgium library
  radar_1h[[i]]=as.geofield(t(apply(pcp,1,rev)),domain=domain_list$qpe)
}


# REGRID 
pcp_reg=regrid(pcp_1h[[i]], radar_1h[[i]])


# NETCDF OUTFILE 
ncfname=paste("ar13_",year ,lz(month),lz(day),lz(starth),"_native_1h.nc", sep="")

  # META DATA 
endlon=-731900
endlat=731900
delx  =-1300.
dely  =1300.

nlon  <- as.array(seq(-0, endlon ,delx ))
nlat  <- as.array(seq( 0, endlat ,dely ))
times <- as.array(dates)  
ntime <- length( times )

# DIMENSIONS  
londim   <- ncdim_def("y","m"  , as.double(nlon)) 
latdim   <- ncdim_def("x","m"  , as.double(nlat))

timedim  <- ncdim_def("time","seconds since 1970-01-01T00:00:00Z",as.double(times))
dimnchar <- ncdim_def("nchar", "", 1:16, create_dimvar=FALSE )
andim    <- ncdim_def("analysis_time","seconds since 1970-01-01T00:00:00Z",as.double(times[1]))

fillvalue <- 1e32

latname="lat"
lonname="lon"
varname="precipitation"
lat_def =ncvar_def("lat","degrees_north"   , list(latdim,londim)         ,fillvalue,latname ,prec="single")
lon_def =ncvar_def("lon","degrees_east"    , list(latdim,londim)         ,fillvalue,lonname ,prec="single")
pcp_def =ncvar_def("precipitation","kg m-2", list(latdim,londim,timedim) ,fillvalue,varname ,prec="single")
proj_def=ncvar_def("Lambert_Conformal", "",dimnchar , prec="char")
an_def  =ncvar_def("analysis_time","", NULL)# list(andim) ,fillvalue,"analysis_t" ,prec="single")

# PUT VAR DIM 
ncout <- nc_create(ncfname,list( lon_def , lat_def ,pcp_def,  proj_def ) ,force_v4=TRUE)

ncvar_put(ncout,lat_def, as.double(yy))
ncvar_put(ncout,lon_def, as.double(xx))

pcp_array <- array(unlist(pcp_1h), c(ntime, ny , nx ))

ncvar_put(ncout,pcp_def,  pcp_array  )


# LAT LON ATTRIBUTES
ncatt_put(ncout,"lat","long_name","latitude")
ncatt_put(ncout,"lon","long_name","longitude")



# x, y  ATTRIBUTES 
ncatt_put(ncout ,"x","long_name","projection_x_coordinate") ; ncatt_put(ncout ,"y","long_name","projection_y_coordinate")
ncatt_put(ncout ,"x","type","uniform")                      ; ncatt_put(ncout ,"y","type","uniform")
ncatt_put(ncout ,"x","axis","X")                            ; ncatt_put(ncout ,"y","axis","Y")
ncatt_put(ncout ,"x","valid_min" ,0. )                      ; ncatt_put(ncout ,"y","valid_min" ,-731900. )
ncatt_put(ncout ,"x","valid_max",731900. )                  ; ncatt_put(ncout ,"y","valid_max",-0. )
ncatt_put(ncout ,"x","spacing",1300. )                       ; ncatt_put(ncout ,"y","spacing",1300. )
# TIME ATTRIBUTES
ncatt_put(ncout,"time","calendar","standard")
ncatt_put(ncout,"time","standard_name","time")                     
ncatt_put(ncout,"time","axis","T")                           
# PRECIP ATTRIBUTES 
ncatt_put(ncout,"precipitation","units"    ,"kg m-2" )
ncatt_put(ncout,"precipitation","long_name"    ,"Accumulation forecast")
ncatt_put(ncout,"precipitation","coordinates"  ,"lat lon")
ncatt_put(ncout,"precipitation","standard_name","precipitation_amount")
# PROJECTION ATTRIBUTES
ncatt_put(ncout,"Lambert_Conformal","grid_mapping_name","lambert_conformal_conic")
ncatt_put(ncout,"Lambert_Conformal","longitude_of_central_meridian","4.55000000000001")
ncatt_put(ncout,"Lambert_Conformal","latitude_of_projection_origin","50.8")
ncatt_put(ncout,"Lambert_Conformal","standard_parallel","50.8")
ncatt_put(ncout,"Lambert_Conformal","false_easting","365950.")
ncatt_put(ncout,"Lambert_Conformal","false_northing","365950.000000001")


gvars=c("interpolation_method","proj4string","reference_longitude", "reference_latitude","Conventions","creation_date","history","NCO" )

gvals=c( "nearest_neighbour",
         "+proj=lcc +lon_0=4.55 +lat_1=50.8 +lat_2=50.8 +a=6371229 +es=0 +lat_0=50.8 +x_0=365950 +y_0=-365950.000000001",
         "-1.04076249768297",
         "53.9648250227956" ,
         "CF-1.7",
         "2021-10-26T11:40:27Z",
         "Tue Oct 26 14:41:10 2021: ncks -L 5 /scratch/ledecruz/ao13_2021071400_native_1h_clipped.nc /scratch/ledecruz/ar13_2021071400_native_1h_clipped_compressed.nc\nTue Oct 26 13:43:33 2021: ncks -F -d time,25,48 ao13_2021071400_native_1h.nc ar13_2021071400_native_1h_clipped.nc" ,"netCDF Operators version 4.7.9 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)")

for (i in seq_along( gvars)) {
ncatt_put( ncout, 0, gvars[i],gvals[i], prec="text" ) 

}
nc_close(ncout)
q()



