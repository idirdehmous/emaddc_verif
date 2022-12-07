# R-helper functions

# - easy leading zeros
lz=function(x,n=2) formatC(x,flag="0",format="d",width=n)

# - reverse irainbow
revrainbow=function(n){
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        Spectral <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  }
  else Spectral <- colors()[c(107, 619, 12, 105, 384, 383,
                               38, 573, 504, 35, 645)]
  if (n > 0)
      rev(colorRampPalette(Spectral)(n))
  else colorRampPalette(Spectral)(-n)
}


# - equivalent for pretty() but with 0 inside an interval (usefull for bias plots)
lbias=function(x,n=15){
  l=pretty(x,n*2)
  l[seq(1,length(l),2)]
}

# - automated colours for bias plots (to be used in tandem with lbias)
cbias=function(l,rev=FALSE){
 zero=which(l==max(l[l<0]))
 if(rev){
   neg=colorRampPalette(colors=c("red","white"),space="Lab")(zero)
   pos=colorRampPalette(colors=c("white","blue"),space="Lab")(length(l)-zero)
 }
 else{
   neg=colorRampPalette(colors=c("blue","white"),space="Lab")(zero)
   pos=colorRampPalette(colors=c("white","red"),space="Lab")(length(l)-zero)
 }
 c(neg,pos[-1])
}

# - color palette for topography 
ocean.pal <- colorRampPalette(c("#000000","#000209","#000413","#00061E",
                                "#000728","#000932","#002650","#00426E",
                                "#005E8C","#007AAA","#0096C8","#22A9C2",
                                "#45BCBB","#67CFB5","#8AE2AE","#ACF6A8",
                                "#BCF8B9","#CBF9CA","#DBFBDC","#EBFDED"))

land.pal <- colorRampPalette(c("#336600","#F3CA89","#D9A627","#A49019",
                               "#9F7B0D","#996600","#B27676","#C2B0B0",
                               "#E5E5E5","#FFFFFF"))

# - wrapper for nice cloud plots
cloudview=function(cloudiness,geopotential,lsm,nlevels=20,legend=F,title=NULL){
  levels=seq(0,1,length.out=nlevels+1)
  colors=sapply(seq(0,255,length.out=nlevels),function(i) 
                 rgb(red=255,green=255,blue=255,alpha=i,maxColorValue=255))
  iview(geopotential/9.80665,color.palette=land.pal,title="",legend=F,mask=lsm,drawmap=F)
  iview(geopotential,title="",legend=F,mask=!lsm,add=TRUE,drawmap=F)
  if (is.null(title)){
    iview(cloudiness,col=colors,add=TRUE,legend=legend,levels=levels)
  } else {
    iview(cloudiness,col=colors,add=TRUE,legend=legend,levels=levels,title=title)
  } 
}


# - plot nice legend
plotLegend=function(x,y=NULL,legend,col=par("col"),lty,lwd,pch,cex=1,horiz=FALSE,bty="o",oldpar=TRUE,...){
  .oldpar <- par()
  par(mar=c(0,0,0,0))
  plot(0,type="n",xlab="",ylab="",axes=F)
  legend(x=x,y=y,legend=legend,col=col,lty=lty,lwd=lwd,pch=pch,cex=cex,horiz=horiz,bty=bty,...)
  if(oldpar) suppressWarnings(par(.oldpar))
}

# - plot nice colorbar
plotColorbar=function(levels,title,palette=rainbow,cex=1,col=NULL){
 .oldpar <- par()
  ncol = length(levels)-1
  if(is.null(col)) col  = palette(ncol)
  par(mar=c(0,5,0,0))
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F)
  legw = 1
  legylim = c(0.4,0.7)
  xbreaks = seq(0,1,length.out = length(levels))
  rect(xleft = xbreaks[1:ncol], xright=xbreaks[2:(ncol+1)],
       ybottom = rep(legylim[1],ncol), ytop =rep(legylim[2],ncol),col=col)
  text.x = c(0,xbreaks[2:ncol],1)
  text.y = legylim[1]/4
  text(x = 0, y = text.y, labels = format(min(levels),
            digits = 5), cex = cex, adj = c(0, 0))
  x.index <- seq(2, ncol, by = 1)
  text(x = text.x[x.index], y = text.y, labels = format(levels[x.index],
            digits = 2), cex = cex, adj = c(0.5,0))
  text(x = 1, y = text.y, labels = format(max(levels),
            digits = 5), cex = cex, adj = c(1,0))
  text.y = legylim[2] + (1 - legylim[2])/4
  text(x = 0, y = text.y, labels = title,cex=cex,adj=c(0,0))
}



# - read in SAFNWC files
readSAF=function(filebase,col=NULL){
  nc_file=nc_open(paste0(filebase,".nc"))
  if(is.null(col)) h5_file=paste0(filebase,".h5")
  
  lon=ncvar_get(nc_file,'lon')
  lat=ncvar_get(nc_file,'lat')
  dom=Make.domain(projtype='latlong',
                  clonlat=c(median(lon),median(lat)),
                  dxdy=c(diff(lon)[1],diff(lat)[1]),
                  nxny=c(length(lon),length(lat)))
  if (is.null(col)){
    rgb=h5read(h5_file,'01-PALETTE',native=TRUE)
    col=sapply(seq(1,nrow(rgb)),function(i) rgb(red=rgb[i,1],green=rgb[i,2],blue=rgb[i,3],maxColorValue=255))
  }
  
  data=as.geofield(ncvar_get(nc_file,"Band1"),dom)
  return(list(data=data,col=col))
}

# Precipitation type
col_preciptype=c("#00d533","#d00d24","#4e68fa","#001dd8",
                 "#009895","#7f0c8e","#fffa43","#ffc133",
                 "#21fd56","#964047","#001173","#aafeaf",
                 "#f47bf0","#d54ecf","#64fbf3","#b128b5")

col_inca=c("#e7daa9","#bcc17f","#649969","#629aca","#8c71c1","#dd154d")

legend_preciptype=c("Rain","Freezing Rain","Dry snow","Wet snow",
                    "Rain snow mixture","Ice pellets","Graupel/small hail","Hail",
                    "Drizzle","Freezing drizzle","Moist snow/sleet","Inter. rain",
                    "Inter. dry snow","Inter wet snow","Inter rain snow mix","Inter. moist snow/sleet")
legend_inca=c("Rain","Wet snow","Snow","Freezing rain","Hail","Severe Hail")


breaks_preciptype=c(0.5,2,4,5.5,6.5,7.5,8.5,9.5,10.5,11.5,100,196,203,205.5,206.5,208,214)
breaks_inca=seq(0.5,6.5,1)

# Radar reflectivity 
col_radar=c("#c8c8c8", #greys
            "#a0ffff","#66c9fc","#2394fa","#0061f9", #blues
            "#a8fd40","#7ecd31","#539b23","#286a15","#003907", #greens
            "#fffd42","#ffbe32","#ff7e22","#ff3c16","#c1000a", #yellow -> red
            "#ff1cfa") #purple
breaks_radar=c(seq(4,64,4),100)
radview=function(fa,title=NULL,legend=F,...){
  requireNamespace('Rfa')
  if(!is.list(fa)) fa=Rfa::FAopen(fa)
  if(is.null(title)){
    title = sprintf("%s\n%s", "Simulated MAX(dBZ)", format(attr(FAdec(fa,1),
          "info")$time$basedate, "%Y/%m/%d %H:%M"))
    title = paste0(title, " +", attr(FAdec(fa,1), "info")$time$leadtime)
  }
  domain=attr(fa,"domain")
  rad_max=matrix(0,nrow=domain$nx,ncol=domain$ny)
  iField=Rfa::FAfind(fa,"REFLEC_DBZ")
  for( i in iField) rad_max=pmax(rad_max,Rfa::FAdec(fa,i,outform="G"))
  iview(as.geofield(rad_max,domain),col=col_radar,levels=breaks_radar,legend=legend,title=title,...)
}

# - Visibility
col_visi=c("#fe00ff","#9100fd","#1600ff","#025cf6","#00d1fb","#4dfab6",
           "#b5ff46","#f4f223","#d4d66f","#c3c3c3","#e1e1e1","#ffffff")
breaks_visi=c(0,50,100,250,500,1000,2000,3000,4000,6000,8000,10000,20000)

# - Precipitation 
col_pcp=c("#ffffff","#edddd3","#dcbfb5","#c3a197","#beeff9","#a3cff8",
           "#6b9ef2","#5177ed","#31a31b","#52d643","#87f977","#bafdac",
           "#fbfeae","#f5c54a","#f16723","#f02820","#9b1111","#a400b8","#e606fd")
breaks_pcp=c(0,0.2,0.5,1,2,3,4,5,7,10,15,20,25,30,35,40,50,65,80,150)
# - Precipitation type SYNOP codes (see google docs for more info)
drizzle=c(20,51,52,53,54,55,122,150,151,152,153)
freezing_drizzle=c(56,57,154,55,156)
rain=c(21,25,58,59,60,61,62,63,64,65,80,81,82,91,92,95,97,
       123,141,142,143,144,160,161,162,163,181,182,183,184,192,195)
freezing_rain=c(24,66,67,125,164,165,166)
snow=c(22,70,71,72,73,74,75,85,86,124,170,171,172,173,177,185,186,187)
rain_snow=c(23,26,68,69,83,84,167,168)
graupel=c(79,87,88,174,175,176)
hail=c(27,89,90,93,94,96,189,193,196)

# - Plot METEOSAT RGB composites
mist_channels=c("C001_METEOSAT_11","C006_METEOSAT_11","C007_METEOSAT_11")

mistview=function(fa,title=NULL){
#requireNamespace('raster')
requireNamespace('Rfa')
requireNamespace('meteogrid')
# Decode the FA-fields
  if (!is.list(fa)) fa=Rfa::FAopen(fa)
  channels=lapply(mist_channels,function(f) Rfa::FAdec(fa,f))
# Combine and rescale the channels
  mist=array(dim=c(prod(dim(channels[[1]])),3))
  mist[,1]=as.vector(channels[[3]]-channels[[2]])
  mist[,2]=as.vector(channels[[2]]-channels[[1]])
  mist[,3]=as.vector(channels[[2]])
  mist[,1][mist[,1]<(-4)]=-4
  mist[,1][mist[,1]>2]=2
  mist[,2][mist[,2]<0]=0
  mist[,2][mist[,2]>10]=10
  mist[,3][mist[,3]<243]=243
  mist[,3][mist[,3]>293]=293
  mist[,1]=(mist[,1]+4)*42.5
  mist[,2]=mist[,2]*25.5
  mist[,3]=(mist[,3]-243)*5.1
# Create the RGB-colors
  z=rgb(mist[, 1], mist[, 2], mist[, 3], max = 255)
  z <- matrix(z, nrow = nrow(channels[[1]]), ncol = ncol(channels[[1]]),byrow=T)
  z=apply(z,2,rev)
# Get some domain-info
  gdomain <- attr(channels[[1]],"domain")
  glimits <- meteogrid::DomainExtent(gdomain)
  meteogrid::.Last.domain(gdomain)
# Start plotting
  plot(NA,NA,xlim=c(glimits$x0,glimits$x1),ylim=c(glimits$y0,glimits$y1),type = "n",
                xaxs = "i", yaxs = "i",axes=F,xlab="",ylab="",asp=1)
  rasterImage(z, glimits$x0,glimits$y0,glimits$x1,glimits$y1,interpolate=FALSE)
  meteogrid::plot.geodomain( add = TRUE, add.dx = TRUE, box = TRUE,
            lwd = 0.5, col = "black", interior = TRUE,  fill = FALSE, map.database = "world")
# Add title
  if(is.null(title)){
    mytitle = sprintf("%s\n%s", "Meteosat-11 Night Microphysics", format(attr(channels[[1]],
          "info")$time$basedate, "%Y/%m/%d %H:%M"))
    mytitle = paste0(mytitle, " +", attr(channels[[1]], "info")$time$leadtime)
  } else {
    mytitle = title
  }
  title(main=mytitle)
}

# Fast loading of vertical levels
find.index=function(fa,field,levels=seq(1,46),pressure=F){
  if(!is.list(fa)) fa=FAopen(fa)
  list=fa$list$name
  if(pressure){
    flds=paste0("P",stringr::str_sub(lz(levels*100,5),-5,-1),field)
  } else {
    flds=paste0("S",lz(levels,3),sprintf("%-12s",field))
  }
  unlist(lapply(flds,function(x) which(list==x)))
}

fastLoadLevels=function(fa,field,levels=seq(1,46),pressure=F){
  if(!is.list(fa)) fa=FAopen(fa)
  frame=attr(fa,"frame")
  index=find.index(fa,field,levels,pressure)
  fields=lapply(index,function(x) FAdec(fa=fa,field=x,faframe=frame,outform="G"))
  return(aperm(array(unlist(fields),dim=c(dim(fields[[1]]),length(fields))),c(3,1,2)))
}












